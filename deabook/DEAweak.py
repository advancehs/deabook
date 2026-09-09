# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint
import numpy as np
import pandas as pd
from .constant import (
    CET_ADDI, RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_DEFAULT, OPT_LOCAL,
    CONTEMPORARY, GLOBAL, SEQUENTIAL, WINDOW,
)
from .utils import tools
from .DEA import _PerDmuEnvelopmentModel, _SBM, _SDDF
from .core import (
    PyomoResultMixin,
    build_weak_direction_flags,
    format_ref_indexed_results,
    format_var_indexed_results,
    parse_sent,
)
from .mqcore import (
    MalmquistModel,
    build_distance_frame,
    two_factor,
    normalize_tech,
    _shift,
)


class _WeakOrientationMixin:
    """Shared orientation helpers for weak-disposability DEA classes."""

    def _set_weak_orientation_flags(self):
        flags = build_weak_direction_flags(self.gx, self.gy, self.gb)
        self.direction_flags = flags
        self.input_oriented = flags.input_oriented
        self.output_oriented = flags.output_oriented
        self.unoutput_oriented = flags.unoutput_oriented
        self.undesirable_oriented = flags.unoutput_oriented
        self.hyper_orientedyx = flags.hyper_orientedyx
        self.hyper_orientedyb = flags.hyper_orientedyb
        self.hyper_orientedxb = flags.hyper_orientedxb
        self.hyper_orientedyxb = flags.hyper_orientedyxb
        self.hyper_oriented = flags.is_hyper
        return flags


class _WeakPerDmuResultMixin:
    """Reusable result display/get/info methods for weak per-DMU models."""

    def _models(self):
        preferred = f"_{self.__class__.__name__}__modeldict"
        if hasattr(self, preferred):
            return getattr(self, preferred)
        for name in dir(self):
            if name.endswith("__modeldict"):
                return getattr(self, name)
        raise AttributeError("This object does not expose a model dictionary.")

    def _format_ref_results(self, values):
        return format_ref_indexed_results(values, self.reference_data_index)

    def _format_var_results(self, values, columns=None):
        return format_var_indexed_results(values, columns)

    def _solve_per_dmu_models(self, email=OPT_LOCAL, solver=OPT_DEFAULT, include_objective=False):
        """Solve one weak-disposability Pyomo model per evaluated DMU.

        ``DEAweak`` and its subclasses share this whole solver/result-collection
        lifecycle.
        Keeping that lifecycle here avoids copy-pasted loops while preserving
        the historical ``self.results`` dictionary and returned column names.
        """
        self.results = {}
        all_optimization_statuses = {}
        rho_values, theta_values, objective_values = {}, {}, {}
        lamda_values, mu_values = {}, {}
        use_neos = tools.set_neos_email(email)

        for actual_index, model in self._models().items():
            try:
                optimization_status = tools.optimize_model2(
                    model, actual_index, use_neos, CET_ADDI, solver=solver
                )
            except Exception as e:
                print(f"Error optimizing DMU {actual_index}: {e}")
                optimization_status = "Error"

            all_optimization_statuses[actual_index] = optimization_status
            if optimization_status in ["ok"]:
                try:
                    rho_values[actual_index] = model.rho.value
                    if self.rts == RTS_VRS1:
                        theta_values[actual_index] = model.theta.value
                    lamda_values[actual_index] = self._collect_ref_component(model, "lamda")
                    if include_objective:
                        objective_values[actual_index] = model.objective()
                    if self.rts == RTS_VRS2:
                        mu_values[actual_index] = self._collect_ref_component(model, "mu")
                except Exception as e:
                    print(
                        "Warning: Could not retrieve results for DMU "
                        f"{actual_index} despite status '{optimization_status}': {e}"
                    )
                    rho_values[actual_index] = np.nan
                    if self.rts == RTS_VRS1:
                        theta_values[actual_index] = np.nan
                    lamda_values[actual_index] = None
                    if include_objective:
                        objective_values[actual_index] = np.nan
                    if self.rts == RTS_VRS2:
                        mu_values[actual_index] = None
            else:
                print(f"Error optimizing DMU {actual_index}, optimization_status :{optimization_status}")
                rho_values[actual_index] = np.nan
                if self.rts == RTS_VRS1:
                    theta_values[actual_index] = np.nan
                lamda_values[actual_index] = None
                if include_objective:
                    objective_values[actual_index] = np.nan
                if self.rts == RTS_VRS2:
                    mu_values[actual_index] = None

        self.results["optimization_status"] = all_optimization_statuses
        self.results["rho"] = rho_values
        if self.rts == RTS_VRS1:
            self.results["theta"] = theta_values
        self.results["lamda_df"] = self._format_ref_results(lamda_values)
        if include_objective:
            self.results["objective_value"] = objective_values
        if self.rts == RTS_VRS2:
            self.results["mu_df"] = self._format_ref_results(mu_values)

        data = {
            "optimization_status": pd.Series(all_optimization_statuses),
            "rho": pd.Series(rho_values),
        }
        if include_objective:
            data["objective_value"] = pd.Series(objective_values)
        results_df = pd.DataFrame(data)
        self._add_efficiency_columns(results_df)
        return results_df

    def _collect_ref_component(self, model, component_name):
        component = getattr(model, component_name)
        return {
            self.reference_data_index[r_range]: component[r_range].value
            for r_range in model.R
        }

    def _safe_inverse_one_plus(self, series):
        return series.apply(lambda x: 1 / (1 + x) if pd.notna(x) else np.nan)

    def _safe_inverse(self, series):
        return series.apply(lambda x: 1 / x if pd.notna(x) and x != 0 else np.nan)

    def _safe_sqrt(self, series):
        return series.apply(lambda x: np.sqrt(x) if pd.notna(x) else np.nan)

    def _one_minus(self, series):
        return series.apply(lambda x: (1 - x) if pd.notna(x) else np.nan)

    def display_status(self):
        if not self.results or 'optimization_status' not in self.results:
            print("Optimization has not been run yet.")
            return
        print(pd.Series(self.results['optimization_status']))

    def display_rho(self):
        tools.assert_optimized(self.results)
        print(self.get_rho())

    def display_theta(self):
        tools.assert_optimized(self.results)
        theta = self.get_theta()
        print(theta if theta is not None else "Theta variable is only available for RTS_VRS1.")

    def display_lamda(self):
        tools.assert_optimized(self.results)
        print(self.get_lamda())

    def display_mu(self):
        tools.assert_optimized(self.results)
        mu = self.get_mu()
        print(mu if mu is not None else "Mu variables are only available for RTS_VRS2.")

    def get_status(self):
        return self.results.get('optimization_status', {})

    def get_rho(self):
        tools.assert_optimized(self.results)
        return pd.Series(self.results.get('rho', {}), dtype='float64')

    def get_theta(self):
        tools.assert_optimized(self.results)
        if self.rts != RTS_VRS1 or 'theta' not in self.results:
            return None
        return pd.Series(self.results.get('theta', {}), dtype='float64')

    def get_lamda(self):
        tools.assert_optimized(self.results)
        return self.results.get('lamda_df')

    def get_mu(self):
        tools.assert_optimized(self.results)
        if self.rts != RTS_VRS2 or 'mu_df' not in self.results:
            return None
        return self.results.get('mu_df')

    def info(self, dmu="all"):
        modeldict = self._models()
        if dmu == "all":
            for ind, problem in modeldict.items():
                print(f"\n--- Model for DMU: {ind} ---")
                problem.pprint()
            return None
        dmu_list = [dmu] if isinstance(dmu, (str, int, float)) else list(dmu)
        for ind in dmu_list:
            key = ind
            if key not in modeldict:
                try:
                    key = int(ind)
                except (TypeError, ValueError):
                    pass
            if key in modeldict:
                print(f"\n--- Model for DMU: {key} ---")
                modeldict[key].pprint()
            else:
                print(f"DMU '{ind}' not found in the evaluated set.")
        return None

class DEAweak(_WeakOrientationMixin, _WeakPerDmuResultMixin, _PerDmuEnvelopmentModel):
    """Data Envelopment Analysis (DEA) - Solves per DMU"""

    def __init__(self,data,sent,gy=[1],gx=[0],gb=[1],rts=RTS_VRS1,evalindex=None, refindex=None):
        """DEA: Envelopment problem, solving for each DMU individually.

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y:CO2"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            evalindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEAweak should return numpy arrays for x, y,b, xref, yref,bref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DEAweak is modified or a similar function returns these:
        try:
             self.y, self.x,self.b, self.yref,self.xref,self.bref, self.gy, self.gx,self.gb, \
                 self.evaluated_data_index, self.reference_data_index = tools.assert_DEAweak_with_indices(
                     data, sent, gy, gx,gb, evalindex, refindex
                 )
        except AttributeError:
             # Fallback if assert_DEA_with_indices doesn't exist, assuming original assert_DEA
             # In this case, evaluated_data_index and reference_data_index will be ranges 0..N-1
             print("Warning: tools.assert_DEAweak_with_indices not found. Using range indices.")
             self.y,self.x,self.b, self.yref,self.xref,self.bref, self.gy,self.gx,self.gb = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, evalindex, refindex)
             self.evaluated_data_index = list(range(len(self.y)))
             self.reference_data_index = list(range(len(self.yref)))

        self.rts = rts

        self.evaluated_indices_range = range(len(self.y)) # Range indices for internal numpy arrays
        self.reference_indices_range = range(len(self.yref)) # Range indices for internal numpy arrays
        self.num_inputs = len(self.x[0])
        self.num_outputs = len(self.y[0])
        self.num_unoutputs = len(self.b[0])
        # print(self.gx,self.gy,")))))))))))))")
        # Determine orientation based on gx/gy/gb vectors.
        self._set_weak_orientation_flags()
        if (self.input_oriented and self.output_oriented)   \
            or (self.input_oriented and self.unoutput_oriented) or (self.output_oriented and self.unoutput_oriented):
             # Original code structure suggests one or the other dominates based on the if/elif.
             # Assuming standard DEA where it's one or the other.
             pass # This case implies a potential issue if both sums >= 1, but following original logic.


        # Dictionary to store individual models for each evaluated DMU
        self.__modeldict = {}
        # Dictionary to store results after optimization
        self.results = {}

        # Loop through each DMU to be evaluated (using range index to access numpy arrays)
        for i_range in self.evaluated_indices_range:
            # Get the actual index label for this DMU
            actual_index = self.evaluated_data_index[i_range]

            # Create a new model for this specific DMU
            model = ConcreteModel()

            # Define sets for this specific model
            # R is the set of reference DMU range indices
            model.R = Set(initialize=self.reference_indices_range)
            # J is the set of input indices
            model.J = Set(initialize=range(self.num_inputs))
            # K is the set of output indices
            model.K = Set(initialize=range(self.num_outputs))
            # L is the set of unoutput indices
            model.L = Set(initialize=range(self.num_unoutputs))
            # Define variables for this specific model
            # rho is the efficiency score for the current DMU (not indexed by DMU within this model)
            # Bounds based on standard DEA orientation (matching DEAt)
            if self.input_oriented:  # Minimize rho, efficiency <= 1
                if self.rts == RTS_CRS:
                    model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')
                elif self.rts == RTS_VRS1:
                    model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')
                    model.theta = Var(bounds=(0.0, 1.0),doc='emission reduction factor')
                elif self.rts == RTS_VRS2:
                    model.rho = Var(bounds=(0, None), doc=f'efficiency for DMU {actual_index}')
                    model.mu = Var(model.R, bounds=(0.0, None),doc='emission reduction factor2')
            elif self.unoutput_oriented:
                if self.rts == RTS_CRS:
                    model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')
                elif self.rts == RTS_VRS1:
                    model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')
                    model.theta = Var(bounds=(0.0, 1.0),doc='emission reduction factor')
                elif self.rts == RTS_VRS2:
                    model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')
                    model.mu = Var(model.R, bounds=(0.0, None),doc='emission reduction factor2')
            elif self.output_oriented: # Maximize rho, efficiency >= 1
                if self.rts == RTS_CRS:
                    model.rho = Var(bounds=(1, None), doc=f'efficiency for DMU {actual_index}')
                elif self.rts == RTS_VRS1:
                    model.rho = Var(bounds=(1, None), doc=f'efficiency for DMU {actual_index}')
                    model.theta = Var(bounds=(0.0, 1.0),doc='emission reduction factor')
                elif self.rts == RTS_VRS2:
                    model.rho = Var(bounds=(1, None), doc=f'efficiency for DMU {actual_index}')
                    model.mu = Var(model.R, bounds=(0.0, None),doc='emission reduction factor2')
            elif self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb: # Hyper orientation, rho can be any value, 但正常<1
                if self.rts == RTS_CRS: # Hyper orientation with CRS, rho efficiency <= 1
                    model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')  
                    model.theta = Var(bounds=(0.0, 1.0),doc='emission reduction factor')
                elif self.rts == RTS_VRS1:
                    model.rho = Var(bounds=(0, None), doc=f'efficiency for DMU {actual_index}')
                    model.theta = Var(bounds=(0.0, 1.0),doc='emission reduction factor')
                elif self.rts == RTS_VRS2:
                    model.rho = Var(bounds=(0, None), doc=f'efficiency for DMU {actual_index}')
                    model.mu = Var(model.R, bounds=(0.0, None),doc='emission reduction factor2')
            else:
                 # This else should not be reached due to check above, but for safety
                 raise ValueError("Invalid orientation configuration.")
            # lamda are intensity variables, indexed by the reference set
            model.lamda = Var(model.R, bounds=(0.0, None), doc='intensity variables')

            # Setup the objective function and constraints for THIS DMU
            # Objective: Minimize/Maximize the SINGLE rho variable in this model
            if self.input_oriented or self.unoutput_oriented:
                model.objective = Objective(rule=lambda m: m.rho, sense=minimize, doc='objective function')
            elif self.output_oriented:
                model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc='objective function')
            elif self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb:
                if self.rts == RTS_CRS:
                    model.objective = Objective(rule=lambda m: m.rho, sense=minimize, doc='objective function')
                elif self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                    # Hyper orientation with VRS, maximize rho (DDF-style beta)
                    model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc='objective function')
                else:
                    raise ValueError("Invalid RTS configuration for hyper orientation.")
            else:
                raise ValueError("Invalid orientation configuration.")

            # Input constraint for THIS DMU, referencing current DMU's data (self.x[i_range])
            # Use a factory function to capture the current DMU index
            input_rule_factory = self.__create_input_rule(i_range)
            model.input = Constraint(model.J, rule=input_rule_factory, doc='input constraint')

            # Output constraint for THIS DMU, referencing current DMU's data (self.y[i_range])
            output_rule_factory = self.__create_output_rule(i_range)
            model.output = Constraint(model.K, rule=output_rule_factory, doc='output constraint')

            # UNOutput constraint for THIS DMU, referencing current DMU's data (self.b[i_range])
            unoutput_rule_factory = self.__create_unoutput_rule(i_range)
            model.unoutput = Constraint(model.L, rule=unoutput_rule_factory, doc='unoutput constraint')
            # VRS constraint for THIS DMU
            if self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                vrs_rule_factory = self.__create_vrs_rule() # VRS rule doesn't need current DMU index
                model.vrs = Constraint(rule=vrs_rule_factory, doc='variable return to scale rule')

            # Store the created model in the dictionary using the actual DMU index as the key
            self.__modeldict[actual_index] = model

        # Optimization status will be stored per DMU after optimize is called
        # self.optimization_status = 0 # Not needed as overall status
        # self.problem_status = 0 # Not needed as overall status


    # --- Rule Factories (Helper methods to create constraint rules) ---
    # These capture the data for the specific DMU being modeled

    def __create_input_rule(self, current_dmu_range_index):
        """Factory for creating the input constraint rule for a specific DMU."""

        if self.input_oriented:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Access current DMU's input data self.x[current_dmu_range_index][j]
                    # Access reference DMUs' input data self.xref[r][j]
                    if self.gx[j] == 1:
                        # Input is scaled by rho
                        return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= model.rho * self.x[current_dmu_range_index][j]
                    else:
                        # Input is not scaled by rho
                        return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Access current DMU's input data self.x[current_dmu_range_index][j]
                    # Access reference DMUs' input data self.xref[r][j]
                    if self.gx[j] == 1:
                        # Input is scaled by rho
                        return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= model.rho * self.x[current_dmu_range_index][j]
                    else:
                        # Input is not scaled by rho
                        return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
        elif self.output_oriented or self.unoutput_oriented:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                def input_rule(model, j):
                    # Input is scaled by theta in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=model.theta*self.x[current_dmu_range_index][j] 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum((model.lamda[r]+model.mu[r])* self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j]
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Hyper orientation: Input is scaled by rho if gx[j] == 1
                    if self.gx[j] == 1:
                        return sum(model.lamda[r]*self.xref[r][j] for r in model.R) <= model.rho*self.x[current_dmu_range_index][j]
                    else:
                        return sum(model.lamda[r]*self.xref[r][j] for r in model.R) <= model.theta*self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rho in input orientation
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rho*self.gx[j]*self.x[current_dmu_range_index][j]
 
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Hyper orientation: Input is scaled by rho if gx[j] == 1
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= model.theta * self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                def input_rule(model, j):
                    # Input is scaled by rho in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= model.theta * self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum((model.lamda[r]+model.mu[r])* self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j]
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Hyper xb: Input is scaled by rho if gx[j]==1, else theta
                    if self.gx[j] == 1:
                        return sum(model.lamda[r]*self.xref[r][j] for r in model.R) <= model.rho*self.x[current_dmu_range_index][j]
                    else:
                        return sum(model.lamda[r]*self.xref[r][j] for r in model.R) <= model.theta*self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperxb orientation.")
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rho via DDF-style
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rho*self.gx[j]*self.x[current_dmu_range_index][j]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Hyper yxb: Input is scaled by rho if gx[j]==1, else theta
                    if self.gx[j] == 1:
                        return sum(model.lamda[r]*self.xref[r][j] for r in model.R) <= model.rho*self.x[current_dmu_range_index][j]
                    else:
                        return sum(model.lamda[r]*self.xref[r][j] for r in model.R) <= model.theta*self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperyxb orientation.")
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rho via DDF-style
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rho*self.gx[j]*self.x[current_dmu_range_index][j]

        else:
            # Should not be reached
            def input_rule(model, j):
                return Constraint.Skip # Or raise error
        # Return the rule function
        return input_rule

    def __create_output_rule(self, current_dmu_range_index):
        """Factory for creating the output constraint rule for a specific DMU."""
        if self.input_oriented:
            if self.rts == RTS_CRS  or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.output_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is scaled by rho in output orientation (original code applies rho to all outputs if sum(gy)>=1)
                    if self.gy[k] == 1:
                        return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                            model.rho * self.y[current_dmu_range_index][k]
                    else:
                        # Output is not scaled by rho in output orientation
                        return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.unoutput_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rho if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 

            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rho if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper xb: Output is not scaled (gy==0)
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperxb orientation.")
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled (gy==0)
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper yxb: Output is not scaled by rho in CRS hyperbolic
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperyxb orientation.")
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is scaled by rho via DDF-style
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
        else:
             # Should not be reached
             def output_rule(model, k):
                 return Constraint.Skip # Or raise error
        return output_rule

    def __create_unoutput_rule(self, current_dmu_range_index):
        """Factory for creating the unoutput constraint rule for a specific DMU."""
        if self.input_oriented :
            if self.rts == RTS_CRS or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # unOutput is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.output_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # unOutput is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l]
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.unoutput_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Output is scaled by rho in output orientation (original code applies rho to all outputs if sum(gy)>=1)
                    if self.gb[l] == 1:
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == model.rho * self.b[current_dmu_range_index][l]
                    else:
                        # Output is not scaled by rho in output orientation
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l]
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rho if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==model.theta* self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:  
                def unoutput_rule(model, l):
                    # Output is not scaled by rho in input orientation
                    return sum((model.lamda[r])*self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l] 
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rho if gb[k] == 1
                    if self.gb[l] == 1:
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==model.rho* self.b[current_dmu_range_index][l]
                    else:
                        # Output is not scaled by rho in output orientation
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == model.theta*self.b[current_dmu_range_index][l] 
            elif self.rts == RTS_VRS1:
                def unoutput_rule(model, l):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r]*self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper xb: Non-desirable output is scaled by rho if gb[l]==1, else theta
                    if self.gb[l] == 1:
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == model.rho * self.b[current_dmu_range_index][l]
                    else:
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == model.theta * self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperxb orientation.")
            elif self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Non-desirable output is scaled by rho via DDF-style
                    return sum(model.lamda[r]*self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper yxb: Non-desirable output is scaled by rho if gb[l]==1, else theta
                    if self.gb[l] == 1:
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == model.rho * self.b[current_dmu_range_index][l]
                    else:
                        return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == model.theta * self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperyxb orientation.")
            elif self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Non-desirable output is scaled by rho via DDF-style
                    return sum(model.lamda[r]*self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
        else:
             # Should not be reached
             def unoutput_rule(model, l):
                 return Constraint.Skip # Or raise error
        return unoutput_rule

    def __create_vrs_rule(self):
        """Factory for creating the VRS constraint rule."""
        if self.rts==RTS_VRS1:
            def vrs_rule(model):
                return sum(model.lamda[r] for r in model.R) == model.theta
            return vrs_rule
        elif self.rts==RTS_VRS2:
            def vrs_rule(model):
                return sum((model.lamda[r]+model.mu[r]) for r in model.R) == 1
            return vrs_rule





    # --- Optimization and Results Methods ---

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize every evaluated DMU and return legacy DEAweak columns."""
        return self._solve_per_dmu_models(email=email, solver=solver, include_objective=False)

    def _add_efficiency_columns(self, results_df):
        rho = results_df['rho']
        if self.input_oriented or self.unoutput_oriented:
            results_df['te'] = rho
        elif self.output_oriented:
            results_df['te'] = self._safe_inverse(rho)
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                results_df['te'] = self._safe_sqrt(rho)
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).")
            elif self.rts == RTS_VRS2:
                results_df['tei'] = self._one_minus(rho)
                results_df['teo'] = self._safe_inverse_one_plus(rho)
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS:
                results_df['te'] = self._safe_sqrt(rho)
            elif self.rts in (RTS_VRS1, RTS_VRS2):
                results_df['teuo'] = self._one_minus(rho)
                results_df['teo'] = self._safe_inverse_one_plus(rho)
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                results_df['te'] = self._safe_sqrt(rho)
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperxb orientation.")
            elif self.rts == RTS_VRS2:
                results_df['tei'] = self._one_minus(rho)
                results_df['teuo'] = self._one_minus(rho)
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                results_df['te'] = self._safe_sqrt(rho)
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for hyperyxb orientation.")
            elif self.rts == RTS_VRS2:
                results_df['tei'] = self._one_minus(rho)
                results_df['teuo'] = self._one_minus(rho)
                results_df['teo'] = self._safe_inverse_one_plus(rho)
        return results_df


class DEAweakRho(DEAweak):
    """Radial weak-disposability DEA that reports only the raw LP value ``rho``
    (plus ``objective_value``), mirroring ``DEA.ddf`` / ``DEAweak.ddf`` /
    ``DEAweak.nddf``: efficiency is derived by the caller from ``rho`` at the
    call site.

    For the default ``ushephard`` direction (``gx=0, gy=0, gb=1``) the radial
    efficiency is ``PEE = rho``.  For the default ``hyper`` direction
    (``gx=0, gy=1, gb=1``): CRS gives ``MEE = sqrt(rho)``; VRS1/VRS2 give
    ``MEE_BO = 1 - rho`` and ``MEE_O = 1/(1 + rho)``.
    """

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        return self._solve_per_dmu_models(email=email, solver=solver,
                                          include_objective=True)

    def _add_efficiency_columns(self, results_df):
        return results_df

# ---------------------------------------------------------------------------
# Module-level convenience interfaces
#
# The ``DEAweak`` class uses a radial distance value.  The functions below expose a
# distance-function-oriented API: each function carries a sensible default
# direction vector so that the canonical orientation no longer needs to be
# specified by hand.  They are reached as ``DEAweak.ushephard(...)`` etc. when
# the module is imported via ``from deabook import DEAweak``.
# ---------------------------------------------------------------------------


def _weak_var_counts(sent):
    """Parse a ``sent`` string and return ``(n_inputs, n_outputs, n_unoutputs)``.

    ``sent`` follows ``"inputs = outputs : unoutputs"``, e.g. ``"K L E=Y:CO2"``.
    """
    parts = sent.split('=')
    if len(parts) != 2:
        raise ValueError(
            "Invalid 'sent' format. Expected 'inputs=outputs:unoutputs'")
    inputvars = [v for v in parts[0].strip().split() if v]
    out_parts = parts[1].strip().split(':')
    if len(out_parts) != 2:
        raise ValueError(
            "Invalid 'sent' format. Expected 'outputs:unoutputs' after '='")
    outputvars = [v for v in out_parts[0].strip().split() if v]
    unoutputvars = [v for v in out_parts[1].strip().split() if v]
    return len(inputvars), len(outputvars), len(unoutputvars)


def ushephard(data, sent, rts=RTS_CRS, evalindex=None, refindex=None):
    """非期望产出方向的谢泼德距离函数（Shephard distance function for undesirable outputs）。

    该接口方向固定为非期望产出导向（``gx=0, gy=0, gb=1``：投入与期望产出保持不变，
    仅压缩非期望产出），对应纯环境效率（PEE），**不支持**设定方向向量。与
    ``DEA.ddf`` 一致，本接口只返回线性规划的目标值 ``rho``（外加 ``objective_value``），
    即纯环境效率 ``PEE = rho``；不生成 ``te`` 等效率分解列。

    Args:
        data (pandas.DataFrame): 输入数据。
        sent (str): ``"投入=期望产出:非期望产出"``，如 ``"K L E=Y:CO2"``。
        rts: 规模报酬与弱可处置性设定，默认 ``RTS_CRS``。
        evalindex, refindex (str, optional): 评价/参考样本筛选条件。
    """
    n_in, n_out, n_bad = _weak_var_counts(sent)
    gx = [0] * n_in
    gy = [0] * n_out
    gb = [1] * n_bad
    return DEAweakRho(data, sent, gy=gy, gx=gx, gb=gb, rts=rts,
                      evalindex=evalindex, refindex=refindex)


def hyper(data, sent, gx=None, gy=None, gb=None, rts=RTS_CRS,
          evalindex=None, refindex=None):
    """包含非期望产出的双曲距离函数（hyperbolic distance function with undesirable outputs）。

    默认方向 ``gx=0, gy=1, gb=1``：同时扩张期望产出并压缩非期望产出，对应混合环境效率。
    与 ``DEA.ddf`` 一致，本接口只返回线性规划的目标值 ``rho``（外加 ``objective_value``），
    不生成 ``te``/``teuo``/``teo`` 等效率分解列；混合环境效率由调用方据方向在调用处换算：
    CRS 下 ``MEE = sqrt(rho)``，VRS1/VRS2 下 ``MEE_BO = 1 - rho``、``MEE_O = 1/(1 + rho)``。

    Args:
        data (pandas.DataFrame): 输入数据。
        sent (str): ``"投入=期望产出:非期望产出"``，如 ``"K L E=Y:CO2"``。
        gx, gy, gb (list, optional): 方向向量，不指定时取上述混合环境效率默认值。
        rts: 规模报酬与弱可处置性设定，默认 ``RTS_CRS``。
        evalindex, refindex (str, optional): 评价/参考样本筛选条件。
    """
    n_in, n_out, n_bad = _weak_var_counts(sent)
    if gx is None:
        gx = [0] * n_in
    if gy is None:
        gy = [1] * n_out
    if gb is None:
        gb = [1] * n_bad
    return DEAweakRho(data, sent, gy=gy, gx=gx, gb=gb, rts=rts,
                      evalindex=evalindex, refindex=refindex)


def ddf(data, sent, gx=None, gy=None, gb=None, rts=RTS_CRS,
        evalindex=None, refindex=None):
    """将方向距离函数（DDF）扩展到包含非期望产出。

    默认方向 ``gx=0, gy=1, gb=1``：沿期望产出扩张、沿非期望产出压缩。与
    ``DEA.ddf`` 一致，本接口只返回线性规划的目标值 ``rho``（即方向距离函数
    值，外加 ``objective_value``），**不**生成 ``teuo``/``teo`` 等效率分解列；
    环境效率 ``EB`` 与产出效率 ``EO`` 由调用方据方向在调用处换算
    （``EB = 1 - rho``、``EO = 1 / (1 + rho)``）。

    Args:
        data (pandas.DataFrame): 输入数据。
        sent (str): ``"投入=期望产出:非期望产出"``，如 ``"K L E=Y:CO2"``。
        gx, gy, gb (list, optional): 方向向量，不指定时取上述期望产出—非期望产出方向默认值。
        rts: 规模报酬与弱可处置性设定，默认 ``RTS_CRS``。
        evalindex, refindex (str, optional): 评价/参考样本筛选条件。
    """
    n_in, n_out, n_bad = _weak_var_counts(sent)
    if gx is None:
        gx = [0] * n_in
    if gy is None:
        gy = [1] * n_out
    if gb is None:
        gb = [1] * n_bad
    return DDFweakRho(data, sent, gy=gy, gx=gx, gb=gb, rts=rts,
                      evalindex=evalindex, refindex=refindex)


def nddf(data, sent, gx=None, gy=None, gb=None, rts=RTS_CRS,
         evalindex=None, refindex=None):
    """非径向方向距离函数（nonradial DDF），各变量独立调整。

    默认方向 ``gx=0, gy=1, gb=1``。与 ``DEA.ddf`` 一致，本接口只返回线性规划
    的目标值 ``objective_value``（即各变量 ``beta`` 的加权和），以及每个变量的
    缩放因子 ``beta``，命名规则为 ``beta_<变量名>``（默认方向下为 ``beta_Y``、
    ``beta_CO2``）。本接口**不**生成 ``teY``/``teCO2``/``teuo``/``teo``/``teuo2o``
    等效率分解列；环境效率 ``EB``、产出效率 ``EO`` 与相对排放强度 ``EBO`` 由调用方
    据方向在调用处换算（标准等权设定 ``wy = wb = 0.5`` 下：
    ``EB = 0.5*(1 - beta_CO2)``、``EO = 0.5/(1 + beta_Y)``、``EBO = EB / EO``）。

    Args:
        data (pandas.DataFrame): 输入数据。
        sent (str): ``"投入=期望产出:非期望产出"``，如 ``"K L E=Y:CO2"``。
        gx, gy, gb (list, optional): 方向向量，不指定时取上述期望产出—非期望产出方向默认值。
        rts: 规模报酬与弱可处置性设定，默认 ``RTS_CRS``。
        evalindex, refindex (str, optional): 评价/参考样本筛选条件。
    """
    n_in, n_out, n_bad = _weak_var_counts(sent)
    if gx is None:
        gx = [0] * n_in
    if gy is None:
        gy = [1] * n_out
    if gb is None:
        gb = [1] * n_bad
    return NDDFweakRho(data, sent, gy=gy, gx=gx, gb=gb, rts=rts,
                       evalindex=evalindex, refindex=refindex)


# ---------------------------------------------------------------------------
# Dynamic (panel) environmental productivity indices — Chapter 7.
#
# These follow the Chapter-4 ``DEA.mq`` engine (deabook.mqcore) but plug in the
# weak-disposability distance primitives above.  They deliberately do NOT depend
# on the legacy ``MQDEAweak`` module.  All three share the same ``tech``
# reference-set semantics (Chapter 3 §3.4) and a common ``.optimize()``
# lifecycle, exactly like ``DEA.mq`` / ``DEA.mqddf``:
#
#   * tech=CONTEMPORARY (default) / SEQUENTIAL / WINDOW / GLOBAL
#   * .optimize(email, solver) returns a DataFrame indexed like the input data
#
#   * dep     -- dynamic environmental performance index (MEPI), built on
#                :func:`ushephard` (undesirable-output-oriented Shephard distance,
#                §7.4.1).  Ratio decomposition: MQ / MEFFCH / MTECHCH.
#   * depddf  -- DDF-based dynamic environmental-performance index (ML, §7.4.2),
#                built on the bad-output-oriented DDF (:func:`ddf` with direction
#                (0, 0, 1)).  Ratio decomposition: MQ / MEFFCH / MTECHCH.
#   * lddf    -- additive Luenberger indicator, built on the single-period
#                Luenberger DDF (:func:`_ldf_distance`).  Additive decomposition:
#                L / LEFFCH / LTECHCH (benchmark 0).  ``orient`` = "oo" / "io" / "both"
#                (§7.2.1 / §7.2.2 / §7.2.3).
# ---------------------------------------------------------------------------

# --- single (eval, ref) distance solvers -> one-column DataFrame["E"] ---------
def _ushephard_distance(data, sent, rts, evalindex, refindex, email, solver):
    """MEPI distance: undesirable-output-oriented Shephard distance = ushephard rho."""
    res = ushephard(data, sent, rts=rts, evalindex=evalindex, refindex=refindex).optimize(email, solver)
    out = pd.DataFrame(index=res.index)
    out["E"] = res["rho"].astype(float)
    return out


def _ddf_distance(data, sent, gx, gy, gb, rts, orient, evalindex, refindex, email, solver):
    """ML distance: D = 1 + rho (oo, §7.2.1) or D = 1 - rho (io, §7.2.2) from ddf."""
    res = ddf(data, sent, gx=gx, gy=gy, gb=gb, rts=rts, evalindex=evalindex, refindex=refindex).optimize(email, solver)
    rho = res["rho"].astype(float)
    out = pd.DataFrame(index=res.index)
    out["E"] = (1.0 + rho) if orient == "oo" else (1.0 - rho)
    return out


def _ddf_bad_distance(data, sent, rts, evalindex, refindex, email, solver):
    """DDF-based dynamic environmental-performance distance (§7.4.2).

    Direction is fixed to the undesirable-output contraction
    ``(gx, gy, gb) = (0, 0, 1)`` -- the directional analogue of
    :func:`ushephard` -- so that the returned term ``E = 1 - rho`` equals
    ``1 - D_b(vec)`` entering the Malmquist-Luenberger index of eq. (7-51).
    """
    n_in, n_out, n_bad = _weak_var_counts(sent)
    gx, gy, gb = [0] * n_in, [0] * n_out, [1] * n_bad
    res = ddf(data, sent, gx=gx, gy=gy, gb=gb, rts=rts, evalindex=evalindex, refindex=refindex).optimize(email, solver)
    rho = res["rho"].astype(float)
    out = pd.DataFrame(index=res.index)
    out["E"] = 1.0 - rho
    return out


def _ldf_distance(data, sent, orient="oo", rts=RTS_CRS,
                 evalindex=None, refindex=None, email=OPT_LOCAL, solver=OPT_DEFAULT):
    """Single-period Luenberger DDF solve: returns a one-column DataFrame["E"]
    whose value is the directional distance :math:`\\overrightarrow{D}` for the
    chosen orientation (§7.2.1 / §7.2.2 / §7.2.3).  This is the internal building
    block for the dynamic :func:`lddf` indicator; the public ``lddf`` name is taken
    by that dynamic indicator, so the primitive lives here under ``_ldf_distance``.
    """
    if orient not in ("io", "oo", "both"):
        raise ValueError(
            "orient must be 'oo' (output, §7.2.1), 'io' (input, §7.2.2) "
            "or 'both' (graph/comprehensive, §7.2.3 / §7.4)."
        )
    n_in, n_out, n_bad = _weak_var_counts(sent)
    if orient == "io":
        gx, gy, gb = [1] * n_in, [0] * n_out, [0] * n_bad
    elif orient == "oo":
        gx, gy, gb = [0] * n_in, [1] * n_out, [1] * n_bad
    else:  # both
        gx, gy, gb = [1] * n_in, [1] * n_out, [1] * n_bad
    if orient in ("io", "both") and rts == RTS_VRS1:
        raise ValueError(
            "rts=RTS_VRS1 is not supported for the input/both orientation of the "
            "Luenberger DDF: the directional distance function with an input direction "
            "is not linearizable under VRS1. Use RTS_CRS or RTS_VRS2 instead."
        )
    res = DDFweakRho(data, sent, gy=gy, gx=gx, gb=gb, rts=rts,
                     evalindex=evalindex, refindex=refindex).optimize(email, solver)
    out = pd.DataFrame(index=res.index)
    out["E"] = res["rho"].astype(float)
    return out


# --- decomposition helpers ----------------------------------------------------
def _ratio_global(glob, cont, id):
    """Pastor & Lovell (2005) global-frontier split for a ratio index.

    Returns ``D11`` (global), ``mqpi``, ``MEFFCH`` (contemporary efficiency
    change) and ``MTECHCH`` (best-practice-gap change), so that
    ``mqpi = MEFFCH * MTECHCH``.
    """
    bpg = pd.DataFrame({id: glob[id], "BPG": glob["D11"] / cont["D11"]}, index=glob.index)
    bpg.attrs = cont.attrs.copy()
    out = pd.DataFrame(index=glob.index)
    out["D11"] = glob["D11"]
    out["mqpi"] = glob["D11"] / _shift(glob, id, "D11")
    out["MEFFCH"] = cont["D11"] / _shift(cont, id, "D11")
    out["MTECHCH"] = bpg["BPG"] / _shift(bpg, id, "BPG")
    return out


def _luenberger_decompose(frame, id):
    """Additive Luenberger decomposition (§7.2): L = LEFFCH + LTECHCH, benchmark 0.

    Uses the same D11 / D12 / D21 layout produced by ``build_distance_frame``;
    under GLOBAL (where D12 == D21 == D11) it collapses to LTECHCH = 0 and
    L = LEFFCH, the additive analogue of MTECHCH = 1.
    """
    d11 = frame["D11"]
    d12 = frame["D12"]
    d21 = frame["D21"]
    d11_prev = _shift(frame, id, "D11")
    d21_prev = _shift(frame, id, "D21")
    out = pd.DataFrame(index=frame.index)
    out["LEFFCH"] = d11_prev - d11
    out["LTECHCH"] = 0.5 * (d11 - d12) + 0.5 * (d21_prev - d11_prev)
    out["L"] = 0.5 * (d11_prev + d21_prev - d11 - d12)
    return out


# --- dynamic indicator factories (return a mqcore.MalmquistModel) ------------
def dep(data, id, year, sent="inputvar=outputvar:unoutputvar",
        rts=RTS_CRS, tech=CONTEMPORARY, window=3, verbose=False):
    """Dynamic environmental-performance index (MEPI, §7.4.1).

    Built on :func:`ushephard` (undesirable-output-oriented Shephard distance);
    corresponds to the MEPI of eq. (7-45).  Returns a
    :class:`~deabook.mqcore.MalmquistModel`; call ``.optimize(email, solver)``
    to solve.  ``tech`` accepts ``CONTEMPORARY`` (default), ``WINDOW``,
    ``SEQUENTIAL`` and ``GLOBAL`` (Chapter 3 §3.4).  Non-global technologies
    return ``D11``/``D12``/``D21`` plus ``MQ``/``MEFFCH``/``MTECHCH``; ``GLOBAL``
    returns ``D11``/``mqpi``/``MEFFCH``/``MTECHCH`` (Pastor & Lovell, 2005).
    """
    if rts not in (RTS_CRS, RTS_VRS1, RTS_VRS2):
        raise ValueError("rts must be RTS_CRS, RTS_VRS1 or RTS_VRS2.")

    def build(work, years, t, w, year_col, email, solver):
        solve = lambda e, r, em, sv: _ushephard_distance(work, sent, rts, e, r, em, sv)
        if t == GLOBAL:
            glob = build_distance_frame(solve, work, years, GLOBAL, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            cont = build_distance_frame(solve, work, years, CONTEMPORARY, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            return _ratio_global(glob, cont, id)
        frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
        return pd.concat([frame[["D11", "D12", "D21"]], two_factor(frame, id)], axis=1)

    return MalmquistModel(data, id, year, tech=normalize_tech(tech), window=int(window),
                          build=build, kind="dep", orient=None, verbose=verbose)


def depddf(data, id, year, sent="inputvar=outputvar:unoutputvar",
           rts=RTS_CRS, tech=CONTEMPORARY, window=3, verbose=False):
    """DDF-based dynamic environmental-performance index (ML, §7.4.2).

    Built on the bad-output-oriented directional distance function (direction
    ``(gx, gy, gb) = (0, 0, 1)``, i.e. pollution abatement only), the directional
    analogue of :func:`ushephard`; corresponds to the ML index of eq. (7-51).
    Unlike the general radial/output-oriented DDF, this indicator keeps the
    same non-desirable-output direction as :func:`dep`, so the two are directly
    comparable.  Returns a :class:`~deabook.mqcore.MalmquistModel`; call
    ``.optimize(email, solver)`` to solve.  ``tech`` accepts ``CONTEMPORARY``
    (default), ``WINDOW``, ``SEQUENTIAL`` and ``GLOBAL`` (Chapter 3 §3.4).
    Non-global technologies return ``D11``/``D12``/``D21`` plus ``MQ``/``MEFFCH``/
    ``MTECHCH``; ``GLOBAL`` returns ``D11``/``mqpi``/``MEFFCH``/``MTECHCH``
    (Pastor & Lovell, 2005).
    """
    if rts not in (RTS_CRS, RTS_VRS1, RTS_VRS2):
        raise ValueError("rts must be RTS_CRS, RTS_VRS1 or RTS_VRS2.")

    def build(work, years, t, w, year_col, email, solver):
        solve = lambda e, r, em, sv: _ddf_bad_distance(work, sent, rts, e, r, em, sv)
        if t == GLOBAL:
            glob = build_distance_frame(solve, work, years, GLOBAL, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            cont = build_distance_frame(solve, work, years, CONTEMPORARY, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            return _ratio_global(glob, cont, id)
        frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
        return pd.concat([frame[["D11", "D12", "D21"]], two_factor(frame, id)], axis=1)

    return MalmquistModel(data, id, year, tech=normalize_tech(tech), window=int(window),
                          build=build, kind="depddf", orient=None, verbose=verbose)


def lddf(data, id, year, sent="inputvar=outputvar:unoutputvar",
         orient="oo", rts=RTS_CRS, tech=CONTEMPORARY, window=3, verbose=False):
    """Additive Luenberger indicator based on the directional distance function.

    Built on the single-period Luenberger DDF (:func:`_ldf_distance`).  ``orient``
    selects ``"oo"`` (§7.2.1), ``"io"`` (§7.2.2) or ``"both"`` (§7.2.3 / §7.4).
    Returns a :class:`~deabook.mqcore.MalmquistModel`; ``.optimize(email,
    solver)`` solves and returns ``D11``/``D12``/``D21`` plus the additive
    decomposition ``L``/``LEFFCH``/``LTECHCH`` (benchmark 0).  ``tech`` accepts
    ``CONTEMPORARY`` (default), ``WINDOW``, ``SEQUENTIAL`` and ``GLOBAL``
    (Chapter 3 §3.4); under ``GLOBAL`` the technical-change component ``LTECHCH``
    is 0 and ``L`` equals the efficiency-change component.
    """
    if orient not in ("oo", "io", "both"):
        raise ValueError("orient must be 'oo' (output, §7.2.1), 'io' (input, §7.2.2) or 'both' (graph, §7.2.3).")
    if rts not in (RTS_CRS, RTS_VRS1, RTS_VRS2):
        raise ValueError("rts must be RTS_CRS, RTS_VRS1 or RTS_VRS2.")
    if orient in ("io", "both") and rts == RTS_VRS1:
        raise ValueError(
            "rts=RTS_VRS1 is not supported for the input/both orientation of the "
            "Luenberger DDF: the directional distance function with an input direction "
            "is not linearizable under VRS1. Use RTS_CRS or RTS_VRS2 instead."
        )

    def build(work, years, t, w, year_col, email, solver):
        solve = lambda e, r, em, sv: _ldf_distance(work, sent, orient, rts, e, r, em, sv)
        if t == GLOBAL:
            # Pastor & Lovell (2005) global Malmquist is a *ratio* decomposition
            # (MEFFCH x MTECHCH). The additive Luenberger needs no analogous
            # contemporary-vs-global split: under GLOBAL build_distance_frame sets
            # D12 = D21 = D11, so _luenberger_decompose collapses to LTECHCH = 0
            # and L = LEFFCH (the additive analogue of MTECHCH = 1). We therefore
            # apply the same additive decomposition used for the other technologies
            # rather than branching the formula, keeping L = LEFFCH + LTECHCH valid.
            frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            return pd.concat([frame[["D11", "D12", "D21"]], _luenberger_decompose(frame, id)], axis=1)
        frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
        return pd.concat([frame[["D11", "D12", "D21"]], _luenberger_decompose(frame, id)], axis=1)

    return MalmquistModel(data, id, year, tech=normalize_tech(tech), window=int(window),
                          build=build, kind="lddf", orient=orient, verbose=verbose)


# --- static slacks-based measure (SBM) factory ---------------------------------
def _sbm_sent_has_bad(sent):
    """Return ``True`` if ``sent`` declares an undesirable-output component.

    The package ``sent`` convention is ``"inputs = outputs"`` or
    ``"inputs = outputs : undesirable_outputs"``; an undesirable part exists iff
    the sub-string after the first ``=`` contains a ``:``.
    """
    return len(parse_sent(sent).b) > 0


def sbmu(data, sent="inputvar=outputvar:unoutputvar", gx=None, gy=None, gb=None,
         rts=RTS_VRS1, evalindex=None, refindex=None):
    """Static slacks-based measure (SBM) *with* undesirable outputs (§8.4).

    Extends :func:`sbm` by letting the ``:``-separated undesirable outputs enter
    the SBM objective as an additional slack-excess dimension (Tone, 2003/2010).
    Returns the (unsolved) :class:`~deabook.DEA._SBM` model object; call
    ``.optimize(email, solver)`` to solve and obtain ``rho`` plus the
    ``slack_x_*`` / ``slack_y_*`` / ``slack_b_*`` columns.  ``sent`` MUST
    declare an undesirable-output component.

    Args:
        data (pandas.DataFrame): Input-output data.
        sent (str): ``"inputs = outputs : undesirable"`` expression, e.g.
            ``"K L E=Y:CO2"``.
        gx (list, optional): Input slack direction.
        gy (list, optional): Desirable-output slack direction.
        gb (list, optional): Undesirable-output slack direction.
        rts (str): ``RTS_CRS`` or ``RTS_VRS1``.
        evalindex (str, optional): Evaluated-sample filter.
        refindex (str, optional): Reference-sample filter.

    Returns:
        deabook.DEA._SBM: An unsolved SBM model; call ``.optimize()`` on it.

    Example:
        >>> m = sbmu(data, sent="K L E=Y:CO2", gy=[1], gx=[1, 1, 1], gb=[1],
        ...          rts=RTS_VRS1)
        >>> res = m.optimize(OPT_LOCAL, "highs")
    """
    if not _sbm_sent_has_bad(sent):
        raise ValueError(
            "DEAweak.sbmu is for SBM *with* undesirable outputs, but the sent "
            f"{sent!r} has no undesirable-output component. Separate it with ':' "
            "(e.g. \"K L E=Y:CO2\"), or use DEAweak.sbm for the standard SBM (§8.3)."
        )
    kwargs = {}
    if gx is not None:
        kwargs["gx"] = gx
    if gy is not None:
        kwargs["gy"] = gy
    if gb is not None:
        kwargs["gb"] = gb
    return _SBM(data, sent, rts=rts, evalindex=evalindex, refindex=refindex, **kwargs)


def sddf(data, sent="inputvar=outputvar:unoutputvar", gx=None, gy=None, gb=None,
         rts=RTS_VRS1, evalindex=None, refindex=None):
    """Slack-based directional distance function with undesirable outputs
    (Pham & Zelenyuk, 2018).

    Extends :func:`DEA.sddf` by letting the ``:``-separated undesirable outputs
    enter the SD objective as an additional slack-excess dimension under joint
    weak disposability of good and bad outputs (§8.5.2).  Returns the (unsolved)
    :class:`~deabook.DEA._SDDF` model object; call ``.optimize(email, solver)``
    to obtain ``objective_value`` (the SD inefficiency) plus the
    ``slack_x_*`` / ``slack_y_*`` / ``slack_b_*`` columns.  ``sent`` MUST
    declare an undesirable-output component.

    Args:
        data (pandas.DataFrame): Input-output data.
        sent (str): ``"inputs = outputs : undesirable"`` expression, e.g.
            ``"K L E=Y:CO2"``.
        gx (list, optional): Input slack direction.
        gy (list, optional): Desirable-output slack direction.
        gb (list, optional): Undesirable-output slack direction.
        rts (str): ``RTS_CRS`` or ``RTS_VRS1``.
        evalindex (str, optional): Evaluated-sample filter.
        refindex (str, optional): Reference-sample filter.

    Returns:
        deabook.DEA._SDDF: An unsolved SDDF model; call ``.optimize()`` on it.

    Example:
        >>> m = sddf(data, sent="K L E=Y:CO2", rts=RTS_VRS1)
        >>> res = m.optimize(OPT_LOCAL, "highs")
    """
    if not _sbm_sent_has_bad(sent):
        raise ValueError(
            "DEAweak.sddf is the Pham-Zelenyuk (2018) SD model *with* "
            f"undesirable outputs, but the sent {sent!r} has no undesirable-"
            "output component. Separate it with ':' (e.g. \"K L E=Y:CO2\"), or "
            "use DEA.sddf for the standard model (§8.5.1)."
        )
    kwargs = {}
    if gx is not None:
        kwargs["gx"] = gx
    if gy is not None:
        kwargs["gy"] = gy
    if gb is not None:
        kwargs["gb"] = gb
    return _SDDF(data, sent, rts=rts, evalindex=evalindex, refindex=refindex, **kwargs)


# ===========================================================================
# Migrated from backup (_deabook_backup_20260905/DEAweak.py) on 2026-09-05.
# These four classes (DDFweak, DDFweakRho, NDDFweak, NDDFweakRho) were
# removed earlier alongside the deleted inner `DEA.DDF`; ch6 uses
# DEAweak.ddf / nddf / dep / depddf / lddf which require them. The deleted
# `DDF` base class is dropped; the classes now inherit the existing weak
# mixins only. No other code in the package is changed.
# ===========================================================================
class DDFweak(_WeakOrientationMixin, _WeakPerDmuResultMixin):

    def __init__(self,data,sent,gy=[1],gx=[0],gb=[1],rts=RTS_VRS1,evalindex=None, refindex=None):
        """ DDF: Envelopment problem, solving for each DMU individually.
        
        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y:CO2"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            evalindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEAweak should return numpy arrays for x, y,b, xref, yref,bref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DEAweak is modified or a similar function returns these:
        try:
             self.y, self.x,self.b, self.yref,self.xref,self.bref, self.gy, self.gx,self.gb, \
                 self.evaluated_data_index, self.reference_data_index = tools.assert_DEAweak_with_indices(
                     data, sent, gy, gx,gb, evalindex, refindex
                 )
        except AttributeError:
             # Fallback if assert_DEA_with_indices doesn't exist, assuming original assert_DEA
             # In this case, evaluated_data_index and reference_data_index will be ranges 0..N-1
             print("Warning: tools.assert_DEAweak_with_indices not found. Using range indices.")
             self.y,self.x,self.b, self.yref,self.xref,self.bref, self.gy,self.gx,self.gb = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, evalindex, refindex)
             self.evaluated_data_index = list(range(len(self.y)))
             self.reference_data_index = list(range(len(self.yref)))

        self.rts = rts

        self.evaluated_indices_range = range(len(self.y)) # Range indices for internal numpy arrays
        self.reference_indices_range = range(len(self.yref)) # Range indices for internal numpy arrays
        self.num_inputs = len(self.x[0])
        self.num_outputs = len(self.y[0])
        self.num_unoutputs = len(self.b[0])
        # print(self.gx,self.gy,")))))))))))))")
        # Determine orientation based on gx/gy/gb vectors.
        self._set_weak_orientation_flags()
        if (self.input_oriented and self.output_oriented)   \
            or (self.input_oriented and self.unoutput_oriented) or (self.output_oriented and self.unoutput_oriented):
             # Original code structure suggests one or the other dominates based on the if/elif.
             # Assuming standard DEA where it's one or the other.
             pass # This case implies a potential issue if both sums >= 1, but following original logic.


        # Dictionary to store individual models for each evaluated DMU
        self.__modeldict = {}
        # Dictionary to store results after optimization
        self.results = {}

        # Loop through each DMU to be evaluated (using range index to access numpy arrays)
        for i_range in self.evaluated_indices_range:
            # Get the actual index label for this DMU
            actual_index = self.evaluated_data_index[i_range]

            # Create a new model for this specific DMU
            model = ConcreteModel()


            # Define sets for this specific model
            # R is the set of reference DMU range indices
            model.R = Set(initialize=self.reference_indices_range)
            # J is the set of input indices
            model.J = Set(initialize=range(self.num_inputs))
            # K is the set of output indices
            model.K = Set(initialize=range(self.num_outputs))
            # L is the set of unoutput indices
            model.L = Set(initialize=range(self.num_unoutputs))
            # Define variables for this specific model
            # rho is the efficiency score for the current DMU (not indexed by DMU within this model)
            # Bounds based on standard DEA orientation (matching DEAt)
            if self.rts == RTS_CRS:
                model.rho = Var(bounds=(0, None), doc=f'beta for DMU {actual_index}')
            elif self.rts == RTS_VRS1:
                model.rho = Var(bounds=(0, None), doc=f'beta for DMU {actual_index}')
                model.theta = Var(bounds=(0.0, 1.0),doc='emission reduction factor')
            elif self.rts == RTS_VRS2:
                model.rho = Var(bounds=(0, None), doc=f'beta for DMU {actual_index}')
                model.mu = Var(model.R, bounds=(0.0, None),doc='emission reduction factor2')
           
            else:
                 # This else should not be reached due to check above, but for safety
                 raise ValueError("Invalid orientation configuration.")
            # lamda are intensity variables, indexed by the reference set
            model.lamda = Var(model.R, bounds=(0.0, None), doc='intensity variables')

            # Setup the objective function and constraints for THIS DMU
            # Objective: Maximize the SINGLE rho variable in this model
            model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc='objective function')
            

            # Input constraint for THIS DMU, referencing current DMU's data (self.x[i_range])
            # Use a factory function to capture the current DMU index
            input_rule_factory = self.__create_input_rule(i_range)
            model.input = Constraint(model.J, rule=input_rule_factory, doc='input constraint')

            # Output constraint for THIS DMU, referencing current DMU's data (self.y[i_range])
            output_rule_factory = self.__create_output_rule(i_range)
            model.output = Constraint(model.K, rule=output_rule_factory, doc='output constraint')

            # UNOutput constraint for THIS DMU, referencing current DMU's data (self.b[i_range])
            unoutput_rule_factory = self.__create_unoutput_rule(i_range)
            model.unoutput = Constraint(model.L, rule=unoutput_rule_factory, doc='unoutput constraint')
            # VRS constraint for THIS DMU
            if self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                vrs_rule_factory = self.__create_vrs_rule() # VRS rule doesn't need current DMU index
                model.vrs = Constraint(rule=vrs_rule_factory, doc='variable return to scale rule')

            # Store the created model in the dictionary using the actual DMU index as the key
            self.__modeldict[actual_index] = model

        # Optimization status will be stored per DMU after optimize is called
        # self.optimization_status = 0 # Not needed as overall status
        # self.problem_status = 0 # Not needed as overall status


    # --- Rule Factories (Helper methods to create constraint rules) ---
    # These capture the data for the specific DMU being modeled

    def __create_input_rule(self, current_dmu_range_index):
        """Factory for creating the input constraint rule for a specific DMU."""

        if self.input_oriented:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Access current DMU's input data self.x[current_dmu_range_index][j]
                    # Access reference DMUs' input data self.xref[r][j]
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rho*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Access current DMU's input data self.x[current_dmu_range_index][j]
                    # Access reference DMUs' input data self.xref[r][j]
                    # Input is scaled by rho
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <=\
                        self.x[current_dmu_range_index][j] - model.rho*self.gx[j]* self.x[current_dmu_range_index][j]
        elif self.output_oriented or self.unoutput_oriented:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                def input_rule(model, j):
                    # Input is scaled by theta in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=model.theta*self.x[current_dmu_range_index][j] 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum((model.lamda[r]+model.mu[r])* self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j]
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rho*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rho in input orientation
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rho*self.gx[j]*self.x[current_dmu_range_index][j]
 
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Hyper orientation: Input is scaled by rho if gx[j] == 1
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                def input_rule(model, j):
                    # Input is scaled by rho in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= model.theta * self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum((model.lamda[r]+model.mu[r])* self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j]
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rho*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rho in input orientation
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rho*self.gx[j]*self.x[current_dmu_range_index][j]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rho*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rho in input orientation
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rho*self.gx[j]*self.x[current_dmu_range_index][j]
 
        else:
            # Should not be reached
            def input_rule(model, j):
                return Constraint.Skip # Or raise error
        # Return the rule function
        return input_rule

    def __create_output_rule(self, current_dmu_range_index):
        """Factory for creating the output constraint rule for a specific DMU."""
        if self.input_oriented:
            if self.rts == RTS_CRS  or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.output_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is scaled by rho in output orientation (original code applies rho to all outputs if sum(gy)>=1)
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.unoutput_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) \
                        >= self.y[current_dmu_range_index][k]
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rho if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]+ model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rho if gy[k] == 1
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]           
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rho if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rho if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]+ model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
        else:
             # Should not be reached
             def output_rule(model, k):
                 return Constraint.Skip # Or raise error
        return output_rule

    def __create_unoutput_rule(self, current_dmu_range_index):
        """Factory for creating the unoutput constraint rule for a specific DMU."""
        if self.input_oriented :
            if self.rts == RTS_CRS or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # unOutput is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.output_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # unOutput is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l]
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.unoutput_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Output is scaled by rho in output orientation (original code applies rho to all outputs if sum(gy)>=1)
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rho if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:  
                def unoutput_rule(model, l):
                    # Output is not scaled by rho in input orientation
                    return sum((model.lamda[r])*self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l] 
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS  or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rho if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rho if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:  
                def unoutput_rule(model, l):
                    # Output is not scaled by rho in input orientation
                    return sum((model.lamda[r])*self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rho if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:  
                def unoutput_rule(model, l):
                    # Output is not scaled by rho in input orientation
                    return sum((model.lamda[r])*self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rho*self.gb[l]*self.b[current_dmu_range_index][l]
        
        else:
             # Should not be reached
             def unoutput_rule(model, l):
                 return Constraint.Skip # Or raise error
        return unoutput_rule

    def __create_vrs_rule(self):
        """Factory for creating the VRS constraint rule."""
        if self.rts==RTS_VRS1:
            def vrs_rule(model):
                return sum(model.lamda[r] for r in model.R) == model.theta
            return vrs_rule
        elif self.rts==RTS_VRS2:
            def vrs_rule(model):
                return sum((model.lamda[r]+model.mu[r]) for r in model.R) == 1
            return vrs_rule



    # --- Optimization and Results Methods ---

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize every evaluated DMU and return legacy DDFweak columns."""
        return self._solve_per_dmu_models(email=email, solver=solver, include_objective=True)

    def _add_efficiency_columns(self, results_df):
        rho = results_df['rho']
        if self.input_oriented:
            results_df['tei'] = self._one_minus(rho)
        elif self.output_oriented:
            results_df['teo'] = self._safe_inverse_one_plus(rho)
        elif self.unoutput_oriented:
            results_df['teuo'] = self._one_minus(rho)
        elif self.hyper_orientedyx:
            results_df['tei'] = self._one_minus(rho)
            results_df['teo'] = self._safe_inverse_one_plus(rho)
        elif self.hyper_orientedyb:
            results_df['teuo'] = self._one_minus(rho)
            results_df['teo'] = self._safe_inverse_one_plus(rho)
        elif self.hyper_orientedxb:
            results_df['tei'] = self._one_minus(rho)
            results_df['teuo'] = self._one_minus(rho)
        elif self.hyper_orientedyxb:
            results_df['tei'] = self._one_minus(rho)
            results_df['teuo'] = self._one_minus(rho)
            results_df['teo'] = self._safe_inverse_one_plus(rho)
        return results_df

class DDFweakRho(DDFweak):
    """Weak-disposability directional distance function that reports only the
    raw LP value ``rho`` (plus ``objective_value``), mirroring ``DEA.ddf``'s
    ``_DdfRho`` behaviour.

    ``DEA.ddf`` leaves the efficiency conversion (``1 - rho`` on an input
    direction, ``1 / (1 + rho)`` on an output direction) to the caller, and so
    does this class: it does **not** generate the ``teuo``/``teo`` decomposition
    columns that :class:`DDFweak` produces.  The chapter derives environmental
    efficiency ``EB`` and output efficiency ``EO`` from ``rho`` at the call
    site, exactly as :func:`DEA.ddf` callers derive ``E_tei``/``E_teo``.
    """

    def _add_efficiency_columns(self, results_df):
        return results_df


class NDDFweak(_WeakOrientationMixin, _WeakPerDmuResultMixin):

    def __init__(self,data,sent,gy=[1],gx=[0],gb=[1],rts=RTS_VRS1,evalindex=None, refindex=None):
        """ DDF: Envelopment problem, solving for each DMU individually.
        
        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y:CO2"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            evalindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEAweak should return numpy arrays for x, y,b, xref, yref,bref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DEAweak is modified or a similar function returns these:
        self.y, self.x,self.b, self.outputvars,self.inputvars,self.unoutputvars, \
            self.yref,self.xref,self.bref, self.gy, self.gx,self.gb, self.wx,self.wy,self.wb, \
            self.evaluated_data_index, self.reference_data_index = tools.assert_NDDFweak_with_indices(
                data, sent, gy, gx,gb, evalindex, refindex
            )

        self.rts = rts

        self.evaluated_indices_range = range(len(self.y)) # Range indices for internal numpy arrays
        self.reference_indices_range = range(len(self.yref)) # Range indices for internal numpy arrays
        self.num_inputs = len(self.x[0])
        self.num_outputs = len(self.y[0])
        self.num_unoutputs = len(self.b[0])
        # print(self.gx,self.gy,"############")

        self.sum_gx = sum(self.gx)
        self.sum_gy = sum(self.gy)
        self.sum_gb = sum(self.gb)
        # Determine orientation based on gx/gy/gb vectors.
        self._set_weak_orientation_flags()
        # Return-shape flags for the nddf interface:
        #  - NDDFweak    (full): emits per-variable beta ("rho<var>") AND the
        #    te*/weight/aggregate efficiency columns.
        #  - NDDFweakRho (minimal): emits only the LP value (objective_value)
        #    plus per-variable beta ("beta_<var>"); efficiency is derived by the
        #    caller.  The prefix lets subclasses rename the beta columns.
        self._minimal = False
        self._beta_prefix = "rho"
        if (self.input_oriented and self.output_oriented)   \
            or (self.input_oriented and self.unoutput_oriented) or (self.output_oriented and self.unoutput_oriented):
             # Original code structure suggests one or the other dominates based on the if/elif.
             # Assuming standard DEA where it's one or the other.
             pass # This case implies a potential issue if both sums >= 1, but following original logic.


        # Dictionary to store individual models for each evaluated DMU
        self.__modeldict = {}
        # Dictionary to store results after optimization
        self.results = {}

        # Loop through each DMU to be evaluated (using range index to access numpy arrays)
        for i_range in self.evaluated_indices_range:
            # Get the actual index label for this DMU
            actual_index = self.evaluated_data_index[i_range]

            # Create a new model for this specific DMU
            model = ConcreteModel()

            def rhox_initialize_rule(model, j):
                return 0.0

            def rhoy_initialize_rule(model, k):
                    # 类似地为 rhoy 定义初始值函数
                return 0.0

            def rhob_initialize_rule(model, l):
                    # 类似地为 rhob 定义初始值函数
                return 0.0


            # Define sets for this specific model
            # R is the set of reference DMU range indices
            model.R = Set(initialize=self.reference_indices_range)
            # J is the set of input indices
            model.J = Set(initialize=range(self.num_inputs))
            # K is the set of output indices
            model.K = Set(initialize=range(self.num_outputs))
            # L is the set of unoutput indices
            model.L = Set(initialize=range(self.num_unoutputs))
            # Define variables for this specific model
            # rho is the efficiency score for the current DMU (not indexed by DMU within this model)
            # Bounds based on standard DEA orientation (matching DEAt)
            if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb:
                model.rhox = Var(model.J,initialize=rhox_initialize_rule,bounds=(0, None), doc=f'beta in input orientation for DMU {actual_index}')
            if self.output_oriented or self.hyper_orientedyb or self.hyper_orientedyx or self.hyper_orientedyxb:
                model.rhoy = Var(model.K,initialize=rhoy_initialize_rule,bounds=(0, None), doc=f'beta in output orientation for DMU {actual_index}')
            if self.unoutput_oriented or self.hyper_orientedxb or self.hyper_orientedyb or self.hyper_orientedyxb:
                model.rhob = Var(model.L,initialize=rhob_initialize_rule,bounds=(0, None), doc=f'beta in unoutput orientation for DMU {actual_index}')
            if not (self.input_oriented or self.output_oriented or self.unoutput_oriented \
                or self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb):
                raise ValueError("Invalid orientation configuration.")

            if self.rts == RTS_VRS1:
                model.theta = Var(bounds=(0.0, 1.0),doc='emission reduction factor')
            elif self.rts == RTS_VRS2:
                model.mu = Var(model.R, bounds=(0.0, None),doc='emission reduction factor2')
           

            # lamda are intensity variables, indexed by the reference set
            model.lamda = Var(model.R, bounds=(0.0, None), doc='intensity variables')

            # Setup the objective function and constraints for THIS DMU
            # Objective: Maximize the weighted sum of rhox, rhoy, rhob
            objective_rule_factory = self.__create_objective_rule()
            model.objective = Objective(rule=objective_rule_factory, sense=maximize, doc='objective function')

            # Input constraint for THIS DMU, referencing current DMU's data (self.x[i_range])
            input_rule_factory = self.__create_input_rule(i_range)
            model.input = Constraint(model.J, rule=input_rule_factory, doc='input constraint')

            # Output constraint for THIS DMU, referencing current DMU's data (self.y[i_range])
            output_rule_factory = self.__create_output_rule(i_range)
            model.output = Constraint(model.K, rule=output_rule_factory, doc='output constraint')

            # UNOutput constraint for THIS DMU, referencing current DMU's data (self.b[i_range])
            unoutput_rule_factory = self.__create_unoutput_rule(i_range)
            model.unoutput = Constraint(model.L, rule=unoutput_rule_factory, doc='unoutput constraint')
            # VRS constraint for THIS DMU
            if self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                vrs_rule_factory = self.__create_vrs_rule() # VRS rule doesn't need current DMU index
                model.vrs = Constraint(rule=vrs_rule_factory, doc='variable return to scale rule')

            # Store the created model in the dictionary using the actual DMU index as the key
            self.__modeldict[actual_index] = model

        # Optimization status will be stored per DMU after optimize is called
        # self.optimization_status = 0 # Not needed as overall status
        # self.problem_status = 0 # Not needed as overall status


    # --- Rule Factories (Helper methods to create constraint rules) ---
    # These capture the data for the specific DMU being modeled


    def __create_objective_rule(self):
            """Return the proper objective function based on variables defined on the model."""
            def objective_rule(model):
                objective_expr = 0

                # Check if rhox was defined on the model instance
                if hasattr(model, 'rhox'):
                    objective_expr += sum(self.wx[j] * model.rhox[j] for j in model.J)

                # Check if rhoy was defined on the model instance
                if hasattr(model, 'rhoy'):
                    objective_expr += sum(self.wy[k] * model.rhoy[k] for k in model.K)

                # Check if rhob was defined on the model instance
                if hasattr(model, 'rhob'):
                    objective_expr += sum(self.wb[l] * model.rhob[l] for l in model.L)

                return objective_expr
            return objective_rule
    
    
    def __create_input_rule(self, current_dmu_range_index):
        """Factory for creating the input constraint rule for a specific DMU."""

        if self.input_oriented:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Access current DMU's input data self.x[current_dmu_range_index][j]
                    # Access reference DMUs' input data self.xref[r][j]
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Access current DMU's input data self.x[current_dmu_range_index][j]
                    # Access reference DMUs' input data self.xref[r][j]
                    # Input is scaled by rhox
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <=\
                        self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]* self.x[current_dmu_range_index][j]
        elif self.output_oriented or self.unoutput_oriented:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                def input_rule(model, j):
                    # Input is scaled by theta in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=model.theta*self.x[current_dmu_range_index][j] 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum((model.lamda[r]+model.mu[r])* self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j]
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rhox in input orientation
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]*self.x[current_dmu_range_index][j]
 
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    # Hyper orientation: Input is scaled by rhox if gx[j] == 1
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                def input_rule(model, j):
                    # Input is scaled by rhox in input orientation
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= model.theta * self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is not scaled by theta in input orientation
                    return sum((model.lamda[r]+model.mu[r])* self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j]
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rhox in input orientation
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]*self.x[current_dmu_range_index][j]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def input_rule(model, j):
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <=\
                          self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]* self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def input_rule(model, j):
                    # Input is scaled by rhox in input orientation
                    return sum((model.lamda[r]+model.mu[r])*self.xref[r][j] for r in model.R) <= \
                        self.x[current_dmu_range_index][j] - model.rhox[j]*self.gx[j]*self.x[current_dmu_range_index][j]
 
        else:
            # Should not be reached
            def input_rule(model, j):
                return Constraint.Skip # Or raise error
        # Return the rule function
        return input_rule

    def __create_output_rule(self, current_dmu_range_index):
        """Factory for creating the output constraint rule for a specific DMU."""
        if self.input_oriented:
            if self.rts == RTS_CRS  or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rhoy in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.output_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is scaled by rhoy in output orientation (original code applies rhoy to all outputs if sum(gy)>=1)
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rhoy[k]*self.gy[k]*self.y[current_dmu_range_index][k]
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.unoutput_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rhoy in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) \
                        >= self.y[current_dmu_range_index][k]
            else:   
                def output_rule(model, k):
                    return Constraint.Skip # Or raise error
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rhoy if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]+ model.rhoy[k]*self.gy[k]*self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rhoy in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rhoy[k]*self.gy[k]*self.y[current_dmu_range_index][k]
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rhoy if gy[k] == 1
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rhoy[k]*self.gy[k]*self.y[current_dmu_range_index][k]           
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rhoy if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rhoy in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rhoy if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k]+ model.rhoy[k]*self.gy[k]*self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:
                def output_rule(model, k):
                    # Output is not scaled by rhoy in input orientation
                    return sum(model.lamda[r]*self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rhoy[k]*self.gy[k]*self.y[current_dmu_range_index][k]
        else:
             # Should not be reached
             def output_rule(model, k):
                 return Constraint.Skip # Or raise error
        return output_rule

    def __create_unoutput_rule(self, current_dmu_range_index):
        """Factory for creating the unoutput constraint rule for a specific DMU."""
        if self.input_oriented :
            if self.rts == RTS_CRS or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # unOutput is not scaled by rhob in input orientation
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for input orientation.")
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.output_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # unOutput is not scaled by rhob in input orientation
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l]
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.unoutput_oriented:
            if self.rts == RTS_CRS or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Output is scaled by rhob in output orientation (original code applies rhob to all outputs if sum(gy)>=1)
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rhob[l]*self.gb[l]*self.b[current_dmu_range_index][l]
            else:   
                def unoutput_rule(model, l):
                    return Constraint.Skip # Or raise error
        elif self.hyper_orientedyx:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rhob if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:  
                def unoutput_rule(model, l):
                    # Output is not scaled by rhob in input orientation
                    return sum((model.lamda[r])*self.bref[r][l] for r in model.R) == self.b[current_dmu_range_index][l] 
        elif self.hyper_orientedyb:
            if self.rts == RTS_CRS  or self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rhob if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rhob[l]*self.gb[l]*self.b[current_dmu_range_index][l]
        elif self.hyper_orientedxb:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rhob if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rhob[l]*self.gb[l]*self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:  
                def unoutput_rule(model, l):
                    # Output is not scaled by rhob in input orientation
                    return sum((model.lamda[r])*self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rhob[l]*self.gb[l]*self.b[current_dmu_range_index][l]
        elif self.hyper_orientedyxb:
            if self.rts == RTS_CRS:
                def unoutput_rule(model, l):
                    # Hyper orientation: Output is not scaled by rhob if gb[k] == 1
                    return sum(model.lamda[r] * self.bref[r][l] for r in model.R) ==\
                        self.b[current_dmu_range_index][l] - model.rhob[l]*self.gb[l]*self.b[current_dmu_range_index][l]
            elif self.rts == RTS_VRS1:
                raise ValueError("RTS_VRS1 not supported for this hyper orientation (input direction + same reduction factor is not linearizable).") 
            elif self.rts == RTS_VRS2:  
                def unoutput_rule(model, l):
                    # Output is not scaled by rhob in input orientation
                    return sum((model.lamda[r])*self.bref[r][l] for r in model.R) == \
                        self.b[current_dmu_range_index][l] - model.rhob[l]*self.gb[l]*self.b[current_dmu_range_index][l]
        
        else:
             # Should not be reached
             def unoutput_rule(model, l):
                 return Constraint.Skip # Or raise error
        return unoutput_rule

    def __create_vrs_rule(self):
        """Factory for creating the VRS constraint rule."""
        if self.rts==RTS_VRS1:
            def vrs_rule(model):
                return sum(model.lamda[r] for r in model.R) == model.theta
            return vrs_rule
        elif self.rts==RTS_VRS2:
            def vrs_rule(model):
                return sum((model.lamda[r]+model.mu[r]) for r in model.R) == 1
            return vrs_rule




    # --- Optimization and Results Methods ---

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the model for each DMU individually.

        Args:
            email (string): The email address for remote optimization (ignored if solver is local).
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.

        Returns:
            pandas.DataFrame: DataFrame containing optimization status and rho (efficiency) for each evaluated DMU.
        """
        # 1. Centralized Initialization
        self.results = {}  # Clear previous results
        all_optimization_statuses = {}
        # Initialize dictionaries to collect raw results
        rhoy_values = {}
        rhox_values = {}
        rhob_values = {}
        theta_values = {}
        objective_values = {}
        lamda_values = {}
        mu_values = {}

        # Loop through each DMU's model and solve it
        # The dictionary keys are the actual data indices

        use_neos = tools.set_neos_email(email)

        # 2. Consolidated Optimization Loop
        for actual_index, model in self.__modeldict.items():
            # print(f"Optimizing for DMU: {actual_index}...") # Optional: print progress
            optimization_status = "Error" # Default status
            try:
                optimization_status = tools.optimize_model2(
                    model, actual_index, use_neos, CET_ADDI, solver=solver)
            except Exception as e:
                print(f"Error optimizing DMU {actual_index}: {e}")
                # Status remains 'Error'

            all_optimization_statuses[actual_index] = optimization_status

            # 3. Dynamic Result Collection (inside loop)
            if optimization_status in ['ok']: # Check standard Pyomo status strings
                try:
                    # Collect rhoy if applicable
                    if self.output_oriented or self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedyxb:
                        rhoy_data_for_sub = {}
                        # model.K is the index set for output variables
                        for k_range in model.K:
                            # Use value() to get the numeric value
                            rhoy_data_for_sub[k_range] = model.rhoy[k_range].value
                        rhoy_values[actual_index] = rhoy_data_for_sub

                    # Collect rhox if applicable
                    if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb:
                        rhox_data_for_sub = {}
                        # model.J is the index set for input variables
                        for j_range in model.J:
                            rhox_data_for_sub[j_range] = model.rhox[j_range].value
                        rhox_values[actual_index] = rhox_data_for_sub

                    # Collect rhob if applicable
                    if self.unoutput_oriented or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb:
                        rhob_data_for_sub = {}
                         # model.L is the index set for undesirable output variables
                        for l_range in model.L:
                            rhob_data_for_sub[l_range] = model.rhob[l_range].value
                        rhob_values[actual_index] = rhob_data_for_sub

                    # Collect theta if applicable (VRS1)
                    if hasattr(model, 'theta') and self.rts == RTS_VRS1:
                        theta_values[actual_index] = model.theta.value
                    objective_values[actual_index] = model.objective()

                    # Collect lamda (always collected in the original code for these cases)
                    lamda_data_for_dmu = {}
                    # model.R is likely the index set for reference DMUs
                    for r_range in model.R:
                         actual_ref_index = self.reference_data_index[r_range]
                         lamda_data_for_dmu[actual_ref_index] = model.lamda[r_range].value
                    lamda_values[actual_index] = lamda_data_for_dmu

                    # Collect mu if applicable (VRS2)
                    if hasattr(model, 'mu') and self.rts == RTS_VRS2:
                         mu_data_for_dmu = {}
                         # model.R is likely the index set for reference DMUs
                         for r_range in model.R:
                             actual_ref_index = self.reference_data_index[r_range]
                             mu_data_for_dmu[actual_ref_index] = model.mu[r_range].value
                         mu_values[actual_index] = mu_data_for_dmu

                except Exception as e:
                    # Handle errors during result retrieval
                    print(f"Warning: Could not retrieve results for DMU {actual_index} despite status '{optimization_status}': {e}")
                    # Indicate failure to retrieve specific values
                    if self.output_oriented or self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedyxb: 
                        rhoy_values[actual_index] = None
                    if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb: 
                        rhox_values[actual_index] = None
                    if self.unoutput_oriented or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb: 
                        rhob_values[actual_index] = None
                    if self.rts == RTS_VRS1: 
                        theta_values[actual_index] = np.nan # Use NaN for scalar theta
                    if self.rts == RTS_VRS2: 
                        mu_values[actual_index] = None # Use None for dict results
                    objective_values[actual_index] = None # Use None for dict results
                    lamda_values[actual_index] = None # Use None for dict results

            else: # Optimization status is not 'ok'
                # Mark all relevant results as failed/NaN for this DMU
                if self.output_oriented or self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedyxb: 
                    rhoy_values[actual_index] = None
                if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb: 
                    rhox_values[actual_index] = None
                if self.unoutput_oriented or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb: 
                    rhob_values[actual_index] = None
                if self.rts == RTS_VRS1: 
                    theta_values[actual_index] = np.nan
                if self.rts == RTS_VRS2: 
                    mu_values[actual_index] = None
                objective_values[actual_index] = None
                lamda_values[actual_index] = None


        # 4. Consolidated Post-Loop Processing and Storage
        self.results['optimization_status'] = all_optimization_statuses

        # Helper function to convert result dictionaries to DataFrames
        # This handles dictionaries where keys are evaluated DMUs and values are
        # dictionaries indexed by reference DMUs (for lamda, mu)
        def _process_ref_indexed_dict(results_dict):
            df_list = []
            # Need a consistent column index for failed rows. This is self.reference_data_index.
            for evaluated_index, data_for_dmu in results_dict.items():
                if data_for_dmu is not None:
                    # data_for_dmu is a dict {actual_ref_index: value}
                    series = pd.Series(data_for_dmu, name=evaluated_index)
                    df_list.append(series)
                else:
                    # For failed DMUs, add a row of NaNs with correct reference index columns
                    nan_series = pd.Series(np.nan, index=self.reference_data_index, name=evaluated_index)
                    df_list.append(nan_series)

            if df_list:
                # Concatenate and transpose to get evaluated DMUs as index
                return pd.concat(df_list, axis=1).T
            else:
                return None # No data

        # Helper function to convert result dictionaries to DataFrames
        # This handles dictionaries where keys are evaluated DMUs and values are
        # dictionaries indexed by variable indices (for rhoy, rhox, rhob)
        def _process_var_indexed_dict(results_dict, var_index_list):
            # print(results_dict,"results_dict")
            # print(var_index_list,"var_index_list")
            try:
                return pd.DataFrame(results_dict).T
            except:
                return None # No data


        # Process reference-indexed results (lamda, mu)
        self.results['lamda_df'] = _process_ref_indexed_dict(lamda_values)
        if self.rts == RTS_VRS2:
            self.results['mu_df'] = _process_ref_indexed_dict(mu_values)

        # Process variable-indexed results (rhoy, rhox, rhob)
        if self.output_oriented or self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedyxb:
             self.results['rhoy_df'] = _process_var_indexed_dict(rhoy_values, self.outputvars)
        if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb:
             self.results['rhox_df'] = _process_var_indexed_dict(rhox_values, self.inputvars)
        if self.unoutput_oriented or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb:
             self.results['rhob_df'] = _process_var_indexed_dict(rhob_values, self.unoutputvars)

        # Store theta directly (it's a dict of scalars)
        if self.rts == RTS_VRS1:
             self.results['theta'] = theta_values

        self.results['objective_values'] = objective_values


        # 5. Abstract Summary DataFrame Creation
        results_df = pd.DataFrame({
            'optimization_status': pd.Series(all_optimization_statuses),
        })

        # Define parameters and process rho based on orientation
        if self.output_oriented or self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedyxb:
            rho_df_summary = self.results.get('rhoy_df')
            vars_list = self.outputvars
            weights_list = self.wy
            te_formula = lambda x: 1 / (1 + x) if pd.notna(x) else np.nan
            aggregate_col_name = 'teo' # Assuming 'teo' for both output and hyper output part

            if rho_df_summary is not None:
                 # Ensure columns are named consistently and in correct order
                 rho_cols = [f"{self._beta_prefix}{i}" for i in vars_list]
                 te_cols = [f"te{i}" for i in vars_list]
                 wei_te_cols = [f"weight_te{i}" for i in vars_list]

                 # Rename columns for clarity in the output DataFrame
                 rho_df_summary.columns = rho_cols

                 # The per-variable scaling factors (beta) are always reported.
                 results_df = pd.concat([results_df, rho_df_summary[rho_cols]], axis=1)

                 # The efficiency decomposition (te*/weight/aggregate) is only
                 # emitted in the full NDDFweak variant, not in the minimal
                 # NDDFweakRho variant.
                 if not self._minimal:
                     for k, v in enumerate(weights_list):
                         rho_df_summary[te_cols[k]] = rho_df_summary[rho_cols[k]].apply(te_formula)
                         rho_df_summary[wei_te_cols[k]] = rho_df_summary[te_cols[k]] * v # Simpler calculation

                     # Calculate aggregate efficiency
                     rho_df_summary[aggregate_col_name] = rho_df_summary[wei_te_cols].sum(axis=1)

                     # Concatenate only the relevant new columns to results_df
                     results_df = pd.concat([results_df, rho_df_summary[te_cols + wei_te_cols + [aggregate_col_name]]], axis=1)

        if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb:
            rho_df_summary = self.results.get('rhox_df')
            vars_list = self.inputvars
            weights_list = self.wx
            te_formula = lambda x: (1 - x) if pd.notna(x) else np.nan # Assuming same formula as unoutput
            aggregate_col_name = 'tei' # Assuming 'tei' for both input and hyper input part

            if rho_df_summary is not None:
                 rho_cols = [f"{self._beta_prefix}{i}" for i in vars_list]
                 te_cols = [f"te{i}" for i in vars_list]
                 wei_te_cols = [f"weight_te{i}" for i in vars_list]

                 rho_df_summary.columns = rho_cols

                 # The per-variable scaling factors (beta) are always reported.
                 results_df = pd.concat([results_df, rho_df_summary[rho_cols]], axis=1)

                 if not self._minimal:
                     for k, v in enumerate(weights_list):
                         rho_df_summary[te_cols[k]] = rho_df_summary[rho_cols[k]].apply(te_formula)
                         rho_df_summary[wei_te_cols[k]] = rho_df_summary[te_cols[k]] * v

                     rho_df_summary[aggregate_col_name] = rho_df_summary[wei_te_cols].sum(axis=1)
                     results_df = pd.concat([results_df, rho_df_summary[te_cols + wei_te_cols + [aggregate_col_name]]], axis=1)


        if self.unoutput_oriented or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb:
            rho_df_summary = self.results.get('rhob_df')
            vars_list = self.unoutputvars
            weights_list = self.wb
            te_formula = lambda x: (1 - x) if pd.notna(x) else np.nan # Assuming same formula as input
            aggregate_col_name = 'teuo' # Assuming this name

            if rho_df_summary is not None:
                 rho_cols = [f"{self._beta_prefix}{i}" for i in vars_list]
                 te_cols = [f"te{i}" for i in vars_list]
                 wei_te_cols = [f"weight_te{i}" for i in vars_list]

                 rho_df_summary.columns = rho_cols

                 # The per-variable scaling factors (beta) are always reported.
                 results_df = pd.concat([results_df, rho_df_summary[rho_cols]], axis=1)

                 if not self._minimal:
                     for k, v in enumerate(weights_list):
                         rho_df_summary[te_cols[k]] = rho_df_summary[rho_cols[k]].apply(te_formula)
                         rho_df_summary[wei_te_cols[k]] = rho_df_summary[te_cols[k]] * v

                     rho_df_summary[aggregate_col_name] = rho_df_summary[wei_te_cols].sum(axis=1)
                     results_df = pd.concat([results_df, rho_df_summary[te_cols + wei_te_cols + [aggregate_col_name]]], axis=1)
        if not self._minimal:
            if self.hyper_orientedyb:
                results_df['teuo2o'] = results_df['teuo'] / results_df['teo']
            if self.hyper_orientedyx:
                results_df['tei2o'] = results_df['tei'] / results_df['teo']
            if self.hyper_orientedyxb:
                results_df['teiuo2o'] = 0.5*(results_df['tei']+results_df['teuo'] ) / results_df['teo']
            if self.hyper_orientedxb:
                results_df['teuo2i'] = results_df['teuo'] / results_df['tei']

        objective_df = self.get_objective_value()
        objective_df.name = 'objective_value'
        results_df = pd.concat([results_df, objective_df], axis=1)
        # results_df['1-beta'] = 1-results_df['objective_value']

        # Return the combined results DataFrame
        return results_df






    # --- Display and Get Methods ---
    # Update these to reflect that 'rho' is now split into rhox, rhoy, rhob

    def display_status(self):
        """Display the optimization status for each DMU."""
        if not self.results or 'optimization_status' not in self.results:
            print("Optimization has not been run yet.")
            return
        print("Optimization Status per DMU:")
        # Use the index of the results_df for consistent ordering and available DMUs
        results_df = self.get_results_df() # Get the summary df
        if results_df is not None:
             print(results_df['optimization_status'])
        else:
             print("No optimization status available.")


    def display_objective_value(self):
        """Display objective value for each DMU."""
        tools.assert_optimized(self.results)
        print("Objective values per DMU:")
        objective_series = pd.Series(self.results.get('objective_values', {}))
        print(objective_series)

    def get_objective_value(self):
        """Return objective values as a pandas Series."""
        tools.assert_optimized(self.results)
        return pd.Series(self.results.get('objective_values', {}))

    def display_theta(self):
        """Display theta value for each DMU (if RTS_VRS1)."""
        tools.assert_optimized(self.results)
        if self.rts != RTS_VRS1 or 'theta' not in self.results:
             print("Theta variable is only applicable and available for RTS_VRS1 after optimization.")
             return
        print("Theta values per DMU:")
        theta_series = pd.Series(self.results.get('theta', {}))
        print(theta_series)

    def get_theta(self):
        """Return theta values as a pandas Series (if RTS_VRS1)."""
        tools.assert_optimized(self.results)
        if self.rts != RTS_VRS1 or 'theta' not in self.results:
             print("Theta variable is only applicable and available for RTS_VRS1 after optimization.")
             return None
        return pd.Series(self.results.get('theta', {}))

    def display_lamda(self):
        """Display lamda values (intensity variables) for each DMU."""
        tools.assert_optimized(self.results)
        print("Lamda values per DMU:")
        lamda_df = self.results.get('lamda_df')
        if lamda_df is not None and not lamda_df.empty:
            print(lamda_df)
        else:
            print("No lamda values available (optimization may have failed for all DMUs).")

    def display_rhoy(self):
        """Display rhoy values (intensity variables) for each DMU."""
        tools.assert_optimized(self.results)
        print("rhoy values per DMU:")
        rhoy_df = self.results.get('rhoy_df')
        if rhoy_df is not None and not rhoy_df.empty:
            print(rhoy_df)
        else:
            print("No rhoy values available (optimization may have failed for all DMUs).")

    def get_rhoy(self):
        """Return rhoy values as a pandas DataFrame."""
        tools.assert_optimized(self.results)
        return self.results.get('rhoy_df')

    def display_rhob(self):
        """Display rhob values (intensity variables) for each DMU."""
        tools.assert_optimized(self.results)
        print("rhob values per DMU:")
        rhob_df = self.results.get('rhob_df')
        if rhob_df is not None and not rhob_df.empty:
            print(rhob_df)
        else:
            print("No rhob values available (optimization may have failed for all DMUs).")

    def get_rhob(self):
        """Return rhob values as a pandas DataFrame."""
        tools.assert_optimized(self.results)
        return self.results.get('rhob_df')

    def display_rhox(self):
        """Display rhox values (intensity variables) for each DMU."""
        tools.assert_optimized(self.results)
        print("rhox values per DMU:")
        rhox_df = self.results.get('rhox_df')
        if rhox_df is not None and not rhox_df.empty:
            print(rhox_df)
        else:
            print("No rhox values available (optimization may have failed for all DMUs).")

    def get_rhox(self):
        """Return rhox values as a pandas DataFrame."""
        tools.assert_optimized(self.results)
        return self.results.get('rhox_df')


    def get_lamda(self):
        """Return lamda values as a pandas DataFrame."""
        tools.assert_optimized(self.results)
        return self.results.get('lamda_df')

    def display_mu(self):
        """Display mu values (intensity variables) for each DMU (if RTS_VRS2)."""
        tools.assert_optimized(self.results)
        if self.rts != RTS_VRS2 or 'mu_df' not in self.results:
             print("Mu variables are only applicable and available for RTS_VRS2 after optimization.")
             return
        print("Mu values per DMU:")
        mu_df = self.results.get('mu_df')
        if mu_df is not None and not mu_df.empty:
            print(mu_df)
        else:
            print("No mu values available (optimization may have failed for all DMUs).")

    def get_mu(self):
        """Return mu values as a pandas DataFrame (if RTS_VRS2)."""
        tools.assert_optimized(self.results)
        if self.rts != RTS_VRS2 or 'mu_df' not in self.results:
             print("Mu variables are only applicable and available for RTS_VRS2 after optimization.")
             return None
        return self.results.get('mu_df')


    def get_results_df(self):
        """Return the primary results DataFrame including expanded rho columns."""
        tools.assert_optimized(self.results)

        if not self.results or 'optimization_status' not in self.results:
             return None # No results yet

        # Start with scalar results
        results_data = {
            'optimization_status': pd.Series(self.results.get('optimization_status', {})),
            'objective_value': pd.Series(self.results.get('objective_value', {})),
        }

        if self.rts == RTS_VRS1 and 'theta' in self.results:
             results_data['theta'] = pd.Series(self.results.get('theta', {}))

        # Create the initial DataFrame with scalar results
        results_df = pd.DataFrame(results_data)

        # Get the DataFrame containing the expanded rho columns
        rho_df = self.get_rho()

        # Concatenate the scalar results DataFrame and the rho DataFrame if rho_df exists
        if rho_df is not None:
            # Ensure both DataFrames have the same index (evaluated_data_index)
            results_df = results_df.reindex(self.evaluated_data_index)
            rho_df = rho_df.reindex(self.evaluated_data_index) # Should already have this index from get_rho

            # Join them side by side
            results_df = pd.concat([results_df, rho_df], axis=1)

        return results_df


    def info(self, dmu="all"):
        """Show the information of the lp model for specified DMU(s).

        Args:
            dmu (str or list): The actual index label of the DMU(s) to display, or "all". Default is "all".
        """
        if not self.__modeldict:
            print("No models have been initialized.")
            return

        if dmu == "all":
            print("Displaying all DMU models:")
            # Iterate through the dictionary using actual index labels
            for ind, problem in self.__modeldict.items():
                print(f"\n--- Model for DMU: {ind} ---")
                try:
                    problem.pprint()
                except Exception as e:
                    print(f"Could not print model for {ind}: {e}")
                print("-" * (len(f"--- Model for DMU: {ind} ---")))
        else:
            if isinstance(dmu, str):
                dmu_list = [dmu]
            else:
                dmu_list = dmu

            for ind in dmu_list:
                # Check if the actual index exists in the model dictionary
                if ind in self.__modeldict:
                    print(f"\n--- Model for DMU: {ind} ---")
                    try:
                        self.__modeldict[ind].pprint()
                    except Exception as e:
                         print(f"Could not print model for {ind}: {e}")
                    print("-" * (len(f"--- Model for DMU: {ind} ---")))
                else:
                    print(f"DMU '{ind}' not found in the evaluated set.")


class NDDFweakRho(NDDFweak):
    """Nonradial DDF that reports only the raw LP value (``objective_value``)
    plus each variable's scaling factor ``beta`` as columns ``beta_<varname>``,
    mirroring ``DEAweak.ddf`` / ``DEA.ddf`` (which leave efficiency decomposition
    to the caller).

    In the default direction (``gx=0, gy=1, gb=1``) the returned beta columns are
    ``beta_Y`` and ``beta_CO2``.  Environmental efficiency ``EB`` and output
    efficiency ``EO`` are derived at the call site.  Under the equal group-weight
    specification used by the standard case (``wy = wb = 0.5``):

        EB = 0.5 * (1 - beta_CO2),   EO = 0.5 / (1 + beta_Y),   EBO = EB / EO.

    These coincide with the legacy ``teuo`` / ``teo`` / ``teuo2o`` columns that
    this class no longer emits.
    """

    def __init__(self, data, sent, gy=[1], gx=[0], gb=[1], rts=RTS_VRS1,
                 evalindex=None, refindex=None):
        super().__init__(data, sent, gy=gy, gx=gx, gb=gb, rts=rts,
                         evalindex=evalindex, refindex=refindex)
        self._minimal = True
        self._beta_prefix = "beta_"
