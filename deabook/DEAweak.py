# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint
import numpy as np
import pandas as pd
from .constant import CET_ADDI, RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from .utils import tools
from .DEA import DEA,DEA2,DDF,DDF2
from .core import (
    PyomoResultMixin,
    build_weak_direction_flags,
    format_ref_indexed_results,
    format_var_indexed_results,
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

        ``DEAweak2`` and ``DDFweak2`` differ in their formulas and efficiency
        columns, but they share the whole solver/result-collection lifecycle.
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

class _WeakDualResultMixin(PyomoResultMixin):
    """Shared display/get accessors for weak-disposability dual models."""

    def display_delta(self):
        tools.assert_optimized(self.optimization_status)
        self._display_component("delta")

    def display_gamma(self):
        tools.assert_optimized(self.optimization_status)
        self._display_component("gamma")

    def display_kappa(self):
        tools.assert_optimized(self.optimization_status)
        self._display_component("kappa")

    def display_alpha(self):
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale_alpha(self.rts)
        self._display_component("alpha")

    def get_delta(self):
        tools.assert_optimized(self.optimization_status)
        return self._component_array("delta")

    def get_gamma(self):
        tools.assert_optimized(self.optimization_status)
        return self._component_array("gamma")

    def get_kappa(self):
        tools.assert_optimized(self.optimization_status)
        return self._component_array("kappa")

    def get_alpha(self):
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale_alpha(self.rts)
        return self._component_series("alpha")



class DEAweak(PyomoResultMixin, DEA):
    """weak dispsbnility of Data Envelopment Analysis (DEA)
    """

    def __init__(self, data, sent, gy, gx,gb , rts, baseindex=None, refindex=None):
        """DEA: Envelopment problem

        Args:
            data
            sent
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.b,  self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, baseindex, refindex)
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()
        self.__model__.R = Set(initialize=range(len(self.yref)))

        # Initialize sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable
        self.__model__.rho = Var(self.__model__.I, doc='efficiency')
        if self.rts == RTS_VRS1:
            self.__model__.theta = Var(self.__model__.I, bounds=(0.0, 1.0),doc='emission reduction factor')
        elif self.rts == RTS_VRS2:
            self.__model__.mu = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None),doc='emission reduction factor2')

        self.__model__.lamda = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None), doc='intensity variables')

        # Setup the objective function and constraints
        if sum(self.gy) >= 1:
            self.__model__.objective = Objective(
                rule=self.__objective_rule(), sense=maximize, doc='objective function')
        else:
            self.__model__.objective = Objective(
                rule=self.__objective_rule(), sense=minimize, doc='objective function')

        self.__model__.input = Constraint(
            self.__model__.I, self.__model__.J, rule=self.__input_rule(), doc='input constraint')
        self.__model__.output = Constraint(
            self.__model__.I, self.__model__.K, rule=self.__output_rule(), doc='output constraint')
        self.__model__.unoutput = Constraint(
            self.__model__.I, self.__model__.L, rule=self.__unoutput_rule(), doc='undesirable output constraint')
        if self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
            self.__model__.vrs = Constraint(self.__model__.I, rule=self.__vrs_rule(), doc='variable return to scale rule')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.rho[i] for i in model.I)
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        if sum(self.gx)>=1:
            if self.rts == RTS_VRS2:
                def input_rule(model, o, j):
                    if self.gx[j]==1:
                        return sum((model.lamda[o, r]+model.mu[o,r])*self.xref[r][j] for r in model.R)<=model.rho[o]*self.x[o][j]
                    else:
                        return sum((model.lamda[o, r]+model.mu[o,r])*self.xref[r][j] for r in model.R)<=self.x[o][j]

                return input_rule
            else:
                def input_rule(model, o, j):
                    if self.gx[j]==1:
                        return sum(model.lamda[o, r]*self.xref[r][j] for r in model.R)<=model.rho[o]*self.x[o][j]
                    else:
                        return sum(model.lamda[o, r]*self.xref[r][j] for r in model.R)<=self.x[o][j]
                return input_rule
        else:
            if self.rts == RTS_VRS2:
                def input_rule(model, o, j):
                    return sum((model.lamda[o, r]+model.mu[o,r])*self.xref[r][j] for r in model.R)<=self.x[o][j]
                return input_rule
            else:
                def input_rule(model, o, j):
                    return sum(model.lamda[o, r] * self.xref[r][j] for r in model.R) <= self.x[o][j]

                return input_rule



    def __output_rule(self):
        """Return the proper output constraint"""
        if sum(self.gy)>=1:

            def output_rule(model, o, k):
                if self.gy[k] == 1:
                    return sum(model.lamda[o, r]*self.yref[r][k] for r in model.R) >= model.rho[o]*self.y[o][k]
                else:
                    return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k]

            return output_rule

        else:
            def output_rule(model, o, k):
                return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k]
            return output_rule

    def __unoutput_rule(self):
        """Return the proper undesirable output constraint"""
        if sum(self.gb)>=1:
            def unoutput_rule(model, o, l):
                if self.gb[l] == 1:
                    return sum(model.lamda[o, r]*self.bref[r][l] for r in model.R) == model.rho[o]*self.b[o][l]
                else:
                    return sum(model.lamda[o, r] * self.bref[r][l] for r in model.R) == self.b[o][l]

            return unoutput_rule
        else:
            def unoutput_rule(model, o, l):
                return sum(model.lamda[o, r] * self.bref[r][l] for r in model.R) == self.b[o][l]
            return unoutput_rule

    def __vrs_rule(self):
        if self.rts==RTS_VRS1:
            def vrs_rule(model, o):
                return sum(model.lamda[o, r] for r in model.R) == model.theta[o]
            return vrs_rule

        elif self.rts==RTS_VRS2:
            def vrs_rule(model, o):
                return sum((model.lamda[o, r]+model.mu[o, r] )for r in model.R) == 1

            return vrs_rule

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization
        self.problem_status, self.optimization_status = tools.optimize_model(
            self._single_model(), email, CET_ADDI, solver)

    def display_status(self):
        """Display the status of problem"""
        print(self.optimization_status)

    def display_theta(self):
        """Display theta value"""
        tools.assert_optimized(self.optimization_status)
        self._display_component("theta")

    def display_rho(self):
        """Display rho value"""
        tools.assert_optimized(self.optimization_status)
        self._display_component("rho")

    def display_lamda(self):
        """Display lamda value"""
        tools.assert_optimized(self.optimization_status)
        self._display_component("lamda")

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_theta(self):
        """Return theta value by array"""
        tools.assert_optimized(self.optimization_status)
        return self._component_series("theta")

    def get_rho(self):
        """Return rho value by array"""
        tools.assert_optimized(self.optimization_status)
        return self._component_series("rho")

    def get_lamda(self):
        """Return lamda value by array"""
        tools.assert_optimized(self.optimization_status)
        return self._component_array("lamda")


class DEAweak2(_WeakOrientationMixin, _WeakPerDmuResultMixin, DEA2):
    """Data Envelopment Analysis (DEA) - Solves per DMU"""

    def __init__(self,data,sent,gy=[1],gx=[0],gb=[1],rts=RTS_VRS1,baseindex=None, refindex=None):
        """DEA: Envelopment problem, solving for each DMU individually.

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y:CO2"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEAweak should return numpy arrays for x, y,b, xref, yref,bref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DEAweak is modified or a similar function returns these:
        try:
             self.y, self.x,self.b, self.yref,self.xref,self.bref, self.gy, self.gx,self.gb, \
                 self.evaluated_data_index, self.reference_data_index = tools.assert_DEAweak_with_indices(
                     data, sent, gy, gx,gb, baseindex, refindex
                 )
        except AttributeError:
             # Fallback if assert_DEA_with_indices doesn't exist, assuming original assert_DEA
             # In this case, evaluated_data_index and reference_data_index will be ranges 0..N-1
             print("Warning: tools.assert_DEAweak_with_indices not found. Using range indices.")
             self.y,self.x,self.b, self.yref,self.xref,self.bref, self.gy,self.gx,self.gb = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, baseindex, refindex)
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
        """Optimize every evaluated DMU and return legacy DEAweak2 columns."""
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



class DDFweak(DEAweak):
    def __init__(self,  data, sent, gy=[1], gx=[1], gb=[1], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DDFweak: Directional distance function with undesirable output
        
        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        self.y, self.x, self.b,  self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DDFweak(data,sent, gy, gx, gb,baseindex,refindex)
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()

        # Initialize sets
        self.__model__.R = Set(initialize=range(len(self.yref)))
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable

        self.__model__.rho = Var(self.__model__.I, doc='directional distance')
        if self.rts == RTS_VRS1:
            self.__model__.theta = Var(self.__model__.I, bounds=(0.0, 1.0),doc='emission reduction factor')
        elif self.rts == RTS_VRS2:
            self.__model__.mu = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None),
                                       doc='emission reduction factor2')

        self.__model__.lamda = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None), doc='intensity variables')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self._DEAweak__objective_rule(), sense=maximize, doc='objective function')
        self.__model__.input = Constraint(
            self.__model__.I, self.__model__.J, rule=self.__input_rule(), doc='input constraint')
        self.__model__.output = Constraint(
            self.__model__.I, self.__model__.K, rule=self.__output_rule(), doc='output constraint')
        self.__model__.undesirable_output = Constraint(
            self.__model__.I, self.__model__.L, rule=self.__undesirable_output_rule(), doc='undesirable output constraint')

        if self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
            self.__model__.vrs = Constraint(self.__model__.I, rule=self.__vrs_rule(), doc='various return to scale rule')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __input_rule(self):
        """Return the proper input constraint"""
        if self.rts == RTS_VRS2:

            def input_rule(model, o, j):
                return sum((model.lamda[o, r]+model.mu[o, r]) * self.xref[r][j] for r in model.R) \
                            <= self.x[o][j] - model.rho[o]*self.gx[j]*self.x[o][j]
            return input_rule
        else:

            def input_rule(model, o, j):
                return sum(model.lamda[o, r] * self.xref[r][j] for r in model.R) \
                            <= self.x[o][j] - model.rho[o]*self.gx[j]*self.x[o][j]
            return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, o, k):
            return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k] + model.rho[o]*self.gy[k]*self.y[o][k]
        return output_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, o, l):
            return sum(model.lamda[o, r] * self.bref[r][l] for r in model.R) == self.b[o][l] - model.rho[o]*self.gb[l]*self.b[o][l]
        return undesirable_output_rule

    def __vrs_rule(self):
        """Return the VRS constraint"""
        if self.rts == RTS_VRS1:

            def vrs_rule(model, o):
                return sum(model.lamda[o, r] for r in model.R) == model.theta[o]
            return vrs_rule
        elif self.rts == RTS_VRS2:

            def vrs_rule(model, o):
                return sum((model.lamda[o, r]+model.mu[o, r] )for r in model.R) == 1
            return vrs_rule


class DDFweak2(_WeakOrientationMixin, _WeakPerDmuResultMixin, DDF2):

    def __init__(self,data,sent,gy=[1],gx=[0],gb=[1],rts=RTS_VRS1,baseindex=None, refindex=None):
        """ DDF: Envelopment problem, solving for each DMU individually.
        
        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y:CO2"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEAweak should return numpy arrays for x, y,b, xref, yref,bref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DEAweak is modified or a similar function returns these:
        try:
             self.y, self.x,self.b, self.yref,self.xref,self.bref, self.gy, self.gx,self.gb, \
                 self.evaluated_data_index, self.reference_data_index = tools.assert_DEAweak_with_indices(
                     data, sent, gy, gx,gb, baseindex, refindex
                 )
        except AttributeError:
             # Fallback if assert_DEA_with_indices doesn't exist, assuming original assert_DEA
             # In this case, evaluated_data_index and reference_data_index will be ranges 0..N-1
             print("Warning: tools.assert_DEAweak_with_indices not found. Using range indices.")
             self.y,self.x,self.b, self.yref,self.xref,self.bref, self.gy,self.gx,self.gb = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, baseindex, refindex)
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
        """Optimize every evaluated DMU and return legacy DDFweak2 columns."""
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



class NDDFweak2(_WeakOrientationMixin, DDF2):

    def __init__(self,data,sent,gy=[1],gx=[0],gb=[1],rts=RTS_VRS1,baseindex=None, refindex=None):
        """ DDF: Envelopment problem, solving for each DMU individually.
        
        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y:CO2"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEAweak should return numpy arrays for x, y,b, xref, yref,bref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DEAweak is modified or a similar function returns these:
        self.y, self.x,self.b, self.outputvars,self.inputvars,self.unoutputvars, \
            self.yref,self.xref,self.bref, self.gy, self.gx,self.gb, self.wx,self.wy,self.wb, \
            self.evaluated_data_index, self.reference_data_index = tools.assert_NDDFweak_with_indices(
                data, sent, gy, gx,gb, baseindex, refindex
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
                 rho_cols = [f"rho{i}" for i in vars_list]
                 te_cols = [f"te{i}" for i in vars_list]
                 wei_te_cols = [f"weight_te{i}" for i in vars_list]

                 # Rename columns for clarity in the output DataFrame
                 rho_df_summary.columns = rho_cols

                 # Calculate derived columns
                 for k, v in enumerate(weights_list):
                     rho_df_summary[te_cols[k]] = rho_df_summary[rho_cols[k]].apply(te_formula)
                     rho_df_summary[wei_te_cols[k]] = rho_df_summary[te_cols[k]] * v # Simpler calculation

                 # Calculate aggregate efficiency
                 rho_df_summary[aggregate_col_name] = rho_df_summary[wei_te_cols].sum(axis=1)

                 # Concatenate only the relevant new columns to results_df
                 results_df = pd.concat([results_df, rho_df_summary[rho_cols + te_cols + wei_te_cols + [aggregate_col_name]]], axis=1)

        if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb:
            rho_df_summary = self.results.get('rhox_df')
            vars_list = self.inputvars
            weights_list = self.wx
            te_formula = lambda x: (1 - x) if pd.notna(x) else np.nan # Assuming same formula as unoutput
            aggregate_col_name = 'tei' # Assuming 'tei' for both input and hyper input part

            if rho_df_summary is not None:
                 rho_cols = [f"rho{i}" for i in vars_list]
                 te_cols = [f"te{i}" for i in vars_list]
                 wei_te_cols = [f"weight_te{i}" for i in vars_list]

                 rho_df_summary.columns = rho_cols

                 for k, v in enumerate(weights_list):
                     rho_df_summary[te_cols[k]] = rho_df_summary[rho_cols[k]].apply(te_formula)
                     rho_df_summary[wei_te_cols[k]] = rho_df_summary[te_cols[k]] * v

                 rho_df_summary[aggregate_col_name] = rho_df_summary[wei_te_cols].sum(axis=1)
                 results_df = pd.concat([results_df, rho_df_summary[rho_cols + te_cols + wei_te_cols + [aggregate_col_name]]], axis=1)


        if self.unoutput_oriented or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb:
            rho_df_summary = self.results.get('rhob_df')
            vars_list = self.unoutputvars
            weights_list = self.wb
            te_formula = lambda x: (1 - x) if pd.notna(x) else np.nan # Assuming same formula as input
            aggregate_col_name = 'teuo' # Assuming this name

            if rho_df_summary is not None:
                 rho_cols = [f"rho{i}" for i in vars_list]
                 te_cols = [f"te{i}" for i in vars_list]
                 wei_te_cols = [f"weight_te{i}" for i in vars_list]

                 rho_df_summary.columns = rho_cols

                 for k, v in enumerate(weights_list):
                     rho_df_summary[te_cols[k]] = rho_df_summary[rho_cols[k]].apply(te_formula)
                     rho_df_summary[wei_te_cols[k]] = rho_df_summary[te_cols[k]] * v

                 rho_df_summary[aggregate_col_name] = rho_df_summary[wei_te_cols].sum(axis=1)
                 results_df = pd.concat([results_df, rho_df_summary[rho_cols + te_cols + wei_te_cols + [aggregate_col_name]]], axis=1)
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





class DEAweakDUAL(_WeakDualResultMixin, DEAweak):

    def __init__(self, data, sent, gy, gx,gb , rts, baseindex=None, refindex=None):
        """DEA: Envelopment problem

        Args:
            data
            sent
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.b,  self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, baseindex, refindex)
        self.rts = rts


        # Initialize DEA model
        self.__model__ = ConcreteModel()

        # Initialize sets
        self.__model__.R = Set(initialize=range(len(self.yref)))
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable
        self.__model__.delta = Var(self.__model__.I, self.__model__.J, bounds=(
            0.0, None), doc='multiplier x')
        self.__model__.gamma = Var(self.__model__.I, self.__model__.K, bounds=(
            0.0, None), doc='multiplier y')
        self.__model__.kappa = Var(self.__model__.I, self.__model__.L,  doc='multiplier b')

        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            0.0, None), doc='variable return to scale')
        elif self.rts == RTS_VRS2:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            None, None), doc='variable return to scale')
        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self.__objective_rule(), sense=minimize, doc='objective function')

        self.__model__.first = Constraint(
            self.__model__.I, self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
        self.__model__.second = Constraint(
            self.__model__.I,                   rule=self.__second_rule(), doc='normalization constraint')
        if self.rts == RTS_VRS2:
            self.__model__.third = Constraint(
                self.__model__.I, rule=self.__third_rule(), doc='weak disposability constraint')
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(sum(model.delta[o, j] * self.x[o][j]* (1-self.gx[j]) for o in model.I) for j in model.J) - \
                sum(sum(model.gamma[o, k] * self.y[o][k]* (1-self.gb[k]) for o in model.I) for k in model.K) +\
                sum(sum(model.kappa[o,l] * self.b[o][l]* (1-self.gb[l]) for o in model.I) for l in model.L)

        return objective_rule


    def __first_rule(self):
        """Return the proper technology constraint"""
        if self.rts == RTS_VRS1:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) - model.alpha[o] >= 0
            return first_rule

        elif self.rts == RTS_VRS2:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) + model.alpha[o] >= 0
            return first_rule

        elif self.rts == RTS_CRS:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) >= 0
            return first_rule

    def __second_rule(self):
        """Return the proper normalization constraint"""
        if sum(self.gx) >= 1:
            def second_rule(model, o):
                return sum(model.delta[o, j] * self.x[o][j]* self.gx[j] for j in model.J) == 1
            return second_rule
        elif sum(self.gy) >= 1:

            def second_rule(model, o):
                return sum(model.gamma[o, k] * self.y[o][k]* self.gy[k] for k in model.K) == 1
            return second_rule
        elif sum(self.gb) >= 1:

            def second_rule(model, o):
                return sum(model.kappa[o, l] * self.b[o][l]* self.gb[l] for l in model.L) == 1
            return second_rule

    def __third_rule(self):
        """Return the proper weak disposability constraint"""

        def third_rule(model, o, r):
            return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) + model.alpha[o] >= 0

        return third_rule

    def get_efficiency(self):
        """Return efficiency value by array"""
        tools.assert_optimized(self.optimization_status)
        if sum(self.gx) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1)
            else:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1) + self.get_alpha().reshape(len(self.y), 1)
        if sum(self.gy) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1)
            else:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1) + self.get_alpha().reshape(len(self.x), 1)
        if sum(self.gb) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1)
            else:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1) + self.get_alpha().reshape(len(self.y), 1)




class DDFweakDUAL(DEAweakDUAL):
    def __init__(self, y, x, b, gy=[1], gx=[1], gb=[1], rts=RTS_VRS1, yref=None, xref=None, bref=None):
        """DDFweak: Directional distance function with undesirable output

        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        if yref is None:
            yref = y
        if xref is None:
            xref = x
        if bref is None:
            bref = b
        self.y, self.x, self.b, self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DDFweak1(y, x, b, gy, gx, gb, yref, xref, bref)
        self.rts = rts


        # Initialize DEA model
        self.__model__ = ConcreteModel()
        self.__model__.R = Set(initialize=range(len(self.yref)))

        # Initialize sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable
        self.__model__.delta = Var(self.__model__.I, self.__model__.J, bounds=(
            0.0, None), doc='multiplier x')
        self.__model__.gamma = Var(self.__model__.I, self.__model__.K, bounds=(
            0.0, None), doc='multiplier y')
        self.__model__.kappa = Var(self.__model__.I, self.__model__.L,  doc='multiplier b')

        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            0.0, None), doc='variable return to scale')
        elif self.rts == RTS_VRS2:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            None, None), doc='variable return to scale')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self.__objective_rule(), sense=minimize, doc='objective function')

        self.__model__.first = Constraint(
            self.__model__.I, self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
        self.__model__.second = Constraint(
            self.__model__.I, rule=self.__second_rule(), doc='normalization constraint')
        if self.rts == RTS_VRS2:
            self.__model__.third = Constraint(
                self.__model__.I, rule=self.__third_rule(), doc='weak disposability constraint')
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        if self.rts == RTS_VRS1:

            def objective_rule(model):
                return sum(sum(model.delta[o, j] * self.x[o][j] for o in model.I) for j in model.J) -\
                    sum(sum(model.gamma[o, k] * self.y[o][k] for o in model.I) for k in model.K) + \
                    sum(sum(model.kappa[o, l] * self.b[o][l] for o in model.I) for l in model.L)
            return objective_rule

        elif self.rts == RTS_VRS2:

            def objective_rule(model):
                return sum(sum(model.delta[o, j] * self.x[o][j] for o in model.I) for j in model.J) -\
                    sum(sum(model.gamma[o, k] * self.y[o][k] for o in model.I) for k in model.K) + \
                    sum(sum(model.kappa[o, l] * self.b[o][l] for o in model.I) for l in model.L) + \
                    sum(model.alpha[o] for o in model.I)
            return objective_rule
    def __first_rule(self):
        """Return the proper technology constraint"""
        if self.rts == RTS_VRS1:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) - model.alpha[o] >= 0
            return first_rule
        elif self.rts == RTS_VRS2:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) + model.alpha[o] >= 0
            return first_rule
        elif self.rts == RTS_CRS:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) >= 0
            return first_rule


    def __second_rule(self):
        """Return the proper normalization constraint"""
        def second_rule(model, o):
            return sum(model.delta[o, j] *self.gx[j]* self.x[o][j] for j in model.J) +\
                sum(model.gamma[o, k] *self.gy[k]* self.y[o][k] for k in model.K) +\
                sum(model.kappa[o, l] *self.gb[l]* self.b[o][l] for l in model.L) == 1
        return second_rule

    def __third_rule(self):
        """Return the proper weak disposability constraint"""

        def third_rule(model, o, r):
            return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) + model.alpha[o] >= 0

        return third_rule

