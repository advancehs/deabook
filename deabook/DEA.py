# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint
import numpy as np
import pandas as pd
from .constant import CET_ADDI, RTS_VRS1, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from .utils import tools


class DEA:
    """Data Envelopment Analysis (DEA)
    """

    def __init__(self, data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Envelopment problem

        Args:
            data
            sent
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.yref, self.xref, self.gy, self.gx ,data_index= tools.assert_DEA(data, sent, gy, gx,baseindex,refindex)
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()

        # Initialize sets
        self.__model__.R = Set(initialize=range(len(self.yref)))
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize variable
        self.__model__.rho = Var(self.__model__.I, doc='efficiency')
        self.__model__.lamda = Var(self.__model__.I, self.__model__.R, bounds=(
            0.0, None), doc='intensity variables')

        # Setup the objective function and constraints
        if sum(self.gx) >= 1:
            self.__model__.objective = Objective(
                rule=self.__objective_rule(), sense=minimize, doc='objective function')
        elif sum(self.gy) >=1:
            self.__model__.objective = Objective(
                rule=self.__objective_rule(), sense=maximize, doc='objective function')
        else:
            raise ValueError("gx and gy must either be 1 or 0.")
        # self.__model__.objective.pprint()

        self.__model__.input = Constraint(
            self.__model__.I, self.__model__.J, rule=self.__input_rule(), doc='input constraint')
        # self.__model__.input.pprint()
        self.__model__.output = Constraint(
            self.__model__.I, self.__model__.K, rule=self.__output_rule(), doc='output constraint')
        # self.__model__.output.pprint()

        if self.rts == RTS_VRS1:
            self.__model__.vrs = Constraint(
                self.__model__.I, rule=self.__vrs_rule(), doc='variable return to scale rule')

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
        if sum(self.gx) >= 1:
            def input_rule(model, o, j):
                if self.gx[j] == 1:

                    return (sum(model.lamda[o, r]*self.xref[r][j] for r in model.R) <= \
                            model.rho[o]*self.x[o][j]   )
                else:
                    return (sum(model.lamda[o, r]*self.xref[r][j] for r in model.R) <= \
                            self.x[o][j] )
            return input_rule
        elif sum(self.gy) >=1:
            def input_rule(model, o, j):
                return sum(model.lamda[o, r] * self.xref[r][j] for r in model.R) <= self.x[o][j]
            return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        if sum(self.gx) >= 1:
            def output_rule(model, o, k):
                return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k]
            return output_rule
        elif sum(self.gy) >= 1:
            def output_rule(model, o, k):
                if sum(self.gy) >= 1:

                    return (sum(model.lamda[o, r]*self.yref[r][k] for r in model.R) >= \
                            model.rho[o]*self.y[o][k]   )
                else:
                    return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k]

            return output_rule

    def __vrs_rule(self):
        def vrs_rule(model, o):
            return sum(model.lamda[o, r] for r in model.R) == 1
        return vrs_rule

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization
        self.problem_status, self.optimization_status = tools.optimize_model(
            self.__model__, email, CET_ADDI, solver)

    def display_status(self):
        """Display the status of problem"""
        print(self.optimization_status)

    def display_rho(self):
        """Display rho value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.rho.display()

    def display_lamda(self):
        """Display lamda value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.lamda.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_rho(self):
        """Return rho value by array"""
        tools.assert_optimized(self.optimization_status)
        rho = list(self.__model__.rho[:].value)
        return np.asarray(rho)

    def get_lamda(self):
        """Return lamda value by array"""
        tools.assert_optimized(self.optimization_status)
        lamda = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.lamda),
                                                          list(self.__model__.lamda[:, :].value))])
        lamda = pd.DataFrame(lamda, columns=['Name', 'Key', 'Value'])
        lamda = lamda.pivot(index='Name', columns='Key', values='Value')
        return lamda.to_numpy()







# ---------------------------------------------------------------------------
# OOP refactor for the public per-DMU models used throughout the book.
# ---------------------------------------------------------------------------

from .core import (
    SolverConfig,
    build_direction_spec,
    prepare_dea_data,
    pyomo_value,
    solve_pyomo_model,
    validate_rts,
)


class _PerDmuEnvelopmentModel:
    """Base class for one-model-per-DMU DEA estimators.

    Subclasses implement only the model-specific variable, objective and
    constraint rules.  Data parsing, sample/reference splitting, solver calls,
    result tables and ``info``/getter methods are shared here.
    """

    allowed_rts = {RTS_CRS, RTS_VRS1}
    result_value_name = "rho"

    def __init__(self, data, sent, gy, gx, rts=RTS_VRS1, baseindex=None, refindex=None):
        self.data = data
        self.sent = sent
        self.rts = validate_rts(rts, self.allowed_rts)
        self.baseindex = baseindex
        self.refindex = refindex
        self.prepared = prepare_dea_data(
            data,
            sent,
            eval_query=baseindex,
            ref_query=refindex,
            input_side="left",
        )
        self.formula = self.prepared.formula
        self.xcol = self.formula.x
        self.ycol = self.formula.y
        self.y = self.prepared.y
        self.x = self.prepared.x
        self.yref = self.prepared.yref
        self.xref = self.prepared.xref
        self.evaluated_data_index = self.prepared.evaluated_index
        self.reference_data_index = self.prepared.reference_index
        self.direction = build_direction_spec(
            gx,
            gy,
            None,
            num_inputs=len(self.xcol),
            num_outputs=len(self.ycol),
        )
        self.gx = self.direction.gx
        self.gy = self.direction.gy
        self.input_oriented = self.direction.orientation == "input"
        self.output_oriented = self.direction.orientation == "output"
        self.hyper_oriented = self.direction.orientation == "hyper_input_output"
        self.evaluated_indices_range = range(len(self.evaluated_data_index))
        self.reference_indices_range = range(len(self.reference_data_index))
        self.num_inputs = len(self.xcol)
        self.num_outputs = len(self.ycol)
        self._models = {}
        self.results = {}
        for position, actual_index in enumerate(self.evaluated_data_index):
            self._models[actual_index] = self._build_model(position, actual_index)

    def _build_model(self, current_dmu_range_index, actual_index):
        model = ConcreteModel()
        model.R = Set(initialize=self.reference_indices_range)
        model.J = Set(initialize=range(self.num_inputs))
        model.K = Set(initialize=range(self.num_outputs))
        self._add_decision_variables(model, actual_index)
        self._add_objective(model)
        model.input = Constraint(
            model.J,
            rule=self._create_input_rule(current_dmu_range_index),
            doc="input constraint",
        )
        model.output = Constraint(
            model.K,
            rule=self._create_output_rule(current_dmu_range_index),
            doc="output constraint",
        )
        if self.rts == RTS_VRS1:
            model.vrs = Constraint(rule=self._create_vrs_rule(), doc="variable return to scale rule")
        return model

    def _add_decision_variables(self, model, actual_index):
        raise NotImplementedError

    def _add_objective(self, model):
        raise NotImplementedError

    def _create_input_rule(self, current_dmu_range_index):
        raise NotImplementedError

    def _create_output_rule(self, current_dmu_range_index):
        raise NotImplementedError

    def _create_vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[r] for r in model.R) == 1
        return vrs_rule

    def _efficiency_columns(self, result_df):
        return result_df

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Solve every evaluated DMU and return a result table.

        Parameters
        ----------
        email : str
            ``OPT_LOCAL`` for local solving; otherwise a NEOS email address.
        solver : str or None
            Pyomo solver name. ``None`` keeps the historical default solver.

        Returns
        -------
        pandas.DataFrame
            Rows use the original evaluated-data index and include
            ``optimization_status``, ``rho`` and model-specific efficiency
            columns.
        """
        status_values = {}
        objective_values = {}
        scalar_values = {}
        lamda_values = {}
        config = SolverConfig(solver=solver, email=email, tee=False)
        for actual_index, model in self._models.items():
            try:
                outcome = solve_pyomo_model(model, config=config, strict=False)
                status_values[actual_index] = outcome.legacy_status
                if outcome.ok:
                    scalar_values[actual_index] = pyomo_value(getattr(model, self.result_value_name))
                    objective_values[actual_index] = pyomo_value(model.objective)
                    lamda_values[actual_index] = {
                        self.reference_data_index[r]: pyomo_value(model.lamda[r])
                        for r in model.R
                    }
                else:
                    scalar_values[actual_index] = np.nan
                    objective_values[actual_index] = np.nan
                    lamda_values[actual_index] = None
            except Exception as exc:
                status_values[actual_index] = f"error: {exc}"
                scalar_values[actual_index] = np.nan
                objective_values[actual_index] = np.nan
                lamda_values[actual_index] = None

        self.results = {
            "optimization_status": status_values,
            self.result_value_name: scalar_values,
            "objective_value": objective_values,
            "lamda_df": self._format_lamda(lamda_values),
        }
        result_df = pd.DataFrame(
            {
                "optimization_status": pd.Series(status_values),
                self.result_value_name: pd.Series(scalar_values, dtype="float64"),
                "objective_value": pd.Series(objective_values, dtype="float64"),
            }
        )
        result_df = self._efficiency_columns(result_df)
        return result_df

    def _format_lamda(self, lamda_values):
        rows = []
        for actual_index, values in lamda_values.items():
            if values is None:
                rows.append(pd.Series(np.nan, index=self.reference_data_index, name=actual_index))
            else:
                rows.append(pd.Series(values, name=actual_index))
        if not rows:
            return pd.DataFrame(index=self.evaluated_data_index, columns=self.reference_data_index)
        return pd.concat(rows, axis=1).T

    def display_status(self):
        """Display solver status for each DMU."""
        print(pd.Series(self.get_status()))

    def display_rho(self):
        """Display rho values for each DMU."""
        print(self.get_rho())

    def display_lamda(self):
        """Display intensity variables for each DMU."""
        print(self.get_lamda())

    def get_status(self):
        """Return solver status by original evaluated-data index."""
        return self.results.get("optimization_status", {})

    def get_rho(self):
        """Return rho as a pandas Series indexed by original DMU labels."""
        if not self.results:
            raise RuntimeError("The model has not been optimized yet.")
        return pd.Series(self.results.get(self.result_value_name, {}), dtype="float64")

    def get_lamda(self):
        """Return intensity variables as a DataFrame."""
        if not self.results:
            raise RuntimeError("The model has not been optimized yet.")
        return self.results.get("lamda_df")

    def info(self, dmu="all"):
        """Print Pyomo model details for selected DMUs."""
        if dmu == "all":
            dmu_list = list(self._models)
        elif isinstance(dmu, (str, int, float)):
            dmu_list = [dmu]
        else:
            dmu_list = list(dmu)
        for actual_index in dmu_list:
            if actual_index not in self._models:
                print(f"DMU '{actual_index}' not found in the evaluated set.")
                continue
            print(f"\n--- Model for DMU: {actual_index} ---")
            self._models[actual_index].pprint()


class DEA2(_PerDmuEnvelopmentModel):
    """Classic radial DEA model solved separately for each evaluated DMU."""

    def __init__(self, data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        super().__init__(data, sent, gy=gy, gx=gx, rts=rts, baseindex=baseindex, refindex=refindex)

    def _add_decision_variables(self, model, actual_index):
        if self.input_oriented:
            bounds = (0.0, 1.0)
        elif self.output_oriented:
            bounds = (1.0, None)
        elif self.hyper_oriented and self.rts == RTS_CRS:
            bounds = (0.0, 1.0)
        else:
            bounds = (0.0, None)
        model.rho = Var(bounds=bounds, doc=f"efficiency for DMU {actual_index}")
        model.lamda = Var(model.R, bounds=(0.0, None), doc="intensity variables")

    def _add_objective(self, model):
        if self.output_oriented or (self.hyper_oriented and self.rts == RTS_VRS1):
            model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc="objective function")
        else:
            model.objective = Objective(rule=lambda m: m.rho, sense=minimize, doc="objective function")

    def _create_input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            lhs = sum(model.lamda[r] * self.xref[r][j] for r in model.R)
            x0 = self.x[current_dmu_range_index][j]
            if self.input_oriented:
                return lhs <= (model.rho * x0 if self.gx[j] > 0 else x0)
            if self.output_oriented:
                return lhs <= x0
            # hyper_oriented: CRS vs VRS
            if self.rts == RTS_CRS:
                return lhs <= (model.rho * x0 if self.gx[j] > 0 else x0)
            return lhs <= x0 - model.rho * self.gx[j] * x0
        return input_rule

    def _create_output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            lhs = sum(model.lamda[r] * self.yref[r][k] for r in model.R)
            y0 = self.y[current_dmu_range_index][k]
            if self.input_oriented:
                return lhs >= y0
            if self.output_oriented:
                return lhs >= (model.rho * y0 if self.gy[k] > 0 else y0)
            # hyper_oriented: CRS vs VRS
            if self.rts == RTS_CRS:
                return lhs >= y0
            return lhs >= y0 + model.rho * self.gy[k] * y0
        return output_rule

    def _efficiency_columns(self, result_df):
        rho = result_df["rho"]
        if self.input_oriented:
            result_df["te"] = rho
        elif self.output_oriented:
            result_df["te"] = rho.apply(lambda item: 1 / item if pd.notna(item) and item != 0 else np.nan)
        elif self.hyper_oriented and self.rts == RTS_CRS:
            # CRS hyperbolic linearization (3.15): rho = eta'^*, te = sqrt(rho) = eta*
            # Valid only when ALL inputs are on the hyperbolic path (gx all > 0).
            # For partial-input CRS hyper, the LP degenerates to input-oriented
            # radial (rho = theta*), and sqrt(rho) has no hyperbolic meaning.
            # Fall back to te = rho (same as input-oriented) for partial-input case.
            if all(g > 0 for g in self.gx):
                result_df["te"] = rho.apply(lambda item: np.sqrt(item) if pd.notna(item) and item >= 0 else np.nan)
            else:
                result_df["te"] = rho
        elif self.hyper_oriented:
            result_df["tei"] = rho.apply(lambda item: 1 - item if pd.notna(item) else np.nan)
            result_df["teo"] = rho.apply(lambda item: 1 / (1 + item) if pd.notna(item) else np.nan)
        return result_df


class DDF2(_PerDmuEnvelopmentModel):
    """Directional distance function solved separately for each evaluated DMU."""

    def __init__(self, data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        super().__init__(data, sent, gy=gy, gx=gx, rts=rts, baseindex=baseindex, refindex=refindex)

    def _add_decision_variables(self, model, actual_index):
        model.rho = Var(bounds=(0.0, None), doc=f"directional distance for DMU {actual_index}")
        model.lamda = Var(model.R, bounds=(0.0, None), doc="intensity variables")

    def _add_objective(self, model):
        model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc="objective function")

    def _create_input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= (
                self.x[current_dmu_range_index][j]
                - model.rho * self.gx[j] * self.x[current_dmu_range_index][j]
            )
        return input_rule

    def _create_output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= (
                self.y[current_dmu_range_index][k]
                + model.rho * self.gy[k] * self.y[current_dmu_range_index][k]
            )
        return output_rule

    def _efficiency_columns(self, result_df):
        rho = result_df["rho"]
        if self.direction.has_input_direction:
            result_df["tei"] = rho.apply(lambda item: 1 - item if pd.notna(item) else np.nan)
        if self.direction.has_output_direction:
            result_df["teo"] = rho.apply(lambda item: 1 / (1 + item) if pd.notna(item) else np.nan)
        return result_df



class DDF(DEA):
    def __init__(self,  data, sent, gy=[1], gx=[1], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Directional distance function

        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        self.y, self.x, self.yref, self.xref, self.gy, self.gx, _data_index = tools.assert_DDF(
            data, sent, gy, gx, baseindex, refindex
        )
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()
        self.__model__.R = Set(initialize=range(len(self.yref)))

        # Initialize sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize variable
        self.__model__.rho = Var(
            self.__model__.I, doc='directional distance')

        self.__model__.lamda = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None), doc='intensity variables')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self._DEA__objective_rule(), sense=maximize, doc='objective function')
        # self.__model__.objective.pprint()
        self.__model__.input = Constraint(
            self.__model__.I, self.__model__.J, rule=self.__input_rule(), doc='input constraint')
        # self.__model__.input.pprint()


        self.__model__.output = Constraint(
            self.__model__.I, self.__model__.K, rule=self.__output_rule(), doc='output constraint')

        # self.__model__.output.pprint()

        if self.rts == RTS_VRS1:
            self.__model__.vrs = Constraint(
                self.__model__.I, rule=self.__vrs_rule(), doc='various return to scale rule')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __input_rule(self):
        """Return the proper input constraint"""
        def input_rule(model, o, j):
            return sum(model.lamda[o, r] * self.xref[r][j] for r in model.R) <= \
                    self.x[o][j] - model.rho[o]*self.gx[j]*self.x[o][j]

        return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, o, k):
            return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= \
                        self.y[o][k] + model.rho[o]*self.gy[k]*self.y[o][k]

        return output_rule



    def __vrs_rule(self):
        """Return the VRS constraint"""
        def vrs_rule(model, o):
            return sum(model.lamda[o, r] for r in model.R) == 1
        return vrs_rule




class DEADUAL(DEA):
    def __init__(self, data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Multiplier problem

        Args:
            data
            sent
            orient (String): ORIENT_IO (input orientation) or ORIENT_OO (output orientation)
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        self.y, self.x, self.yref, self.xref, self.gy, self.gx, _data_index = tools.assert_DEA(
            data, sent, gy, gx, baseindex, refindex
        )
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()
        # Initialize sets
        self.__model__.R = Set(initialize=range(len(self.yref)))
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize variable
        self.__model__.delta = Var(self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='multiplier x')
        self.__model__.gamma = Var(self.__model__.I, self.__model__.K, bounds=(0.0, None), doc='multiplier y')
        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I, doc='variable return to scale')

        # Setup the objective function and constraints
        if sum(self.gx)>=1:
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
        elif sum(self.gy)>=1:
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
        else:
            raise ValueError("gx and gy must either be 1 or 0.")
        # self.__model__.objective.pprint()
        self.__model__.first = Constraint(
            self.__model__.I, self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
        # self.__model__.first.pprint()

        self.__model__.second = Constraint(
            self.__model__.I,                   rule=self.__second_rule(), doc='normalization constraint')
        # self.__model__.second.pprint()
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        if sum(self.gx) >= 1:

            def objective_rule(model):
                if self.rts == RTS_VRS1:
                    return sum(sum(model.delta[o, j] * self.x[o][j] * (1-self.gx[j]) for o in model.I) for j in model.J) - \
                        sum(sum(model.gamma[o, k] * self.y[o][k] * (1-self.gy[k]) for o in model.I) for k in model.K) +\
                        sum(model.alpha[o] for o in model.I)
                elif self.rts == RTS_CRS:
                    return sum(sum(model.delta[o, j] * self.x[o][j] * (1-self.gx[j]) for o in model.I) for j in model.J) - \
                        sum(sum(model.gamma[o, k] * self.y[o][k] * (1-self.gy[k]) for o in model.I) for k in model.K)
            return objective_rule

        elif sum(self.gy) >= 1:

            def objective_rule(model):
                if self.rts == RTS_VRS1:
                    return sum(sum(model.delta[o, j] * self.x[o][j] * (1-self.gx[j]) for o in model.I) for j in model.J) - \
                        sum(sum(model.gamma[o, k] * self.y[o][k] * (1-self.gy[k]) for o in model.I) for k in model.K) +\
                        sum(model.alpha[o] for o in model.I)
                elif self.rts == RTS_CRS:
                    return sum(sum(model.delta[o, j] * self.x[o][j] * (1-self.gx[j]) for o in model.I) for j in model.J) - \
                        sum(sum(model.gamma[o, k] * self.y[o][k] * (1-self.gy[k]) for o in model.I) for k in model.K)

            return objective_rule

    def __first_rule(self):
        """Return the proper technology constraint"""
        if sum(self.gx) >= 1:
            if self.rts == RTS_VRS1:
                def first_rule(model, o, r):
                    return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) \
                        - sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + model.alpha[o] >= 0
                return first_rule
            elif self.rts == RTS_CRS:
                def first_rule(model, o, r):
                    return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) \
                        - sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) >= 0
                return first_rule
        elif sum(self.gy) >= 1:
            if self.rts == RTS_VRS1:
                def first_rule(model, o, r):
                    return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) - \
                        sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + model.alpha[o] >= 0
                return first_rule
            elif self.rts == RTS_CRS:
                def first_rule(model, o, r):
                    return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) - \
                        sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) >= 0
                return first_rule

    def __second_rule(self):
        """Return the proper normalization constraint"""
        if sum(self.gx) >= 1:
            def second_rule(model, o):
                return sum(model.delta[o, j] * self.x[o][j]*self.gx[j] for j in model.J) <= 1
            return second_rule
        elif sum(self.gy) >= 1:

            def second_rule(model, o):
                return sum(model.gamma[o, k] * self.y[o][k]*self.gy[k] for k in model.K) == 1
            return second_rule

    def display_delta(self):
        """Display delta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.delta.display()

    def display_gamma(self):
        """Display gamma value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.gamma.display()

    def display_alpha(self):
         """Display omega value"""
         tools.assert_optimized(self.optimization_status)
         tools.assert_various_return_to_scale_alpha(self.rts)
         self.__model__.alpha.display()

    def get_delta(self):
        """Return delta value by array"""
        tools.assert_optimized(self.optimization_status)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                          list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy()

    def get_gamma(self):
        """Return nu value by array"""
        tools.assert_optimized(self.optimization_status)
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma),
                                                          list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value'])
        gamma = gamma.pivot(index='Name', columns='Key', values='Value')
        return gamma.to_numpy()

    def get_alpha(self):
        """Return omega value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale_alpha(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_efficiency(self):
        """Return efficiency value by array"""
        tools.assert_optimized(self.optimization_status)
        if sum(self.gx) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1)
            elif self.rts == RTS_VRS1:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1) + self.get_alpha().reshape(len(self.y), 1)
        elif sum(self.gy) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1)
            elif self.rts == RTS_VRS1:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1) + self.get_alpha().reshape(len(self.x), 1)



class DDFDUAL(DEADUAL):

    def __init__(self,  data, sent, gy=[1], gx=[1], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Directional distance function

        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        self.y, self.x, self.yref, self.xref, self.gy, self.gx, _data_index = tools.assert_DDF(
            data, sent, gy, gx, baseindex, refindex
        )
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()
        self.__model__.R = Set(initialize=range(len(self.yref)))

        # Initialize sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize variable
        self.__model__.delta = Var(self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='multiplier x')
        self.__model__.gamma = Var(self.__model__.I, self.__model__.K, bounds=(0.0, None), doc='multiplier y')

        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I, doc='variable return to scale')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self.__objective_rule(), sense=minimize, doc='objective function')
        self.__model__.first = Constraint(
            self.__model__.I, self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
        self.__model__.second = Constraint(
            self.__model__.I, rule=self.__second_rule(), doc='normalization constraint')


        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            if self.rts == RTS_VRS1:
                return sum(sum(model.delta[o, j] * self.x[o][j] for o in model.I) for j in model.J) - \
                    sum(sum(model.gamma[o, k] * self.y[o][k] for o in model.I) for k in model.K) + \
                    sum(model.alpha[o] for o in model.I)
            elif self.rts == RTS_CRS:
                return sum(sum(model.delta[o, j] * self.x[o][j] for o in model.I) for j in model.J) - \
                    sum(sum(model.gamma[o, k] * self.y[o][k] for o in model.I) for k in model.K)
        return objective_rule


    def __first_rule(self):
        """Return the proper technology constraint"""
        if self.rts == RTS_VRS1:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) - \
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + model.alpha[o] >= 0
            return first_rule
        elif self.rts == RTS_CRS:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) - \
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) >= 0
            return first_rule


    def __second_rule(self):
        """Return the proper normalization constraint"""
        def second_rule(model, o):
            return sum(model.delta[o, j] *self.gx[j]*self.x[o][j] for j in model.J) + \
                sum(model.gamma[o, k] *self.gy[k]* self.y[o][k] for k in model.K) == 1
        return second_rule


    def display_delta(self):
        """Display delta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.delta.display()

    def display_gamma(self):
        """Display gamma value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.gamma.display()

    def display_alpha(self):
         """Display omega value"""
         tools.assert_optimized(self.optimization_status)
         tools.assert_various_return_to_scale_alpha(self.rts)
         self.__model__.alpha.display()

    def get_delta(self):
        """Return delta value by array"""
        tools.assert_optimized(self.optimization_status)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                          list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy()

    def get_gamma(self):
        """Return nu value by array"""
        tools.assert_optimized(self.optimization_status)
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma),
                                                          list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value'])
        gamma = gamma.pivot(index='Name', columns='Key', values='Value')
        return gamma.to_numpy()

    def get_alpha(self):
        """Return omega value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale_alpha(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_efficiency(self):
        """Return efficiency value by array"""
        tools.assert_optimized(self.optimization_status)
        if sum(self.gx) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1)
            elif self.rts == RTS_VRS1:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1) + self.get_alpha().reshape(len(self.y), 1)
        elif sum(self.gy) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1)
            elif self.rts == RTS_VRS1:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1) + self.get_alpha().reshape(len(self.x), 1)

