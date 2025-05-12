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






class DEA2:
    """Data Envelopment Analysis (DEA) - Solves per DMU"""

    def __init__(self, data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Envelopment problem, solving for each DMU individually.

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEA should return numpy arrays for x, y, xref, yref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DEA is modified or a similar function returns these:
        try:
             self.y, self.x, self.yref, self.xref, self.gy, self.gx, \
                 self.evaluated_data_index, self.reference_data_index = tools.assert_DEA_with_indices(
                     data, sent, gy, gx, baseindex, refindex
                 )
        except AttributeError:
             # Fallback if assert_DEA_with_indices doesn't exist, assuming original assert_DEA
             # In this case, evaluated_data_index and reference_data_index will be ranges 0..N-1
             print("Warning: tools.assert_DEA_with_indices not found. Using range indices.")
             self.y, self.x, self.yref, self.xref, self.gy, self.gx, _ = tools.assert_DEA(
                 data, sent, gy, gx, baseindex, refindex
             )
             self.evaluated_data_index = list(range(len(self.y)))
             self.reference_data_index = list(range(len(self.yref)))


        self.rts = rts

        self.evaluated_indices_range = range(len(self.y)) # Range indices for internal numpy arrays
        self.reference_indices_range = range(len(self.yref)) # Range indices for internal numpy arrays
        self.num_inputs = len(self.x[0])
        self.num_outputs = len(self.y[0])
        # print(self.gx,self.gy,")))))))))))))")
        # Determine orientation based on gx/gy vectors
        self.input_oriented = sum(self.gx) >= 1 and sum(self.gy) == 0
        self.output_oriented = sum(self.gy) >= 1 and sum(self.gx) == 0
        self.hyper_oriented = sum(self.gx) >= 1 and sum(self.gy) >= 1

        if not (self.input_oriented or self.output_oriented or self.hyper_oriented):
             raise ValueError("gx and gy must represent either input or output or hyper orientation.")
        if self.input_oriented and self.output_oriented:
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

            # Define variables for this specific model
            # rho is the efficiency score for the current DMU (not indexed by DMU within this model)
            # Bounds based on standard DEA orientation (matching DEAt)
            if self.input_oriented: # Minimize rho, efficiency <= 1
                model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')
            elif self.output_oriented: # Maximize rho, efficiency >= 1
                model.rho = Var(bounds=(1, None), doc=f'efficiency for DMU {actual_index}')
            elif self.hyper_oriented: # Hyper orientation, rho can be any value, 但正常<1
                if self.rts == RTS_VRS1:
                    model.rho = Var(bounds=(0, None), doc=f'efficiency for DMU {actual_index}')  
                elif self.rts == RTS_CRS:
                    # Hyper orientation with CRS, rho efficiency <= 1
                    model.rho = Var(bounds=(0, 1), doc=f'efficiency for DMU {actual_index}')
            else:
                 # This else should not be reached due to check above, but for safety
                 model.rho = Var(bounds=(0, None), doc=f'efficiency for DMU {actual_index}')


            # lamda are intensity variables, indexed by the reference set
            model.lamda = Var(model.R, bounds=(0.0, None), doc='intensity variables')

            # Setup the objective function and constraints for THIS DMU
            # Objective: Minimize/Maximize the SINGLE rho variable in this model
            if self.input_oriented:
                model.objective = Objective(rule=lambda m: m.rho, sense=minimize, doc='objective function')
            elif self.output_oriented:
                model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc='objective function')
            elif self.hyper_oriented:
                if self.rts == RTS_VRS1:
                    model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc='objective function')
                elif self.rts == RTS_CRS:
                    # Hyper orientation with CRS, minimize rho
                    model.objective = Objective(rule=lambda m: m.rho, sense=minimize, doc='objective function')
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

            # VRS constraint for THIS DMU
            if self.rts == RTS_VRS1:
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
            def input_rule(model, j):
                # Access current DMU's input data self.x[current_dmu_range_index][j]
                # Access reference DMUs' input data self.xref[r][j]
                if self.gx[j] == 1:
                    # Input is scaled by rho
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= model.rho * self.x[current_dmu_range_index][j]
                else:
                    # Input is not scaled by rho
                    return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
        elif self.output_oriented:
            def input_rule(model, j):
                # Input is not scaled by rho in output orientation
                return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
        elif self.hyper_oriented:
            if self.rts == RTS_CRS:

                def input_rule(model, j):
                    # Hyper orientation: Input is scaled by rho if gx[j] == 1
                    if self.gx[j] == 1:
                        return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= model.rho * self.x[current_dmu_range_index][j]
                    else:
                        return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
            elif self.rts == RTS_VRS1:
                def input_rule(model, j):
                    return sum(model.lamda[ r] * self.xref[r][j] for r in model.R) <= \
                            self.x[current_dmu_range_index][j] - model.rho*self.gx[j]*self.x[current_dmu_range_index][j]
        else:
             # Should not be reached
             def input_rule(model, j):
                 return Constraint.Skip # Or raise error
        return input_rule

    def __create_output_rule(self, current_dmu_range_index):
        """Factory for creating the output constraint rule for a specific DMU."""
        if self.input_oriented:
            def output_rule(model, k):
                # Output is not scaled by rho in input orientation
                return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
        elif self.output_oriented:
            def output_rule(model, k):
                # Output is scaled by rho in output orientation (original code applies rho to all outputs if sum(gy)>=1)
                # Note: The original code's output rule didn't check gy[k]==1 like the input rule checked gx[j]==1.
                # Replicating original logic:
                return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= model.rho * self.y[current_dmu_range_index][k]
        elif self.hyper_oriented:
            if self.rts == RTS_CRS:

                def output_rule(model, k):
                    # Hyper orientation: Output is not scaled by rho if gy[k] == 1
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
            elif self.rts == RTS_VRS1:
                def output_rule(model, k):
                    # Output is not scaled by rho in input orientation
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                        self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
        else:
             # Should not be reached
             def output_rule(model, k):
                 return Constraint.Skip # Or raise error
        return output_rule


    def __create_vrs_rule(self):
        """Factory for creating the VRS constraint rule."""
        def vrs_rule(model):
            return sum(model.lamda[r] for r in model.R) == 1
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
        self.results = {} # Clear previous results
        all_problem_statuses = {}
        all_optimization_statuses = {}
        rho_values = {}
        lamda_values = {} # Store lamda results temporarily as dicts

        # Loop through each DMU's model and solve it
        # The dictionary keys are the actual data indices

        use_neos = tools.set_neos_email(email)

        for actual_index, model in self.__modeldict.items():
            # print(f"Optimizing for DMU: {actual_index}...") # Optional: print progress
            try:
                # Assuming tools.optimize_model2 exists and works like DEAt
                # It should return (problem_status, optimization_status)

                optimization_status = tools.optimize_model2(
                   model, actual_index, use_neos, CET_ADDI, solver=solver)

            except Exception as e:
                print(f"sasaa Error optimizing DMU {actual_index}: {e}")
                optimization_status = "Error"

            all_optimization_statuses[actual_index] = optimization_status

            # Store results if optimization was successful (check standard Pyomo status strings)
            # Common successful statuses include 'optimal', 'feasible'
            if optimization_status in ['ok']:
                 try:
                     rho_values[actual_index] = model.rho.value
                     # Collect lamda values (indexed by reference range index in the model)
                     lamda_data_for_dmu = {}
                     for r_range in model.R:
                         # Map the reference range index back to the actual reference index label
                         actual_ref_index = self.reference_data_index[r_range]
                         lamda_data_for_dmu[actual_ref_index] = model.lamda[r_range].value
                     lamda_values[actual_index] = lamda_data_for_dmu # Store as dictionary for easier conversion
                 except Exception as e:
                      print(f"Warning: Could not retrieve results for DMU {actual_index} despite status '{optimization_status}': {e}")
                      rho_values[actual_index] = np.nan # Indicate failure to retrieve value
                      lamda_values[actual_index] = None
            else:
                 print(f"Error optimizing DMU {actual_index}, optimization_status :{optimization_status}")

                 rho_values[actual_index] = np.nan # Use NaN for failed optimizations
                 lamda_values[actual_index] = None
            



        # Store collected results in self.results
        self.results['optimization_status'] = all_optimization_statuses
        self.results['rho'] = rho_values

        # Process lamda values into a DataFrame
        lamda_df_list = []
        # Keys of lamda_values are actual evaluated DMU indices
        for actual_index, lam_data_for_dmu in lamda_values.items():
            if lam_data_for_dmu is not None:
                 # Create a Series for this DMU's lamda values
                 # The keys of lam_data_for_dmu are actual reference DMU indices
                 lam_series = pd.Series(lam_data_for_dmu, name=actual_index)
                 lamda_df_list.append(lam_series)
            else:
                 # For failed DMUs, add a row of NaNs with correct reference index columns
                 nan_series = pd.Series(np.nan, index=self.reference_data_index, name=actual_index)
                 lamda_df_list.append(nan_series)


        if lamda_df_list:
             # Concatenate the Series into a DataFrame. Transpose to get evaluated DMUs as index.
             # Columns will be the actual reference DMU indices.
             self.results['lamda_df'] = pd.concat(lamda_df_list, axis=1).T
        else:
             self.results['lamda_df'] = None # No successful optimizations or no DMUs evaluated

        # Create a summary DataFrame to return (similar to DEAt)
        results_df = pd.DataFrame({
            'optimization_status': pd.Series(all_optimization_statuses),
            'rho': pd.Series(rho_values),
            # Could add problem_status here too if needed
        })

        # Add 'te' column if applicable, based on rho interpretation
        # Standard DEA rho is efficiency, often 0-1 or >=1.
        # If rho is the efficiency score directly, maybe no 'te' calculation needed?
        # DEAt calculates 'te' from 'beta'. Let's assume rho *is* the efficiency score desired.
        # If input-oriented (minimize rho, 0-1), rho is efficiency.
        # If output-oriented (maximize rho, >=1), 1/rho is efficiency.
        # Let's add a 'te' column that aligns with 0-1 efficiency measure.
        # If input_oriented, te = rho
        # If output_oriented, te = 1/rho
        if self.input_oriented:
            results_df['te'] = results_df['rho']
        elif self.output_oriented:
             # Avoid division by zero or NaN
             results_df['te'] = results_df['rho'].apply(lambda x: 1/x if pd.notna(x) and x != 0 else np.nan)
        elif self.hyper_oriented:
             if self.rts == RTS_CRS:
                 results_df['te'] = results_df['rho'].apply(lambda x: np.sqrt(x) if pd.notna(x) else np.nan)
             elif self.rts == RTS_VRS1:
                 # Avoid  NaN    
                results_df['tei'] = results_df['rho'].apply(lambda x: (1-x) if pd.notna(x) else np.nan)
                results_df['teo'] = results_df['rho'].apply(lambda x: 1/(1+x) if pd.notna(x) else np.nan)

        return results_df


    # --- Display and Get Methods ---

    def display_status(self):
        """Display the optimization status for each DMU."""
        if not self.results:
            print("Optimization has not been run yet.")
            return
        print("Optimization Status per DMU:")
        for dmu, status in self.results.get('optimization_status', {}).items():
            print(f"  {dmu}: {status}")

    def display_rho(self):
        """Display rho value for each DMU."""
        # Use assert_optimized to check if results exist
        tools.assert_optimized(self.results)
        print("Rho values per DMU:")
        rho_series = pd.Series(self.results.get('rho', {}))
        print(rho_series)


    def display_lamda(self):
        """Display lamda values (intensity variables) for each DMU."""
        tools.assert_optimized(self.results)
        print("Lamda values per DMU:")
        lamda_df = self.results.get('lamda_df')
        if lamda_df is not None:
            print(lamda_df)
        else:
            print("No lamda values available (optimization may have failed for all DMUs).")


    def get_status(self):
        """Return optimization status dictionary."""
        if not self.results:
             return {}
        return self.results.get('optimization_status', {})

    def get_rho(self):
        """Return rho values as a pandas Series."""
        tools.assert_optimized(self.results)
        return pd.Series(self.results.get('rho', {}))

    def get_lamda(self):
        """Return lamda values as a pandas DataFrame."""
        tools.assert_optimized(self.results)
        return self.results.get('lamda_df')

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
            for ind, problem in self.__modeldict.items():
                print(f"\n--- Model for DMU: {ind} ---")
                problem.pprint()
                print("-" * (len(f"--- Model for DMU: {ind} ---")))
        else:
            if isinstance(dmu, str):
                dmu_list = [dmu]
            else:
                dmu_list = dmu

            for ind in dmu_list:
                if ind in self.__modeldict:
                    print(f"\n--- Model for DMU: {ind} ---")
                    self.__modeldict[ind].pprint()
                    print("-" * (len(f"--- Model for DMU: {ind} ---")))
                else:
                    print(f"DMU '{ind}' not found in the evaluated set.")









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
        self.y, self.x, self.yref, self.xref, self.gy, self.gx = tools.assert_DDF(data, sent, gy, gx,baseindex,refindex)
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


class DDF2(DEA2):
    """Data Envelopment Analysis (DEA) - Solves per DMU"""

    def __init__(self, data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Envelopment problem, solving for each DMU individually.

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # assert_DEA should return numpy arrays for x, y, xref, yref, and the actual indices for evaluated and reference DMUs
        # Let's assume tools.assert_DDF is modified or a similar function returns these:
        try:
             self.y, self.x, self.yref, self.xref, self.gy, self.gx, \
                 self.evaluated_data_index, self.reference_data_index = tools.assert_DDF_with_indices(
                     data, sent, gy, gx, baseindex, refindex
                 )
             
        except AttributeError:
             # Fallback if assert_DEA_with_indices doesn't exist, assuming original assert_DEA
             # In this case, evaluated_data_index and reference_data_index will be ranges 0..N-1
             print("Warning: tools.assert_DEA_with_indices not found. Using range indices.")
             self.y, self.x, self.yref, self.xref, self.gy, self.gx, _ = tools.assert_DDF(
                 data, sent, gy, gx, baseindex, refindex
             )
             self.evaluated_data_index = list(range(len(self.y)))
             self.reference_data_index = list(range(len(self.yref)))


        self.rts = rts

        self.evaluated_indices_range = range(len(self.y)) # Range indices for internal numpy arrays
        self.reference_indices_range = range(len(self.yref)) # Range indices for internal numpy arrays
        self.num_inputs = len(self.x[0])
        self.num_outputs = len(self.y[0])
        # print(self.gx,self.gy,"############")



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

            # Define variables for this specific model
            # rho is the efficiency score for the current DMU (not indexed by DMU within this model)
            model.rho = Var(bounds=(0, None), doc=f'beta for DMU {actual_index}')

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

            # VRS constraint for THIS DMU
            if self.rts == RTS_VRS1:
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
        def input_rule(model, j):
            # Access current DMU's input data self.x[current_dmu_range_index][j]
            # Access reference DMUs' input data self.xref[r][j]
            return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= \
                self.x[current_dmu_range_index][j]- model.rho*self.gx[j]*self.x[current_dmu_range_index][j]
        return input_rule

    def __create_output_rule(self, current_dmu_range_index):
        """Factory for creating the output constraint rule for a specific DMU."""
        def output_rule(model, k):
            # Output is not scaled by rho in input orientation
            return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= \
                self.y[current_dmu_range_index][k] + model.rho*self.gy[k]*self.y[current_dmu_range_index][k]
        return output_rule


    def __create_vrs_rule(self):
        """Factory for creating the VRS constraint rule."""
        def vrs_rule(model):
            return sum(model.lamda[r] for r in model.R) == 1
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
        self.results = {} # Clear previous results
        all_problem_statuses = {}
        all_optimization_statuses = {}
        rho_values = {}
        lamda_values = {} # Store lamda results temporarily as dicts

        # Loop through each DMU's model and solve it
        # The dictionary keys are the actual data indices

        use_neos = tools.set_neos_email(email)

        for actual_index, model in self.__modeldict.items():
            # print(f"Optimizing for DMU: {actual_index}...") # Optional: print progress
            try:
                # Assuming tools.optimize_model2 exists and works like DEAt
                # It should return (problem_status, optimization_status)

                optimization_status = tools.optimize_model2(
                   model, actual_index, use_neos, CET_ADDI, solver=solver)

            except Exception as e:
                print(f"Error optimizing DMU {actual_index}: {e}")
                optimization_status = "Error"

            all_optimization_statuses[actual_index] = optimization_status

            # Store results if optimization was successful (check standard Pyomo status strings)
            # Common successful statuses include 'optimal', 'feasible'
            if optimization_status in ['ok']:
                 try:
                     rho_values[actual_index] = model.rho.value
                     # Collect lamda values (indexed by reference range index in the model)
                     lamda_data_for_dmu = {}
                     for r_range in model.R:
                         # Map the reference range index back to the actual reference index label
                         actual_ref_index = self.reference_data_index[r_range]
                         lamda_data_for_dmu[actual_ref_index] = model.lamda[r_range].value
                     lamda_values[actual_index] = lamda_data_for_dmu # Store as dictionary for easier conversion
                 except Exception as e:
                      print(f"Warning: Could not retrieve results for DMU {actual_index} despite status '{optimization_status}': {e}")
                      rho_values[actual_index] = np.nan # Indicate failure to retrieve value
                      lamda_values[actual_index] = None
            else:
                 print(f"Error optimizing DMU {actual_index}, optimization_status :{optimization_status}")

                 rho_values[actual_index] = np.nan # Use NaN for failed optimizations
                 lamda_values[actual_index] = None


        # Store collected results in self.results
        self.results['optimization_status'] = all_optimization_statuses
        self.results['rho'] = rho_values

        # Process lamda values into a DataFrame
        lamda_df_list = []
        # Keys of lamda_values are actual evaluated DMU indices
        for actual_index, lam_data_for_dmu in lamda_values.items():
            if lam_data_for_dmu is not None:
                 # Create a Series for this DMU's lamda values
                 # The keys of lam_data_for_dmu are actual reference DMU indices
                 lam_series = pd.Series(lam_data_for_dmu, name=actual_index)
                 lamda_df_list.append(lam_series)
            else:
                 # For failed DMUs, add a row of NaNs with correct reference index columns
                 nan_series = pd.Series(np.nan, index=self.reference_data_index, name=actual_index)
                 lamda_df_list.append(nan_series)


        if lamda_df_list:
             # Concatenate the Series into a DataFrame. Transpose to get evaluated DMUs as index.
             # Columns will be the actual reference DMU indices.
             self.results['lamda_df'] = pd.concat(lamda_df_list, axis=1).T
        else:
             self.results['lamda_df'] = None # No successful optimizations or no DMUs evaluated

        # Create a summary DataFrame to return (similar to DEAt)
        results_df = pd.DataFrame({
            'optimization_status': pd.Series(all_optimization_statuses),
            'rho': pd.Series(rho_values),
            # Could add problem_status here too if needed
        })

        # Add 'te' column if applicable, based on rho interpretation
        # Standard DEA rho is efficiency, often 0-1 or >=1.
        # If rho is the efficiency score directly, maybe no 'te' calculation needed?
        # DEAt calculates 'te' from 'beta'. Let's assume rho *is* the efficiency score desired.
        # If input-oriented (minimize rho, 0-1), rho is efficiency.
        # If output-oriented (maximize rho, >=1), 1/rho is efficiency.
        # Let's add a 'te' column that aligns with 0-1 efficiency measure.

        if sum(self.gx) >= 1 and sum(self.gy) == 0:
            results_df['tei'] = results_df['rho'].apply(lambda x: 1-x if pd.notna(x) else np.nan)
        elif sum(self.gy) >= 1 and sum(self.gx) == 0:
             # Avoid division by zero or NaN
             results_df['teo'] = results_df['rho'].apply(lambda x: 1/(1+x) if pd.notna(x) else np.nan)
        elif sum(self.gx) >= 1 and sum(self.gy) >= 1:
                results_df['tei'] = results_df['rho'].apply(lambda x: 1-x if pd.notna(x) else np.nan)
                results_df['teo'] = results_df['rho'].apply(lambda x: 1/(1+x) if pd.notna(x) else np.nan)

        return results_df
    
    # --- Display and Get Methods ---

    def display_status(self):
        """Display the optimization status for each DMU."""
        if not self.results:
            print("Optimization has not been run yet.")
            return
        print("Optimization Status per DMU:")
        for dmu, status in self.results.get('optimization_status', {}).items():
            print(f"  {dmu}: {status}")

    def display_rho(self):
        """Display rho value for each DMU."""
        # Use assert_optimized to check if results exist
        tools.assert_optimized(self.results)
        print("Rho values per DMU:")
        rho_series = pd.Series(self.results.get('rho', {}))
        print(rho_series)


    def display_lamda(self):
        """Display lamda values (intensity variables) for each DMU."""
        tools.assert_optimized(self.results)
        print("Lamda values per DMU:")
        lamda_df = self.results.get('lamda_df')
        if lamda_df is not None:
            print(lamda_df)
        else:
            print("No lamda values available (optimization may have failed for all DMUs).")


    def get_status(self):
        """Return optimization status dictionary."""
        if not self.results:
             return {}
        return self.results.get('optimization_status', {})

    def get_rho(self):
        """Return rho values as a pandas Series."""
        tools.assert_optimized(self.results)
        return pd.Series(self.results.get('rho', {}))

    def get_lamda(self):
        """Return lamda values as a pandas DataFrame."""
        tools.assert_optimized(self.results)
        return self.results.get('lamda_df')

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
            for ind, problem in self.__modeldict.items():
                print(f"\n--- Model for DMU: {ind} ---")
                problem.pprint()
                print("-" * (len(f"--- Model for DMU: {ind} ---")))
        else:
            if isinstance(dmu, str):
                dmu_list = [dmu]
            else:
                dmu_list = dmu

            for ind in dmu_list:
                if ind in self.__modeldict:
                    print(f"\n--- Model for DMU: {ind} ---")
                    self.__modeldict[ind].pprint()
                    print("-" * (len(f"--- Model for DMU: {ind} ---")))
                else:
                    print(f"DMU '{ind}' not found in the evaluated set.")







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

        self.y, self.x, self.yref, self.xref, self.gy, self.gx = tools.assert_DEA(data, sent, gy, gx,baseindex,refindex)
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

        self.y, self.x, self.yref, self.xref, self.gy, self.gx = tools.assert_DDF(data, sent, gy, gx, baseindex,
                                                                                  refindex)
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
