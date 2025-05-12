"""Main module."""
# import dependencies

import numpy as np
import pandas as pd
from .constant import CET_ADDI, RTS_VRS1,RTS_VRS2, RTS_CRS, OPT_DEFAULT, OPT_LOCAL,TOTAL,CONTEMPORARY,LUE,MAL
from .utils import tools
from .DEAweak import DEAweak2,DDFweak2,NDDFweak2




class MQDEAweak:
    """Malmquist production index (MQPI)
    """

    def __init__(self, data,id,year,sent = "inputvar=outputvar:unoutputvar",  gy=[1], gx=[0], gb=[0], rts=RTS_VRS1, \
                 tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """MQDEAt: Calculates Malmquist index using DEA2 for underlying efficiency scores.

        Args:
            data (pandas.DataFrame): input pandas.
            id (str): column name to specify id.
            year (str): column name to specify time.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y:CO2"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output distance vector. Defaults to [0].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Note: DEA2 uses RTS_VRS1.
            tech (str): TOTAL or CONTEMPORARY.
            solver (str): The solver to use (e.g., "mosek", "cbc").
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # Initialize MQDEAt model

        # Ensure year column exists and is sortable
        if year not in data.columns:
            raise ValueError(f"Year column '{year}' not found in data.")
        self.tlt = pd.Series(data[year]).drop_duplicates().sort_values()  # 生成时间的列表

        # Parse input/output variables

        self.gy, self.gx, self.gb, self.inputvars,self.outputvars, self.unoutputvars = tools.assert_MQDEAweak(
                        data, sent, gy, gx, gb
                    )

        self.xcol = list(self.inputvars)  # Ensure it's a list for indexing
        self.ycol = list(self.outputvars)  # Ensure it's a list for indexing
        self.bcol = list(self.unoutputvars)  # Ensure it's a list for indexing


        self.tech = tech
        self.rts = rts
        self.email = email
        self.solver = solver

        # Determine orientation based on gx/gy vectors
        self.input_oriented = sum(self.gx) >= 1 and sum(self.gy) == 0 and sum(self.gb) == 0
        self.output_oriented = sum(self.gy) >= 1 and sum(self.gx) == 0 and sum(self.gb) == 0
        self.undesirable_oriented = sum(self.gb) >= 1 and sum(self.gx) == 0 and sum(self.gy) == 0
        self.hyper_orientedyx = sum(self.gx) >= 1 and sum(self.gy) >= 1 and sum(self.gb) == 0
        self.hyper_orientedyb = sum(self.gb) >= 1 and sum(self.gy) >= 1 and sum(self.gx) == 0

        # Create a copy of the original data to add results columns
        self.datazz = data.copy()

        # --- Perform DEA calculations using DEA2 based on the chosen technology ---

        if self.tech == TOTAL:
            print("Calculating D11 (Total frontier) for all periods...")

            self.get_total(data,sent,id,year)

            print("TOTAL tech calculation finished.")

        elif self.tech == CONTEMPORARY:
            print("Calculating CONTEMPORARY tech components (D11, D12, D21)...")

            self.get_contemp(data,sent,id,year)

        else:
            raise ValueError(f"Unsupported technology type '{self.tech}'. Must be '{TOTAL}' or '{CONTEMPORARY}'.")



    def optimize(self):
        """Returns the calculated Malmquist index and components DataFrame."""
        # In this implementation, optimize() just returns the pre-calculated results
        # from the __init__ method.
        if not hasattr(self, 'datazz'):
             raise RuntimeError("Malmquist index calculation failed during initialization.")
        return self.datazz





    def get_total(self,data,sent,id,year):
        """Calculate the total efficiency scores for all years."""
        dataz11_list = []  # List to store D11 results (or components) for each year
        # For Total frontier, evaluate each DMU in each year against the frontier of ALL years
        # The baseindex selects the DMU(s) for evaluation in a specific year.
        # The refindex should select ALL DMUs in ALL years for the reference set.
        all_years_ref_index = f"{year}=[{','.join(map(str, self.tlt.tolist()))}]" # Reference set includes all years
        # Loop through each year in the time list

        # Determine which columns to expect and how to process based on orientation and RTS
        if self.input_oriented or self.output_oriented or self.undesirable_oriented \
            or (self.hyper_orientedyx and self.rts == RTS_CRS) or (self.hyper_orientedyb and self.rts == RTS_CRS):
            # Standard case: Expect 'te' and calculate a single 'D11'
            expected_cols = ['te']
            output_cols = ['D11']
            process_type = 'single_d11'

        elif (self.hyper_orientedyb and self.rts == RTS_VRS1) or (self.hyper_orientedyb and self.rts == RTS_VRS2) :
            # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
            expected_cols = ['teuo', 'teo']
            output_cols = ['D11_teuo', 'D11_teo'] # Renaming for clarity
            process_type = 'separate_teuo_teo'

        elif  (self.hyper_orientedyx and self.rts == RTS_VRS2)  :
            # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
            expected_cols = ['tei', 'teo']
            output_cols = ['D11_tei', 'D11_teo'] # Renaming for clarity
            process_type = 'separate_tei_teo'

        else:
            raise ValueError(f"Unsupported orientation/RTS combination: input={self.input_oriented}, output={self.output_oriented}, hyper={self.hyper_oriented}, rts={self.rts}")

        # --- Loop through years and perform DEA ---
        for tindex in self.tlt.index:
            current_year = self.tlt.iloc[tindex]
            print(f"  Evaluating year {current_year} against Total frontier...")

            # Call DEA2 instead of DEAt
            # Use the calculated gx and gy, and the mapped RTS
            model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                            rts=self.rts, baseindex=f"{year}=[{current_year}]",
                            refindex=all_years_ref_index) # Reference set is all years

            # model.optimize() should return a DataFrame with DMU index and result columns
            data11_results = model.optimize(self.email,self.solver)

            # --- Extract/Select the relevant efficiency column(s) ---
            if not all(col in data11_results.columns for col in expected_cols):
                # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                # Consider adding more specific error messages or handling based on model.optimize() status
                raise KeyError(f"DEA2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

            if process_type == 'single_d11':
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            elif process_type == 'separate_tei_teo':
                # Select both 'tei' and 'teo' and rename them
                # Assuming the order in expected_cols is ['tei', 'teo'] if process_type is 'separate_tei_teo'
                data11_component = data11_results[expected_cols].rename(columns={'tei': 'D11_tei', 'teo': 'D11_teo'})
            elif process_type == 'separate_teuo_teo':
                data11_component = data11_results[expected_cols].rename(columns={'teuo': 'D11_teuo', 'teo': 'D11_teo'})
            # Ensure the index matches the actual DMU index from DEA2 results
            # Assuming data11_results' index is the actual DMU identifier
            data11_component.index = data11_results.index

            dataz11_list.append(data11_component)

        # --- Concatenate results for all years ---
        # pd.concat handles the case where dataz11_list might be empty
        # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
        # If process_type is 'single_d11', it will have a 'D11' column
        dataz11 = pd.concat(dataz11_list)

        # --- Join results with the main datazz DataFrame ---
        # Assumes self.datazz is initialized before this method is called and has the correct index structure
        # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
        self.datazz = self.datazz.join(dataz11, how='left')

        # --- Calculate Malmquist components based on the processed D11 values ---
        # This calculation relies on the data being sorted correctly (e.g., by DMU index then by year)
        # for the shift(1) operation to compare consecutive years for the same DMU.
        # It's highly recommended that the input 'data' DataFrame is sorted this way
        # before being passed to MQDEAt.
        # Assumes 'id' variable holds the name of the DMU identifier column/index level used for grouping.

        if process_type == 'single_d11':
            # Calculate single MQPI based on D11 (likely EC_total)
            # Check if the D11 column exists and has at least some non-null data after join
            if self.datazz.empty or "D11" not in self.datazz.columns or self.datazz["D11"].isnull().all():
                print("Warning: D11 calculation resulted in no valid data. Cannot compute MQPI.")
                self.datazz["mqpi"] = np.nan
            else:
                self.datazz["D11"] = pd.to_numeric(self.datazz["D11"], errors='coerce')
                # Calculate previous period value for D11 within each DMU group
                # Assumes 'id' is the grouping key (column name or index level name)
                self.datazz['D11_prev'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                # Compute the ratio (Current D11 / Previous D11)
                # This is the Efficiency Change component relative to the Total frontier (EC_total)
                self.datazz["mqpi"] = self.datazz["D11"] / self.datazz["D11_prev"]
                self.datazz.drop(columns = ['D11_prev'], inplace=True) # Drop the intermediate column

        elif process_type == 'separate_tei_teo':
            # Calculate two separate components based on D11_tei and D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_tei', 'D11_teo']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_tei and EC_total_teo components
            result_cols = ['mqpi_tei', 'mqpi_teo']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')

                # Calculate previous period values for both D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev'], inplace=True)


        elif process_type == 'separate_teuo_teo':
            # Calculate two separate components based on D11_tei and D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_teuo', 'D11_teo']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_tei and EC_total_teo components
            result_cols = ['mqpi_teuo', 'mqpi_teo']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')

                # Calculate previous period values for both D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev'], inplace=True)

        # The function modifies self.datazz in place. It doesn't explicitly return datazz.
        # If the calling code expects a return value, uncomment the line below:
        # return self.datazz
        pass # Or return self.datazz if needed by the calling code









    def get_contemp(self,data,sent,id,year):

        # For Total frontier, evaluate each DMU in each year against the frontier of ALL years
        # The baseindex selects the DMU(s) for evaluation in a specific year.
        # The refindex should select ALL DMUs in ALL years for the reference set.
        all_years_ref_index = f"{year}=[{','.join(map(str, self.tlt.tolist()))}]" # Reference set includes all years
        # Loop through each year in the time list

        # Determine which columns to expect and how to process based on orientation and RTS
        if self.input_oriented or self.output_oriented or self.undesirable_oriented \
            or (self.hyper_orientedyx and self.rts == RTS_CRS) or (self.hyper_orientedyb and self.rts == RTS_CRS):
            # D11: Current period tech, current period frontier (Efficiency change component)
            dataz11_list = []  # List to store D11 results (or components) for each year
            expected_cols = ['te']
            output_cols = ['D11']
            # --- Loop through years and perform DEA ---
            for tindex in self.tlt.index:
                current_year = self.tlt.iloc[tindex]
                print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{current_year}]")
                data11_results = model.optimize(self.email,self.solver)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data11_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

    
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data11_results' index is the actual DMU identifier
                data11_component.index = data11_results.index

                dataz11_list.append(data11_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz11_list might be empty
            # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
            # If process_type is 'single_d11', it will have a 'D11' column
            dataz11 = pd.concat(dataz11_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz11, how='left')

            # D12: Current period tech, previous period frontier (Used in Tech change component)
            dataz12_list = []
            expected_cols = ['te']
            output_cols = ['D12']
            # Loop starts from the second year (index 1)
            for tindex in self.tlt.index[1:]:
                current_year = self.tlt.iloc[tindex]
                previous_year = self.tlt.iloc[tindex - 1]
                print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                data12_results = model.optimize(self.email,self.solver)
                # print(data12_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data12_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                # Select the single efficiency column and rename it
                data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            
                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data12_results' index is the actual DMU identifier
                data12_component.index = data12_results.index

                dataz12_list.append(data12_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz12_list might be empty
            dataz12 = pd.concat(dataz12_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz12, how='left')





            # D21: Previous period tech, current period frontier (Used in Tech change component)
            dataz21_list = []
            expected_cols = ['te']
            output_cols = ['D21']
            # Loop goes up to the second to last year (index -1)
            for tindex in self.tlt.index[:-1]:
                current_year = self.tlt.iloc[tindex]
                next_year = self.tlt.iloc[tindex + 1]
                print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{next_year}]") # Reference set is next year
                data21_results = model.optimize(self.email,self.solver)
                # print(data21_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data21_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                # Select the single efficiency column and rename it
                data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            
                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data21_results' index is the actual DMU identifier
                data21_component.index = data21_results.index

                dataz21_list.append(data21_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz21_list might be empty
            dataz21 = pd.concat(dataz21_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz21, how='left')



            # --- Calculate Malmquist Indices and components ---
            # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
            for col in ["D11", "D12", "D21"]:
                    if col in self.datazz.columns:
                        self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


            # Calculate ratios, handling potential division by zero or NaN
            # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
            self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
            self.datazz['D12_除以_D11上一期'] = (self.datazz['D12'] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

            self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
            self.datazz['D11_除以_D21上一期'] = (self.datazz['D11'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
            # Malmquist Index (MQ)
            # Handle cases where either ratio is NaN or negative (sqrt of negative)
            self.datazz['product_of_ratios'] = self.datazz['D12_除以_D11上一期'] * self.datazz['D11_除以_D21上一期']
            self.datazz["MQ"] = np.sqrt(self.datazz['product_of_ratios'].clip(lower=0)) # Ensure non-negative before sqrt

            # Technical Efficiency Change (MEFFCH)
            self.datazz["MEFFCH"] = (self.datazz["D11"] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

            # Technical Change (MTECHCH)
            # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
            # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
            ratio3 = (self.datazz["D12"] / self.datazz["D11"]).replace([np.inf, -np.inf], np.nan)
            ratio4 = ( self.datazz['D11_上一期'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
            product_of_ratios_tech = ratio3 * ratio4
            self.datazz["MTECHCH"] = np.sqrt(product_of_ratios_tech.clip(lower=0)) # Ensure non-negative before sqrt

            self.datazz.drop(columns = ['D11_上一期','D12_除以_D11上一期','D21_上一期','D11_除以_D21上一期','product_of_ratios'], inplace=True) # Optional: drop intermediate columns


            print("CONTEMPORARY tech calculation finished.")






        elif (self.hyper_orientedyx and self.rts == RTS_VRS2):
            # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
            # D11: Current period tech, current period frontier (Efficiency change component)
            dataz11_list = []  # List to store D11 results (or components) for each year
            expected_cols =['tei', 'teo']
            output_cols = ['D11_tei', 'D11_teo']
            # --- Loop through years and perform DEA ---
            for tindex in self.tlt.index:
                current_year = self.tlt.iloc[tindex]
                print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{current_year}]")
                data11_results = model.optimize(self.email,self.solver)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data11_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

    
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].\
                    rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data11_results' index is the actual DMU identifier
                data11_component.index = data11_results.index

                dataz11_list.append(data11_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz11_list might be empty
            # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
            # If process_type is 'single_d11', it will have a 'D11' column
            dataz11 = pd.concat(dataz11_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz11, how='left')

            # D12: Current period tech, previous period frontier (Used in Tech change component)
            dataz12_list = []
            expected_cols =['tei', 'teo']
            output_cols = ['D12_tei', 'D12_teo']
            # Loop starts from the second year (index 1)
            for tindex in self.tlt.index[1:]:
                current_year = self.tlt.iloc[tindex]
                previous_year = self.tlt.iloc[tindex - 1]
                print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                data12_results = model.optimize(self.email,self.solver)
                # print(data12_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data12_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                # Select the single efficiency column and rename it
                data12_component = data12_results[expected_cols].\
                    rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data12_results' index is the actual DMU identifier
                data12_component.index = data12_results.index

                dataz12_list.append(data12_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz12_list might be empty
            dataz12 = pd.concat(dataz12_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz12, how='left')





            # D21: Previous period tech, current period frontier (Used in Tech change component)
            dataz21_list = []
            expected_cols =['tei', 'teo']
            output_cols = ['D21_tei', 'D21_teo']
            # Loop goes up to the second to last year (index -1)
            for tindex in self.tlt.index[:-1]:
                current_year = self.tlt.iloc[tindex]
                next_year = self.tlt.iloc[tindex + 1]
                print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{next_year}]") # Reference set is next year
                data21_results = model.optimize(self.email,self.solver)
                # print(data21_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data21_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                # Select the single efficiency column and rename it
                data21_component = data21_results[expected_cols].\
                    rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
            
                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data21_results' index is the actual DMU identifier
                data21_component.index = data21_results.index

                dataz21_list.append(data21_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz21_list might be empty
            dataz21 = pd.concat(dataz21_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz21, how='left')



 
            # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
            # Ensure D11_tei, D11_teo, D12_tei, D12_teo, D21_tei, D21_teo are numeric
            cols_to_numeric = ['D11_tei', 'D11_teo', 'D12_tei', 'D12_teo', 'D21_tei', 'D21_teo']
            for col in cols_to_numeric:
                    if col in self.datazz.columns:
                        self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

            # Calculate previous period's D11 values for both tei and teo
            # Using transform to keep the original DataFrame structure and align by id
            self.datazz['D11_tei_上一期'] = self.datazz.groupby(id)['D11_tei'].transform(lambda x: x.shift(1))
            self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
            self.datazz['D21_tei_上一期'] = self.datazz.groupby(id)['D21_tei'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
            self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change


            # Calculate ratios for tei (Input-oriented Malmquist)
            ratio1_tei = (self.datazz['D12_tei'] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
            ratio2_tei = (self.datazz['D11_tei'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)

            # Calculate ratios for teo (Output-oriented Malmquist)
            ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
            ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

            # Malmquist Index (MQ) - Separate for tei and teo
            self.datazz['product_of_ratios_tei'] = ratio1_tei * ratio2_tei
            self.datazz["MQ_tei"] = np.sqrt(self.datazz['product_of_ratios_tei'].clip(lower=0)) # Ensure non-negative before sqrt

            self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
            self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

            # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
            self.datazz["MEFFCH_tei"] = (self.datazz["D11_tei"] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
            self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)


            # Technical Change (MTECHCH) - Separate for tei and teo
            # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
            # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
            ratio3_tei = (self.datazz["D12_tei"] / self.datazz["D11_tei"]).replace([np.inf, -np.inf], np.nan)
            ratio4_tei = (self.datazz['D11_tei_上一期'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)
            product_of_ratios_tech_tei = ratio3_tei * ratio4_tei
            self.datazz["MTECHCH_tei"] = np.sqrt(product_of_ratios_tech_tei.clip(lower=0)) # Ensure non-negative before sqrt

            ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
            ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
            product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
            self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt


            # Optional: drop intermediate columns
            intermediate_cols_to_drop = [
                'D11_tei_上一期', 'D11_teo_上一期', 'D21_tei_上一期', 'D21_teo_上一期',
                'product_of_ratios_tei', 'product_of_ratios_teo'
            ]
            self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


            print("CONTEMPORARY tech (Hyperbolic yx) calculation finished.")

        elif (self.hyper_orientedyb and self.rts == RTS_VRS1) or (self.hyper_orientedyb and self.rts == RTS_VRS2):
            # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
            # D11: Current period tech, current period frontier (Efficiency change component)
            dataz11_list = []  # List to store D11 results (or components) for each year
            expected_cols =['teuo', 'teo']
            output_cols = ['D11_teuo', 'D11_teo']
            # --- Loop through years and perform DEA ---
            for tindex in self.tlt.index:
                current_year = self.tlt.iloc[tindex]
                print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{current_year}]")
                data11_results = model.optimize(self.email,self.solver)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data11_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

    
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].\
                    rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data11_results' index is the actual DMU identifier
                data11_component.index = data11_results.index

                dataz11_list.append(data11_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz11_list might be empty
            # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
            # If process_type is 'single_d11', it will have a 'D11' column
            dataz11 = pd.concat(dataz11_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz11, how='left')

            # D12: Current period tech, previous period frontier (Used in Tech change component)
            dataz12_list = []
            expected_cols =['teuo', 'teo']
            output_cols = ['D12_teuo', 'D12_teo']
            # Loop starts from the second year (index 1)
            for tindex in self.tlt.index[1:]:
                current_year = self.tlt.iloc[tindex]
                previous_year = self.tlt.iloc[tindex - 1]
                print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                data12_results = model.optimize(self.email,self.solver)
                # print(data12_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data12_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                # Select the single efficiency column and rename it
                data12_component = data12_results[expected_cols].\
                    rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data12_results' index is the actual DMU identifier
                data12_component.index = data12_results.index

                dataz12_list.append(data12_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz12_list might be empty
            dataz12 = pd.concat(dataz12_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz12, how='left')





            # D21: Previous period tech, current period frontier (Used in Tech change component)
            dataz21_list = []
            expected_cols =['teuo', 'teo']
            output_cols = ['D21_teuo', 'D21_teo']
            # Loop goes up to the second to last year (index -1)
            for tindex in self.tlt.index[:-1]:
                current_year = self.tlt.iloc[tindex]
                next_year = self.tlt.iloc[tindex + 1]
                print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                model = DEAweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{next_year}]") # Reference set is next year
                data21_results = model.optimize(self.email,self.solver)
                # print(data21_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data21_results.columns for col in expected_cols):
                    # This check is crucial. If DEA2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DEAweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                # Select the single efficiency column and rename it
                data21_component = data21_results[expected_cols].\
                    rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
            
                # Ensure the index matches the actual DMU index from DEA2 results
                # Assuming data21_results' index is the actual DMU identifier
                data21_component.index = data21_results.index

                dataz21_list.append(data21_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz21_list might be empty
            dataz21 = pd.concat(dataz21_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz21, how='left')



 
            # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
            # Ensure D11_teuo, D11_teo, D12_teuo, D12_teo, D21_teuo, D21_teo are numeric
            cols_to_numeric = ['D11_teuo', 'D11_teo', 'D12_teuo', 'D12_teo', 'D21_teuo', 'D21_teo']
            for col in cols_to_numeric:
                    if col in self.datazz.columns:
                        self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

            # Calculate previous period's D11 values for both tei and teo
            # Using transform to keep the original DataFrame structure and align by id
            self.datazz['D11_teuo_上一期'] = self.datazz.groupby(id)['D11_teuo'].transform(lambda x: x.shift(1))
            self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
            self.datazz['D21_teuo_上一期'] = self.datazz.groupby(id)['D21_teuo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
            self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change


            # Calculate ratios for tei (Input-oriented Malmquist)
            ratio1_teuo = (self.datazz['D12_teuo'] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
            ratio2_teuo = (self.datazz['D11_teuo'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)

            # Calculate ratios for teo (Output-oriented Malmquist)
            ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
            ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

            # Malmquist Index (MQ) - Separate for teuo and teo
            self.datazz['product_of_ratios_teuo'] = ratio1_teuo * ratio2_teuo
            self.datazz["MQ_teuo"] = np.sqrt(self.datazz['product_of_ratios_teuo'].clip(lower=0)) # Ensure non-negative before sqrt

            self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
            self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

            # Technical Efficiency Change (MEFFCH) - Separate for teuo and teo
            self.datazz["MEFFCH_teuo"] = (self.datazz["D11_teuo"] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
            self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)


            # Technical Change (MTECHCH) - Separate for teuo and teo
            # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
            # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
            ratio3_teuo = (self.datazz["D12_teuo"] / self.datazz["D11_teuo"]).replace([np.inf, -np.inf], np.nan)
            ratio4_teuo = (self.datazz['D11_teuo_上一期'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
            product_of_ratios_tech_teuo = ratio3_teuo * ratio4_teuo
            self.datazz["MTECHCH_teuo"] = np.sqrt(product_of_ratios_tech_teuo.clip(lower=0)) # Ensure non-negative before sqrt

            ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
            ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
            product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
            self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt


            # Optional: drop intermediate columns
            intermediate_cols_to_drop = [
                'D11_teuo_上一期', 'D11_teo_上一期', 'D21_teuo_上一期', 'D21_teo_上一期',
                'product_of_ratios_teuo', 'product_of_ratios_teo'
            ]
            self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


            print("CONTEMPORARY tech (Hyperbolic  yb) calculation finished.")


        else:
            raise ValueError(f"Unsupported orientation/RTS combination: input={self.input_oriented}, output={self.output_oriented}, hyper={self.hyper_oriented}, rts={self.rts}")
        
        






class MQDDFweak(MQDEAweak):
    """Malmquist production index (MQPI)
    """

    def __init__(self, data,id,year,sent = "inputvar=outputvar:unoutputvar",  gy=[1], gx=[0], gb=[0], rts=RTS_VRS1, \
                 tech=TOTAL, dynamic = MAL, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """MQDEAt: Calculates Malmquist index using DDF2 for underlying efficiency scores.

        Args:
            data (pandas.DataFrame): input pandas.
            id (str): column name to specify id.
            year (str): column name to specify time.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output distance vector. Defaults to [0].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Note: DDF2 uses RTS_VRS1.
            dynamic (String): MAL (malmquist index) or LUE (luenberger index)
            tech (str): TOTAL or CONTEMPORARY.
            solver (str): The solver to use (e.g., "mosek", "cbc").
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # Initialize MQDEAt model

        # Ensure year column exists and is sortable
        if year not in data.columns:
            raise ValueError(f"Year column '{year}' not found in data.")
        self.tlt = pd.Series(data[year]).drop_duplicates().sort_values()  # 生成时间的列表

        # Parse input/output variables

        self.gy, self.gx, self.gb, self.inputvars,self.outputvars,self.unoutputvars = tools.assert_MQDEAweak(
                        data, sent, gy, gx, gb
                    )

        self.xcol = list(self.inputvars)  # Ensure it's a list for indexing
        self.ycol = list(self.outputvars)  # Ensure it's a list for indexing
        self.bcol = list(self.unoutputvars)  # Ensure it's a list for indexing


        self.tech = tech
        self.rts = rts
        self.dynamic = dynamic
        self.email = email
        self.solver = solver

        # Determine orientation based on gx/gy vectors
        self.input_oriented = sum(self.gx) >= 1 and sum(self.gy) == 0 and sum(self.gb) == 0
        self.output_oriented = sum(self.gy) >= 1 and sum(self.gx) == 0 and sum(self.gb) == 0
        self.unoutput_oriented = sum(self.gb) >= 1 and sum(self.gx) == 0 and sum(self.gy) == 0
        self.hyper_orientedyx = sum(self.gx) >= 1 and sum(self.gy) >= 1 and sum(self.gb) == 0
        self.hyper_orientedyb = sum(self.gb) >= 1 and sum(self.gy) >= 1 and sum(self.gx) == 0
        self.hyper_orientedxb = sum(self.gb) >= 1 and sum(self.gx) >= 1 and sum(self.gy) == 0
        self.hyper_orientedyxb = sum(self.gb) >= 1 and sum(self.gx) >= 1 and sum(self.gy) >= 1 

        # Create a copy of the original data to add results columns
        self.datazz = data.copy()

        # --- Perform DDF calculations using DDF2 based on the chosen technology ---

        if self.tech == TOTAL:
            print("Calculating D11 (Total frontier) for all periods...")

            self.get_total(data,sent,id,year)
            print("TOTAL tech calculation finished.")

        elif self.tech == CONTEMPORARY:
            print("Calculating CONTEMPORARY tech components (D11, D12, D21)...")

            self.get_contemp(data,sent,id,year)
            print("TOTAL tech calculation finished.")

        else:
            raise ValueError(f"Unsupported technology type '{self.tech}'. Must be '{TOTAL}' or '{CONTEMPORARY}'.")



    def optimize(self):
        """Returns the calculated Malmquist index and components DataFrame."""
        # In this implementation, optimize() just returns the pre-calculated results
        # from the __init__ method.
        if not hasattr(self, 'datazz'):
             raise RuntimeError("Malmquist index calculation failed during initialization.")
        return self.datazz





    def get_total(self,data,sent,id,year):
        """Calculate the total efficiency scores for all years."""
        dataz11_list = []  # List to store D11 results (or components) for each year
        # For Total frontier, evaluate each DMU in each year against the frontier of ALL years
        # The baseindex selects the DMU(s) for evaluation in a specific year.
        # The refindex should select ALL DMUs in ALL years for the reference set.
        all_years_ref_index = f"{year}=[{','.join(map(str, self.tlt.tolist()))}]" # Reference set includes all years
        # Loop through each year in the time list
        if self.dynamic == MAL:
            # Determine which columns to expect and how to process based on orientation and RTS
            if self.input_oriented  :
                # Standard case: Expect 'te' and calculate a single 'D11'
                expected_cols = ['tei']
                output_cols = ['D11']
                process_type = 'single_d11'

            elif self.output_oriented:
                # Output-oriented case: Expect 'te' and calculate a single 'D11'
                expected_cols = ['teo']
                output_cols = ['D11']
                process_type = 'single_d11'
            elif self.unoutput_oriented:
                # Undesirable Output-oriented case: Expect 'te' and calculate a single 'D11'
                expected_cols = ['teuo']
                output_cols = ['D11']
                process_type = 'single_d11'
            elif (self.hyper_orientedyx  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['tei', 'teo']
                output_cols = ['D11_tei', 'D11_teo'] # Renaming for clarity
                process_type = 'separate_tei_teo'
            elif (self.hyper_orientedyb  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['teuo', 'teo']
                output_cols = ['D11_teuo', 'D11_teo'] # Renaming for clarity
                process_type = 'separate_teuo_teo'
            elif (self.hyper_orientedxb  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['teuo', 'tei']
                output_cols = ['D11_teuo', 'D11_tei'] # Renaming for clarity
                process_type = 'separate_teuo_tei'
            elif (self.hyper_orientedyxb  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['teuo', 'tei','teo']
                output_cols = ['D11_teuo', 'D11_tei', 'D11_teo'  ] # Renaming for clarity
                process_type = 'separate_teuo_tei_teo'
            else:
                raise ValueError(f"Unsupported orientation/RTS combination: input={self.input_oriented}, output={self.output_oriented}, hyper={self.hyper_oriented}, rts={self.rts}")
        elif self.dynamic == LUE:
            expected_cols = ['objective_value']
            output_cols = ['D11']
            process_type = 'single_d11'
        else:
            raise ValueError(f"Unsupported dynamic type '{self.dynamic}'. Must be '{MAL}' or '{LUE}'.")
        # --- Loop through years and perform DEA ---
        for tindex in self.tlt.index:
            current_year = self.tlt.iloc[tindex]
            print(f"  Evaluating year {current_year} against Total frontier...")

            # Call DDF2 instead of DDFt
            # Use the calculated gx and gy, and the mapped RTS
            model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                            rts=self.rts, baseindex=f"{year}=[{current_year}]",
                            refindex=all_years_ref_index) # Reference set is all years

            # model.optimize() should return a DataFrame with DMU index and result columns
            data11_results = model.optimize(self.email,self.solver)

            # --- Extract/Select the relevant efficiency column(s) ---
            if not all(col in data11_results.columns for col in expected_cols):
                # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                # Consider adding more specific error messages or handling based on model.optimize() status
                raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

            if process_type == 'single_d11':
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            elif process_type == 'separate_tei_teo':
                # Select both 'tei' and 'teo' and rename them
                # Assuming the order in expected_cols is ['tei', 'teo'] if process_type is 'separate_tei_teo'
                data11_component = data11_results[expected_cols].rename(columns={'tei': 'D11_tei', 'teo': 'D11_teo'})
            elif process_type == 'separate_teuo_teo':
                # Select both 'teuo' and 'teo' and rename them
                # Assuming the order in expected_cols is ['tei', 'teo'] if process_type is 'separate_teuo_teo'
                data11_component = data11_results[expected_cols].rename(columns={'teuo': 'D11_teuo', 'teo': 'D11_teo'})
            elif process_type == 'separate_teuo_tei':
                # Select both 'teuo' and 'tei' and rename them
                # Assuming the order in expected_cols is ['teuo', 'tei'] if process_type is 'separate_teuo_tei'
                data11_component = data11_results[expected_cols].rename(columns={'teuo': 'D11_teuo', 'tei': 'D11_tei'})
            elif process_type == 'separate_teuo_tei_teo':
                # Select both ''teuo' 'tei' and 'teo' and rename them
                # Assuming the order in expected_cols is ['teuo','tei', 'teo'] if process_type is 'separate_teuo_tei_teo'
                data11_component = data11_results[expected_cols].rename(columns={'teuo': 'D11_teuo', 'tei': 'D11_tei', 'teo': 'D11_teo'})


            # Ensure the index matches the actual DMU index from DDF2 results
            # Assuming data11_results' index is the actual DMU identifier
            data11_component.index = data11_results.index

            dataz11_list.append(data11_component)

        # --- Concatenate results for all years ---
        # pd.concat handles the case where dataz11_list might be empty
        # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
        # If process_type is 'single_d11', it will have a 'D11' column
        dataz11 = pd.concat(dataz11_list)

        # --- Join results with the main datazz DataFrame ---
        # Assumes self.datazz is initialized before this method is called and has the correct index structure
        # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
        self.datazz = self.datazz.join(dataz11, how='left')

        # --- Calculate Malmquist components based on the processed D11 values ---
        # This calculation relies on the data being sorted correctly (e.g., by DMU index then by year)
        # for the shift(1) operation to compare consecutive years for the same DMU.
        # It's highly recommended that the input 'data' DataFrame is sorted this way
        # before being passed to MQDEAt.
        # Assumes 'id' variable holds the name of the DMU identifier column/index level used for grouping.

        if process_type == 'single_d11':
            # Calculate single MQPI based on D11 (likely EC_total)
            # Check if the D11 column exists and has at least some non-null data after join
            if self.datazz.empty or "D11" not in self.datazz.columns or self.datazz["D11"].isnull().all():
                print("Warning: D11 calculation resulted in no valid data. Cannot compute MQPI.")
                self.datazz["mqpi"] = np.nan
            else:
                if self.dynamic == MAL:
                    self.datazz["D11"] = pd.to_numeric(self.datazz["D11"], errors='coerce')
                    # Calculate previous period value for D11 within each DMU group
                    # Assumes 'id' is the grouping key (column name or index level name)
                    self.datazz['D11_prev'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                    # Compute the ratio (Current D11 / Previous D11)
                    # This is the Efficiency Change component relative to the Total frontier (EC_total)
                    self.datazz["mqpi"] = self.datazz["D11"] / self.datazz["D11_prev"]
                    self.datazz.drop(columns = ['D11_prev'], inplace=True) # Drop the intermediate column
                elif self.dynamic == LUE:
                    self.datazz["D11"] = pd.to_numeric(self.datazz["D11"], errors='coerce')
                    # Calculate previous period value for D11 within each DMU group
                    # Assumes 'id' is the grouping key (column name or index level name)
                    self.datazz['D11_prev'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                    # Compute the ratio (Current D11 / Previous D11)
                    # This is the Efficiency Change component relative to the Total frontier (EC_total)
                    self.datazz["lueni"] = self.datazz["D11_prev"] - self.datazz["D11"]
                    self.datazz.drop(columns = ['D11_prev'], inplace=True) # Drop the intermediate column

        elif process_type == 'separate_tei_teo':
            # Calculate two separate components based on D11_tei and D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_tei', 'D11_teo']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_tei and EC_total_teo components
            result_cols = ['mqpi_tei', 'mqpi_teo']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')

                # Calculate previous period values for both D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev'], inplace=True)


        elif process_type == 'separate_teuo_teo':
            # Calculate two separate components based on D11_teuo and D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_teuo', 'D11_teo']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_teuo and EC_total_teo components
            result_cols = ['mqpi_teuo', 'mqpi_teo']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')

                # Calculate previous period values for both D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev'], inplace=True)

        elif process_type == 'separate_teuo_tei':
            # Calculate two separate components based on D11_teuo and D11_tei
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_teuo', 'D11_tei']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_teuo and EC_total_tei components
            result_cols = ['mqpi_teuo', 'mqpi_tei']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')

                # Calculate previous period values for both D11_teuo and D11_tei within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev'], inplace=True)

        elif process_type == 'separate_teuo_tei_teo':
            # Calculate two separate components based on D11_teuo and D11_tei D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_teuo', 'D11_tei', 'D11_teo']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_teuo and EC_total_tei components
            result_cols = ['mqpi_teuo', 'mqpi_tei', 'mqpi_teo']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all() and self.datazz[calc_cols[2]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
                self.datazz[result_cols[2]] = np.nan

            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')
                self.datazz[calc_cols[2]] = pd.to_numeric(self.datazz[calc_cols[2]], errors='coerce')

                # Calculate previous period values for both D11_teuo and D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[2]}_prev'] = self.datazz.groupby(id)[calc_cols[2]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both teuo and tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                self.datazz[result_cols[2]] = self.datazz[calc_cols[2]] / self.datazz[f'{calc_cols[2]}_prev']

                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev',f'{calc_cols[2]}_prev'], inplace=True)



        # The function modifies self.datazz in place. It doesn't explicitly return datazz.
        # If the calling code expects a return value, uncomment the line below:
        # return self.datazz
        pass # Or return self.datazz if needed by the calling code









    def get_contemp(self,data,sent,id,year):

        # For Total frontier, evaluate each DMU in each year against the frontier of ALL years
        # The baseindex selects the DMU(s) for evaluation in a specific year.
        # The refindex should select ALL DMUs in ALL years for the reference set.
        all_years_ref_index = f"{year}=[{','.join(map(str, self.tlt.tolist()))}]" # Reference set includes all years
        # Loop through each year in the time list

        # Determine which columns to expect and how to process based on orientation and RTS
        if self.dynamic == MAL:
            if self.input_oriented:
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols = ['tei']
                output_cols = ['D11']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDFweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols = ['tei']
                output_cols = ['D12']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols = ['tei']
                output_cols = ['D21']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



                # --- Calculate Malmquist Indices and components ---
                # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
                for col in ["D11", "D12", "D21"]:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


                # Calculate ratios, handling potential division by zero or NaN
                # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
                self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                self.datazz['D12_除以_D11上一期'] = (self.datazz['D12'] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
                self.datazz['D11_除以_D21上一期'] = (self.datazz['D11'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                # Malmquist Index (MQ)
                # Handle cases where either ratio is NaN or negative (sqrt of negative)
                self.datazz['product_of_ratios'] = self.datazz['D12_除以_D11上一期'] * self.datazz['D11_除以_D21上一期']
                self.datazz["MQ"] = np.sqrt(self.datazz['product_of_ratios'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH)
                self.datazz["MEFFCH"] = (self.datazz["D11"] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                # Technical Change (MTECHCH)
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3 = (self.datazz["D12"] / self.datazz["D11"]).replace([np.inf, -np.inf], np.nan)
                ratio4 = ( self.datazz['D11_上一期'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech = ratio3 * ratio4
                self.datazz["MTECHCH"] = np.sqrt(product_of_ratios_tech.clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz.drop(columns = ['D11_上一期','D12_除以_D11上一期','D21_上一期','D11_除以_D21上一期','product_of_ratios'], inplace=True) # Optional: drop intermediate columns


                print("CONTEMPORARY tech calculation finished.")


            elif self.output_oriented:
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols = ['teo']
                output_cols = ['D11']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols = ['teo']
                output_cols = ['D12']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols = ['teo']
                output_cols = ['D21']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



                # --- Calculate Malmquist Indices and components ---
                # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
                for col in ["D11", "D12", "D21"]:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


                # Calculate ratios, handling potential division by zero or NaN
                # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
                self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                self.datazz['D12_除以_D11上一期'] = (self.datazz['D12'] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
                self.datazz['D11_除以_D21上一期'] = (self.datazz['D11'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                # Malmquist Index (MQ)
                # Handle cases where either ratio is NaN or negative (sqrt of negative)
                self.datazz['product_of_ratios'] = self.datazz['D12_除以_D11上一期'] * self.datazz['D11_除以_D21上一期']
                self.datazz["MQ"] = np.sqrt(self.datazz['product_of_ratios'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH)
                self.datazz["MEFFCH"] = (self.datazz["D11"] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                # Technical Change (MTECHCH)
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3 = (self.datazz["D12"] / self.datazz["D11"]).replace([np.inf, -np.inf], np.nan)
                ratio4 = ( self.datazz['D11_上一期'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech = ratio3 * ratio4
                self.datazz["MTECHCH"] = np.sqrt(product_of_ratios_tech.clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz.drop(columns = ['D11_上一期','D12_除以_D11上一期','D21_上一期','D11_除以_D21上一期','product_of_ratios'], inplace=True) # Optional: drop intermediate columns


                print("CONTEMPORARY tech calculation finished.")


            elif self.unoutput_oriented:
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols = ['teuo']
                output_cols = ['D11']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols = ['teuo']
                output_cols = ['D12']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols = ['teuo']
                output_cols = ['D21']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



                # --- Calculate Malmquist Indices and components ---
                # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
                for col in ["D11", "D12", "D21"]:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


                # Calculate ratios, handling potential division by zero or NaN
                # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
                self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                self.datazz['D12_除以_D11上一期'] = (self.datazz['D12'] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
                self.datazz['D11_除以_D21上一期'] = (self.datazz['D11'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                # Malmquist Index (MQ)
                # Handle cases where either ratio is NaN or negative (sqrt of negative)
                self.datazz['product_of_ratios'] = self.datazz['D12_除以_D11上一期'] * self.datazz['D11_除以_D21上一期']
                self.datazz["MQ"] = np.sqrt(self.datazz['product_of_ratios'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH)
                self.datazz["MEFFCH"] = (self.datazz["D11"] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                # Technical Change (MTECHCH)
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3 = (self.datazz["D12"] / self.datazz["D11"]).replace([np.inf, -np.inf], np.nan)
                ratio4 = ( self.datazz['D11_上一期'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech = ratio3 * ratio4
                self.datazz["MTECHCH"] = np.sqrt(product_of_ratios_tech.clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz.drop(columns = ['D11_上一期','D12_除以_D11上一期','D21_上一期','D11_除以_D21上一期','product_of_ratios'], inplace=True) # Optional: drop intermediate columns


                print("CONTEMPORARY tech calculation finished.")


            elif (self.hyper_orientedyx):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['tei', 'teo']
                output_cols = ['D11_tei', 'D11_teo']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['tei', 'teo']
                output_cols = ['D12_tei', 'D12_teo']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['tei', 'teo']
                output_cols = ['D21_tei', 'D21_teo']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_tei, D11_teo, D12_tei, D12_teo, D21_tei, D21_teo are numeric
                cols_to_numeric = ['D11_tei', 'D11_teo', 'D12_tei', 'D12_teo', 'D21_tei', 'D21_teo']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_tei_上一期'] = self.datazz.groupby(id)['D11_tei'].transform(lambda x: x.shift(1))
                self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
                self.datazz['D21_tei_上一期'] = self.datazz.groupby(id)['D21_tei'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change


                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_tei = (self.datazz['D12_tei'] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_tei = (self.datazz['D11_tei'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for tei and teo
                self.datazz['product_of_ratios_tei'] = ratio1_tei * ratio2_tei
                self.datazz["MQ_tei"] = np.sqrt(self.datazz['product_of_ratios_tei'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
                self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_tei"] = (self.datazz["D11_tei"] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)


                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_tei = (self.datazz["D12_tei"] / self.datazz["D11_tei"]).replace([np.inf, -np.inf], np.nan)
                ratio4_tei = (self.datazz['D11_tei_上一期'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_tei = ratio3_tei * ratio4_tei
                self.datazz["MTECHCH_tei"] = np.sqrt(product_of_ratios_tech_tei.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
                self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt


                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_tei_上一期', 'D11_teo_上一期', 'D21_tei_上一期', 'D21_teo_上一期',
                    'product_of_ratios_tei', 'product_of_ratios_teo'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic yx) calculation finished.")



            elif (self.hyper_orientedyb):
                # Hyper + VRS case: Expect 'teuo' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['teuo', 'teo']
                output_cols = ['D11_teuo', 'D11_teo']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['teuo', 'teo']
                output_cols = ['D12_teuo', 'D12_teo']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['teuo', 'teo']
                output_cols = ['D21_teuo', 'D21_teo']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_teuo, D11_teo, D12_teuo, D12_teo, D21_teuo, D21_teo are numeric
                cols_to_numeric = ['D11_teuo', 'D11_teo', 'D12_teuo', 'D12_teo', 'D21_teuo', 'D21_teo']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_teuo_上一期'] = self.datazz.groupby(id)['D11_teuo'].transform(lambda x: x.shift(1))
                self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
                self.datazz['D21_teuo_上一期'] = self.datazz.groupby(id)['D21_teuo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change


                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_teuo = (self.datazz['D12_teuo'] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo = (self.datazz['D11_teuo'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for tei and teo
                self.datazz['product_of_ratios_teuo'] = ratio1_teuo * ratio2_teuo
                self.datazz["MQ_teuo"] = np.sqrt(self.datazz['product_of_ratios_teuo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
                self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_teuo"] = (self.datazz["D11_teuo"] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)


                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_teuo = (self.datazz["D12_teuo"] / self.datazz["D11_teuo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo = (self.datazz['D11_teuo_上一期'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo = ratio3_teuo * ratio4_teuo
                self.datazz["MTECHCH_teuo"] = np.sqrt(product_of_ratios_tech_teuo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
                self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt


                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_teuo_上一期', 'D11_teo_上一期', 'D21_teuo_上一期', 'D21_teo_上一期',
                    'product_of_ratios_teuo', 'product_of_ratios_teo'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic yb) calculation finished.")


            elif (self.hyper_orientedxb):
                # Hyper + VRS case: Expect 'teuo' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['teuo', 'tei']
                output_cols = ['D11_teuo', 'D11_tei']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['teuo', 'tei']
                output_cols = ['D12_teuo', 'D12_tei']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['teuo', 'tei']
                output_cols = ['D21_teuo', 'D21_tei']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_teuo, D11_tei, D12_teuo, D12_tei, D21_teuo, D21_tei are numeric
                cols_to_numeric = ['D11_teuo', 'D11_tei', 'D12_teuo', 'D12_tei', 'D21_teuo', 'D21_tei']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_teuo_上一期'] = self.datazz.groupby(id)['D11_teuo'].transform(lambda x: x.shift(1))
                self.datazz['D11_tei_上一期'] = self.datazz.groupby(id)['D11_tei'].transform(lambda x: x.shift(1))
                self.datazz['D21_teuo_上一期'] = self.datazz.groupby(id)['D21_teuo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_tei_上一期'] = self.datazz.groupby(id)['D21_tei'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change


                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_teuo = (self.datazz['D12_teuo'] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo = (self.datazz['D11_teuo'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_tei = (self.datazz['D12_tei'] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_tei = (self.datazz['D11_tei'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for tei and teo
                self.datazz['product_of_ratios_teuo'] = ratio1_teuo * ratio2_teuo
                self.datazz["MQ_teuo"] = np.sqrt(self.datazz['product_of_ratios_teuo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_tei'] = ratio1_tei * ratio2_tei
                self.datazz["MQ_tei"] = np.sqrt(self.datazz['product_of_ratios_tei'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_teuo"] = (self.datazz["D11_teuo"] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_tei"] = (self.datazz["D11_tei"] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)


                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_teuo = (self.datazz["D12_teuo"] / self.datazz["D11_teuo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo = (self.datazz['D11_teuo_上一期'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo = ratio3_teuo * ratio4_teuo
                self.datazz["MTECHCH_teuo"] = np.sqrt(product_of_ratios_tech_teuo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_tei = (self.datazz["D12_tei"] / self.datazz["D11_tei"]).replace([np.inf, -np.inf], np.nan)
                ratio4_tei = (self.datazz['D11_tei_上一期'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_tei = ratio3_tei * ratio4_tei
                self.datazz["MTECHCH_tei"] = np.sqrt(product_of_ratios_tech_tei.clip(lower=0)) # Ensure non-negative before sqrt


                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_teuo_上一期', 'D11_tei_上一期', 'D21_teuo_上一期', 'D21_tei_上一期',
                    'product_of_ratios_teuo', 'product_of_ratios_tei'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic xb) calculation finished.")


            elif (self.hyper_orientedyxb):
                # Hyper + VRS case: Expect 'tei' 'teuo' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['teuo', 'tei', 'teo']
                output_cols = ['D11_teuo', 'D11_tei', 'D11_teo']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['teuo', 'tei', 'teo']
                output_cols = ['D12_teuo', 'D12_tei', 'D12_teo']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['teuo', 'tei',  'teo']
                output_cols = ['D21_teuo', 'D21_tei', 'D21_teo']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_teuo, D11_tei, D12_teuo, D12_tei, D21_teuo, D21_tei are numeric
                cols_to_numeric = ['D11_teuo', 'D11_tei', 'D11_teo','D12_teuo', 'D12_tei','D12_teo', 'D21_teuo', 'D21_tei','D21_teo']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_teuo_上一期'] = self.datazz.groupby(id)['D11_teuo'].transform(lambda x: x.shift(1))
                self.datazz['D11_tei_上一期'] = self.datazz.groupby(id)['D11_tei'].transform(lambda x: x.shift(1))
                self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
                self.datazz['D21_teuo_上一期'] = self.datazz.groupby(id)['D21_teuo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_tei_上一期'] = self.datazz.groupby(id)['D21_tei'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1))

                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_teuo = (self.datazz['D12_teuo'] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo = (self.datazz['D11_teuo'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for tei (Output-oriented Malmquist)
                ratio1_tei = (self.datazz['D12_tei'] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_tei = (self.datazz['D11_tei'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for tei and teo
                self.datazz['product_of_ratios_teuo'] = ratio1_teuo * ratio2_teuo
                self.datazz["MQ_teuo"] = np.sqrt(self.datazz['product_of_ratios_teuo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_tei'] = ratio1_tei * ratio2_tei
                self.datazz["MQ_tei"] = np.sqrt(self.datazz['product_of_ratios_tei'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
                self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_teuo"] = (self.datazz["D11_teuo"] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_tei"] = (self.datazz["D11_tei"] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)


                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_teuo = (self.datazz["D12_teuo"] / self.datazz["D11_teuo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo = (self.datazz['D11_teuo_上一期'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo = ratio3_teuo * ratio4_teuo
                self.datazz["MTECHCH_teuo"] = np.sqrt(product_of_ratios_tech_teuo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_tei = (self.datazz["D12_tei"] / self.datazz["D11_tei"]).replace([np.inf, -np.inf], np.nan)
                ratio4_tei = (self.datazz['D11_tei_上一期'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_tei = ratio3_tei * ratio4_tei
                self.datazz["MTECHCH_tei"] = np.sqrt(product_of_ratios_tech_tei.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
                self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt


                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_teuo_上一期', 'D11_tei_上一期', 'D11_teo_上一期', 'D21_teuo_上一期', 'D21_tei_上一期','D21_teo_上一期',
                    'product_of_ratios_teuo', 'product_of_ratios_tei' ,'product_of_ratios_teo'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic yxb) calculation finished.")

            else:
                raise ValueError(f"Unsupported orientation/RTS combination: input={self.input_oriented}, output={self.output_oriented}, hyper={self.hyper_oriented}, rts={self.rts}")
        elif self.dynamic == LUE:
            dataz11_list = []  # List to store D11 results (or components) for each year
            expected_cols = ['objective_value']
            output_cols = ['D11']
            # --- Loop through years and perform DEA ---
            for tindex in self.tlt.index:
                current_year = self.tlt.iloc[tindex]
                print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{current_year}]")
                data11_results = model.optimize(self.email,self.solver)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data11_results.columns for col in expected_cols):
                    # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DDFweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

    
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                # Ensure the index matches the actual DMU index from DDF2 results
                # Assuming data11_results' index is the actual DMU identifier
                data11_component.index = data11_results.index

                dataz11_list.append(data11_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz11_list might be empty
            # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
            # If process_type is 'single_d11', it will have a 'D11' column
            dataz11 = pd.concat(dataz11_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz11, how='left')

            # D12: Current period tech, previous period frontier (Used in Tech change component)
            dataz12_list = []
            expected_cols = ['objective_value']
            output_cols = ['D12']
            # Loop starts from the second year (index 1)
            for tindex in self.tlt.index[1:]:
                current_year = self.tlt.iloc[tindex]
                previous_year = self.tlt.iloc[tindex - 1]
                print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                data12_results = model.optimize(self.email,self.solver)
                # print(data12_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data12_results.columns for col in expected_cols):
                    # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                # Select the single efficiency column and rename it
                data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            
                # Ensure the index matches the actual DMU index from DDF2 results
                # Assuming data12_results' index is the actual DMU identifier
                data12_component.index = data12_results.index

                dataz12_list.append(data12_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz12_list might be empty
            dataz12 = pd.concat(dataz12_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz12, how='left')





            # D21: Previous period tech, current period frontier (Used in Tech change component)
            dataz21_list = []
            expected_cols = ['objective_value']
            output_cols = ['D21']
            # Loop goes up to the second to last year (index -1)
            for tindex in self.tlt.index[:-1]:
                current_year = self.tlt.iloc[tindex]
                next_year = self.tlt.iloc[tindex + 1]
                print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                model = DDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{next_year}]") # Reference set is next year
                data21_results = model.optimize(self.email,self.solver)
                # print(data21_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data21_results.columns for col in expected_cols):
                    # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                # Select the single efficiency column and rename it
                data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            
                # Ensure the index matches the actual DMU index from DDF2 results
                # Assuming data21_results' index is the actual DMU identifier
                data21_component.index = data21_results.index

                dataz21_list.append(data21_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz21_list might be empty
            dataz21 = pd.concat(dataz21_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz21, how='left')



            # --- Calculate Malmquist Indices and components ---
            # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
            for col in ["D11", "D12", "D21"]:
                    if col in self.datazz.columns:
                        self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


            # Calculate ratios, handling potential division by zero or NaN
            # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
            self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
            self.datazz['D11上一期_减_D12'] = (self.datazz['D11_上一期'] - self.datazz['D12']).replace([np.inf, -np.inf], np.nan)

            self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
            self.datazz['D21上一期_减_D11'] = (self.datazz['D21_上一期'] / self.datazz['D11']).replace([np.inf, -np.inf], np.nan)
            # Malmquist Index (LQ)
            # Handle cases where either ratio is NaN or negative (sqrt of negative)
            self.datazz['LQ'] = 1/2*(self.datazz['D11上一期_减_D12'] + self.datazz['D21上一期_减_D11'])

            # Technical Efficiency Change (LEFFCH)
            self.datazz["LEFFCH"] = (self.datazz["D11_上一期"] - self.datazz['D11']).replace([np.inf, -np.inf], np.nan)

            # Technical Change (LTECHCH)
            # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
            # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
            diff3 = (self.datazz["D11"] - self.datazz["D12"]).replace([np.inf, -np.inf], np.nan)
            diff4 = ( self.datazz['D21_上一期'] - self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)
            self.datazz["LTECHCH"] = 1/2*(diff3 + diff4)
            self.datazz.drop(columns = ['D11_上一期','D11上一期_减_D12','D21_上一期','D21上一期_减_D11'], inplace=True) # Optional: drop intermediate columns


            print("CONTEMPORARY tech calculation finished.")   
    



class MQNDDFweak(MQDEAweak):
    """Malmquist production index (MQPI)
    """

    def __init__(self, data,id,year,sent = "inputvar=outputvar:unoutputvar",  gy=[1], gx=[0], gb=[0], rts=RTS_VRS1, \
                 tech=TOTAL, dynamic = MAL, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """MQDEAt: Calculates Malmquist index using DDF2 for underlying efficiency scores.

        Args:
            data (pandas.DataFrame): input pandas.
            id (str): column name to specify id.
            year (str): column name to specify time.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L= Y"
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output distance vector. Defaults to [0].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Note: DDF2 uses RTS_VRS1.
            dynamic (String): MAL (malmquist index) or LUE (luenberger index)
            tech (str): TOTAL or CONTEMPORARY.
            solver (str): The solver to use (e.g., "mosek", "cbc").
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        # Initialize MQDEAt model

        # Ensure year column exists and is sortable
        if year not in data.columns:
            raise ValueError(f"Year column '{year}' not found in data.")
        self.tlt = pd.Series(data[year]).drop_duplicates().sort_values()  # 生成时间的列表

        # Parse input/output variables

        self.gy, self.gx, self.gb, self.inputvars,self.outputvars,self.unoutputvars = tools.assert_MQDEAweak(
                        data, sent, gy, gx, gb
                    )

        self.xcol = list(self.inputvars)  # Ensure it's a list for indexing
        self.ycol = list(self.outputvars)  # Ensure it's a list for indexing
        self.bcol = list(self.unoutputvars)  # Ensure it's a list for indexing


        self.tech = tech
        self.rts = rts
        self.dynamic = dynamic
        self.email = email
        self.solver = solver

        # Determine orientation based on gx/gy vectors
        self.input_oriented = sum(self.gx) >= 1 and sum(self.gy) == 0 and sum(self.gb) == 0
        self.output_oriented = sum(self.gy) >= 1 and sum(self.gx) == 0 and sum(self.gb) == 0
        self.unoutput_oriented = sum(self.gb) >= 1 and sum(self.gx) == 0 and sum(self.gy) == 0
        self.hyper_orientedyx = sum(self.gx) >= 1 and sum(self.gy) >= 1 and sum(self.gb) == 0
        self.hyper_orientedyb = sum(self.gb) >= 1 and sum(self.gy) >= 1 and sum(self.gx) == 0
        self.hyper_orientedxb = sum(self.gb) >= 1 and sum(self.gx) >= 1 and sum(self.gy) == 0
        self.hyper_orientedyxb = sum(self.gb) >= 1 and sum(self.gx) >= 1 and sum(self.gy) >= 1 

        # Create a copy of the original data to add results columns
        self.datazz = data.copy()

        # --- Perform DDF calculations using DDF2 based on the chosen technology ---

        if self.tech == TOTAL:
            print("Calculating D11 (Total frontier) for all periods...")

            self.get_total(data,sent,id,year)
            print("TOTAL tech calculation finished.")

        elif self.tech == CONTEMPORARY:
            print("Calculating CONTEMPORARY tech components (D11, D12, D21)...")

            self.get_contemp(data,sent,id,year)
            print("TOTAL tech calculation finished.")

        else:
            raise ValueError(f"Unsupported technology type '{self.tech}'. Must be '{TOTAL}' or '{CONTEMPORARY}'.")



    def optimize(self):
        """Returns the calculated Malmquist index and components DataFrame."""
        # In this implementation, optimize() just returns the pre-calculated results
        # from the __init__ method.
        if not hasattr(self, 'datazz'):
             raise RuntimeError("Malmquist index calculation failed during initialization.")
        return self.datazz




    def get_total(self,data,sent,id,year):
        """Calculate the total efficiency scores for all years."""
        dataz11_list = []  # List to store D11 results (or components) for each year
        # For Total frontier, evaluate each DMU in each year against the frontier of ALL years
        # The baseindex selects the DMU(s) for evaluation in a specific year.
        # The refindex should select ALL DMUs in ALL years for the reference set.
        all_years_ref_index = f"{year}=[{','.join(map(str, self.tlt.tolist()))}]" # Reference set includes all years
        # Loop through each year in the time list

        # Determine which columns to expect and how to process based on orientation and RTS
        if self.dynamic == MAL:
            if self.input_oriented  :
                # Standard case: Expect 'te' and calculate a single 'D11'
                expected_cols = ['tei']
                output_cols = ['D11']
                process_type = 'single_d11'

            elif self.output_oriented:
                # Output-oriented case: Expect 'te' and calculate a single 'D11'
                expected_cols = ['teo']
                output_cols = ['D11']
                process_type = 'single_d11'
            elif self.unoutput_oriented:
                # Undesirable Output-oriented case: Expect 'te' and calculate a single 'D11'
                expected_cols = ['teuo']
                output_cols = ['D11']
                process_type = 'single_d11'
            elif (self.hyper_orientedyx  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['tei', 'teo','tei2o']
                output_cols = ['D11_tei', 'D11_teo', 'D11_tei2o'] # Renaming for clarity
                process_type = 'separate_tei_teo'
            elif (self.hyper_orientedyb  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['teuo', 'teo','teuo2o']
                output_cols = ['D11_teuo', 'D11_teo', 'D11_teuo2o'] # Renaming for clarity
                process_type = 'separate_teuo_teo'
            elif (self.hyper_orientedxb  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['teuo', 'tei', 'teuo2i']
                output_cols = ['D11_teuo', 'D11_tei','D11_teuo2i'] # Renaming for clarity
                process_type = 'separate_teuo_tei'
            elif (self.hyper_orientedyxb  ):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                expected_cols = ['teuo', 'tei','teo','teiuo2o']
                output_cols = ['D11_teuo', 'D11_tei', 'D11_teo' , 'D11_teiuo2o' ] # Renaming for clarity
                process_type = 'separate_teuo_tei_teo'
            else:
                raise ValueError(f"Unsupported orientation/RTS combination: input={self.input_oriented}, output={self.output_oriented}, hyper={self.hyper_oriented}, rts={self.rts}")
        elif self.dynamic == LUE:
            # Luenberger case: Expect 'te' and calculate a single 'D11'
            expected_cols = ['objective_value']
            output_cols = ['D11']
            process_type = 'single_d11'
        else:
            raise ValueError(f"Unsupported dynamic type '{self.dynamic}'. Must be '{MAL}' or '{LUE}'.")
        # --- Loop through years and perform DEA ---
        for tindex in self.tlt.index:
            current_year = self.tlt.iloc[tindex]
            print(f"  Evaluating year {current_year} against Total frontier...")

            # Call DDF2 instead of DDFt
            # Use the calculated gx and gy, and the mapped RTS
            model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                            rts=self.rts, baseindex=f"{year}=[{current_year}]",
                            refindex=all_years_ref_index) # Reference set is all years

            # model.optimize() should return a DataFrame with DMU index and result columns
            data11_results = model.optimize(self.email,self.solver)

            # --- Extract/Select the relevant efficiency column(s) ---
            if not all(col in data11_results.columns for col in expected_cols):
                # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                # Consider adding more specific error messages or handling based on model.optimize() status
                raise KeyError(f"NDDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

            if process_type == 'single_d11':
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            elif process_type == 'separate_tei_teo':
                # Select both 'tei' and 'teo' and rename them
                # Assuming the order in expected_cols is ['tei', 'teo'] if process_type is 'separate_tei_teo'
                data11_component = data11_results[expected_cols].rename(columns={'tei': 'D11_tei', 'teo': 'D11_teo', 'tei2o': 'D11_tei2o'})
            elif process_type == 'separate_teuo_teo':
                # Select both 'teuo' and 'teo' and rename them
                # Assuming the order in expected_cols is ['tei', 'teo'] if process_type is 'separate_teuo_teo'
                data11_component = data11_results[expected_cols].rename(columns={'teuo': 'D11_teuo', 'teo': 'D11_teo', 'teuo2o': 'D11_teuo2o'})
            elif process_type == 'separate_teuo_tei':
                # Select both 'teuo' and 'tei' and rename them
                # Assuming the order in expected_cols is ['teuo', 'tei'] if process_type is 'separate_teuo_tei'
                data11_component = data11_results[expected_cols].rename(columns={'teuo': 'D11_teuo', 'tei': 'D11_tei', 'teuo2i': 'D11_teuo2i'})
            elif process_type == 'separate_teuo_tei_teo':
                # Select both ''teuo' 'tei' and 'teo' and rename them
                # Assuming the order in expected_cols is ['teuo','tei', 'teo'] if process_type is 'separate_teuo_tei_teo'
                data11_component = data11_results[expected_cols].rename(columns={'teuo': 'D11_teuo', 'tei': 'D11_tei', 'teo': 'D11_teo', 'teiuo2o': 'D11_teiuo2o'})


            # Ensure the index matches the actual DMU index from DDF2 results
            # Assuming data11_results' index is the actual DMU identifier
            data11_component.index = data11_results.index

            dataz11_list.append(data11_component)

        # --- Concatenate results for all years ---
        # pd.concat handles the case where dataz11_list might be empty
        # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
        # If process_type is 'single_d11', it will have a 'D11' column
        dataz11 = pd.concat(dataz11_list)

        # --- Join results with the main datazz DataFrame ---
        # Assumes self.datazz is initialized before this method is called and has the correct index structure
        # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
        self.datazz = self.datazz.join(dataz11, how='left')

        # --- Calculate Malmquist components based on the processed D11 values ---
        # This calculation relies on the data being sorted correctly (e.g., by DMU index then by year)
        # for the shift(1) operation to compare consecutive years for the same DMU.
        # It's highly recommended that the input 'data' DataFrame is sorted this way
        # before being passed to MQDEAt.
        # Assumes 'id' variable holds the name of the DMU identifier column/index level used for grouping.

        if process_type == 'single_d11':
            # Calculate single MQPI based on D11 (likely EC_total)
            # Check if the D11 column exists and has at least some non-null data after join
            if self.datazz.empty or "D11" not in self.datazz.columns or self.datazz["D11"].isnull().all():
                print("Warning: D11 calculation resulted in no valid data. Cannot compute MQPI.")
                self.datazz["mqpi"] = np.nan
            else:
                if self.dynamic == MAL:
                    self.datazz["D11"] = pd.to_numeric(self.datazz["D11"], errors='coerce')
                    # Calculate previous period value for D11 within each DMU group
                    # Assumes 'id' is the grouping key (column name or index level name)
                    self.datazz['D11_prev'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                    # Compute the ratio (Current D11 / Previous D11)
                    # This is the Efficiency Change component relative to the Total frontier (EC_total)
                    self.datazz["mqpi"] = self.datazz["D11"] / self.datazz["D11_prev"]
                    self.datazz.drop(columns = ['D11_prev'], inplace=True) # Drop the intermediate column
                elif self.dynamic == LUE:
                    self.datazz["D11"] = pd.to_numeric(self.datazz["D11"], errors='coerce')
                    # Calculate previous period value for D11 within each DMU group
                    # Assumes 'id' is the grouping key (column name or index level name)
                    self.datazz['D11_prev'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                    # Compute the ratio (Current D11 / Previous D11)
                    # This is the Efficiency Change component relative to the Total frontier (EC_total)
                    self.datazz["lueni"] = self.datazz["D11_prev"] - self.datazz["D11"]
                    self.datazz.drop(columns = ['D11_prev'], inplace=True) # Drop the intermediate column
        elif process_type == 'separate_tei_teo':
            # Calculate two separate components based on D11_tei and D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_tei', 'D11_teo']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_tei and EC_total_teo components
            result_cols = ['mqpi_tei', 'mqpi_teo']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all() and self.datazz[calc_cols[2]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
                self.datazz[result_cols[2]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')
                self.datazz[calc_cols[2]] = pd.to_numeric(self.datazz[calc_cols[2]], errors='coerce')

                # Calculate previous period values for both D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[2]}_prev'] = self.datazz.groupby(id)[calc_cols[2]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                self.datazz[result_cols[2]] = self.datazz[calc_cols[2]] / self.datazz[f'{calc_cols[2]}_prev']

                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev', f'{calc_cols[2]}_prev'], inplace=True)


        elif process_type == 'separate_teuo_teo':
            # Calculate two separate components based on D11_teuo and D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_teuo', 'D11_teo', 'D11_teuo2o']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_teuo and EC_total_teo components
            result_cols = ['mqpi_teuo', 'mqpi_teo', 'mqpi_teuo2o']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all() and self.datazz[calc_cols[2]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
                self.datazz[result_cols[2]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')
                self.datazz[calc_cols[2]] = pd.to_numeric(self.datazz[calc_cols[2]], errors='coerce')
                # Calculate previous period values for both D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[2]}_prev'] = self.datazz.groupby(id)[calc_cols[2]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                self.datazz[result_cols[2]] = self.datazz[calc_cols[2]] / self.datazz[f'{calc_cols[2]}_prev']

                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev', f'{calc_cols[2]}_prev'], inplace=True)

        elif process_type == 'separate_teuo_tei':
            # Calculate two separate components based on D11_teuo and D11_tei
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_teuo', 'D11_tei', 'D11_teuo2i']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_teuo and EC_total_tei components
            result_cols = ['mqpi_teuo', 'mqpi_tei', 'mqpi_teuo2i']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all() and self.datazz[calc_cols[2]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
                self.datazz[result_cols[2]] = np.nan
            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')
                self.datazz[calc_cols[2]] = pd.to_numeric(self.datazz[calc_cols[2]], errors='coerce')

                # Calculate previous period values for both D11_teuo and D11_tei within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[2]}_prev'] = self.datazz.groupby(id)[calc_cols[2]].transform(lambda x: x.shift(1))

                # Compute the ratios (Current / Previous) for both tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                self.datazz[result_cols[2]] = self.datazz[calc_cols[2]] / self.datazz[f'{calc_cols[2]}_prev']
                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev', f'{calc_cols[2]}_prev'], inplace=True)

        elif process_type == 'separate_teuo_tei_teo':
            # Calculate two separate components based on D11_teuo and D11_tei D11_teo
            # These are likely the Efficiency Change components relative to the Total frontier
            # in the input and output directions for the hyper-oriented case.
            calc_cols = ['D11_teuo', 'D11_tei', 'D11_teo', 'D11_teiuo2o']
            # Naming the output columns based on user request to modify the 'mqpi' part
            # These represent the EC_total_teuo and EC_total_tei components
            result_cols = ['mqpi_teuo', 'mqpi_tei', 'mqpi_teo', 'mqpi_teiuo2o']

            # Check if the necessary columns exist and have at least some non-null data after join
            if (
                self.datazz.empty or
                not all(col in self.datazz.columns for col in calc_cols) or
                (self.datazz[calc_cols[0]].isnull().all() and self.datazz[calc_cols[1]].isnull().all() \
                 and self.datazz[calc_cols[2]].isnull().all() and self.datazz[calc_cols[3]].isnull().all())
            ):
                print(f"Warning: {', '.join(calc_cols)} calculation resulted in no valid data. Cannot compute {', '.join(result_cols)}.")
                # Create empty columns with NaN values if calculation is not possible
                self.datazz[result_cols[0]] = np.nan
                self.datazz[result_cols[1]] = np.nan
                self.datazz[result_cols[2]] = np.nan
                self.datazz[result_cols[3]] = np.nan

            else:
                # Ensure columns are numeric, coercing errors to NaN
                self.datazz[calc_cols[0]] = pd.to_numeric(self.datazz[calc_cols[0]], errors='coerce')
                self.datazz[calc_cols[1]] = pd.to_numeric(self.datazz[calc_cols[1]], errors='coerce')
                self.datazz[calc_cols[2]] = pd.to_numeric(self.datazz[calc_cols[2]], errors='coerce')
                self.datazz[calc_cols[3]] = pd.to_numeric(self.datazz[calc_cols[3]], errors='coerce')
                # Calculate previous period values for both D11_teuo and D11_tei and D11_teo within each DMU group
                # Assumes 'id' is the grouping key
                self.datazz[f'{calc_cols[0]}_prev'] = self.datazz.groupby(id)[calc_cols[0]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[1]}_prev'] = self.datazz.groupby(id)[calc_cols[1]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[2]}_prev'] = self.datazz.groupby(id)[calc_cols[2]].transform(lambda x: x.shift(1))
                self.datazz[f'{calc_cols[3]}_prev'] = self.datazz.groupby(id)[calc_cols[3]].transform(lambda x: x.shift(1))
                # Compute the ratios (Current / Previous) for both teuo and tei and teo
                self.datazz[result_cols[0]] = self.datazz[calc_cols[0]] / self.datazz[f'{calc_cols[0]}_prev']
                self.datazz[result_cols[1]] = self.datazz[calc_cols[1]] / self.datazz[f'{calc_cols[1]}_prev']
                self.datazz[result_cols[2]] = self.datazz[calc_cols[2]] / self.datazz[f'{calc_cols[2]}_prev']
                self.datazz[result_cols[3]] = self.datazz[calc_cols[3]] / self.datazz[f'{calc_cols[3]}_prev']
                # print(self.datazz)
                # Drop the intermediate columns
                self.datazz.drop(columns = [f'{calc_cols[0]}_prev', f'{calc_cols[1]}_prev',f'{calc_cols[2]}_prev',f'{calc_cols[3]}_prev'], inplace=True)



        # The function modifies self.datazz in place. It doesn't explicitly return datazz.
        # If the calling code expects a return value, uncomment the line below:
        # return self.datazz
        pass # Or return self.datazz if needed by the calling code









    def get_contemp(self,data,sent,id,year):

        # For Total frontier, evaluate each DMU in each year against the frontier of ALL years
        # The baseindex selects the DMU(s) for evaluation in a specific year.
        # The refindex should select ALL DMUs in ALL years for the reference set.
        all_years_ref_index = f"{year}=[{','.join(map(str, self.tlt.tolist()))}]" # Reference set includes all years
        # Loop through each year in the time list

        # Determine which columns to expect and how to process based on orientation and RTS
        if self.dynamic == MAL:
            if self.input_oriented:
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols = ['tei']
                output_cols = ['D11']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDFweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols = ['tei']
                output_cols = ['D12']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols = ['tei']
                output_cols = ['D21']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



                # --- Calculate Malmquist Indices and components ---
                # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
                for col in ["D11", "D12", "D21"]:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


                # Calculate ratios, handling potential division by zero or NaN
                # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
                self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                self.datazz['D12_除以_D11上一期'] = (self.datazz['D12'] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
                self.datazz['D11_除以_D21上一期'] = (self.datazz['D11'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                # Malmquist Index (MQ)
                # Handle cases where either ratio is NaN or negative (sqrt of negative)
                self.datazz['product_of_ratios'] = self.datazz['D12_除以_D11上一期'] * self.datazz['D11_除以_D21上一期']
                self.datazz["MQ"] = np.sqrt(self.datazz['product_of_ratios'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH)
                self.datazz["MEFFCH"] = (self.datazz["D11"] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                # Technical Change (MTECHCH)
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3 = (self.datazz["D12"] / self.datazz["D11"]).replace([np.inf, -np.inf], np.nan)
                ratio4 = ( self.datazz['D11_上一期'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech = ratio3 * ratio4
                self.datazz["MTECHCH"] = np.sqrt(product_of_ratios_tech.clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz.drop(columns = ['D11_上一期','D12_除以_D11上一期','D21_上一期','D11_除以_D21上一期','product_of_ratios'], inplace=True) # Optional: drop intermediate columns


                print("CONTEMPORARY tech calculation finished.")


            elif self.output_oriented:
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols = ['teo']
                output_cols = ['D11']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols = ['teo']
                output_cols = ['D12']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols = ['teo']
                output_cols = ['D21']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



                # --- Calculate Malmquist Indices and components ---
                # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
                for col in ["D11", "D12", "D21"]:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


                # Calculate ratios, handling potential division by zero or NaN
                # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
                self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                self.datazz['D12_除以_D11上一期'] = (self.datazz['D12'] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
                self.datazz['D11_除以_D21上一期'] = (self.datazz['D11'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                # Malmquist Index (MQ)
                # Handle cases where either ratio is NaN or negative (sqrt of negative)
                self.datazz['product_of_ratios'] = self.datazz['D12_除以_D11上一期'] * self.datazz['D11_除以_D21上一期']
                self.datazz["MQ"] = np.sqrt(self.datazz['product_of_ratios'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH)
                self.datazz["MEFFCH"] = (self.datazz["D11"] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                # Technical Change (MTECHCH)
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3 = (self.datazz["D12"] / self.datazz["D11"]).replace([np.inf, -np.inf], np.nan)
                ratio4 = ( self.datazz['D11_上一期'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech = ratio3 * ratio4
                self.datazz["MTECHCH"] = np.sqrt(product_of_ratios_tech.clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz.drop(columns = ['D11_上一期','D12_除以_D11上一期','D21_上一期','D11_除以_D21上一期','product_of_ratios'], inplace=True) # Optional: drop intermediate columns


                print("CONTEMPORARY tech calculation finished.")


            elif self.unoutput_oriented:
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols = ['teuo']
                output_cols = ['D11']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols = ['teuo']
                output_cols = ['D12']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols = ['teuo']
                output_cols = ['D21']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



                # --- Calculate Malmquist Indices and components ---
                # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
                for col in ["D11", "D12", "D21"]:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


                # Calculate ratios, handling potential division by zero or NaN
                # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
                self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
                self.datazz['D12_除以_D11上一期'] = (self.datazz['D12'] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
                self.datazz['D11_除以_D21上一期'] = (self.datazz['D11'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                # Malmquist Index (MQ)
                # Handle cases where either ratio is NaN or negative (sqrt of negative)
                self.datazz['product_of_ratios'] = self.datazz['D12_除以_D11上一期'] * self.datazz['D11_除以_D21上一期']
                self.datazz["MQ"] = np.sqrt(self.datazz['product_of_ratios'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH)
                self.datazz["MEFFCH"] = (self.datazz["D11"] / self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)

                # Technical Change (MTECHCH)
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3 = (self.datazz["D12"] / self.datazz["D11"]).replace([np.inf, -np.inf], np.nan)
                ratio4 = ( self.datazz['D11_上一期'] / self.datazz['D21_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech = ratio3 * ratio4
                self.datazz["MTECHCH"] = np.sqrt(product_of_ratios_tech.clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz.drop(columns = ['D11_上一期','D12_除以_D11上一期','D21_上一期','D11_除以_D21上一期','product_of_ratios'], inplace=True) # Optional: drop intermediate columns


                print("CONTEMPORARY tech calculation finished.")


            elif (self.hyper_orientedyx):
                # Hyper + VRS case: Expect 'tei' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['tei', 'teo', 'tei2o']
                output_cols = ['D11_tei', 'D11_teo', 'D11_tei2o']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['tei', 'teo', 'tei2o']
                output_cols = ['D12_tei', 'D12_teo', 'D12_tei2o']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['tei', 'teo', 'tei2o']
                output_cols = ['D21_tei', 'D21_teo', 'D21_tei2o']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_tei, D11_teo, D12_tei, D12_teo, D21_tei, D21_teo are numeric
                cols_to_numeric = ['D11_tei', 'D11_teo','D11_tei2o', 'D12_tei', 'D12_teo','D12_tei2o', 'D21_tei', 'D21_teo','D21_tei2o']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_tei_上一期'] = self.datazz.groupby(id)['D11_tei'].transform(lambda x: x.shift(1))
                self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
                self.datazz['D11_tei2o_上一期'] = self.datazz.groupby(id)['D11_tei2o'].transform(lambda x: x.shift(1))

                self.datazz['D21_tei_上一期'] = self.datazz.groupby(id)['D21_tei'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_tei2o_上一期'] = self.datazz.groupby(id)['D21_tei2o'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change


                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_tei = (self.datazz['D12_tei'] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_tei = (self.datazz['D11_tei'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for tei2o (Output-oriented Malmquist)
                ratio1_tei2o = (self.datazz['D12_tei2o'] / self.datazz['D11_tei2o_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_tei2o = (self.datazz['D11_tei2o'] / self.datazz['D21_tei2o_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for tei and teo
                self.datazz['product_of_ratios_tei'] = ratio1_tei * ratio2_tei
                self.datazz["MQ_tei"] = np.sqrt(self.datazz['product_of_ratios_tei'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
                self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_tei2o'] = ratio1_tei2o * ratio2_tei2o
                self.datazz["MQ_tei2o"] = np.sqrt(self.datazz['product_of_ratios_tei2o'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_tei"] = (self.datazz["D11_tei"] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_tei2o"] = (self.datazz["D11_tei2o"] / self.datazz['D11_tei2o_上一期']).replace([np.inf, -np.inf], np.nan)


                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_tei = (self.datazz["D12_tei"] / self.datazz["D11_tei"]).replace([np.inf, -np.inf], np.nan)
                ratio4_tei = (self.datazz['D11_tei_上一期'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_tei = ratio3_tei * ratio4_tei
                self.datazz["MTECHCH_tei"] = np.sqrt(product_of_ratios_tech_tei.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
                self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_tei2o = (self.datazz["D12_tei2o"] / self.datazz["D11_tei2o"]).replace([np.inf, -np.inf], np.nan)
                ratio4_tei2o = (self.datazz['D11_tei2o_上一期'] / self.datazz['D21_tei2o_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_tei2o = ratio3_tei2o * ratio4_tei2o
                self.datazz["MTECHCH_tei2o"] = np.sqrt(product_of_ratios_tech_tei2o.clip(lower=0)) # Ensure non-negative before sqrt

                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_tei_上一期', 'D11_teo_上一期','D11_tei2o_上一期', 'D21_tei_上一期', 'D21_teo_上一期', 'D21_tei2o_上一期',
                    'product_of_ratios_tei', 'product_of_ratios_teo','product_of_ratios_tei2o'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic yx) calculation finished.")



            elif (self.hyper_orientedyb):
                # Hyper + VRS case: Expect 'teuo' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['teuo', 'teo','teuo2o']
                output_cols = ['D11_teuo', 'D11_teo', 'D11_teuo2o']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['teuo', 'teo', 'teuo2o']
                output_cols = ['D12_teuo', 'D12_teo', 'D12_teuo2o']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['teuo', 'teo', 'teuo2o']
                output_cols = ['D21_teuo', 'D21_teo', 'D21_teuo2o']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_teuo, D11_teo, D12_teuo, D12_teo, D21_teuo, D21_teo are numeric
                cols_to_numeric = ['D11_teuo', 'D11_teo',  'D11_teuo2o','D12_teuo', 'D12_teo', 'D12_teuo2o','D21_teuo', 'D21_teo','D21_teuo2o']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_teuo_上一期'] = self.datazz.groupby(id)['D11_teuo'].transform(lambda x: x.shift(1))
                self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
                self.datazz['D11_teuo2o_上一期'] = self.datazz.groupby(id)['D11_teuo2o'].transform(lambda x: x.shift(1))
                self.datazz['D21_teuo_上一期'] = self.datazz.groupby(id)['D21_teuo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teuo2o_上一期'] = self.datazz.groupby(id)['D21_teuo2o'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change

                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_teuo = (self.datazz['D12_teuo'] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo = (self.datazz['D11_teuo'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teuo2o (Output-oriented Malmquist)
                ratio1_teuo2o = (self.datazz['D12_teuo2o'] / self.datazz['D11_teuo2o_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo2o = (self.datazz['D11_teuo2o'] / self.datazz['D21_teuo2o_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for teuo and teo  
                self.datazz['product_of_ratios_teuo'] = ratio1_teuo * ratio2_teuo
                self.datazz["MQ_teuo"] = np.sqrt(self.datazz['product_of_ratios_teuo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
                self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teuo2o'] = ratio1_teuo2o * ratio2_teuo2o
                self.datazz["MQ_teuo2o"] = np.sqrt(self.datazz['product_of_ratios_teuo2o'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_teuo"] = (self.datazz["D11_teuo"] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teuo2o"] = (self.datazz["D11_teuo2o"] / self.datazz['D11_teuo2o_上一期']).replace([np.inf, -np.inf], np.nan)


                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_teuo = (self.datazz["D12_teuo"] / self.datazz["D11_teuo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo = (self.datazz['D11_teuo_上一期'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo = ratio3_teuo * ratio4_teuo
                self.datazz["MTECHCH_teuo"] = np.sqrt(product_of_ratios_tech_teuo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
                self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teuo2o = (self.datazz["D12_teuo2o"] / self.datazz["D11_teuo2o"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo2o = (self.datazz['D11_teuo2o_上一期'] / self.datazz['D21_teuo2o_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo2o = ratio3_teuo2o * ratio4_teuo2o
                self.datazz["MTECHCH_teuo2o"] = np.sqrt(product_of_ratios_tech_teuo2o.clip(lower=0)) # Ensure non-negative before sqrt
                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_teuo_上一期', 'D11_teo_上一期', 'D11_teuo2o_上一期','D21_teuo_上一期', 'D21_teo_上一期','D21_teuo2o_上一期',
                    'product_of_ratios_teuo', 'product_of_ratios_teo', 'product_of_ratios_teuo2o'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic yb) calculation finished.")


            elif (self.hyper_orientedxb):
                # Hyper + VRS case: Expect 'teuo' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['teuo', 'tei','teuo2i']
                output_cols = ['D11_teuo', 'D11_tei', 'D11_teuo2i']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['teuo', 'tei','teuo2i']
                output_cols = ['D12_teuo', 'D12_tei', 'D12_teuo2i']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['teuo', 'tei','teuo2i']
                output_cols = ['D21_teuo', 'D21_tei', 'D21_teuo2i']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_teuo, D11_tei, D12_teuo, D12_tei, D21_teuo, D21_tei are numeric
                cols_to_numeric = ['D11_teuo', 'D11_tei','D11_teuo2i', 'D12_teuo', 'D12_tei','D12_teuo2i', 'D21_teuo', 'D21_tei','D21_teuo2i']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_teuo_上一期'] = self.datazz.groupby(id)['D11_teuo'].transform(lambda x: x.shift(1))
                self.datazz['D11_tei_上一期'] = self.datazz.groupby(id)['D11_tei'].transform(lambda x: x.shift(1))
                self.datazz['D11_teuo2i_上一期'] = self.datazz.groupby(id)['D11_teuo2i'].transform(lambda x: x.shift(1))
                self.datazz['D21_teuo_上一期'] = self.datazz.groupby(id)['D21_teuo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_tei_上一期'] = self.datazz.groupby(id)['D21_tei'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teuo2i_上一期'] = self.datazz.groupby(id)['D21_teuo2i'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change

                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_teuo = (self.datazz['D12_teuo'] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo = (self.datazz['D11_teuo'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_tei = (self.datazz['D12_tei'] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_tei = (self.datazz['D11_tei'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teuo2i = (self.datazz['D12_teuo2i'] / self.datazz['D11_teuo2i_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo2i = (self.datazz['D11_teuo2i'] / self.datazz['D21_teuo2i_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for tei and teo
                self.datazz['product_of_ratios_teuo'] = ratio1_teuo * ratio2_teuo
                self.datazz["MQ_teuo"] = np.sqrt(self.datazz['product_of_ratios_teuo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_tei'] = ratio1_tei * ratio2_tei
                self.datazz["MQ_tei"] = np.sqrt(self.datazz['product_of_ratios_tei'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teuo2i'] = ratio1_teuo2i * ratio2_teuo2i
                self.datazz["MQ_teuo2i"] = np.sqrt(self.datazz['product_of_ratios_teuo2i'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_teuo"] = (self.datazz["D11_teuo"] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_tei"] = (self.datazz["D11_tei"] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teuo2i"] = (self.datazz["D11_teuo2i"] / self.datazz['D11_teuo2i_上一期']).replace([np.inf, -np.inf], np.nan)

                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_teuo = (self.datazz["D12_teuo"] / self.datazz["D11_teuo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo = (self.datazz['D11_teuo_上一期'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo = ratio3_teuo * ratio4_teuo
                self.datazz["MTECHCH_teuo"] = np.sqrt(product_of_ratios_tech_teuo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_tei = (self.datazz["D12_tei"] / self.datazz["D11_tei"]).replace([np.inf, -np.inf], np.nan)
                ratio4_tei = (self.datazz['D11_tei_上一期'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_tei = ratio3_tei * ratio4_tei
                self.datazz["MTECHCH_tei"] = np.sqrt(product_of_ratios_tech_tei.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teuo2i = (self.datazz["D12_teuo2i"] / self.datazz["D11_teuo2i"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo2i = (self.datazz['D11_teuo2i_上一期'] / self.datazz['D21_teuo2i_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo2i = ratio3_teuo2i * ratio4_teuo2i
                self.datazz["MTECHCH_teuo2i"] = np.sqrt(product_of_ratios_tech_teuo2i.clip(lower=0)) # Ensure non-negative before sqrt

                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_teuo_上一期', 'D11_tei_上一期', 'D11_teuo2i_上一期', 'D21_teuo_上一期', 'D21_tei_上一期', 'D21_teuo2i_上一期',
                    'product_of_ratios_teuo', 'product_of_ratios_tei', 'product_of_ratios_teuo2i'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic xb) calculation finished.")


            elif (self.hyper_orientedyxb):
                # Hyper + VRS case: Expect 'tei' 'teuo' and 'teo' and keep them separate
                # D11: Current period tech, current period frontier (Efficiency change component)
                dataz11_list = []  # List to store D11 results (or components) for each year
                expected_cols =['teuo', 'tei', 'teo', 'teiuo2o']
                output_cols = ['D11_teuo', 'D11_tei', 'D11_teo', 'D11_teiuo2o']
                # --- Loop through years and perform DEA ---
                for tindex in self.tlt.index:
                    current_year = self.tlt.iloc[tindex]
                    print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{current_year}]")
                    data11_results = model.optimize(self.email,self.solver)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data11_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

        
                    # Select the single efficiency column and rename it
                    data11_component = data11_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})


                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data11_results' index is the actual DMU identifier
                    data11_component.index = data11_results.index

                    dataz11_list.append(data11_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz11_list might be empty
                # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
                # If process_type is 'single_d11', it will have a 'D11' column
                dataz11 = pd.concat(dataz11_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz11, how='left')

                # D12: Current period tech, previous period frontier (Used in Tech change component)
                dataz12_list = []
                expected_cols =['teuo', 'tei', 'teo', 'teiuo2o']
                output_cols = ['D12_teuo', 'D12_tei', 'D12_teo', 'D12_teiuo2o']
                # Loop starts from the second year (index 1)
                for tindex in self.tlt.index[1:]:
                    current_year = self.tlt.iloc[tindex]
                    previous_year = self.tlt.iloc[tindex - 1]
                    print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                    data12_results = model.optimize(self.email,self.solver)
                    # print(data12_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data12_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                    # Select the single efficiency column and rename it
                    data12_component = data12_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})

                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data12_results' index is the actual DMU identifier
                    data12_component.index = data12_results.index

                    dataz12_list.append(data12_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz12_list might be empty
                dataz12 = pd.concat(dataz12_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz12, how='left')





                # D21: Previous period tech, current period frontier (Used in Tech change component)
                dataz21_list = []
                expected_cols =['teuo', 'tei',  'teo', 'teiuo2o']
                output_cols = ['D21_teuo', 'D21_tei', 'D21_teo', 'D21_teiuo2o']
                # Loop goes up to the second to last year (index -1)
                for tindex in self.tlt.index[:-1]:
                    current_year = self.tlt.iloc[tindex]
                    next_year = self.tlt.iloc[tindex + 1]
                    print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                    model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx, gb=self.gb,
                                    rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                    refindex=f"{year}=[{next_year}]") # Reference set is next year
                    data21_results = model.optimize(self.email,self.solver)
                    # print(data21_results)

                    # --- Extract/Select the relevant efficiency column(s) ---
                    if not all(col in data21_results.columns for col in expected_cols):
                        # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                        # Consider adding more specific error messages or handling based on model.optimize() status
                        raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                    # Select the single efficiency column and rename it
                    data21_component = data21_results[expected_cols].\
                        rename(columns={expected_cols[0]: output_cols[0],expected_cols[1]: output_cols[1]})
                
                    # Ensure the index matches the actual DMU index from DDF2 results
                    # Assuming data21_results' index is the actual DMU identifier
                    data21_component.index = data21_results.index

                    dataz21_list.append(data21_component)

                # --- Concatenate results for all years ---
                # pd.concat handles the case where dataz21_list might be empty
                dataz21 = pd.concat(dataz21_list)

                # --- Join results with the main datazz DataFrame ---
                # Assumes self.datazz is initialized before this method is called and has the correct index structure
                # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
                self.datazz = self.datazz.join(dataz21, how='left')



    
                # --- Calculate Malmquist Indices and components for Hyperbolic VRS ---
                # Ensure D11_teuo, D11_tei, D12_teuo, D12_tei, D21_teuo, D21_tei are numeric
                cols_to_numeric = ['D11_teuo', 'D11_tei', 'D11_teo','D11_teiuo2o', 'D12_teuo', 'D12_tei','D12_teo','D12_teiuo2o',  'D21_teuo', 'D21_tei','D21_teo','D21_teiuo2o']
                for col in cols_to_numeric:
                        if col in self.datazz.columns:
                            self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')

                # Calculate previous period's D11 values for both tei and teo
                # Using transform to keep the original DataFrame structure and align by id
                self.datazz['D11_teuo_上一期'] = self.datazz.groupby(id)['D11_teuo'].transform(lambda x: x.shift(1))
                self.datazz['D11_tei_上一期'] = self.datazz.groupby(id)['D11_tei'].transform(lambda x: x.shift(1))
                self.datazz['D11_teo_上一期'] = self.datazz.groupby(id)['D11_teo'].transform(lambda x: x.shift(1))
                self.datazz['D11_teiuo2o_上一期'] = self.datazz.groupby(id)['D11_teiuo2o'].transform(lambda x: x.shift(1))
                self.datazz['D21_teuo_上一期'] = self.datazz.groupby(id)['D21_teuo'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_tei_上一期'] = self.datazz.groupby(id)['D21_tei'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                self.datazz['D21_teo_上一期'] = self.datazz.groupby(id)['D21_teo'].transform(lambda x: x.shift(1))
                self.datazz['D21_teiuo2o_上一期'] = self.datazz.groupby(id)['D21_teiuo2o'].transform(lambda x: x.shift(1)) # Need previous D21 for tech change
                # Calculate ratios for tei (Input-oriented Malmquist)
                ratio1_teuo = (self.datazz['D12_teuo'] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teuo = (self.datazz['D11_teuo'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for tei (Output-oriented Malmquist)
                ratio1_tei = (self.datazz['D12_tei'] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_tei = (self.datazz['D11_tei'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teo = (self.datazz['D12_teo'] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teo = (self.datazz['D11_teo'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)

                # Calculate ratios for teo (Output-oriented Malmquist)
                ratio1_teiuo2o = (self.datazz['D12_teiuo2o'] / self.datazz['D11_teiuo2o_上一期']).replace([np.inf, -np.inf], np.nan)
                ratio2_teiuo2o = (self.datazz['D11_teiuo2o'] / self.datazz['D21_teiuo2o_上一期']).replace([np.inf, -np.inf], np.nan)

                # Malmquist Index (MQ) - Separate for tei and teo
                self.datazz['product_of_ratios_teuo'] = ratio1_teuo * ratio2_teuo
                self.datazz["MQ_teuo"] = np.sqrt(self.datazz['product_of_ratios_teuo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_tei'] = ratio1_tei * ratio2_tei
                self.datazz["MQ_tei"] = np.sqrt(self.datazz['product_of_ratios_tei'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teo'] = ratio1_teo * ratio2_teo
                self.datazz["MQ_teo"] = np.sqrt(self.datazz['product_of_ratios_teo'].clip(lower=0)) # Ensure non-negative before sqrt

                self.datazz['product_of_ratios_teiuo2o'] = ratio1_teiuo2o * ratio2_teiuo2o
                self.datazz["MQ_teiuo2o"] = np.sqrt(self.datazz['product_of_ratios_teiuo2o'].clip(lower=0)) # Ensure non-negative before sqrt

                # Technical Efficiency Change (MEFFCH) - Separate for tei and teo
                self.datazz["MEFFCH_teuo"] = (self.datazz["D11_teuo"] / self.datazz['D11_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_tei"] = (self.datazz["D11_tei"] / self.datazz['D11_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teo"] = (self.datazz["D11_teo"] / self.datazz['D11_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                self.datazz["MEFFCH_teiuo2o"] = (self.datazz["D11_teiuo2o"] / self.datazz['D11_teiuo2o_上一期']).replace([np.inf, -np.inf], np.nan)


                # Technical Change (MTECHCH) - Separate for tei and teo
                # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
                # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
                ratio3_teuo = (self.datazz["D12_teuo"] / self.datazz["D11_teuo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teuo = (self.datazz['D11_teuo_上一期'] / self.datazz['D21_teuo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teuo = ratio3_teuo * ratio4_teuo
                self.datazz["MTECHCH_teuo"] = np.sqrt(product_of_ratios_tech_teuo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_tei = (self.datazz["D12_tei"] / self.datazz["D11_tei"]).replace([np.inf, -np.inf], np.nan)
                ratio4_tei = (self.datazz['D11_tei_上一期'] / self.datazz['D21_tei_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_tei = ratio3_tei * ratio4_tei
                self.datazz["MTECHCH_tei"] = np.sqrt(product_of_ratios_tech_tei.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teo = (self.datazz["D12_teo"] / self.datazz["D11_teo"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teo = (self.datazz['D11_teo_上一期'] / self.datazz['D21_teo_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teo = ratio3_teo * ratio4_teo
                self.datazz["MTECHCH_teo"] = np.sqrt(product_of_ratios_tech_teo.clip(lower=0)) # Ensure non-negative before sqrt

                ratio3_teiuo2o = (self.datazz["D12_teiuo2o"] / self.datazz["D11_teiuo2o"]).replace([np.inf, -np.inf], np.nan)
                ratio4_teiuo2o = (self.datazz['D11_teiuo2o_上一期'] / self.datazz['D21_teiuo2o_上一期']).replace([np.inf, -np.inf], np.nan)
                product_of_ratios_tech_teiuo2o = ratio3_teiuo2o * ratio4_teiuo2o
                self.datazz["MTECHCH_teiuo2o"] = np.sqrt(product_of_ratios_tech_teiuo2o.clip(lower=0)) # Ensure non-negative before sqrt

                # Optional: drop intermediate columns
                intermediate_cols_to_drop = [
                    'D11_teuo_上一期', 'D11_tei_上一期', 'D11_teo_上一期','D11_teiuo2o_上一期', 'D21_teuo_上一期', 'D21_tei_上一期','D21_teo_上一期','D21_teiuo2o_上一期',
                    'product_of_ratios_teuo', 'product_of_ratios_tei' ,'product_of_ratios_teo','product_of_ratios_teiuo2o'
                ]
                self.datazz.drop(columns = intermediate_cols_to_drop, inplace=True, errors='ignore')


                print("CONTEMPORARY tech (Hyperbolic yxb) calculation finished.")

            else:
                raise ValueError(f"Unsupported orientation/RTS combination: input={self.input_oriented}, output={self.output_oriented}, hyper={self.hyper_oriented}, rts={self.rts}")
        elif self.dynamic==LUE:

            # D11: Current period tech, current period frontier (Efficiency change component)
            dataz11_list = []  # List to store D11 results (or components) for each year
            expected_cols = ['objective_value']
            output_cols = ['D11']
            # --- Loop through years and perform DEA ---
            for tindex in self.tlt.index:
                current_year = self.tlt.iloc[tindex]
                print(f"  Evaluating year {current_year} against year {current_year} (D11 component)...") # Verbose print

                model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{current_year}]")
                data11_results = model.optimize(self.email,self.solver)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data11_results.columns for col in expected_cols):
                    # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DDFweak2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data11_results.columns)}")

    
                # Select the single efficiency column and rename it
                data11_component = data11_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})


                # Ensure the index matches the actual DMU index from DDF2 results
                # Assuming data11_results' index is the actual DMU identifier
                data11_component.index = data11_results.index

                dataz11_list.append(data11_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz11_list might be empty
            # If process_type is 'hyper_vrs', the resulting DataFrame will have 'D11_tei' and 'D11_teo' columns
            # If process_type is 'single_d11', it will have a 'D11' column
            dataz11 = pd.concat(dataz11_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz11 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz11, how='left')

            # D12: Current period tech, previous period frontier (Used in Tech change component)
            dataz12_list = []
            expected_cols = ['objective_value']
            output_cols = ['D12']
            # Loop starts from the second year (index 1)
            for tindex in self.tlt.index[1:]:
                current_year = self.tlt.iloc[tindex]
                previous_year = self.tlt.iloc[tindex - 1]
                print(f"  Evaluating year {current_year} against year {previous_year} (D12 component)...") # Verbose print

                model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{previous_year}]") # Reference set is previous year
                data12_results = model.optimize(self.email,self.solver)
                # print(data12_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data12_results.columns for col in expected_cols):
                    # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data12_results.columns)}")


                # Select the single efficiency column and rename it
                data12_component = data12_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            
                # Ensure the index matches the actual DMU index from DDF2 results
                # Assuming data12_results' index is the actual DMU identifier
                data12_component.index = data12_results.index

                dataz12_list.append(data12_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz12_list might be empty
            dataz12 = pd.concat(dataz12_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz12 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz12, how='left')





            # D21: Previous period tech, current period frontier (Used in Tech change component)
            dataz21_list = []
            expected_cols = ['objective_value']
            output_cols = ['D21']
            # Loop goes up to the second to last year (index -1)
            for tindex in self.tlt.index[:-1]:
                current_year = self.tlt.iloc[tindex]
                next_year = self.tlt.iloc[tindex + 1]
                print(f"  Evaluating year {current_year} against year {next_year} (D21 component)...") # Verbose print

                model = NDDFweak2(data=data, sent=sent, gy=self.gy, gx=self.gx,gb=self.gb,
                                rts=self.rts, baseindex=f"{year}=[{current_year}]",
                                refindex=f"{year}=[{next_year}]") # Reference set is next year
                data21_results = model.optimize(self.email,self.solver)
                # print(data21_results)

                # --- Extract/Select the relevant efficiency column(s) ---
                if not all(col in data21_results.columns for col in expected_cols):
                    # This check is crucial. If DDF2 doesn't return the expected columns, stop.
                    # Consider adding more specific error messages or handling based on model.optimize() status
                    raise KeyError(f"DDF2 results DataFrame for year {current_year} does not contain required columns '{expected_cols}' for the chosen orientation/RTS combination. Available columns: {list(data21_results.columns)}")


                # Select the single efficiency column and rename it
                data21_component = data21_results[expected_cols].rename(columns={expected_cols[0]: output_cols[0]})
            
                # Ensure the index matches the actual DMU index from DDF2 results
                # Assuming data21_results' index is the actual DMU identifier
                data21_component.index = data21_results.index

                dataz21_list.append(data21_component)

            # --- Concatenate results for all years ---
            # pd.concat handles the case where dataz21_list might be empty
            dataz21 = pd.concat(dataz21_list)

            # --- Join results with the main datazz DataFrame ---
            # Assumes self.datazz is initialized before this method is called and has the correct index structure
            # Assumes the index of self.datazz and dataz21 should be aligned for joining (e.g., MultiIndex of (DMU, Year))
            self.datazz = self.datazz.join(dataz21, how='left')




            # --- Calculate Malmquist Indices and components ---
            # Ensure D11, D12, D21 are numeric and handle potential NaNs or Infs from division
            for col in ["D11", "D12", "D21"]:
                    if col in self.datazz.columns:
                        self.datazz[col] = pd.to_numeric(self.datazz[col], errors='coerce')


            # Calculate ratios, handling potential division by zero or NaN
            # 使用 transform 是因为我们希望结果的长度与原DataFrame相同，并且索引对齐
            self.datazz['D11_上一期'] = self.datazz.groupby(id)['D11'].transform(lambda x: x.shift(1))
            self.datazz['D11上一期_减_D12'] = (self.datazz['D11_上一期'] - self.datazz['D12']).replace([np.inf, -np.inf], np.nan)

            self.datazz['D21_上一期'] = self.datazz.groupby(id)['D21'].transform(lambda x: x.shift(1))
            self.datazz['D21上一期_减_D11'] = (self.datazz['D21_上一期'] / self.datazz['D11']).replace([np.inf, -np.inf], np.nan)
            # Malmquist Index (LQ)
            # Handle cases where either ratio is NaN or negative (sqrt of negative)
            self.datazz['LQ'] = 1/2*(self.datazz['D11上一期_减_D12'] + self.datazz['D21上一期_减_D11'])

            # Technical Efficiency Change (LEFFCH)
            self.datazz["LEFFCH"] = (self.datazz["D11_上一期"] - self.datazz['D11']).replace([np.inf, -np.inf], np.nan)

            # Technical Change (LTECHCH)
            # Formula: sqrt((D_{t-1}(x_t, y_t) / D_t(x_t, y_t)) * (D_{t-1}(x_{t-1}, y_{t-1}) / D_t(x_{t-1}, y_{t-1})))
            # Using D notation from code: sqrt((D12 / D11) * (D11.shift(1) / D21.shift(1)))
            diff3 = (self.datazz["D11"] - self.datazz["D12"]).replace([np.inf, -np.inf], np.nan)
            diff4 = ( self.datazz['D21_上一期'] - self.datazz['D11_上一期']).replace([np.inf, -np.inf], np.nan)
            self.datazz["LTECHCH"] = 1/2*(diff3 + diff4)
            self.datazz.drop(columns = ['D11_上一期','D11上一期_减_D12','D21_上一期','D21上一期_减_D11'], inplace=True) # Optional: drop intermediate columns


            print("CONTEMPORARY tech calculation finished.")

