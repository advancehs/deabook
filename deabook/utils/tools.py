# import dependencies
from re import compile
import ast
from os import environ
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverFactory, SolverManagerFactory, check_available_solvers
from ..constant import CET_ADDI, CET_MULT, CET_Model_Categories, OPT_LOCAL, OPT_DEFAULT, RTS_CRS
__email_re = compile(r'([^@]+@[^@]+\.[a-zA-Z0-9]+)$')
from deabook.constant import FUN_PROD, OPT_LOCAL,RTS_VRS1, RTS_VRS2, CET_ADDI, CET_MULT

def get_remote_solvers():
    import pyomo.neos.kestrel
    kestrel = pyomo.neos.kestrel.kestrelAMPL()
    return list(
        set(
            [
                name.split(":")[0].lower()
                for name in kestrel.solvers()
            ]
        )
    )

def check_remote_solver(solver="mosek"):
    solver_list = get_remote_solvers()
    return bool(solver in solver_list)

def check_local_solver(solver="mosek"):
    return bool(check_available_solvers(solver))


def set_neos_email(address):
    """pass email address to NEOS server

    Args:
        address (String): your own vaild email address.
    """
    if address == OPT_LOCAL:
        print("Optimizing locally.")
        return False
    if not __email_re.match(address):
        raise ValueError("Invalid email address.")
    environ['NEOS_EMAIL'] = address
    return True

def optimize_model4(model, ind, solver=OPT_DEFAULT):
    if solver is not OPT_DEFAULT:
        assert_solver_available_locally(solver)

    solver_instance = SolverFactory(solver)
    # print("Estimating dmu{} locally with {} solver.".format(
    #     ind, solver), flush=True)
    return str(solver_instance.solve(model, tee=False)), 1


def optimize_model2(model, actual_index, use_neos, cet, solver=OPT_DEFAULT):
    """
    Optimizes a single Pyomo model for a specific DMU, handling solver selection
    based on CET type and NEOS configuration.

    Args:
        model (ConcreteModel): The Pyomo model for the current DMU.
        actual_index: The original data index of the current DMU being optimized.
                      Used for logging/messages.
        use_neos (bool): False for locally.
        cet (str): Model category ('additive' or 'multiplicative').
        solver (string): The name of the solver to use. Defaults to OPT_DEFAULT.

    Returns:
        tuple: (problem_status_str, optimization_status_str)
               problem_status_str (str): String representation of the problem status
                                         (e.g., 'optimal', 'infeasible', 'error').
               optimization_status_str (str): String representation of the solver status
                                            (e.g., 'optimal', 'feasible', 'aborted', 'error').
    """
    chosen_solver_name = solver

    # Call set_neos_email first, as in the original function
    # This function's return value determines the local/remote branch

    try:
        if not use_neos: # Local solving branch
            if chosen_solver_name is not OPT_DEFAULT:
                 # Assert local availability only if a specific solver is requested
                 assert_solver_available_locally(chosen_solver_name)
            elif cet == CET_ADDI:
                chosen_solver_name = "mosek"
            elif cet == CET_MULT:
                # Original function raised an error here for local multiplicative
                raise ValueError(
                    "Please specify the solver for optimizing multiplicative model locally.")

            # print(f"DMU {actual_index}: Estimating the {CET_Model_Categories.get(cet, 'unknown model type')} locally with {chosen_solver_name} solver.", flush=True)

            # Get the solver instance for local execution
            solver_instance = SolverFactory(chosen_solver_name)

            # Check if SolverFactory successfully loaded a solver instance
            if solver_instance is None:
                 raise RuntimeError(f"Local solver '{chosen_solver_name}' not found or not available.")

            # Solve the model locally
            results = solver_instance.solve(model, tee=False) # No 'opt' argument for local solve generally

        else: # Remote (NEOS) solving branch
            if chosen_solver_name is OPT_DEFAULT:
                if cet is CET_ADDI:
                    chosen_solver_name = "mosek"
                elif cet == CET_MULT:
                    chosen_solver_name = "knitro"
                # If solver is not OPT_DEFAULT, the provided solver name is used for NEOS too.

            # print(f"DMU {actual_index}: Estimating the {CET_Model_Categories.get(cet, 'unknown model type')} remotely with {chosen_solver_name} solver via NEOS.", flush=True)

            # Get the solver instance for NEOS execution
            solver_instance = SolverFactory(chosen_solver_name)

            # Check if SolverFactory successfully loaded a solver instance
            if solver_instance is None:
                 raise RuntimeError(f"Remote solver '{chosen_solver_name}' not found or not available via NEOS.")


            # Solve the model remotely via NEOS
            # Pass the solver name via 'opt' argument for NEOS
            results = solver_instance.solve(model, tee=False, opt=chosen_solver_name)

        # --- Common logic after solving ---

        # Check if the solver returned a valid results object
        if results is None or not hasattr(results, 'Problem') or not hasattr(results, 'Solver'):
             raise RuntimeError("Solver did not return a valid results object or results format is unexpected.")

        optimization_status = results.Solver.status

        # Convert Pyomo status enums to string names
        # Use .name attribute if the status object is not None
        optimization_status_str = optimization_status.name if optimization_status else "Unknown"

        # print(f"DMU {actual_index} status:  Solver='{optimization_status_str}'", flush=True)

        return optimization_status_str

    except ValueError as ve:
        # Catch specific ValueErrors raised (like the multiplicative local solver requirement)
        print(f"DMU {actual_index}: Configuration Error - {ve}", flush=True)
        return "config_error", "config_error" # Custom status for configuration problems
    except Exception as e:
        # Catch any other exceptions during the solving process
        print(f"DMU {actual_index}: Error during optimization - {e}", flush=True)
        # Return specific error strings to indicate failure
        return "solve_error", "solve_error" # Using custom status for solve errors




def optimize_model(model, email, cet, solver=OPT_DEFAULT):
    optimization_status = 0
    if not set_neos_email(email):
        if solver is not OPT_DEFAULT:
            # assert_solver_available_locally(solver)
            pass
        elif cet == CET_ADDI:
            solver = "mosek"
        elif cet == CET_MULT:
            raise ValueError(
                "Please specify the solver for optimizing multiplicative model locally.")
        try:    
            solver_instance = SolverFactory(solver)
            print("Estimating the {} locally with {} solver.".format(
            CET_Model_Categories[cet], solver), flush=True)
            return solver_instance.solve(model, tee=True), 1
        except:
            from amplpy import modules
            solver_instance = SolverFactory(solver+"nl", executable=modules.find(solver), 
                                                solve_io="nl",
                                                options = {'soltype':1}
                                                )
            print("Estimating the {} locally with {} solver.".format(
            CET_Model_Categories[cet], solver), flush=True)
            return solver_instance.solve(model, tee=True), 1
    else:
        if solver is OPT_DEFAULT and cet is CET_ADDI:
            solvers = ["mosek"]
        elif solver is OPT_DEFAULT and cet == CET_MULT:
            solvers = ["knitro"]
        else:
            solvers = [solver]
        for solver in solvers:
            model, optimization_status = __try_remote_solver(
                model, cet, solver)
            if optimization_status == 1:
                return model, optimization_status
        raise Exception("Remote solvers are temporarily not available.")

def __try_remote_solver(model, cet, solver):
    solver_instance = SolverManagerFactory('neos')
    try:
        print("Estimating the {} remotely with {} solver.".format(
            CET_Model_Categories[cet], solver), flush=True)
        return solver_instance.solve(model, tee=True, opt=solver), 1
    except:
        print("Remote {} solver is not available now.".format(solver))
        return model, 0





def optimize_model3(model, solver=OPT_DEFAULT):
    if solver is not OPT_DEFAULT:
        assert_solver_available_locally(solver)

    solver_instance = SolverFactory(solver)
    # print("Estimating dmu{} locally with {} solver.".format(
    #     ind, solver), flush=True)
    return str(solver_instance.solve(model, tee=False)), 1

# def optimize_model(model, email, cet, solver=OPT_DEFAULT):
#     if not set_neos_email(email):
#         if solver is not OPT_DEFAULT:
#             assert_solver_available_locally(solver)
#         elif cet == CET_ADDI:
#             solver = "mosek"
#         elif cet == CET_MULT:
#             raise ValueError(
#                 "Please specify the solver for optimizing multiplicative model locally.")
#         solver_instance = SolverFactory(solver )
#         print("Estimating the {} locally with {} solver.".format(
#             CET_Model_Categories[cet], solver), flush=True)
#         return solver_instance.solve(model, tee=False ), 1
#     else:
#         if solver is OPT_DEFAULT and cet is CET_ADDI:
#             solver = "mosek"
#         elif solver is OPT_DEFAULT and cet == CET_MULT:
#             solver = "knitro"
#         solver_instance = SolverFactory(solver )
#         print("Estimating the {} remotely with {} solver.".format(
#             CET_Model_Categories[cet], solver), flush=True)
#         return solver_instance.solve(model, tee=False , opt=solver), 1


def trans_list(li):
    if type(li) == list:
        return li
    return li.tolist()


def to_1d_list(li):
    if type(li) == int or type(li) == float:
        return [li]
    if type(li[0]) == list:
        rl = []
        for i in range(len(li)):
            rl.append(li[i][0])
        return rl
    return li


def to_2d_list(li):
    if type(li[0]) != list:
        rl = []
        for value in li:
            rl.append([value])
        return rl
    return li


def assert_valid_basic_data(y, x, z=None):
    y = trans_list(y)
    x = trans_list(x)

    y = to_1d_list(y)
    x = to_2d_list(x)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    if len(y_shape) == 2 and y_shape[1] != 1:
        raise ValueError(
            "The multidimensional output data is supported by direciontal based models.")

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    return y, x, z

def assert_DEA(data, sent, gy, gx,baseindex,refindex):

    if type(baseindex) != type(None):
        varname1 = baseindex.split('=')[0]
        print(baseindex)
        varvalue1 = ast.literal_eval(baseindex.split('=')[1])
        data_base= data.loc[data[varname1].isin(varvalue1)]
    else:
        data_base= data
    data_index = data_base.index
    if type(refindex) != type(None):
        varname = refindex.split('=')[0]
        varvalue = ast.literal_eval(refindex.split('=')[1])

        data_ref = data.loc[data[varname].isin(varvalue)]
    else:
        data_ref = data

    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')

    x = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in outputvars])
    xref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in inputvars])
    yref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in outputvars])

    y, x, yref, xref, gy, gx = assert_DEA1(y, x, yref, xref, gy, gx)
    return y, x, yref, xref, gy, gx,data_index

def assert_DEA1(y, x, yref, xref, gy, gx):
    y = trans_list(y)
    x = trans_list(x)

    y = to_2d_list(y)
    x = to_2d_list(x)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    print(gx,"#############")
    if sum(gx)>=1 and sum(gy)>=1:
        raise ValueError(
            "gy and gx can not be bigger than 1 together.")

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    yref = trans_list(yref)
    xref = trans_list(xref)

    yref = to_2d_list(yref)
    xref = to_2d_list(xref)

    yref_shape = np.asarray(yref).shape
    xref_shape = np.asarray(xref).shape

    if yref_shape[0] != xref_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in xref and yref.")
    if yref_shape[1] != np.asarray(y).shape[1]:
        raise ValueError(
            "Number of outputs must be the same in y and yref.")
    if xref_shape[1] != np.asarray(x).shape[1]:
        raise ValueError(
            "Number of inputs must be the same in x and xref.")

    return y, x, yref, xref, gy, gx


def assert_DEA_with_indices(data, sent, gy, gx, baseindex, refindex):
    """
    Validates DEA input data and returns processed numpy arrays and original indices.

    Args:
        data (pandas.DataFrame): The input data containing all variables and DMUs.
        sent (str): String specifying input and output variables. e.g., "K L = Y"
        gy (list): Output distance vector.
        gx (list): Input distance vector.
        baseindex (str, optional): Filter for evaluated DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).
        refindex (str, optional): Filter for reference DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).

    Returns:
        tuple: (y, x, yref, xref, gy, gx, evaluated_data_index, reference_data_index)
               y, x: numpy arrays for evaluated outputs and inputs.
               yref, xref: numpy arrays for reference outputs and inputs.
               gy, gx: Processed gy and gx lists.
               evaluated_data_index: List of original index labels for evaluated DMUs.
               reference_data_index: List of original index labels for reference DMUs.

    Raises:
        ValueError: If data format or dimensions are inconsistent.
    """
    # 1. Filter data based on baseindex to get evaluated DMUs
    if baseindex is not None:
        parts = baseindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'baseindex' format. Expected 'varname=value'")
        varname1 = parts[0].strip()
        try:
            varvalue1 = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'baseindex': {parts[1].strip()}")

        # Ensure varvalue1 is iterable for isin
        if not isinstance(varvalue1, (list, tuple, set)):
             varvalue1 = [varvalue1]

        if varname1 not in data.columns:
             raise ValueError(f"Variable '{varname1}' specified in 'baseindex' not found in data columns.")

        # Use .loc for filtering and .copy() to avoid potential SettingWithCopyWarning
        data_base = data.loc[data[varname1].isin(varvalue1)].copy()
    else:
        data_base = data.copy() # Use .copy() for default case too

    if data_base.empty:
         raise ValueError("Filtering with 'baseindex' resulted in an empty dataset for evaluated DMUs.")

    # Store evaluated DMU indices (original index labels)
    evaluated_data_index = data_base.index.tolist()


    # 2. Filter data based on refindex to get reference DMUs
    if refindex is not None:
        parts = refindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'refindex' format. Expected 'varname=value'")
        varname = parts[0].strip()
        try:
            varvalue = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'refindex': {parts[1].strip()}")

        # Ensure varvalue is iterable for isin
        if not isinstance(varvalue, (list, tuple, set)):
             varvalue = [varvalue]

        if varname not in data.columns:
             raise ValueError(f"Variable '{varname}' specified in 'refindex' not found in data columns.")

        # Use .loc for filtering and .copy()
        data_ref = data.loc[data[varname].isin(varvalue)].copy()
    else:
        data_ref = data.copy() # Use .copy() for default case too

    if data_ref.empty:
         raise ValueError("Filtering with 'refindex' resulted in an empty dataset for reference DMUs.")

    # Store reference DMU indices (original index labels)
    reference_data_index = data_ref.index.tolist()


    # 3. Parse sent string to identify input and output variables
    parts = sent.split('=')
    if len(parts) != 2:
        raise ValueError("Invalid 'sent' format. Expected 'inputs = outputs[:unwanted_outputs]'")

    input_part = parts[0].strip()
    # Split output part by ':', take the first part (outputs)
    output_part = parts[1].split(':')[0].strip()

    inputvars = input_part.split(' ')
    outputvars = output_part.split(' ')

    # Basic check if specified variables exist in the original data columns
    all_vars = inputvars + outputvars
    missing_vars = [v for v in all_vars if v not in data.columns]
    if missing_vars:
        raise ValueError(f"Specified variables not found in data columns: {missing_vars}")

    # 4. Extract data into numpy arrays
    # Use .values to get underlying numpy array from pandas Series/DataFrame column
    try:
        x = np.column_stack([data_base[selected].values for selected in inputvars])
        y = np.column_stack([data_base[selected].values for selected in outputvars])
        xref = np.column_stack([data_ref[selected].values for selected in inputvars])
        yref = np.column_stack([data_ref[selected].values for selected in outputvars])
    except KeyError as e:
         # This should ideally be caught by the check above, but as a safeguard
         raise ValueError(f"Error extracting data for variables: {e}")
    except Exception as e:
         # Catch other potential errors during numpy array creation
         raise ValueError(f"An error occurred during data extraction: {e}")


    # 5. Process gy, gx and perform validation checks (replicating assert_DEA1 logic)
    # Apply list transformations as in assert_DEA1 before checks
    # Ensure tools module and these functions are available
    try:
        y_list = to_2d_list(trans_list(y))
        x_list = to_2d_list(trans_list(x))
        yref_list = to_2d_list(trans_list(yref))
        xref_list = to_2d_list(trans_list(xref))
        gy_list = to_1d_list(gy)
        gx_list = to_1d_list(gx)
    except NameError:
        raise NameError("Helper functions (trans_list, to_2d_list, to_1d_list) from tools module are required but not found.")
    except Exception as e:
        raise RuntimeError(f"Error during data list transformation: {e}")


    # Now perform the checks using the list representations (as in original assert_DEA1)
    y_shape = np.asarray(y_list).shape # Use np.asarray to get shape like assert_DEA1
    x_shape = np.asarray(x_list).shape
    yref_shape = np.asarray(yref_list).shape
    xref_shape = np.asarray(xref_list).shape


    # Check number of DMUs match for evaluated set
    if y_shape[0] != x_shape[0]:
        raise ValueError(f"Number of evaluated DMUs differs between outputs ({y_shape[0]}) and inputs ({x_shape[0]}).")

    # Check number of DMUs match for reference set
    if yref_shape[0] != xref_shape[0]:
        raise ValueError(f"Number of reference DMUs differs between outputs ({yref_shape[0]}) and inputs ({xref_shape[0]}).")

    # Check number of variables match between evaluated and reference sets
    if yref_shape[1] != y_shape[1]:
        raise ValueError(f"Number of outputs differs between evaluated ({y_shape[1]}) and reference ({yref_shape[1]}) sets.")
    if xref_shape[1] != x_shape[1]:
        raise ValueError(f"Number of inputs differs between evaluated ({x_shape[1]}) and reference ({xref_shape[1]}) sets.")

    # Optional: Check if gx/gy lengths match variable counts if orientation is active
    # Note: Original assert_DEA1 didn't strictly enforce this, but it's good practice.
    # The code in the first DEA class *assumes* gx[j] and gy[k] exist.
    # Let's add these checks.
    if sum(gx_list) >= 1 and len(gx_list) != x.shape[1]:
         raise ValueError(f"Length of gx ({len(gx_list)}) must match the number of inputs ({x.shape[1]}) when input orientation is used.")
    if sum(gy_list) >= 1 and len(gy_list) != y.shape[1]:
         raise ValueError(f"Length of gy ({len(gy_list)}) must match the number of outputs ({y.shape[1]}) when output orientation is used.")


    # 6. Return processed numpy arrays and the original index lists
    # We return the numpy arrays (y, x, yref, xref) as they are suitable for numerical
    # processing in Pyomo, and the shapes have been validated against the list shapes.
    # We return the processed lists for gy and gx.
    return y, x, yref, xref, gy_list, gx_list, evaluated_data_index, reference_data_index




def assert_MQDEA(data, sent, gy, gx):
    """
    Validates DEA input data and returns processed numpy arrays and original indices.

    Args:
        data (pandas.DataFrame): The input data containing all variables and DMUs.
        sent (str): String specifying input and output variables. e.g., "K L = Y"
        gy (list): Output distance vector.
        gx (list): Input distance vector.

    Returns:
        tuple: (gy, gx)
               y, x: numpy arrays for evaluated outputs and inputs.
               yref, xref: numpy arrays for reference outputs and inputs.
               gy, gx: Processed gy and gx lists.

    Raises:
        ValueError: If data format or dimensions are inconsistent.
    """

    gy_list = to_1d_list(gy)
    gx_list = to_1d_list(gx)

    # 3. Parse sent string to identify input and output variables
    parts = sent.split('=')
    if len(parts) != 2:
        raise ValueError("Invalid 'sent' format. Expected 'inputs = outputs[:unwanted_outputs]'")

    input_part = parts[0].strip()
    # Split output part by ':', take the first part (outputs)
    output_part = parts[1].split(':')[0].strip()

    inputvars = input_part.split(' ')
    outputvars = output_part.split(' ')

    # Basic check if specified variables exist in the original data columns
    all_vars = inputvars + outputvars
    missing_vars = [v for v in all_vars if v not in data.columns]
    if missing_vars:
        raise ValueError(f"Specified variables not found in data columns: {missing_vars}")


    return gy_list, gx_list, inputvars,outputvars



def assert_MQDEAweak(data, sent, gy, gx, gb):
    """
    Validates DEA input data and returns processed numpy arrays and original indices.

    Args:
        data (pandas.DataFrame): The input data containing all variables and DMUs.
        sent (str): String specifying input and output variables. e.g., "K L = Y"
        gy (list): Output distance vector.
        gx (list): Input distance vector.
        gb (list): Undesirable output distance vector.

    Returns:
        tuple: (gy, gx, gb)
               y, x, b: numpy arrays for evaluated outputs, inputs and undesirable outputs.
               yref, xref, bref: numpy arrays for reference outputs, inputs and undesirable outputs.
               gy, gx, gb: Processed gy, gx and gb lists.

    Raises:
        ValueError: If data format or dimensions are inconsistent.
    """

    gy_list = to_1d_list(gy)
    gx_list = to_1d_list(gx)
    gb_list = to_1d_list(gb)

    # 3. Parse sent string to identify input and output variables
    parts = sent.split('=')
    if len(parts) != 2:
        raise ValueError("Invalid 'sent' format. Expected 'inputs = outputs[:unwanted_outputs]'")

    input_part = parts[0].strip()
    # Split output part by ':', take the first part (outputs)
    output_part = parts[1].split(':')[0].strip()
    unoutput_part = parts[1].split(':')[1].strip()

    inputvars = input_part.split(' ')
    outputvars = output_part.split(' ')
    unoutputvars = unoutput_part.split(' ')

    # Basic check if specified variables exist in the original data columns
    all_vars = inputvars + outputvars + unoutputvars
    missing_vars = [v for v in all_vars if v not in data.columns]
    if missing_vars:
        raise ValueError(f"Specified variables not found in data columns: {missing_vars}")


    return gy_list, gx_list, gb_list, inputvars,outputvars,unoutputvars


def assert_valid_mupltiple_y_data(y, x):
    y = trans_list(y)
    x = trans_list(x)

    y = to_2d_list(y)
    x = to_2d_list(x)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")
    return y, x


def assert_DEAweak(data, sent, gy, gx, gb, baseindex, refindex):
    if type(baseindex) != type(None):
        varname1 = baseindex.split('=')[0]
        print(baseindex)
        varvalue1 = ast.literal_eval(baseindex.split('=')[1])
        data_base= data.loc[data[varname1].isin(varvalue1)]
    else:
        data_base= data

    if type(refindex) != type(None):
        varname = refindex.split('=')[0]
        varvalue = ast.literal_eval(refindex.split('=')[1])

        data_ref = data.loc[data[varname].isin(varvalue)]
    else:
        data_ref = data

    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')

    x = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in outputvars])
    b = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in unoutputvars])

    xref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in inputvars])
    yref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in outputvars])
    bref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in unoutputvars])

    y, x, b,  gy, gx, gb, yref, xref, bref = assert_DEAweak1(y, x, b, gy, gx, gb, yref, xref, bref)

    return y, x, b,  gy, gx, gb, yref, xref, bref

def assert_DEAweak1(y, x, b, gy, gx, gb, yref, xref, bref):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    gb = to_1d_list(gb)
    print(gx,"#############")
    if (sum(gx)>=1 and sum(gy)>=1) or (sum(gx)>=1 and sum(gb)>=1) or (sum(gy)>=1 and sum(gb)>=1):
        raise ValueError(
            "gy, gx and gb can not be bigger than 1 together.")


    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")
    if y_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in b and y.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if b_shape[1] != len(gb):
        raise ValueError("Number of inputs must be the same in b and gb.")

    yref = trans_list(yref)
    xref = trans_list(xref)
    bref = trans_list(bref)

    yref = to_2d_list(yref)
    xref = to_2d_list(xref)
    bref = to_2d_list(bref)

    yref_shape = np.asarray(yref).shape
    xref_shape = np.asarray(xref).shape
    bref_shape = np.asarray(bref).shape

    if yref_shape[0] != xref_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in xref and yref.")
    if yref_shape[1] != np.asarray(y).shape[1]:
        raise ValueError(
            "Number of outputs must be the same in y and yref.")
    if xref_shape[1] != np.asarray(x).shape[1]:
        raise ValueError(
            "Number of inputs must be the same in x and xref.")


    return y, x, b,yref, xref, bref,gy, gx, gb


def assert_DEAweak_with_indices(data, sent, gy, gx,gb, baseindex, refindex):
    """
    Validates DEA input data and returns processed numpy arrays and original indices.

    Args:
        data (pandas.DataFrame): The input data containing all variables and DMUs.
        sent (str): String specifying input and output variables. e.g., "K L = Y:CO2"
        gy (list): Output distance vector.
        gx (list): Input distance vector.
        gb (list): Input distance vector.
        baseindex (str, optional): Filter for evaluated DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).
        refindex (str, optional): Filter for reference DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).

    Returns:
        tuple: (y, x, yref, xref, gy, gx, evaluated_data_index, reference_data_index)
               y, x: numpy arrays for evaluated outputs and inputs.
               yref, xref: numpy arrays for reference outputs and inputs.
               gy, gx: Processed gy and gx lists.
               evaluated_data_index: List of original index labels for evaluated DMUs.
               reference_data_index: List of original index labels for reference DMUs.

    Raises:
        ValueError: If data format or dimensions are inconsistent.
    """
    # 1. Filter data based on baseindex to get evaluated DMUs
    if baseindex is not None:
        parts = baseindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'baseindex' format. Expected 'varname=value'")
        varname1 = parts[0].strip()
        try:
            varvalue1 = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'baseindex': {parts[1].strip()}")

        # Ensure varvalue1 is iterable for isin
        if not isinstance(varvalue1, (list, tuple, set)):
             varvalue1 = [varvalue1]

        if varname1 not in data.columns:
             raise ValueError(f"Variable '{varname1}' specified in 'baseindex' not found in data columns.")

        # Use .loc for filtering and .copy() to avoid potential SettingWithCopyWarning
        data_base = data.loc[data[varname1].isin(varvalue1)].copy()
    else:
        data_base = data.copy() # Use .copy() for default case too

    if data_base.empty:
         raise ValueError("Filtering with 'baseindex' resulted in an empty dataset for evaluated DMUs.")

    # Store evaluated DMU indices (original index labels)
    evaluated_data_index = data_base.index.tolist()


    # 2. Filter data based on refindex to get reference DMUs
    if refindex is not None:
        parts = refindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'refindex' format. Expected 'varname=value'")
        varname = parts[0].strip()
        try:
            varvalue = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'refindex': {parts[1].strip()}")

        # Ensure varvalue is iterable for isin
        if not isinstance(varvalue, (list, tuple, set)):
             varvalue = [varvalue]

        if varname not in data.columns:
             raise ValueError(f"Variable '{varname}' specified in 'refindex' not found in data columns.")

        # Use .loc for filtering and .copy()
        data_ref = data.loc[data[varname].isin(varvalue)].copy()
    else:
        data_ref = data.copy() # Use .copy() for default case too

    if data_ref.empty:
         raise ValueError("Filtering with 'refindex' resulted in an empty dataset for reference DMUs.")

    # Store reference DMU indices (original index labels)
    reference_data_index = data_ref.index.tolist()


    # 3. Parse sent string to identify input and output variables
    parts = sent.split('=')
    if len(parts) != 2:
        raise ValueError("Invalid 'sent' format. Expected 'inputs = outputs[:unwanted_outputs]'")

    input_part = parts[0].strip()
    # Split output part by ':', take the first part (outputs)
    output_part = parts[1].strip().split(':')[0].strip()
    unoutput_part = parts[1].strip().split(':')[1].strip()

    inputvars = input_part.split(' ')
    outputvars = output_part.split(' ')
    unoutputvars = unoutput_part.split(' ')

    # Basic check if specified variables exist in the original data columns
    all_vars = inputvars + outputvars + unoutputvars
    missing_vars = [v for v in all_vars if v not in data.columns]
    if missing_vars:
        raise ValueError(f"Specified variables not found in data columns: {missing_vars}")

    # 4. Extract data into numpy arrays
    # Use .values to get underlying numpy array from pandas Series/DataFrame column
    try:
        x = np.column_stack([data_base[selected].values for selected in inputvars])
        y = np.column_stack([data_base[selected].values for selected in outputvars])
        b = np.column_stack([data_base[selected].values for selected in unoutputvars])
        xref = np.column_stack([data_ref[selected].values for selected in inputvars])
        yref = np.column_stack([data_ref[selected].values for selected in outputvars])
        bref = np.column_stack([data_ref[selected].values for selected in unoutputvars])
    except KeyError as e:
         # This should ideally be caught by the check above, but as a safeguard
         raise ValueError(f"Error extracting data for variables: {e}")
    except Exception as e:
         # Catch other potential errors during numpy array creation
         raise ValueError(f"An error occurred during data extraction: {e}")


    # 5. Process gy, gx and perform validation checks (replicating assert_DEA1 logic)
    # Apply list transformations as in assert_DEA1 before checks
    # Ensure tools module and these functions are available
    try:
        y_list = to_2d_list(trans_list(y))
        x_list = to_2d_list(trans_list(x))
        b_list = to_2d_list(trans_list(b))
        yref_list = to_2d_list(trans_list(yref))
        xref_list = to_2d_list(trans_list(xref))
        bref_list = to_2d_list(trans_list(bref))
        gy_list = to_1d_list(gy)
        gx_list = to_1d_list(gx)
        gb_list = to_1d_list(gb)

    except NameError:
        raise NameError("Helper functions (trans_list, to_2d_list, to_1d_list) from tools module are required but not found.")
    except Exception as e:
        raise RuntimeError(f"Error during data list transformation: {e}")


    # Now perform the checks using the list representations (as in original assert_DEA1)
    y_shape = np.asarray(y_list).shape # Use np.asarray to get shape like assert_DEA1
    x_shape = np.asarray(x_list).shape
    b_shape = np.asarray(b_list).shape
    yref_shape = np.asarray(yref_list).shape
    xref_shape = np.asarray(xref_list).shape
    bref_shape = np.asarray(bref_list).shape

    # Check number of DMUs match for evaluated set
    if y_shape[0] != x_shape[0]:
        raise ValueError(f"Number of evaluated DMUs differs between outputs ({y_shape[0]}) and inputs ({x_shape[0]}).")
    if y_shape[0] != b_shape[0]:
        raise ValueError(f"Number of evaluated DMUs differs between outputs ({y_shape[0]}) and unoutputs ({b_shape[0]}).")

    # Check number of DMUs match for reference set
    if yref_shape[0] != xref_shape[0]:
        raise ValueError(f"Number of reference DMUs differs between outputs ({yref_shape[0]}) and inputs ({xref_shape[0]}).")
    if yref_shape[0] != bref_shape[0]:
        raise ValueError(f"Number of reference DMUs differs between outputs ({yref_shape[0]}) and unoutputs ({bref_shape[0]}).")

    # Check number of variables match between evaluated and reference sets
    if yref_shape[1] != y_shape[1]:
        raise ValueError(f"Number of outputs differs between evaluated ({y_shape[1]}) and reference ({yref_shape[1]}) sets.")
    if xref_shape[1] != x_shape[1]:
        raise ValueError(f"Number of inputs differs between evaluated ({x_shape[1]}) and reference ({xref_shape[1]}) sets.")
    if bref_shape[1] != b_shape[1]:
        raise ValueError(f"Number of inputs differs between evaluated ({b_shape[1]}) and reference ({bref_shape[1]}) sets.")

    # Optional: Check if gx/gy lengths match variable counts if orientation is active
    # Note: Original assert_DEA1 didn't strictly enforce this, but it's good practice.
    # The code in the first DEA class *assumes* gx[j] and gy[k] exist.
    # Let's add these checks.
    if sum(gx_list) >= 1 and len(gx_list) != x.shape[1]:
         raise ValueError(f"Length of gx ({len(gx_list)}) must match the number of inputs ({x.shape[1]}) when input orientation is used.")
    if sum(gy_list) >= 1 and len(gy_list) != y.shape[1]:
         raise ValueError(f"Length of gy ({len(gy_list)}) must match the number of outputs ({y.shape[1]}) when output orientation is used.")
    if sum(gb_list) >= 1 and len(gb_list) != b.shape[1]:
         raise ValueError(f"Length of gb ({len(gb_list)}) must match the number of outputs ({b.shape[1]}) when unoutput orientation is used.")


    # 6. Return processed numpy arrays and the original index lists
    # We return the numpy arrays (y, x, yref, xref) as they are suitable for numerical
    # processing in Pyomo, and the shapes have been validated against the list shapes.
    # We return the processed lists for gy and gx.
    return y, x,b, yref, xref,bref, gy_list, gx_list,gb_list, evaluated_data_index, reference_data_index


def assert_valid_reference_data1(y, x, yref, xref):
    yref = trans_list(yref)
    xref = trans_list(xref)

    yref = to_2d_list(yref)
    xref = to_2d_list(xref)

    yref_shape = np.asarray(yref).shape
    xref_shape = np.asarray(xref).shape

    if yref_shape[0] != xref_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in xref and yref.")
    if yref_shape[1] != np.asarray(y).shape[1]:
        raise ValueError(
            "Number of outputs must be the same in y and yref.")
    if xref_shape[1] != np.asarray(x).shape[1]:
        raise ValueError(
            "Number of inputs must be the same in x and xref.")
    return yref, xref

def assert_DEAweakref(y, x, b, yref, xref, bref):
    yref, xref = assert_valid_reference_data(y, x, yref, xref)

    if type(b) == type(None):
        return yref, xref, None

    bref = to_2d_list(bref)
    bref_shape = np.asarray(bref).shape

    if bref_shape[0] != np.asarray(yref).shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in yref and bref.")
    if bref_shape[1] != np.asarray(b).shape[1]:
        raise ValueError(
            "Number of undesirable outputs must be the same in b and bref.")

    return yref, xref, bref

def assert_DDFref(y, x, yref, xref):
    yref, xref = assert_valid_reference_data(y, x, yref, xref)
    return yref, xref

def assert_DDFweakref(y, x, b, yref, xref, bref):
    yref, xref = assert_valid_reference_data(y, x, yref, xref)

    if type(b) == type(None):
        return yref, xref, None

    bref = to_2d_list(bref)
    bref_shape = np.asarray(bref).shape

    if bref_shape[0] != np.asarray(yref).shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in yref and bref.")
    if bref_shape[1] != np.asarray(b).shape[1]:
        raise ValueError(
            "Number of undesirable outputs must be the same in b and bref.")

    return yref, xref, bref

def assert_CNLSSD(data, sent, z, gy=[1], gx=[0]):
    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')

    # try:
    #     unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')
    # except:
    #     outputvars = sent.split('=')[1].strip(' ').split(' ')
    #     unoutputvars = None
    zvars = z.strip(' ').split(' ') if type(z) != type(None) else None

    x = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in outputvars])
    # if unoutputvars != None:
    #     b = np.column_stack(
    #         [np.asanyarray(data[selected]).T for selected in unoutputvars])
    if zvars != None:
        z = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in zvars])
    # print(z)
    y, x, z, gy, gx, basexy = assert_CNLSSD1(y, x, z, gy, gx)


    return y, x, z, gy, gx, basexy

def assert_CNLSSD1(y, x, z=None, gy=[1], gx=[1]):

    y = trans_list(y)
    x = trans_list(x)

    y = to_2d_list(y)
    x = to_2d_list(x)
    # print(x,"#############")

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    # print(gx,"#############")
    if sum(gx)>=1 and sum(gy)>=1:
        raise ValueError(
            "gy and gx can not be bigger than 1 together.")

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    basexy = []
    if sum(gx) >= 1:
        basexy = [sublist[i] for sublist in x for i in range(len(gx)) if gx[i] == 1]

    elif sum(gy) >= 1:
        basexy = [sublist[i] for sublist in y for i in range(len(gy)) if gy[i] == 1]

    # print(basexy,"xxxxxxxx")
    x = [
        [
            sublist[i] / sublist[i] if gx[i] == 1 else sublist[i]  # 根据 gx 的值决定是否进行除法
            for i in range(len(sublist))
        ]
        for sublist in x
    ]

    y = [
        [
            sublist[i] / sublist[i] if gy[i] == 1 else sublist[i]  # 根据 gy 的值决定是否进行除法
            for i in range(len(sublist))
        ]
        for sublist in y
    ]
    # print(x,"#############")

    # if sum(gx) >= 1:
    #     gx = [1 for _ in gx]
    # elif sum(gy) >= 1:
    #     gy = [1 for _ in gy]
    # print(gx,"#############")

    return y, x, z, gy, gx, basexy

def assert_CNLSSDFweak(data, sent, z, gy=[1], gx=[0], gb=[0]):
    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')

    zvars = z.strip(' ').split(' ') if type(z) != type(None) else None

    x = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in outputvars])
    if unoutputvars != None:
        b = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in unoutputvars])
    if zvars != None:
        z = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in zvars])
    # print(z)
    y, x, b, z, gy, gx, gb, basexy = assert_CNLSSDFweak1(y, x, b, z, gy, gx, gb)
    return y, x, b, z, gy, gx, gb, basexy


def assert_CNLSSDFweak1(y, x, b, z=None, gy=[1], gx=[1], gb=[1]):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)

    # print(x,"#############")

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    gb = to_1d_list(gb)
    # print(gx,"#############")
    if (sum(gx)>=1 and sum(gy)>=1) or (sum(gx)>=1 and sum(gb)>=1) or (sum(gy)>=1 and sum(gb)>=1):
        raise ValueError(
            "gy, gx and gb can not be bigger than 1 together.")

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if y_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in b and y.")

    if x_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in b and x.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")
    if b_shape[1] != len(gb):
        raise ValueError("Number of inputs must be the same in b and gb.")
    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    basexyb = []
    if sum(gx) >= 1:
        basexyb = [sublist[i] for sublist in x for i in range(len(gx)) if gx[i] == 1]

    elif sum(gy) >= 1:
        basexyb = [sublist[i] for sublist in y for i in range(len(gy)) if gy[i] == 1]

    elif sum(gb) >= 1:
        basexyb = [sublist[i] for sublist in b for i in range(len(gb)) if gb[i] == 1]
    # print(basexyb,"xxxxxxxx")
    x = [
        [
            sublist[i] / sublist[i] if gx[i] == 1 else sublist[i]  # 根据 gx 的值决定是否进行除法
            for i in range(len(sublist))
        ]
        for sublist in x
    ]

    y = [
        [
            sublist[i] / sublist[i] if gy[i] == 1 else sublist[i]  # 根据 gy 的值决定是否进行除法
            for i in range(len(sublist))
        ]
        for sublist in y
    ]
    b = [
        [
            sublist[i] / sublist[i] if gb[i] == 1 else sublist[i]  # 根据 gb 的值决定是否进行除法
            for i in range(len(sublist))
        ]
        for sublist in b
    ]
    # print(x,"#############")
    # print(y,"#############")
    # print(b,"#############")

    # if sum(gx) >= 1:
    #     gx = [1 for _ in gx]
    # elif sum(gy) >= 1:
    #     gy = [1 for _ in gy]
    # print(gx,"#############")

    return y, x, b, z, gy, gx, gb, basexyb

def assert_CNLSSDFweakmeta(data, sent, z, gy=[1], gx=[0], gb=[0]):
    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')

    zvars = z.strip(' ').split(' ') if type(z) != type(None) else None

    x = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in outputvars])
    if unoutputvars != None:
        b = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in unoutputvars])
    if zvars != None:
        z = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in zvars])
    # print(z)
    y, x, b, z, gy, gx, gb, basexy = assert_CNLSSDFweak1(y, x, b, z, gy, gx, gb)
    return y, x, b, z, gy, gx, gb, basexy






def assert_DDF(data, sent, gy, gx,baseindex, refindex):
    if type(baseindex) != type(None):
        varname1 = baseindex.split('=')[0]
        print(baseindex)
        varvalue1 = ast.literal_eval(baseindex.split('=')[1])
        data_base= data.loc[data[varname1].isin(varvalue1)]
    else:
        data_base= data
    data_index = data_base.index

    if type(refindex) != type(None):
        varname = refindex.split('=')[0]
        varvalue = ast.literal_eval(refindex.split('=')[1])

        data_ref = data.loc[data[varname].isin(varvalue)]
    else:
        data_ref = data

    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')

    x = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in outputvars])
    xref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in inputvars])
    yref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in outputvars])

    y, x, yref, xref, gy, gx = assert_DDF1(y, x, yref, xref, gy, gx)

    return y, x, yref, xref, gy, gx,data_index

def assert_DDF1(y, x, yref, xref, gy, gx):
    y = trans_list(y)
    x = trans_list(x)

    y = to_2d_list(y)
    x = to_2d_list(x)
    # print(x,"#############")

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    # print(gx,"#############")

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    yref = trans_list(yref)
    xref = trans_list(xref)

    yref = to_2d_list(yref)
    xref = to_2d_list(xref)

    yref_shape = np.asarray(yref).shape
    xref_shape = np.asarray(xref).shape

    if yref_shape[0] != xref_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in xref and yref.")
    if yref_shape[1] != np.asarray(y).shape[1]:
        raise ValueError(
            "Number of outputs must be the same in y and yref.")
    if xref_shape[1] != np.asarray(x).shape[1]:
        raise ValueError(
            "Number of inputs must be the same in x and xref.")

    return y, x, yref, xref, gy, gx




def assert_DDF_with_indices(data, sent, gy, gx, baseindex, refindex):
    """
    Validates DEA input data and returns processed numpy arrays and original indices.

    Args:
        data (pandas.DataFrame): The input data containing all variables and DMUs.
        sent (str): String specifying input and output variables. e.g., "K L = Y"
        gy (list): Output distance vector.
        gx (list): Input distance vector.
        baseindex (str, optional): Filter for evaluated DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).
        refindex (str, optional): Filter for reference DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).

    Returns:
        tuple: (y, x, yref, xref, gy, gx, evaluated_data_index, reference_data_index)
               y, x: numpy arrays for evaluated outputs and inputs.
               yref, xref: numpy arrays for reference outputs and inputs.
               gy, gx: Processed gy and gx lists.
               evaluated_data_index: List of original index labels for evaluated DMUs.
               reference_data_index: List of original index labels for reference DMUs.

    Raises:
        ValueError: If data format or dimensions are inconsistent.
    """
    # 1. Filter data based on baseindex to get evaluated DMUs
    if baseindex is not None:
        parts = baseindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'baseindex' format. Expected 'varname=value'")
        varname1 = parts[0].strip()
        try:
            varvalue1 = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'baseindex': {parts[1].strip()}")

        # Ensure varvalue1 is iterable for isin
        if not isinstance(varvalue1, (list, tuple, set)):
             varvalue1 = [varvalue1]

        if varname1 not in data.columns:
             raise ValueError(f"Variable '{varname1}' specified in 'baseindex' not found in data columns.")

        # Use .loc for filtering and .copy() to avoid potential SettingWithCopyWarning
        data_base = data.loc[data[varname1].isin(varvalue1)].copy()
    else:
        data_base = data.copy() # Use .copy() for default case too

    if data_base.empty:
         raise ValueError("Filtering with 'baseindex' resulted in an empty dataset for evaluated DMUs.")

    # Store evaluated DMU indices (original index labels)
    evaluated_data_index = data_base.index.tolist()


    # 2. Filter data based on refindex to get reference DMUs
    if refindex is not None:
        parts = refindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'refindex' format. Expected 'varname=value'")
        varname = parts[0].strip()
        try:
            varvalue = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'refindex': {parts[1].strip()}")

        # Ensure varvalue is iterable for isin
        if not isinstance(varvalue, (list, tuple, set)):
             varvalue = [varvalue]

        if varname not in data.columns:
             raise ValueError(f"Variable '{varname}' specified in 'refindex' not found in data columns.")

        # Use .loc for filtering and .copy()
        data_ref = data.loc[data[varname].isin(varvalue)].copy()
    else:
        data_ref = data.copy() # Use .copy() for default case too

    if data_ref.empty:
         raise ValueError("Filtering with 'refindex' resulted in an empty dataset for reference DMUs.")

    # Store reference DMU indices (original index labels)
    reference_data_index = data_ref.index.tolist()


    # 3. Parse sent string to identify input and output variables
    parts = sent.split('=')
    if len(parts) != 2:
        raise ValueError("Invalid 'sent' format. Expected 'inputs = outputs[:unwanted_outputs]'")

    input_part = parts[0].strip()
    # Split output part by ':', take the first part (outputs)
    output_part = parts[1].split(':')[0].strip()

    inputvars = input_part.split(' ')
    outputvars = output_part.split(' ')

    # Basic check if specified variables exist in the original data columns
    all_vars = inputvars + outputvars
    missing_vars = [v for v in all_vars if v not in data.columns]
    if missing_vars:
        raise ValueError(f"Specified variables not found in data columns: {missing_vars}")

    # 4. Extract data into numpy arrays
    # Use .values to get underlying numpy array from pandas Series/DataFrame column
    try:
        x = np.column_stack([data_base[selected].values for selected in inputvars])
        y = np.column_stack([data_base[selected].values for selected in outputvars])
        xref = np.column_stack([data_ref[selected].values for selected in inputvars])
        yref = np.column_stack([data_ref[selected].values for selected in outputvars])
    except KeyError as e:
         # This should ideally be caught by the check above, but as a safeguard
         raise ValueError(f"Error extracting data for variables: {e}")
    except Exception as e:
         # Catch other potential errors during numpy array creation
         raise ValueError(f"An error occurred during data extraction: {e}")


    # 5. Process gy, gx and perform validation checks (replicating assert_DEA1 logic)
    # Apply list transformations as in assert_DEA1 before checks
    # Ensure tools module and these functions are available
    try:
        y_list = to_2d_list(trans_list(y))
        x_list = to_2d_list(trans_list(x))
        yref_list = to_2d_list(trans_list(yref))
        xref_list = to_2d_list(trans_list(xref))
        gy_list = to_1d_list(gy)
        gx_list = to_1d_list(gx)
    except NameError:
        raise NameError("Helper functions (trans_list, to_2d_list, to_1d_list) from tools module are required but not found.")
    except Exception as e:
        raise RuntimeError(f"Error during data list transformation: {e}")


    # Now perform the checks using the list representations (as in original assert_DEA1)
    y_shape = np.asarray(y_list).shape # Use np.asarray to get shape like assert_DEA1
    x_shape = np.asarray(x_list).shape
    yref_shape = np.asarray(yref_list).shape
    xref_shape = np.asarray(xref_list).shape


    # Check number of DMUs match for evaluated set
    if y_shape[0] != x_shape[0]:
        raise ValueError(f"Number of evaluated DMUs differs between outputs ({y_shape[0]}) and inputs ({x_shape[0]}).")

    # Check number of DMUs match for reference set
    if yref_shape[0] != xref_shape[0]:
        raise ValueError(f"Number of reference DMUs differs between outputs ({yref_shape[0]}) and inputs ({xref_shape[0]}).")

    # Check number of variables match between evaluated and reference sets
    if yref_shape[1] != y_shape[1]:
        raise ValueError(f"Number of outputs differs between evaluated ({y_shape[1]}) and reference ({yref_shape[1]}) sets.")
    if xref_shape[1] != x_shape[1]:
        raise ValueError(f"Number of inputs differs between evaluated ({x_shape[1]}) and reference ({xref_shape[1]}) sets.")

    # Optional: Check if gx/gy lengths match variable counts if orientation is active
    # Note: Original assert_DEA1 didn't strictly enforce this, but it's good practice.
    # The code in the first DEA class *assumes* gx[j] and gy[k] exist.
    # Let's add these checks.
    if sum(gx_list) >= 1 and len(gx_list) != x.shape[1]:
         raise ValueError(f"Length of gx ({len(gx_list)}) must match the number of inputs ({x.shape[1]}) when input orientation is used.")
    if sum(gy_list) >= 1 and len(gy_list) != y.shape[1]:
         raise ValueError(f"Length of gy ({len(gy_list)}) must match the number of outputs ({y.shape[1]}) when output orientation is used.")


    # 6. Return processed numpy arrays and the original index lists
    # We return the numpy arrays (y, x, yref, xref) as they are suitable for numerical
    # processing in Pyomo, and the shapes have been validated against the list shapes.
    # We return the processed lists for gy and gx.
    return y, x, yref, xref, gy_list, gx_list, evaluated_data_index, reference_data_index

# Note: Ensure the 'tools' module is correctly imported and contains
# the functions trans_list, to_2d_list, and to_1d_list with appropriate logic.
# If you don't have these specific functions, you might need to adapt
# the validation logic to work directly with numpy arrays or pandas objects,
# or implement equivalent list transformation functions.











def assert_DDFweak(data,sent, gy, gx, gb,baseindex,refindex):
    if type(baseindex) != type(None):
        varname1 = baseindex.split('=')[0]
        print(baseindex)
        varvalue1 = ast.literal_eval(baseindex.split('=')[1])
        data_base= data.loc[data[varname1].isin(varvalue1)]
    else:
        data_base= data

    if type(refindex) != type(None):
        varname = refindex.split('=')[0]
        varvalue = ast.literal_eval(refindex.split('=')[1])

        data_ref = data.loc[data[varname].isin(varvalue)]
    else:
        data_ref = data

    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')

    x = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in outputvars])
    b = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in unoutputvars])

    xref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in inputvars])
    yref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in outputvars])
    bref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in unoutputvars])

    y, x, b,  gy, gx, gb, yref, xref, bref = assert_DDFweak1(y, x, b, gy, gx, gb, yref, xref, bref)

    return y, x, b,  gy, gx, gb, yref, xref, bref



def assert_DDFweak1(y, x, b, gy, gx, gb, yref, xref, bref):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    gb = to_1d_list(gb)
    print(gx,"#############")
    if (sum(gx)>=1 and sum(gy)>=1) or (sum(gx)>=1 and sum(gb)>=1) or (sum(gy)>=1 and sum(gb)>=1):
        raise ValueError(
            "gy, gx and gb can not be bigger than 1 together.")


    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")
    if y_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in b and y.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if b_shape[1] != len(gb):
        raise ValueError("Number of inputs must be the same in b and gb.")

    yref = trans_list(yref)
    xref = trans_list(xref)
    bref = trans_list(bref)

    yref = to_2d_list(yref)
    xref = to_2d_list(xref)
    bref = to_2d_list(bref)

    yref_shape = np.asarray(yref).shape
    xref_shape = np.asarray(xref).shape
    bref_shape = np.asarray(bref).shape

    if yref_shape[0] != xref_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in xref and yref.")
    if yref_shape[1] != np.asarray(y).shape[1]:
        raise ValueError(
            "Number of outputs must be the same in y and yref.")
    if xref_shape[1] != np.asarray(x).shape[1]:
        raise ValueError(
            "Number of inputs must be the same in x and xref.")

    return y, x, b,gy, gx, gb, yref, xref, bref




def assert_NDDFweak(data,sent, gy, gx, gb,baseindex,refindex):
    if type(baseindex) != type(None):
        varname1 = baseindex.split('=')[0]
        print(baseindex)
        varvalue1 = ast.literal_eval(baseindex.split('=')[1])
        data_base= data.loc[data[varname1].isin(varvalue1)]
    else:
        data_base= data

    if type(refindex) != type(None):
        varname = refindex.split('=')[0]
        varvalue = ast.literal_eval(refindex.split('=')[1])

        data_ref = data.loc[data[varname].isin(varvalue)]
    else:
        data_ref = data

    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')

    x = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in outputvars])
    b = np.column_stack(
        [np.asanyarray(data_base[selected]).T for selected in unoutputvars])

    xref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in inputvars])
    yref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in outputvars])
    bref = np.column_stack(
        [np.asanyarray(data_ref[selected]).T for selected in unoutputvars])

    y, x, b,  gy, gx, gb, yref, xref, bref = assert_NDDFweak1(y, x, b, gy, gx, gb, yref, xref, bref)

    return y, x, b,  gy, gx, gb, yref, xref, bref



def assert_NDDFweak1(y, x, b, gy, gx, gb, yref, xref, bref):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    gb = to_1d_list(gb)
    print(gx,"#############")
    if (sum(gx)>=1 and sum(gy)>=1) or (sum(gx)>=1 and sum(gb)>=1) or (sum(gy)>=1 and sum(gb)>=1):
        raise ValueError(
            "gy, gx and gb can not be bigger than 1 together.")


    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")
    if y_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in b and y.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if b_shape[1] != len(gb):
        raise ValueError("Number of inputs must be the same in b and gb.")

    yref = trans_list(yref)
    xref = trans_list(xref)
    bref = trans_list(bref)

    yref = to_2d_list(yref)
    xref = to_2d_list(xref)
    bref = to_2d_list(bref)

    yref_shape = np.asarray(yref).shape
    xref_shape = np.asarray(xref).shape
    bref_shape = np.asarray(bref).shape

    if yref_shape[0] != xref_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in xref and yref.")
    if yref_shape[1] != np.asarray(y).shape[1]:
        raise ValueError(
            "Number of outputs must be the same in y and yref.")
    if xref_shape[1] != np.asarray(x).shape[1]:
        raise ValueError(
            "Number of inputs must be the same in x and xref.")

    return y, x, b,gy, gx, gb, yref, xref, bref






def assert_NDDFweak_with_indices(data, sent, gy, gx,gb, baseindex, refindex):
    """
    Validates DEA input data and returns processed numpy arrays and original indices.

    Args:
        data (pandas.DataFrame): The input data containing all variables and DMUs.
        sent (str): String specifying input and output variables. e.g., "K L = Y:CO2"
        gy (list): Output distance vector.
        gx (list): Input distance vector.
        gb (list): Input distance vector.
        baseindex (str, optional): Filter for evaluated DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).
        refindex (str, optional): Filter for reference DMUs. e.g., "Year=[2010]". Defaults to None (all DMUs).

    Returns:
        tuple: (y, x, yref, xref, gy, gx, evaluated_data_index, reference_data_index)
               y, x: numpy arrays for evaluated outputs and inputs.
               yref, xref: numpy arrays for reference outputs and inputs.
               gy, gx: Processed gy and gx lists.
               evaluated_data_index: List of original index labels for evaluated DMUs.
               reference_data_index: List of original index labels for reference DMUs.

    Raises:
        ValueError: If data format or dimensions are inconsistent.
    """
    # 1. Filter data based on baseindex to get evaluated DMUs
    if baseindex is not None:
        parts = baseindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'baseindex' format. Expected 'varname=value'")
        varname1 = parts[0].strip()
        try:
            varvalue1 = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'baseindex': {parts[1].strip()}")

        # Ensure varvalue1 is iterable for isin
        if not isinstance(varvalue1, (list, tuple, set)):
             varvalue1 = [varvalue1]

        if varname1 not in data.columns:
             raise ValueError(f"Variable '{varname1}' specified in 'baseindex' not found in data columns.")

        # Use .loc for filtering and .copy() to avoid potential SettingWithCopyWarning
        data_base = data.loc[data[varname1].isin(varvalue1)].copy()
    else:
        data_base = data.copy() # Use .copy() for default case too

    if data_base.empty:
         raise ValueError("Filtering with 'baseindex' resulted in an empty dataset for evaluated DMUs.")

    # Store evaluated DMU indices (original index labels)
    evaluated_data_index = data_base.index.tolist()


    # 2. Filter data based on refindex to get reference DMUs
    if refindex is not None:
        parts = refindex.split('=')
        if len(parts) != 2:
            raise ValueError("Invalid 'refindex' format. Expected 'varname=value'")
        varname = parts[0].strip()
        try:
            varvalue = ast.literal_eval(parts[1].strip())
        except (ValueError, SyntaxError):
             raise ValueError(f"Invalid value format in 'refindex': {parts[1].strip()}")

        # Ensure varvalue is iterable for isin
        if not isinstance(varvalue, (list, tuple, set)):
             varvalue = [varvalue]

        if varname not in data.columns:
             raise ValueError(f"Variable '{varname}' specified in 'refindex' not found in data columns.")

        # Use .loc for filtering and .copy()
        data_ref = data.loc[data[varname].isin(varvalue)].copy()
    else:
        data_ref = data.copy() # Use .copy() for default case too

    if data_ref.empty:
         raise ValueError("Filtering with 'refindex' resulted in an empty dataset for reference DMUs.")

    # Store reference DMU indices (original index labels)
    reference_data_index = data_ref.index.tolist()


    # 3. Parse sent string to identify input and output variables
    parts = sent.split('=')
    if len(parts) != 2:
        raise ValueError("Invalid 'sent' format. Expected 'inputs = outputs[:unwanted_outputs]'")

    input_part = parts[0].strip()
    # Split output part by ':', take the first part (outputs)
    output_part = parts[1].strip().split(':')[0].strip()
    unoutput_part = parts[1].strip().split(':')[1].strip()

    inputvars = input_part.split(' ')
    outputvars = output_part.split(' ')
    unoutputvars = unoutput_part.split(' ')

    # Basic check if specified variables exist in the original data columns
    all_vars = inputvars + outputvars + unoutputvars
    missing_vars = [v for v in all_vars if v not in data.columns]
    if missing_vars:
        raise ValueError(f"Specified variables not found in data columns: {missing_vars}")

    # 4. Extract data into numpy arrays
    # Use .values to get underlying numpy array from pandas Series/DataFrame column
    try:
        x = np.column_stack([data_base[selected].values for selected in inputvars])
        y = np.column_stack([data_base[selected].values for selected in outputvars])
        b = np.column_stack([data_base[selected].values for selected in unoutputvars])
        xref = np.column_stack([data_ref[selected].values for selected in inputvars])
        yref = np.column_stack([data_ref[selected].values for selected in outputvars])
        bref = np.column_stack([data_ref[selected].values for selected in unoutputvars])
    except KeyError as e:
         # This should ideally be caught by the check above, but as a safeguard
         raise ValueError(f"Error extracting data for variables: {e}")
    except Exception as e:
         # Catch other potential errors during numpy array creation
         raise ValueError(f"An error occurred during data extraction: {e}")


    # 5. Process gy, gx and perform validation checks (replicating assert_DEA1 logic)
    # Apply list transformations as in assert_DEA1 before checks
    # Ensure tools module and these functions are available
    try:
        y_list = to_2d_list(trans_list(y))
        x_list = to_2d_list(trans_list(x))
        b_list = to_2d_list(trans_list(b))
        yref_list = to_2d_list(trans_list(yref))
        xref_list = to_2d_list(trans_list(xref))
        bref_list = to_2d_list(trans_list(bref))
        gy_list = to_1d_list(gy)
        gx_list = to_1d_list(gx)
        gb_list = to_1d_list(gb)

    except NameError:
        raise NameError("Helper functions (trans_list, to_2d_list, to_1d_list) from tools module are required but not found.")
    except Exception as e:
        raise RuntimeError(f"Error during data list transformation: {e}")


    # Now perform the checks using the list representations (as in original assert_DEA1)
    y_shape = np.asarray(y_list).shape # Use np.asarray to get shape like assert_DEA1
    x_shape = np.asarray(x_list).shape
    b_shape = np.asarray(b_list).shape
    yref_shape = np.asarray(yref_list).shape
    xref_shape = np.asarray(xref_list).shape
    bref_shape = np.asarray(bref_list).shape

    # Check number of DMUs match for evaluated set
    if y_shape[0] != x_shape[0]:
        raise ValueError(f"Number of evaluated DMUs differs between outputs ({y_shape[0]}) and inputs ({x_shape[0]}).")
    if y_shape[0] != b_shape[0]:
        raise ValueError(f"Number of evaluated DMUs differs between outputs ({y_shape[0]}) and unoutputs ({b_shape[0]}).")

    # Check number of DMUs match for reference set
    if yref_shape[0] != xref_shape[0]:
        raise ValueError(f"Number of reference DMUs differs between outputs ({yref_shape[0]}) and inputs ({xref_shape[0]}).")
    if yref_shape[0] != bref_shape[0]:
        raise ValueError(f"Number of reference DMUs differs between outputs ({yref_shape[0]}) and unoutputs ({bref_shape[0]}).")

    # Check number of variables match between evaluated and reference sets
    if yref_shape[1] != y_shape[1]:
        raise ValueError(f"Number of outputs differs between evaluated ({y_shape[1]}) and reference ({yref_shape[1]}) sets.")
    if xref_shape[1] != x_shape[1]:
        raise ValueError(f"Number of inputs differs between evaluated ({x_shape[1]}) and reference ({xref_shape[1]}) sets.")
    if bref_shape[1] != b_shape[1]:
        raise ValueError(f"Number of inputs differs between evaluated ({b_shape[1]}) and reference ({bref_shape[1]}) sets.")

    # Optional: Check if gx/gy lengths match variable counts if orientation is active
    # Note: Original assert_DEA1 didn't strictly enforce this, but it's good practice.
    # The code in the first DEA class *assumes* gx[j] and gy[k] exist.
    # Let's add these checks.
    if sum(gx_list) >= 1 and len(gx_list) != x.shape[1]:
         raise ValueError(f"Length of gx ({len(gx_list)}) must match the number of inputs ({x.shape[1]}) when input orientation is used.")
    if sum(gy_list) >= 1 and len(gy_list) != y.shape[1]:
         raise ValueError(f"Length of gy ({len(gy_list)}) must match the number of outputs ({y.shape[1]}) when output orientation is used.")
    if sum(gb_list) >= 1 and len(gb_list) != b.shape[1]:
         raise ValueError(f"Length of gb ({len(gb_list)}) must match the number of outputs ({b.shape[1]}) when unoutput orientation is used.")


    # 6. Calculate the weight vector 'w'
    # Count oriented dimensions in each group
    oriented_inputs = sum(1 for gxi in gx_list if gxi == 1)
    oriented_good_outputs = sum(1 for gyi in gy_list if gyi == 1)
    oriented_bad_outputs = sum(1 for gbi in gb_list if gbi == 1)

    # Determine which groups are oriented
    is_input_group_oriented = oriented_inputs > 0
    is_good_output_group_oriented = oriented_good_outputs > 0
    is_bad_output_group_oriented = oriented_bad_outputs > 0

    # Count oriented groups
    num_oriented_groups = int(is_input_group_oriented) + int(is_good_output_group_oriented) + int(is_bad_output_group_oriented)
    
    # Validate dimensions of direction vectors against variable counts
    num_inputs = len(inputvars)
    num_good_outputs = len(outputvars)
    num_bad_outputs = len(unoutputvars)

    # Initialize weight vector w
    total_vars = num_inputs + num_good_outputs + num_bad_outputs
    w = [0.0] * total_vars

    if num_oriented_groups > 0:
        # Calculate total weight assigned to each oriented group
        group_weight = 1.0 / num_oriented_groups

        # Distribute weight within each oriented group
        if is_input_group_oriented:
            weight_per_oriented_input = group_weight / oriented_inputs # oriented_inputs > 0 here
            for i in range(num_inputs):
                if gx_list[i] == 1:
                    w[i] = weight_per_oriented_input

        if is_good_output_group_oriented:
            weight_per_oriented_good_output = group_weight / oriented_good_outputs # oriented_good_outputs > 0 here
            for j in range(num_good_outputs):
                 if gy_list[j] == 1:
                    w[num_inputs + j] = weight_per_oriented_good_output

        if is_bad_output_group_oriented:
            weight_per_oriented_bad_output = group_weight / oriented_bad_outputs # oriented_bad_outputs > 0 here
            for k in range(num_bad_outputs):
                 if gb_list[k] == 1:
                    w[num_inputs + num_good_outputs + k] = weight_per_oriented_bad_output

    # Now, slice w into wx, wy, wb based on the number of variables in each group
    wx = w[:num_inputs]
    wy = w[num_inputs : num_inputs + num_good_outputs]
    wb = w[num_inputs + num_good_outputs :] # Slices to the end of w

    # 7. Return processed numpy arrays and the original index lists
    # We return the numpy arrays (y, x, yref, xref) as they are suitable for numerical
    # processing in Pyomo, and the shapes have been validated against the list shapes.
    # We return the processed lists for gy and gx.
    return y, x,b, outputvars,inputvars,unoutputvars,yref, xref,bref, gy_list, gx_list,gb_list,wx,wy,wb, evaluated_data_index, reference_data_index


def assert_optimized(optimization_status):
    if optimization_status == 0:
        raise Exception(
            "Model isn't optimized. Use optimize() method to estimate the model.")


def assert_contextual_variable(z):
    if type(z) == type(None):
        raise Exception(
            "Estimated coefficient (lamda) cannot be retrieved due to no contextual variable (z variable) included in the model.")

def assert_desirable_output(y):
    if type(y) == type(None):
        raise Exception(
            "Estimated coefficient (gamma) cannot be retrieved due to no desirable output (y variable) included in the model.")

def assert_undesirable_output(b):
    if type(b) == type(None):
        raise Exception(
            "Estimated coefficient (delta) cannot be retrieved due to no undesirable output (b variable) included in the model.")


def assert_various_return_to_scale(rts):
    if rts == RTS_CRS:
        raise Exception(
            "Estimated intercept (alpha) cannot be retrieved due to the constant returns-to-scale assumption.")


def assert_various_return_to_scale_alpha(rts):
    if rts == RTS_CRS:
        raise Exception(
            "Omega cannot be retrieved due to the constant returns-to-scale assumption.")


def assert_solver_available_locally(solver):
    if not SolverFactory(solver).available():
        raise ValueError("Solver {} is not available locally.".format(solver))


def assert_CNLSDDF(data, sent, z, gy=[1], gx=[0]):
    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')

    # try:
    #     unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')
    # except:
    #     outputvars = sent.split('=')[1].strip(' ').split(' ')
    #     unoutputvars = None
    zvars = z.strip(' ').split(' ') if type(z) != type(None) else None

    x = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in outputvars])
    # if unoutputvars != None:
    #     b = np.column_stack(
    #         [np.asanyarray(data[selected]).T for selected in unoutputvars])
    if zvars != None:
        z = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in zvars])

    y, x, z, gy, gx, basexy = assert_CNLSDDF1(y, x, z, gy, gx)
    return y, x, z, gy, gx, basexy

def assert_CNLSDDF1(y, x, z=None, gy=[1], gx=[1]):
    y = trans_list(y)
    x = trans_list(x)

    y = to_2d_list(y)
    x = to_2d_list(x)

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    if sum(gy) >= 1:
        print(y,"#########")
        print(x,"#########")
        # 找到第一个为 1 的索引
        index = gy.index(1)
        # 提取 aa 中对应索引的元素
        basexy = [sublist[index] for sublist in y]
        y = [[elem - basexy[i]*gy[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem + basexy[i]*gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]

        print(basexy)
        print(y,"#########")
        print(x,"#########")

    elif sum(gx) >= 1:
        print(y,"#########")
        print(x,"#########")
        index = gx.index(1)
        basexy = [-sublist[index] for sublist in x]
        print(basexy)

        y = [[elem - basexy[i]*gy[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem + basexy[i]*gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]

        print(y,"#########")
        print(x,"#########")
    else:
        raise ValueError(
            "gx and gy must either be 1")

    return y, x, z, gy, gx, basexy

def assert_CNLSDDFweak(data, sent, z, gy=[1], gx=[0], gb=[0]):
    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')

    zvars = z.strip(' ').split(' ') if type(z) != type(None) else None
    # print(zvars,"ssssssssssss")
    x = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in outputvars])
    if unoutputvars != None:
        b = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in unoutputvars])
    if zvars != None:
        z = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in zvars])
    # print(z)
    y, x, b, z, gy, gx, gb, basexy = assert_CNLSDDFweak1(y, x, b, z, gy, gx, gb)
    return y, x, b, z, gy, gx, gb, basexy

def assert_CNLSDDFweak1(y, x, b, z=None, gy=[1], gx=[1], gb=[1]):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    gb = to_1d_list(gb)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")
    if y_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in b and y.")
    if x_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and b.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if b_shape[1] != len(gb):
        raise ValueError("Number of inputs must be the same in b and gb.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    if sum(gy) >= 1:
        print("y#########",y)
        print("x#########",x)
        print("b#########",b)
        # 找到第一个为 1 的索引
        index = gy.index(1)
        # 提取 aa 中对应索引的元素
        basexyb = [sublist[index] for sublist in y]
        y = [[elem - basexyb[i]*gy[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem + basexyb[i]*gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
        b = [[elem + basexyb[i]*gb[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]

        print("#########",basexyb)
        print("y#########",y)
        print("x#########",x)
        print("b#########",b)

    elif sum(gx) >= 1:
        # print(y,"y#########")
        # print(x,"x#########")
        # print(b,"b#########")
        index = gx.index(1)
        basexyb = [-sublist[index] for sublist in x]
        # print(basexyb)

        y = [[elem - basexyb[i]*gy[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem + basexyb[i]*gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
        b = [[elem + basexyb[i]*gb[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]

        # print('basexybq',basexyb)
        # print(y,"y#########")
        # print(x,"x#########")
        # print(b,"b#########")

    elif sum(gb) >= 1:
        print(y, "y#########")
        print(x, "x#########")
        print(b, "b#########")
        index = gb.index(1)
        basexyb = [-sublist[index] for sublist in b]
        print(basexyb)

        y = [[elem - basexyb[i] * gy[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem + basexyb[i] * gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
        b = [[elem + basexyb[i] * gb[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]

        print(basexyb)
        print(y, "y#########")
        print(x, "x#########")
        print(b, "b#########")
    else:
        raise ValueError(
            "gx and gy and gb must either be 1")

    return y, x, b, z, gy, gx, gb, basexyb

def assert_CNLSDDFweakmeta(data, sent, z, gddf, gy=[1], gx=[0], gb=[0]):
    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')

    zvars = z.strip(' ').split(' ') if type(z) != type(None) else None

    x = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in inputvars])
    y = np.column_stack(
        [np.asanyarray(data[selected]).T for selected in outputvars])
    if unoutputvars != None:
        b = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in unoutputvars])
    if zvars != None:
        z = np.column_stack(
            [np.asanyarray(data[selected]).T for selected in zvars])
    # print(z)
    y, x, b, z, gy, gx, gb, basexy, basexy_old = assert_CNLSDDFweakmeta1(y, x, b, z, gddf, gy, gx, gb)
    return y, x, b, z, gy, gx, gb, basexy, basexy_old

def assert_CNLSDDFweakmeta1(y, x, b, z, gddf, gy=[1], gx=[1], gb=[1]):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)
    # gddf = to_2d_list(gddf)
    print("1dgddf",gddf)


    gy = to_1d_list(gy)
    gx = to_1d_list(gx)
    gb = to_1d_list(gb)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape
    gddf_shape = np.asarray(gddf).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")
    if y_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in b and y.")
    if x_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and b.")
    if y_shape[0] != gddf_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in y and gddf.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if b_shape[1] != len(gb):
        raise ValueError("Number of inputs must be the same in b and gb.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    # y2 = [[elem + gddf[i]*gy[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
    # x2 = [[elem - gddf[i]*gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
    # b2 = [[elem - gddf[i]*gb[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]
    # for i, row in enumerate(x2):
    #     for j, value in enumerate(row):
    #         if value < 0:
    #             print(f"x2[{i}][{j}] = {value} 小于0，将其替换为0")
    #             x2[i][j] = 1


    if sum(gy) >= 1:
        # print("y#########",y2)
        # print("x#########",x2)
        # print("b#########",b2)
        print("gddf#########",gddf)
        # 找到第一个为 1 的索引
        index = gy.index(1)
        # 提取 aa 中对应索引的元素
        basexyb = [sublist[index] for sublist in y]
        y = [[elem - basexyb[i] * gy[j]  for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem + basexyb[i] * gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
        b = [[elem + basexyb[i] * gb[j]  for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]
        # for i, row in enumerate(x):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             # print(f"x[{i}][{j}] = {value} 小于0，将其替换为0")
        #             x[i][j] = 1
        # for i, row in enumerate(y):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             # print(f"y[{i}][{j}] = {value} 小于0，将其替换为0")
        #             y[i][j] = 1
        # for i, row in enumerate(b):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             # print(f"b[{i}][{j}] = {value} 小于0，将其替换为0")
        #             b[i][j] = 1
        print("#########",basexyb)
        print("y#########",y)
        print("x#########",x)
        print("b#########",b)

    elif sum(gx) >= 1:

        index = gx.index(1)

        basexyb_old = [-sublist[index] for sublist in x]

        basexyb = [-sublist[index] for sublist in x]
        print(basexyb)

        # y = [[elem  + gddf[i]*gy[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        # x = [[elem  - gddf[i]*gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
        # b = [[elem  - gddf[i]*gb[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]

        y = [[elem  for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem  for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
        b = [[elem  for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]
        # for i, row in enumerate(x2):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             print(f"x2[{i}][{j}] = {value} 小于0，将其替换为0")
        #             x2[i][j] = 1
        # for i, row in enumerate(y2):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             print(f"y2[{i}][{j}] = {value} 小于0，将其替换为0")
        #             y2[i][j] = 1
        # for i, row in enumerate(b2):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             print(f"b2[{i}][{j}] = {value} 小于0，将其替换为0")
        #             b2[i][j] = 1


        print('basexyb',basexyb)
        print(y,"y#########")
        print(x,"x#########wwwwwwwwwwwwwwwwwwww")
        print(b,"b#########")

    elif sum(gb) >= 1:
        print(y, "y#########")
        print(x, "x#########")
        print(b, "b#########")
        index = gb.index(1)
        basexyb = [-sublist[index] for sublist in b]
        print(basexyb)

        y = [[elem - basexyb[i] * gy[j]  for j, elem in enumerate(sublist)] for i, sublist in enumerate(y)]
        x = [[elem + basexyb[i] * gx[j] for j, elem in enumerate(sublist)] for i, sublist in enumerate(x)]
        b = [[elem + basexyb[i] * gb[j]  for j, elem in enumerate(sublist)] for i, sublist in enumerate(b)]
        # for i, row in enumerate(x):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             # print(f"x[{i}][{j}] = {value} 小于0，将其替换为0")
        #             x[i][j] = 1
        # for i, row in enumerate(y):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             # print(f"y[{i}][{j}] = {value} 小于0，将其替换为0")
        #             y[i][j] = 1
        # for i, row in enumerate(b):
        #     for j, value in enumerate(row):
        #         if value < 0:
        #             # print(f"b[{i}][{j}] = {value} 小于0，将其替换为0")
        #             b[i][j] = 1
        print('basexyb',basexyb)
        print("y#########",y)
        print(x, "x#########")
        print(b, "b#########")
    else:
        raise ValueError(
            "gx and gy and gb must either be 1")

    return y, x, b, z, gy, gx, gb, basexyb,basexyb_old


def assert_valid_direciontal_data_with_z(y, x, b=None,z=None, gy=[1], gx=[1], gb=None):
    y = trans_list(y)
    x = trans_list(x)

    y = to_2d_list(y)
    x = to_2d_list(x)

    gy = to_1d_list(gy)
    gx = to_1d_list(gx)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if y_shape[1] != len(gy):
        raise ValueError("Number of outputs must be the same in y and gy.")

    if x_shape[1] != len(gx):
        raise ValueError("Number of inputs must be the same in x and gx.")

    if type(b) != type(None):
        b = trans_list(b)
        b = to_2d_list(b)
        gb = to_1d_list(gb)
        b_shape = np.asarray(b).shape
        if b_shape[0] != b_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and b.")
        if b_shape[1] != len(gb):
            raise ValueError(
                "Number of undesirable outputs must be the same in b and gb.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")
    return y, x, b,z, gy, gx, gb

def assert_valid_wp_data_x(y, x, b, z=None):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)
    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    if len(y_shape) == 2 and y_shape[1] != 1:
        raise ValueError(
            "The multidimensional output data is supported by direciontal based models.")

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if x_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and b.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    return y, x, b, z

def assert_valid_wp_data_b(y, x, b, z=None):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_2d_list(y)
    x = to_2d_list(x)
    b = to_1d_list(b)
    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    if len(y_shape) == 2 and y_shape[1] != 1:
        raise ValueError(
            "The multidimensional output data is supported by direciontal based models.")

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if x_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and b.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    return y, x, b, z

def assert_valid_wp_data(y, x, b, z=None):
    y = trans_list(y)
    x = trans_list(x)
    b = trans_list(b)

    y = to_1d_list(y)
    x = to_2d_list(x)
    b = to_2d_list(b)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape
    b_shape = np.asarray(b).shape

    if len(y_shape) == 2 and y_shape[1] != 1:
        raise ValueError(
            "The multidimensional output data is supported by direciontal based models.")

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if x_shape[0] != b_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and b.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    return y, x, b, z


def assert_valid_mupltiple_x_y_data(y, x, z=None):
    y = trans_list(y)
    x = trans_list(x)

    y = to_2d_list(y)
    x = to_2d_list(x)

    y_shape = np.asarray(y).shape
    x_shape = np.asarray(x).shape

    if y_shape[0] != x_shape[0]:
        raise ValueError(
            "Number of DMUs must be the same in x and y.")

    if type(z) != type(None):
        z = trans_list(z)
        z = to_2d_list(z)
        z_shape = np.asarray(z).shape
        if y_shape[0] != z_shape[0]:
            raise ValueError(
                "Number of DMUs must be the same in y and z.")

    return y, x, z

def assert_valid_yxbz_nog(sent,z):
    inputvars = sent.split('=')[0].strip(' ').split(' ')

    try:
        outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
        unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')
    except:
        outputvars = sent.split('=')[1].strip(' ').split(' ')
        unoutputvars = None
    zvars=z.strip(' ').split(' ') if type(z)!=type(None) else None

    return outputvars,inputvars,unoutputvars,zvars

def assert_valid_yxb(sent,gy,gx,gb):
    inputvars = sent.split('=')[0].strip(' ').split(' ')

    try:
        outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
        unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')
    except:
        outputvars = sent.split('=')[1].strip(' ').split(' ')
        unoutputvars = None

    if len(outputvars) != gy.shape[1]:
        raise ValueError("Number of outputs must be the same in y and gy.")
    if len(inputvars) != gx.shape[1]:
        raise ValueError("Number of inputs must be the same in x and gx.")

    if type(gb) != type(None):
        if len(unoutputvars) != gb.shape[1]:
            raise ValueError(
                "Number of undesirable outputs must be the same in b and gb.")
        gb.columns = unoutputvars
    gy.columns,gx.columns ,=outputvars,inputvars,
    return outputvars,inputvars,unoutputvars,gy,gx,gb

def assert_valid_yxb2(baseindex,refindex,data,outputvars,inputvars,unoutputvars):

    if type(baseindex) != type(None):
        varname1 = baseindex.split('=')[0]
        varvalue1 = ast.literal_eval(baseindex.split('=')[1])
        y, x, b = data.loc[data[varname1].isin(varvalue1), outputvars
        ], data.loc[data[varname1].isin(varvalue1), inputvars
        ], data.loc[data[varname1].isin(varvalue1), unoutputvars
        ]

    else:
        y, x, b = data.loc[:, outputvars], data.loc[:, inputvars], data.loc[:, unoutputvars ]

    if type(refindex) != type(None):
        varname=refindex.split('=')[0]
        varvalue=ast.literal_eval(refindex.split('=')[1])

        yref, xref, bref = data.loc[data[varname].isin(varvalue), outputvars
                                            ], data.loc[data[varname].isin(varvalue), inputvars
                                            ], data.loc[data[varname].isin(varvalue), unoutputvars
                                            ]
    else:
        yref, xref, bref = data.loc[:, outputvars], data.loc[:, inputvars], data.loc[:, unoutputvars ]
    return y,x,b,yref,xref,bref




def assert_valid_yxbz2(baseindex,refindex,data,outputvars,inputvars,unoutputvars,zvars):


    if type(baseindex) != type(None):
        varname = baseindex.split('=')[0]
        yr = ast.literal_eval(baseindex.split('=')[1])
        y, x, b,z = data.loc[data[varname].isin(yr), outputvars], \
                    data.loc[data[varname].isin(yr), inputvars], \
                    data.loc[data[varname].isin(yr), unoutputvars], \
                    data.loc[data[varname].isin(yr), zvars] if type(zvars) != type(None) else None
        if type(refindex) != type(None):
            yrref = ast.literal_eval(refindex.split('=')[1])

            if len(set(yr) - set(yrref)) > 0:
                print("ssssssssssssss1111111")
                raise ValueError(
                    "You must specify basic data smaller than reference data.")
            else:
                print("ssssssssssssss22222222")
                yrref2 = list(set(yrref) - set(yr))
                try:
                    print(yrref2[0])
                    yref, xref, bref, zref = data.loc[data[varname].isin(yrref2), outputvars], \
                        data.loc[data[varname].isin(yrref2), inputvars], \
                        data.loc[data[varname].isin(yrref2), unoutputvars], \
                        data.loc[data[varname].isin(yrref2), zvars] if type(zvars) != type(None) else None
                except:
                    yref, xref, bref, zref = None, \
                        None, \
                        None, \
                        None
        elif type(refindex) == type(None):
            yrref = list(data[varname].unique())
            if len(set(yr) - set(yrref)) > 0:
                print("ssssssssssssss1111111")
                raise ValueError(
                    "You must specify basic data smaller than reference data.")
            else:

                print("ssssssssssssss22222222")
                yrref2 = list(set(yrref) - set(yr))
                try:
                    print(yrref2[0])
                    yref, xref, bref, zref = data.loc[data[varname].isin(yrref2), outputvars], \
                        data.loc[data[varname].isin(yrref2), inputvars], \
                        data.loc[data[varname].isin(yrref2), unoutputvars], \
                        data.loc[data[varname].isin(yrref2), zvars] if type(zvars) != type(None) else None
                except:
                    yref, xref, bref, zref = None, \
                        None, \
                        None, \
                        None

    else:
        y, x, b,z = data.loc[:, outputvars], data.loc[:, inputvars], data.loc[:, unoutputvars],\
                    data.loc[:, zvars] if type(zvars) != type(None) else None

        if type(refindex) != type(None):
            varname = refindex.split('=')[0]
            yrref = ast.literal_eval(refindex.split('=')[1])
            yr = list(data[varname].unique())
            if len(set(yr) - set(yrref)) > 0:
                print("ssssssssssssss1111111")
                raise ValueError(
                    "You must specify basic data smaller than reference data.")
            else:
                print("ssssssssssssss22222222")
                yrref2 = list(set(yrref) - set(yr))
                try:
                    print(yrref2[0])
                    yref, xref, bref, zref = data.loc[data[varname].isin(yrref2), outputvars], \
                        data.loc[data[varname].isin(yrref2), inputvars], \
                        data.loc[data[varname].isin(yrref2), unoutputvars], \
                        data.loc[data[varname].isin(yrref2), zvars] if type(zvars) != type(None) else None
                except:
                    yref, xref, bref, zref = None, \
                        None, \
                        None, \
                        None
        elif type(refindex) == type(None):
            yref, xref, bref, zref = None, \
                None, \
                None, \
                None


    if type(yref) != type(None):
        referenceflag = True
    else:
        referenceflag = False

    # print("1",y)
    # print("2",yref)
    # print("3",referenceflag)
    return y,x,b,z,yref,xref,bref,zref,referenceflag

def assert_valid_yxbz(sent,gy,gx,gb,z=None):
    inputvars = sent.split('=')[0].strip(' ').split(' ')

    try:
        outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
        unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')
    except:
        outputvars = sent.split('=')[1].strip(' ').split(' ')
        unoutputvars = None

    if type(z)!=type(None):
        zvars = z.strip(' ').split(" ")
    else:
        zvars = None
    if len(outputvars) !=  gy.shape[1]:
        raise ValueError("Number of outputs must be the same in y and gy.")
    if len(inputvars) != gx.shape[1]:
        raise ValueError("Number of inputs must be the same in x and gx.")

    if type(gb) != type(None):
        if len(unoutputvars) != gb.shape[1]:
            raise ValueError(
                "Number of undesirable outputs must be the same in b and gb.")
        gb.columns = unoutputvars
    gy.columns,gx.columns ,=outputvars,inputvars,
    return outputvars,inputvars,unoutputvars,zvars,gy,gx,gb





def assert_valid_yxb_drf(sent,fenmu,fenzi):
    inputvars = sent.split('=')[0].strip(' ').split(' ')
    outputvars = sent.split('=')[1].split(':')[0].strip(' ').split(' ')
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')
    vars=inputvars +outputvars+unoutputvars
    if fenmu not in vars:
        raise ValueError("fenmu must be in sent.")
    if fenzi not in vars:
        raise ValueError("fenzi must be in sent.")

    varslt = {"inputvars": inputvars,
              "outputvars": outputvars,
              "unoutputvars": unoutputvars,
              }
    obj_coeflt = {"xobj_coef": len(inputvars) * [0],
                  "yobj_coef": len(outputvars) * [0],
                  "bobj_coef": len(unoutputvars) * [0]
                  }

    rule4_coeflt = {"xrule4_coef": len(inputvars) * [0],
                    "yrule4_coef": len(outputvars) * [0],
                    "brule4_coef": len(unoutputvars) * [0]
                    }

    for i, j in enumerate(varslt["inputvars"]):
        if fenzi == j:
            obj_coeflt["xobj_coef"][i] = 1
        if fenmu == j:
            rule4_coeflt["xrule4_coef"][i] = 1

    for i, j in enumerate(varslt["outputvars"]):
        if fenzi == j:
            obj_coeflt["yobj_coef"][i] = 1
        if fenmu == j:
            rule4_coeflt["yrule4_coef"][i] = 1
    for i, j in enumerate(varslt["unoutputvars"]):
        if fenzi == j:
            obj_coeflt["bobj_coef"][i] = 1
        if fenmu == j:
            rule4_coeflt["brule4_coef"][i] = 1

    ## 判断分母是x，b or y，是x，b的，目标要加负号。
    if (fenmu in inputvars) or (fenmu in unoutputvars):
        neg_obj = True
    elif fenmu in outputvars:
        neg_obj = False

    return outputvars, inputvars, unoutputvars, obj_coeflt, rule4_coeflt,neg_obj


def split_MB(sent, sx, sy,level):
    inputvars = sent.split('=')[0].strip(' ')
    inputvars_np = inputvars.split('+')[0].strip(' ').split(' ')  ## 假设一定有不含污染的投入，为了简单点
    inputvars_p = inputvars.split('+')[1].strip(' ').split(' ')  ## 一定有含污染的投入

    outputvars = sent.split('=')[1].split(':')[0].strip(' ')
    try:  ## 期望产出中，给了加号
        outputvars_np = outputvars.split('+')[0].strip(' ').split(' ')
        outputvars_p = outputvars.split('+')[1].strip(' ').split(' ')
        if outputvars_np[0] == "":  ## 给了加号后，前面（含污染）是空的
            outputvars_np = None
        if outputvars_p[0] == "":  ## 给了加号后，后面（含污染）是空的
            outputvars_p = None

    except:  ## 期望产出中没有加号
        outputvars_np = outputvars.strip(' ').split(' ')  ## 默认所有都是不含污染
        if outputvars_np[0] == "":  ## 没有加号后，前面是空的，后面（含污染）是空的
            outputvars_np = None


        outputvars_p = None
    unoutputvars = sent.split('=')[1].split(':')[1].strip(' ').split(' ')  ## 一定有非期望产出

    if type(outputvars_np) == type(None):
        if type(outputvars_p) == type(None):
            n1, n2, n3, n4, n5 = len(inputvars_np), len(inputvars_p), 0, 0, len(unoutputvars)
        elif type(outputvars_p) != type(None):
            n1, n2, n3, n4, n5 = len(inputvars_np), len(inputvars_p), 0, len(outputvars_p), len(unoutputvars)

    elif type(outputvars_np) != type(None):
        if type(outputvars_p) == type(None):
            n1, n2, n3, n4, n5 = len(inputvars_np), len(inputvars_p), len(outputvars_np), 0, len(unoutputvars)
        elif type(outputvars_p) != type(None):
            n1, n2, n3, n4, n5 = len(inputvars_np), len(inputvars_p), len(outputvars_np), len(outputvars_p), len(
                unoutputvars)
    # print(np.array(sx).shape[0])

    if np.array(sx).shape[0] != n5:
        raise ValueError(
            "Number of lists in sx must be the same in length of b")
    # print(n1,np.array(sx)[0,0:n1],np.array(sx)[0,0:n1].all(0))

    if not np.array(sx)[0, n1:n1 + n2].any(0):
        raise ValueError(
            "Number of polluted input must be the same in the position of sx")

    if type(outputvars_np) != type(None):
        if type(outputvars_p) != type(None):
           if level >5:
               raise ValueError(
                   "There are input_np, input_p, output_np, output_p in your statement of sent, \n"
                   "so you can state at most 5 level in this model")

        elif type(outputvars_p) == type(None):
            if level > 4:
                raise ValueError(
                    "There are input_np, input_p, output_np in your statement of sent, \n"
                    "so you can state at most 4 level in this model")
    elif type(outputvars_np) == type(None):
        if type(outputvars_p) != type(None):
            if level > 4:
                raise ValueError(
                    "There are input_np, input_p, output_p in your statement of sent, \n"
                    "so you can state at most 4 level in this model")
        elif type(outputvars_p) == type(None):
            if level > 3:
                raise ValueError(
                    "There are input_np, input_p, output_p in your statement of sent, \n"
                    "so you can state at most 3 level in this model")
    return inputvars_np, inputvars_p, outputvars_np, outputvars_p, unoutputvars, sx, sy,level
