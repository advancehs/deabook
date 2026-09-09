# 本文件从教材原文提取；在项目根目录运行。

# Markdown line 51
from deabook import data_path
from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import pandas as pd
from deabook.example.book_cases.loaders import china
import numpy as np
import re
import warnings

warnings.filterwarnings("ignore")


def sbm(formula, dataframe, eval_query=None, ref_query=None, solver="highs"):
    """
    计算SBM（Slacks-Based Measure）效率值。
    
    基于Tone(2001)的SBM模型，通过Charnes-Cooper线性转换将非线性
    规划问题转化为线性规划问题求解。
    
    Parameters
    ----------
    formula : str
        公式字符串，格式为 "产出 = 投入1 投入2 ..."
        例如："Y = K L"
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框
    eval_query : str or None, default None
        筛选待评价决策单元的query字符串。None表示评价全部
    ref_query : str or None, default None
        筛选参照决策单元的query字符串。None表示使用全部
    solver : str, default "highs"
        线性规划求解器名称
    
    Returns
    -------
    pd.DataFrame
        包含效率值(obj)和各投入产出松弛变量的DataFrame
    
    Examples
    --------
    >>> import pandas as pd
from deabook.example.book_cases.loaders import china
    >>> ex4 = china(scaled=True)
    >>> result = sbm("Y = K L", ex4, "t==[1,2,3]", "t==[1,2,3]")
    """
    # === Step 1: 解析公式字符串 ===
    input_vars = formula.split("=")[1].strip(" ")
    xcol = re.compile(" +").sub(" ", input_vars).split(" ")
    output_vars = formula.split("=")[0].strip(" ")
    ycol = re.compile(" +").sub(" ", output_vars).split(" ")
    
    # === Step 2: 根据查询条件准备数据 ===
    if eval_query is None:
        index_eval = dataframe.index
    else:
        index_eval = dataframe.query(eval_query).index
    
    if ref_query is None:
        index_ref = dataframe.index
    else:
        index_ref = dataframe.query(ref_query).index
    
    data = dataframe.loc[index_eval, xcol + ycol]
    dataref = dataframe.loc[index_ref, xcol + ycol]
    
    xref = dataref.loc[:, xcol]   # 参照单元投入
    yref = dataref.loc[:, ycol]   # 参照单元产出
    
    # === Step 3: 数据验证 ===
    assert data[xcol + ycol].notna().all().all(), "数据包含缺失值"
    assert (data[xcol] > 0).all().all(), "投入数据必须严格为正"
    assert (data[ycol] > 0).all().all(), "产出数据必须严格为正"
    
    # 初始化结果存储字典
    slack_x_results = {}   # 投入松弛变量
    slack_y_results = {}   # 产出松弛变量
    obj_results = {}       # 目标函数值
    
    # === Step 4: 对每个决策单元构建并求解线性规划 ===
    for j in data.index:
        x = data.loc[j, xcol]   # 被评价单元投入
        y = data.loc[j, ycol]   # 被评价单元产出
        
        # 构建Pyomo模型
        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)          # 参照决策单元集合
        model.K = Set(initialize=range(len(xcol)))       # 投入变量集合
        model.L = Set(initialize=range(len(ycol)))       # 产出变量集合
        
        # 定义决策变量
        model.t = Var(domain=PositiveReals, doc="CC transformation variable")
        model.slack_x = Var(model.K, bounds=(0.0, None), doc="input slack")
        model.slack_y = Var(model.L, bounds=(0.0, None), doc="output slack")
        model.intensity = Var(model.I, bounds=(0.0, None), doc="intensity variables")
        
        # 目标函数：min t - (1/P) * sum(slack_x[p] / x[p])
        def objective_rule(model):
            return model.t - sum(
                model.slack_x[k] / x.iloc[k] for k in model.K
            ) / len(xcol)
        model.obj = Objective(rule=objective_rule, sense=minimize)
        
        # Charnes-Cooper转换约束
        def cc_trans_rule(model):
            return 1 == model.t + sum(
                model.slack_y[l] / y.iloc[l] for l in model.L
            ) / len(ycol)
        model.cc_trans = Constraint(rule=cc_trans_rule)
        
        # 投入约束：t * x[k] = sum(intensity[i] * xref[i,k]) + slack_x[k]
        def input_rule(model, k):
            return sum(
                model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I
            ) == x.iloc[k] * model.t - model.slack_x[k]
        model.input = Constraint(model.K, rule=input_rule)
        
        # 产出约束：t * y[l] = sum(intensity[i] * yref[i,l]) - slack_y[l]
        def output_rule(model, l):
            return sum(
                model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I
            ) == y.iloc[l] * model.t + model.slack_y[l]
        model.output = Constraint(model.L, rule=output_rule)
        
        # === Step 5: 求解模型 ===
        opt = SolverFactory(solver)
        if not opt.available(exception_flag=False):
            raise RuntimeError(f"求解器 '{solver}' 未安装。请安装后重试。")
        
        solution = opt.solve(model)
        
        if solution.solver.termination_condition == "optimal":
            # 提取结果（线性变换后松弛变量 slack = Slack / t）
            t_val = value(model.t)
            slack_x_results[j] = np.asarray(
                list(model.slack_x[:].value)
            ) / t_val
            slack_y_results[j] = np.asarray(
                list(model.slack_y[:].value)
            ) / t_val
            obj_results[j] = value(model.obj)
    
    # === Step 6: 格式化输出结果 ===
    obj_df = pd.DataFrame(obj_results, index=["obj"]).T
    
    slack_x_df = pd.DataFrame(slack_x_results, index=xcol).T
    slack_x_df.columns = slack_x_df.columns.map(lambda x: "slack " + str(x))
    
    slack_y_df = pd.DataFrame(slack_y_results, index=ycol).T
    slack_y_df.columns = slack_y_df.columns.map(lambda y: "slack " + str(y))
    
    result = pd.concat([obj_df, slack_x_df, slack_y_df], axis=1)
    return result


# Markdown line 214
import pandas as pd
from deabook.example.book_cases.loaders import china

# 读取数据
ex4 = china(scaled=True)

# 计算SBM效率（评价第1-3期数据，参照集为第1-3期）
sbm_result = sbm("Y = K L", ex4, "t==[1,2,3]", "t==[1,2,3]")
print(sbm_result)


# Markdown line 266
def sbm_undesirable(formula, dataframe, eval_query=None, ref_query=None,
                    solver="highs"):
    """
    计算含非期望产出的SBM效率值（Tone, 2003）。
    
    Parameters
    ----------
    formula : str
        公式字符串，格式为 "产出:非期望产出 = 投入1 投入2 ..."
        例如："Y:CO2 = K L"
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框
    eval_query : str or None, default None
        筛选待评价决策单元的query字符串
    ref_query : str or None, default None
        筛选参照决策单元的query字符串
    solver : str, default "highs"
        线性规划求解器名称
    
    Returns
    -------
    pd.DataFrame
        包含效率值(obj)、投入松弛、期望产出松弛、非期望产出松弛的DataFrame
    
    Examples
    --------
    >>> result = sbm_undesirable("Y:CO2 = K L", ex4, "t==[1,2,3]", "t==[1,2,3]")
    """
    # === Step 1: 解析公式字符串（投入:非期望产出 = 投入变量）===
    input_vars = formula.split("=")[1].strip(" ")
    xcol = re.compile(" +").sub(" ", input_vars).split(" ")
    
    output_vars = formula.split("=")[0].split(":")[0].strip(" ")
    ycol = re.compile(" +").sub(" ", output_vars).split(" ")
    
    undesirable_vars = formula.split("=")[0].split(":")[1].strip(" ")
    bcol = re.compile(" +").sub(" ", undesirable_vars).split(" ")
    
    # === Step 2: 根据查询条件准备数据 ===
    if eval_query is None:
        index_eval = dataframe.index
    else:
        index_eval = dataframe.query(eval_query).index
    
    if ref_query is None:
        index_ref = dataframe.index
    else:
        index_ref = dataframe.query(ref_query).index
    
    data = dataframe.loc[index_eval, xcol + ycol + bcol]
    dataref = dataframe.loc[index_ref, xcol + ycol + bcol]
    
    xref = dataref.loc[:, xcol]   # 参照单元投入
    yref = dataref.loc[:, ycol]   # 参照单元期望产出
    bref = dataref.loc[:, bcol]   # 参照单元非期望产出（新增）
    
    # === Step 3: 数据验证 ===
    assert data[xcol + ycol + bcol].notna().all().all(), "数据包含缺失值"
    assert (data[xcol] > 0).all().all(), "投入数据必须严格为正"
    assert (data[ycol] > 0).all().all(), "期望产出数据必须严格为正"
    assert (data[bcol] >= 0).all().all(), "非期望产出不能为负"
    
    # 初始化结果存储
    slack_x_results = {}
    slack_y_results = {}
    slack_b_results = {}   # 非期望产出松弛（新增）
    obj_results = {}
    
    # === Step 4: 对每个决策单元构建并求解线性规划 ===
    for j in data.index:
        x = data.loc[j, xcol]   # 被评价单元投入
        y = data.loc[j, ycol]   # 被评价单元期望产出
        b = data.loc[j, bcol]   # 被评价单元非期望产出（新增）
        
        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)
        model.K = Set(initialize=range(len(xcol)))
        model.L = Set(initialize=range(len(ycol)))
        model.M = Set(initialize=range(len(bcol)))   # 非期望产出集合（新增）
        
        model.t = Var(domain=PositiveReals, doc="CC transformation variable")
        model.slack_x = Var(model.K, bounds=(0.0, None), doc="input slack")
        model.slack_y = Var(model.L, bounds=(0.0, None), doc="desirable output slack")
        model.slack_b = Var(model.M, bounds=(0.0, None), doc="undesirable output slack")
        model.intensity = Var(model.I, bounds=(0.0, None), doc="intensity variables")
        
        # 目标函数（与基础SBM相同）
        def objective_rule(model):
            return model.t - sum(
                model.slack_x[k] / x.iloc[k] for k in model.K
            ) / len(xcol)
        model.obj = Objective(rule=objective_rule, sense=minimize)
        
        # CC转换约束（分母变为 Q+R，增加非期望产出项）
        def cc_trans_rule(model):
            return 1 == model.t + (
                sum(model.slack_y[l] / y.iloc[l] for l in model.L)
                + sum(model.slack_b[m] / b.iloc[m] for m in model.M)
            ) / (len(ycol) + len(bcol))
        model.cc_trans = Constraint(rule=cc_trans_rule)
        
        # 投入约束（与基础SBM相同）
        def input_rule(model, k):
            return sum(
                model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I
            ) == x.iloc[k] * model.t - model.slack_x[k]
        model.input = Constraint(model.K, rule=input_rule)
        
        # 期望产出约束（与基础SBM相同）
        def output_rule(model, l):
            return sum(
                model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I
            ) == y.iloc[l] * model.t + model.slack_y[l]
        model.output = Constraint(model.L, rule=output_rule)
        
        # 非期望产出约束（新增）：t*b = sum(lambda*bref) + slack_b
        def undesirable_output_rule(model, m):
            return sum(
                model.intensity[i] * bref.loc[i, bcol[m]] for i in model.I
            ) == b.iloc[m] * model.t - model.slack_b[m]
        model.undesirable_output = Constraint(model.M, rule=undesirable_output_rule)
        
        # === Step 5: 求解模型 ===
        opt = SolverFactory(solver)
        if not opt.available(exception_flag=False):
            raise RuntimeError(f"求解器 '{solver}' 未安装。请安装后重试。")
        
        solution = opt.solve(model)
        
        if solution.solver.termination_condition == "optimal":
            t_val = value(model.t)
            slack_x_results[j] = np.asarray(list(model.slack_x[:].value)) / t_val
            slack_y_results[j] = np.asarray(list(model.slack_y[:].value)) / t_val
            slack_b_results[j] = np.asarray(list(model.slack_b[:].value)) / t_val
            obj_results[j] = value(model.obj)
    
    # === Step 6: 格式化输出结果 ===
    obj_df = pd.DataFrame(obj_results, index=["obj"]).T
    
    slack_x_df = pd.DataFrame(slack_x_results, index=xcol).T
    slack_x_df.columns = slack_x_df.columns.map(lambda x: "slack " + str(x))
    
    slack_y_df = pd.DataFrame(slack_y_results, index=ycol).T
    slack_y_df.columns = slack_y_df.columns.map(lambda y: "slack " + str(y))
    
    result = pd.concat([obj_df, slack_x_df, slack_y_df], axis=1)
    
    slack_b_df = pd.DataFrame(slack_b_results, index=bcol).T
    slack_b_df.columns = slack_b_df.columns.map(lambda b: "slack " + str(b))
    result = pd.concat([result, slack_b_df], axis=1)
    
    return result


# Markdown line 425
import pandas as pd
from deabook.example.book_cases.loaders import china

ex4 = china(scaled=True)
result = sbm_undesirable("Y:CO2 = K L", ex4, "t==[1,2,3]", "t==[1,2,3]")
print(result)


# Markdown line 470
"""DEA实战演练——公共导入模块"""
from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import pandas as pd
from deabook.example.book_cases.loaders import china
import numpy as np
import re
import warnings

warnings.filterwarnings("ignore")


# Markdown line 487
def normalize_formula(formula_str):
    """
    规范化公式字符串：统一分隔符、去除多余空格。
    
    Parameters
    ----------
    formula_str : str
        原始公式字符串，如 " Y  = K     L " 或 "Y:CO2 = K L"
    
    Returns
    -------
    str
        规范化后的公式字符串
    """
    # 统一全角符号为半角符号
    formula_str = formula_str.replace("：", ":").replace("＝", "=")
    # 将多个连续空格替换为单个空格
    formula_str = re.sub(r'\s+', ' ', formula_str).strip()
    return formula_str


def parse_formula(formula_str):
    """
    解析DEA公式字符串，提取投入变量名和产出变量名。
    
    支持格式：
    - 基础格式："产出 = 投入1 投入2"（如 "Y = K L"）
    - 含非期望产出格式："产出:非期望产出 = 投入1 投入2"（如 "Y:CO2 = K L"）
    
    Parameters
    ----------
    formula_str : str
        公式字符串，等号左边为产出变量，右边为投入变量
    
    Returns
    -------
    tuple
        (xcol, ycol, bcol) —— 分别为投入变量列表、期望产出变量列表、
        非期望产出变量列表。若无非期望产出，bcol为空列表。
    """
    formula_str = normalize_formula(formula_str)
    
    parts = formula_str.split("=")
    if len(parts) != 2:
        raise ValueError(
            f"公式格式错误：'{formula_str}'。应包含一个 '=' 分隔符"
        )
    
    left_side = parts[0].strip()
    right_side = parts[1].strip()
    
    xcol = [v.strip() for v in right_side.split(" ") if v.strip()]
    
    if ":" in left_side:
        output_parts = left_side.split(":")
        ycol = [v.strip() for v in output_parts[0].split(" ") if v.strip()]
        bcol = [v.strip() for v in output_parts[1].split(" ") if v.strip()]
    else:
        ycol = [v.strip() for v in left_side.split(" ") if v.strip()]
        bcol = []
    
    return xcol, ycol, bcol


# Markdown line 554
def validate_data(df, xcol, ycol, bcol=None):
    """
    验证投入产出数据是否符合DEA模型要求。
    
    Parameters
    ----------
    df : pd.DataFrame
        包含投入产出数据的DataFrame
    xcol : list of str
        投入变量列名列表
    ycol : list of str
        期望产出变量列名列表
    bcol : list of str or None, default None
        非期望产出变量列名列表
    
    Raises
    ------
    AssertionError
        数据不符合DEA模型要求时抛出
    """
    cols = xcol + ycol
    if bcol:
        cols = cols + bcol
    
    assert df[cols].notna().all().all(), \
        f"数据包含缺失值(NA)，请检查以下变量: {cols}"
    assert (df[xcol] > 0).all().all(), \
        f"投入变量 {xcol} 中存在非正值。SBM和DDF模型要求投入数据严格为正。"
    assert (df[ycol] > 0).all().all(), \
        f"期望产出变量 {ycol} 中存在非正值。"
    if bcol:
        assert (df[bcol] >= 0).all().all(), \
            f"非期望产出变量 {bcol} 中存在负值。"


def prepare_data(formula, dataframe, eval_query=None, ref_query=None):
    """
    根据公式和查询条件准备评价数据和参照数据。
    
    Parameters
    ----------
    formula : str
        公式字符串
    dataframe : pd.DataFrame
        包含投入产出数据的DataFrame
    eval_query : str or None, default None
        筛选待评价决策单元的query字符串
    ref_query : str or None, default None
        筛选参照决策单元的query字符串
    
    Returns
    -------
    tuple
        (data, dataref, xcol, ycol, bcol)
    """
    xcol, ycol, bcol = parse_formula(formula)
    
    if eval_query is None:
        index_eval = dataframe.index
    else:
        index_eval = dataframe.query(eval_query).index
    
    if ref_query is None:
        index_ref = dataframe.index
    else:
        index_ref = dataframe.query(ref_query).index
    
    all_cols = xcol + ycol + (bcol if bcol else [])
    data = dataframe.loc[index_eval, all_cols]
    dataref = dataframe.loc[index_ref, all_cols]
    
    return data, dataref, xcol, ycol, bcol


# Markdown line 631
def solve_dea_model(model, solver_name="highs", timeout=60):
    """
    求解DEA模型并返回求解结果，包含完整的异常处理和状态检查。
    
    Parameters
    ----------
    model : pyomo.environ.ConcreteModel
        已构建的Pyomo模型对象
    solver_name : str, default "highs"
        求解器名称。推荐 "highs"；备选 "glpk"、"ipopt"
    timeout : int, default 60
        求解超时时间（秒）
    
    Returns
    -------
    pyomo.opt.results.SolverResults
        求解器返回的结果对象
    
    Raises
    ------
    RuntimeError
        求解器未安装、求解失败或未达到最优时抛出
    """
    solver = SolverFactory(solver_name)
    
    if not solver.available(exception_flag=False):
        raise RuntimeError(
            f"求解器 '{solver_name}' 未安装或不可用。\n"
            f"安装建议：conda install {solver_name}"
        )
    
    results = solver.solve(model, timelimit=timeout, load_solutions=False)
    
    if results.solver.status != SolverStatus.ok:
        raise RuntimeError(f"求解器状态异常: {results.solver.status}")
    
    if results.solver.termination_condition != TerminationCondition.optimal:
        if (results.solver.termination_condition ==
            TerminationCondition.infeasible):
            raise RuntimeError("线性规划无可行解。请检查数据是否存在异常值。")
        else:
            raise RuntimeError(
                f"求解未达到最优: {results.solver.termination_condition}"
            )
    
    model.solutions.load_from(results)
    return results


def format_sbm_results(obj_dict, slack_x_dict, slack_y_dict,
                       xcol, ycol, bcol=None, slack_b_dict=None):
    """
    将SBM模型的求解结果格式化为标准DataFrame输出。
    
    Parameters
    ----------
    obj_dict : dict
        目标函数值字典
    slack_x_dict, slack_y_dict : dict
        松弛变量字典
    xcol, ycol : list of str
        变量名列表
    bcol : list of str or None, default None
        非期望产出变量名列表
    slack_b_dict : dict or None, default None
        非期望产出松弛变量字典
    
    Returns
    -------
    pd.DataFrame
        包含效率值和各松弛变量的结果DataFrame
    """
    obj_df = pd.DataFrame(obj_dict, index=["obj"]).T
    
    slack_x_df = pd.DataFrame(slack_x_dict, index=xcol).T
    slack_x_df.columns = slack_x_df.columns.map(lambda x: "slack " + str(x))
    
    slack_y_df = pd.DataFrame(slack_y_dict, index=ycol).T
    slack_y_df.columns = slack_y_df.columns.map(lambda y: "slack " + str(y))
    
    result = pd.concat([obj_df, slack_x_df, slack_y_df], axis=1)
    
    if bcol and slack_b_dict:
        slack_b_df = pd.DataFrame(slack_b_dict, index=bcol).T
        slack_b_df.columns = slack_b_df.columns.map(lambda b: "slack " + str(
            b))
        result = pd.concat([result, slack_b_df], axis=1)
    
    return result


def format_ddf_results(obj_dict, index_list):
    """
    将DDF/NDDF模型的求解结果格式化为标准DataFrame输出。
    
    Parameters
    ----------
    obj_dict : dict
        目标函数值字典
    index_list : list
        DMU索引列表
    
    Returns
    -------
    pd.DataFrame
        包含方向距离值的结果DataFrame
    """
    return pd.DataFrame(obj_dict, index=["te"]).T


# Markdown line 748
# 验证求解器是否可用
from pyomo.environ import SolverFactory
for solver in ["highs", "glpk", "ipopt"]:
    opt = SolverFactory(solver)
    print(f"{solver}: {'可用' if opt.available(exception_flag=False) else '不可用'}")


# Markdown line 810
def ddf(formula, dataframe, gx=None, gy=None, gb=None,
        eval_query=None, ref_query=None, solver="highs"):
    """
    计算方向距离函数（Directional Distance Function, DDF）效率值。
    
    基于Chambers, Chung and Färe (1996, 1998)的DDF模型（含非期望产出的
    应用见Chung, Färe and Grosskopf, 1997），支持投入、期望产出和非期望产出
    的自定义方向向量。
    
    Parameters
    ----------
    formula : str
        公式字符串，格式为 "产出:非期望产出 = 投入1 投入2 ..."
        例如："Y:CO2 = K L"
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框
    gx : list or None, default None
        投入方向向量。None表示默认 [-1, -1, ...]（投入缩减）
    gy : list or None, default None
        期望产出方向向量。None表示默认 [1, 1, ...]（产出扩张）
    gb : list or None, default None
        非期望产出方向向量。None表示默认 [-1, -1, ...]（非期望产出缩减）
    eval_query : str or None, default None
        筛选待评价决策单元的query字符串
    ref_query : str or None, default None
        筛选参照决策单元的query字符串
    solver : str, default "highs"
        线性规划求解器名称
    
    Returns
    -------
    pd.DataFrame
        包含方向距离值（te列）的结果DataFrame
    
    Notes
    -----
    方向向量默认取：
    - gx = [-1, -1, ...]：投入按比例缩减
    - gy = [1, 1, ...]：期望产出按比例扩张
    - gb = [-1, -1, ...]：非期望产出按比例缩减
    """
    # 复用公共工具：解析公式并准备数据
    data, dataref, xcol, ycol, bcol = prepare_data(
        formula, dataframe, eval_query, ref_query
    )
    validate_data(dataframe, xcol, ycol, bcol)
    
    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    bref = dataref.loc[:, bcol]
    
    # 设置默认方向向量（若未提供）
    if gx is None:
        gx = [-1] * len(xcol)
    if gy is None:
        gy = [1] * len(ycol)
    if gb is None:
        gb = [-1] * len(bcol)
    
    obj_results = {}
    
    for j in data.index:
        x = data.loc[j, xcol]
        y = data.loc[j, ycol]
        b = data.loc[j, bcol]
        
        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)
        model.K = Set(initialize=range(len(xcol)))
        model.L = Set(initialize=range(len(ycol)))
        model.M = Set(initialize=range(len(bcol)))
        
        # DDF特有变量：方向距离 beta（无约束）
        model.beta = Var(bounds=(None, None), within=Reals,
                         doc='directional distance')
        model.intensity = Var(model.I, bounds=(0.0, None),
                              doc='intensity variables')
        
        # 目标函数：maximize beta
        def objective_rule(model):
            return model.beta
        model.obj = Objective(rule=objective_rule, sense=maximize)
        
        # 投入约束：sum(lambda * xref) <= x + beta * gx * x
        def input_rule(model, k):
            return (sum(model.intensity[i] * xref.loc[i, xcol[k]]
                        for i in model.I)
                    - model.beta * gx[k] * x.iloc[k] <= x.iloc[k])
        model.input = Constraint(model.K, rule=input_rule)
        
        # 期望产出约束：sum(lambda * yref) >= y + beta * gy * y
        def output_rule(model, l):
            return (-sum(model.intensity[i] * yref.loc[i, ycol[l]]
                         for i in model.I)
                    + model.beta * gy[l] * y.iloc[l] <= -y.iloc[l])
        model.output = Constraint(model.L, rule=output_rule)
        
        # 非期望产出约束：sum(lambda * bref) = b + beta * gb * b
        def undesirable_output_rule(model, m):
            return (sum(model.intensity[i] * bref.loc[i, bcol[m]]
                        for i in model.I)
                    - model.beta * gb[m] * b.iloc[m] == b.iloc[m])
        model.undesirable_output = Constraint(model.M,
            rule=undesirable_output_rule)
        
        # 求解并提取结果
        try:
            solution = solve_dea_model(model, solver)
            obj_results[j] = value(model.obj)
        except RuntimeError as e:
            print(f"DMU {j} 求解失败: {e}")
            obj_results[j] = np.nan
    
    # 复用公共工具：格式化输出
    return format_ddf_results(obj_results, data.index)


# Markdown line 944
import pandas as pd
from deabook.example.book_cases.loaders import china
ex4 = china(scaled=True)
ddf_result = ddf("Y:CO2 = K L", ex4)
print(ddf_result)


# Markdown line 998
def nddf(formula, dataframe, gx=None, gy=None, gb=None,
         weight=None, eval_query=None, ref_query=None, solver="highs"):
    """
    计算非径向方向距离函数（NDDF）效率值（Zhou et al., 2012）。
    
    NDDF允许投入和产出进行不同比例的调整，通过权重向量w对各调整因子
    进行加权求和作为目标函数。
    
    Parameters
    ----------
    formula : str
        公式字符串，格式为 "产出:非期望产出 = 投入1 投入2 ..."
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框
    gx, gy, gb : list or None, default None
        投入、期望产出、非期望产出的方向向量。None时使用默认值
    weight : list or None, default None
        权重向量 [w_X, w_Y, w_B]。None时自动计算等权重
    eval_query : str or None, default None
        筛选待评价决策单元的query字符串
    ref_query : str or None, default None
        筛选参照决策单元的query字符串
    solver : str, default "highs"
        线性规划求解器名称
    
    Returns
    -------
    pd.DataFrame
        包含目标函数值(obj)和各变量调整因子(scale)的结果DataFrame
    """
    # 复用公共工具：解析公式并准备数据
    data, dataref, xcol, ycol, bcol = prepare_data(
        formula, dataframe, eval_query, ref_query
    )
    validate_data(dataframe, xcol, ycol, bcol)
    
    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    bref = dataref.loc[:, bcol]
    
    # 设置默认方向向量
    if gx is None:
        gx = [-1] * len(xcol)
    if gy is None:
        gy = [1] * len(ycol)
    if gb is None:
        gb = [-1] * len(bcol)
    
    # 设置默认权重（等权重）
    if weight is None:
        weight = []
        # 计算非零方向向量的种类数
        num_directions = (
            int(gx[0] != 0) + int(gy[0] != 0) + int(gb[0] != 0)
        )
        # 每种变量类型内部等权重
        for _ in range(len(xcol)):
            weight.append(1.0 / num_directions / len(xcol))
        for _ in range(len(ycol)):
            weight.append(1.0 / num_directions / len(ycol))
        for _ in range(len(bcol)):
            weight.append(1.0 / num_directions / len(bcol))
    
    # 分离各变量类型的权重
    w_x = weight[0:len(xcol)]
    w_y = weight[len(xcol):len(xcol) + len(ycol)]
    w_b = weight[len(xcol) + len(ycol):len(xcol) + len(ycol) + len(bcol)]
    
    # 初始化结果存储
    beta_x_results = {}
    beta_y_results = {}
    beta_b_results = {}
    obj_results = {}
    
    for j in data.index:
        x = data.loc[j, xcol]
        y = data.loc[j, ycol]
        b = data.loc[j, bcol]
        
        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)
        model.K = Set(initialize=range(len(xcol)))
        model.L = Set(initialize=range(len(ycol)))
        model.M = Set(initialize=range(len(bcol)))
        
        # NDDF特有变量：各变量类型的独立调整因子
        model.beta_x = Var(model.K, bounds=(0.0, None), doc='scale factor x')
        model.beta_y = Var(model.L, bounds=(0.0, None), doc='scale factor y')
        model.beta_b = Var(model.M, bounds=(0.0, None), doc='scale factor b')
        model.intensity = Var(model.I, bounds=(0.0, None),
                              doc='intensity variables')
        
        # NDDF特有目标函数：加权求和 w' * beta
        def objective_rule(model):
            return (
                -sum(w_x[k] * gx[k] * model.beta_x[k] for k in model.K)
                + sum(w_y[l] * gy[l] * model.beta_y[l] for l in model.L)
                - sum(w_b[m] * gb[m] * model.beta_b[m] for m in model.M)
            )
        model.obj = Objective(rule=objective_rule, sense=maximize)
        
        # 投入约束
        def input_rule(model, k):
            return (sum(model.intensity[i] * xref.loc[i, xcol[k]]
                        for i in model.I)
                    - gx[k] * x.iloc[k] * model.beta_x[k] <= x.iloc[k])
        model.input = Constraint(model.K, rule=input_rule)
        
        # 期望产出约束
        def output_rule(model, l):
            return (-sum(model.intensity[i] * yref.loc[i, ycol[l]]
                         for i in model.I)
                    + gy[l] * y.iloc[l] * model.beta_y[l] <= -y.iloc[l])
        model.output = Constraint(model.L, rule=output_rule)
        
        # 非期望产出约束
        def undesirable_output_rule(model, m):
            return (sum(model.intensity[i] * bref.loc[i, bcol[m]]
                        for i in model.I)
                    - gb[m] * b.iloc[m] * model.beta_b[m] == b.iloc[m])
        model.undesirable_output = Constraint(model.M,
            rule=undesirable_output_rule)
        
        # 求解并提取结果
        try:
            solution = solve_dea_model(model, solver)
            
            beta_x_results[j] = np.asarray(list(model.beta_x[:].value))
            beta_y_results[j] = np.asarray(list(model.beta_y[:].value))
            beta_b_results[j] = np.asarray(list(model.beta_b[:].value))
            obj_results[j] = value(model.obj)
            
        except RuntimeError as e:
            print(f"DMU {j} 求解失败: {e}")
            beta_x_results[j] = np.full(len(xcol), np.nan)
            beta_y_results[j] = np.full(len(ycol), np.nan)
            beta_b_results[j] = np.full(len(bcol), np.nan)
            obj_results[j] = np.nan
    
    # 格式化输出
    obj_df = pd.DataFrame(obj_results, index=["obj"]).T
    
    beta_x_df = pd.DataFrame(beta_x_results, index=xcol).T
    beta_x_df.columns = beta_x_df.columns.map(lambda x: "scale " + str(x))
    
    beta_y_df = pd.DataFrame(beta_y_results, index=ycol).T
    beta_y_df.columns = beta_y_df.columns.map(lambda y: "scale " + str(y))
    
    beta_b_df = pd.DataFrame(beta_b_results, index=bcol).T
    beta_b_df.columns = beta_b_df.columns.map(lambda b: "scale " + str(b))
    
    result = pd.concat([obj_df, beta_x_df, beta_y_df, beta_b_df], axis=1)
    return result


# Markdown line 1157
import pandas as pd
from deabook.example.book_cases.loaders import china

ex4 = china(scaled=True)
nddf_result = nddf("Y:CO2 = K L", ex4)
print(nddf_result)


# Markdown line 1235
def dea8(formula, dataframe, rts, orient,
         eval_query=None, ref_query=None, solver="highs"):
    """
    计算谢泼德距离函数效率值（Malmquist计算的基础函数）。
    
    支持CRS/VRS规模假设和产出导向(oo)/投入导向(io)两种方向。
    
    Parameters
    ----------
    formula : str
        公式字符串，如 "Y = K L"
    dataframe : pd.DataFrame
        投入产出数据框
    rts : str
        规模报酬假设："crs"（不变）或 "vrs"（可变）
    orient : str
        导向类型："oo"（产出导向）或 "io"（投入导向）
    eval_query : str or None, default None
        筛选待评价决策单元的query字符串
    ref_query : str or None, default None
        筛选参照决策单元的query字符串
    solver : str, default "highs"
        求解器名称
    
    Returns
    -------
    pd.DataFrame
        包含theta值和技术效率值(te)的DataFrame
    """
    # 复用公共工具：解析公式并准备数据
    data, dataref, xcol, ycol, _ = prepare_data(
        formula, dataframe, eval_query, ref_query
    )
    
    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    
    theta_results = {}   # 存储各DMU的theta值（原thetalt）
    
    for j in data.index:
        x = data.loc[j, xcol]
        y = data.loc[j, ycol]
        
        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)
        model.K = Set(initialize=range(len(xcol)))
        model.L = Set(initialize=range(len(ycol)))
        
        model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
        model.intensity = Var(model.I, bounds=(0.0, None),
            doc='intensity variables')
        
        # 目标函数：根据导向确定最大/最小化
        def obj_rule(model):
            return model.theta
        model.obj = Objective(
            rule=obj_rule,
            sense=(maximize if orient == "oo" else minimize)
        )
        
        # 投入约束：根据导向变化
        def input_rule(model, k):
            if orient == "oo":       # 产出导向：投入不超过被评价单元
                return sum(
                    model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I
                ) <= x.iloc[k]
            elif orient == "io":     # 投入导向：投入不超过theta倍被评价单元
                return sum(
                    model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I
                ) <= model.theta * x.iloc[k]
        
        # 产出约束：根据导向变化
        def output_rule(model, l):
            if orient == "oo":       # 产出导向：产出不低于theta倍被评价单元
                return sum(
                    model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I
                ) >= model.theta * y.iloc[l]
            elif orient == "io":     # 投入导向：产出不低于被评价单元
                return sum(
                    model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I
                ) >= y.iloc[l]
        
        # VRS约束
        def vrs_rule(model):
            return sum(model.intensity[i] for i in model.I) == 1
        
        model.input = Constraint(model.K, rule=input_rule)
        model.output = Constraint(model.L, rule=output_rule)
        
        if rts == "vrs":
            model.vrs = Constraint(rule=vrs_rule)
        
        # 求解
        try:
            solution = solve_dea_model(model, solver)
            theta = np.asarray(list(model.theta[:].value))
            theta_results[j] = theta[0]
        except RuntimeError as e:
            print(f"DMU {j} 求解失败: {e}")
            theta_results[j] = np.nan
    
    # 格式化输出
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    # 技术效率：产出导向取倒数，投入导向直接取theta
    thetadf["te"] = (
        1.0 / thetadf["theta"] if orient == "oo" else thetadf["theta"]
    )
    return thetadf


# Markdown line 1365
def mpi(formula, data, id_var, t_var, orient="oo", tech=None):
    """
    计算Malmquist生产率指数（支持多种生产技术设定）。
    
    基于Fare et al. (1994)的Malmquist指数，支持当期(com)、时序(seq)、
    视窗(window)和全局(global)四种生产技术设定。
    
    Parameters
    ----------
    formula : str
        公式字符串，如 "Y = K L"
    data : pd.DataFrame
        面板数据（包含id_var和t_var）
    id_var : str
        个体标识变量名
    t_var : str
        时间标识变量名
    orient : str, default "oo"
        导向类型："oo"（产出导向）或 "io"（投入导向）
    tech : str or None, default None
        生产技术设定：
        - None / "com"：当期生产技术（默认）
        - "seq"：时序生产技术
        - "window h"：视窗生产技术（h为窗宽）
        - "global"：全局生产技术
    
    Returns
    -------
    pd.DataFrame
        原始数据附加Malmquist指数及其分解（effch, techch, mpi）
    
    Examples
    --------
    >>> ex4 = china(scaled=True)
    >>> result = mpi("Y = K L", ex4, "id", "t", "oo", "com")
    >>> result_seq = mpi("Y = K L", ex4, "id", "t", "oo", "seq")
    >>> result_win = mpi("Y = K L", ex4, "id", "t", "oo", "window 2")
    """
    # 获取时间序列（去重并排序）
    t_series = pd.Series(data[t_var]).drop_duplicates().sort_values()
    t_values = t_series.tolist()
    
    # ---------- 计算 D^t(X_t, Y_t)：同期效率 ----------
    D11_list = []
    for t_val in t_values:
        eval_q = "{}=={}".format(t_var, t_val)
        
        # 根据生产技术设定确定参照集
        if tech is None or tech == "com":
            ref_q = "{}=={}".format(t_var, t_val)
        elif tech == "seq":
            ref_q = "{}<={}".format(t_var, t_val)
        elif tech == "global":
            ref_q = None
        elif "window " in str(tech):
            h = int(tech.split(" ")[1])
            t_idx = t_values.index(t_val)
            t_min = t_values[max(0, t_idx - h)]
            t_max = t_values[min(len(t_values) - 1, t_idx + h)]
            ref_q = "{}<={}<={}".format(t_min, t_var, t_max)
        
        D11 = dea8(formula, data, "crs", orient, eval_q, ref_q)[["te"]]
        D11_list.append(D11)
    
    D11_df = pd.concat(D11_list)
    D11_df.rename(columns={"te": "D11"}, inplace=True)
    
    # 全局技术只需D11即可计算
    if tech == "global":

        result = pd.merge(data, D11_df, left_index=True, right_index=True,
                          how="left")
        result["mpi"] = (
            result["D11"] / result.groupby(id_var)["D11"].shift(1)
        )
        result.drop(columns=["D11"], inplace=True)
        return result
    
    # ---------- 计算 D^t(X_{t+1}, Y_{t+1})：混合效率 ----------
    D12_list = []
    for t_val in t_values[1:]:   # 从第2个时期开始
        eval_q = "{}=={}".format(t_var, t_val)
        
        if tech is None or tech == "com":
            ref_q = "{}=={}".format(t_var, t_values[t_values.index(t_val) - 1])
        elif tech == "seq":
            ref_q = "{}<={}".format(t_var, t_values[t_values.index(t_val) - 1])
        elif "window " in str(tech):
            h = int(tech.split(" ")[1])
            t_idx = t_values.index(t_val)
            t_min = t_values[max(0, t_idx - 1 - h)]
            t_max = t_values[min(len(t_values) - 1, t_idx - 1 + h)]
            ref_q = "{}<={}<={}".format(t_min, t_var, t_max)
        
        D12 = dea8(formula, data, "crs", orient, eval_q, ref_q)[["te"]]
        D12_list.append(D12)
    
    D12_df = pd.concat(D12_list)
    D12_df.rename(columns={"te": "D12"}, inplace=True)
    
    # ---------- 计算 D^{t+1}(X_t, Y_t)：混合效率 ----------
    D21_list = []
    for t_val in t_values[1:]:
        eval_q = "{}=={}".format(t_var, t_values[t_values.index(t_val) - 1])
        
        if tech is None or tech == "com":
            ref_q = "{}=={}".format(t_var, t_val)
        elif tech == "seq":
            ref_q = "{}<={}".format(t_var, t_val)
        elif "window " in str(tech):
            h = int(tech.split(" ")[1])
            t_idx = t_values.index(t_val)
            t_min = t_values[max(0, t_idx - h)]
            t_max = t_values[min(len(t_values) - 1, t_idx + h)]
            ref_q = "{}<={}<={}".format(t_min, t_var, t_max)
        
        D21 = dea8(formula, data, "crs", orient, eval_q, ref_q)[["te"]]
        D21_list.append(D21)
    
    D21_df = pd.concat(D21_list)
    D21_df.rename(columns={"te": "D21"}, inplace=True)
    
    # ---------- 合并并计算Malmquist指数及其分解 ----------
    result_df = pd.concat([D11_df, D12_df, D21_df], axis=1)
    result = pd.merge(data, result_df, left_index=True, right_index=True,
                      how="left")
    
    # Malmquist指数
    # dea8 在两种导向下返回的 te 都是取值于 (0,1] 的效率：
    #   产出导向 te = 1/theta = 谢泼德产出距离 D_O；
    #   投入导向 te = theta = Farrell 投入效率 = 1/D_I。
    # 因此两种导向共用同一套比率公式，效率改善时 mpi>1，与 effch 方向一致。
    # MPI = [D^t(X_{t+1},Y_{t+1}) / D^t(X_t,Y_t)
    #        * D^{t+1}(X_{t+1},Y_{t+1}) / D^{t+1}(X_t,Y_t)]^{1/2}
    # 即 sqrt(D12/D11.shift(1) * D11/D21.shift(1))
    result["mpi"] = np.sqrt(result["D12"] / result.groupby(id_var)["D11"].shift(1)
                            * result["D11"] / result.groupby(id_var)["D21"].shift(1))
    
    # Malmquist分解
    # 效率变化 EC = D_{t+1}(X_{t+1},Y_{t+1}) / D_t(X_t,Y_t)
    # 即 同一DMU相邻期 D11 之比（按个体分组 shift，首期为 NaN）
    result["effch"] = result["D11"] / result.groupby(id_var)["D11"].shift(1)
    
    # 技术变化 TC = MPI / EC
    result["techch"] = result["mpi"] / result["effch"]
    
    # 清理中间结果列
    result.drop(columns=["D11", "D12", "D21"], inplace=True)
    
    return result


# Markdown line 1566
import pandas as pd
from deabook.example.book_cases.loaders import china

ex4 = china(scaled=True)

# 产出导向，当期生产技术
result_com = mpi("Y = K L", ex4, "id", "t", orient="oo", tech="com")

# 产出导向，时序生产技术
result_seq = mpi("Y = K L", ex4, "id", "t", orient="oo", tech="seq")

# 产出导向，视窗生产技术（窗宽=2）
result_win = mpi("Y = K L", ex4, "id", "t", orient="oo", tech="window 2")

# 产出导向，全局生产技术
result_glb = mpi("Y = K L", ex4, "id", "t", orient="oo", tech="global")

# 投入导向
result_io = mpi("Y = K L", ex4, "id", "t", orient="io", tech="com")


# Markdown line 1611
def mlpi(formula, data, id_var, t_var, tech=None):
    """
    计算Malmquist-Luenberger生产率指数（含非期望产出）。
    
    复用 mpi() 函数的框架，将距离函数从谢泼德距离函数替换为
    方向距离函数（DDF）。
    
    Parameters
    ----------
    formula : str
        公式字符串，如 "Y:CO2 = K L"
    data : pd.DataFrame
        面板数据
    id_var : str
        个体标识变量名
    t_var : str
        时间标识变量名
    tech : str or None, default None
        生产技术设定（同mpi函数）
    
    Returns
    -------
    pd.DataFrame
        原始数据附加ML指数的结果
    """
    t_series = pd.Series(data[t_var]).drop_duplicates().sort_values()
    t_values = t_series.tolist()
    
    # D^t(x_t, y_t, b_t)
    D11_list = []
    for t_val in t_values:
        eval_q = "{}=={}".format(t_var, t_val)
        
        if tech is None or tech == "com":
            ref_q = "{}=={}".format(t_var, t_val)
        elif tech == "seq":
            ref_q = "{}<={}".format(t_var, t_val)
        elif tech == "global":
            ref_q = None
        elif "window " in str(tech):
            h = int(tech.split(" ")[1])
            t_idx = t_values.index(t_val)
            t_min = t_values[max(0, t_idx - h)]
            t_max = t_values[min(len(t_values) - 1, t_idx + h)]
            ref_q = "{}<={}<={}".format(t_min, t_var, t_max)
        
        D11 = ddf(formula, data, eval_query=eval_q, ref_query=ref_q)[["te"]]
        D11_list.append(D11)
    
    D11_df = pd.concat(D11_list)
    D11_df.rename(columns={"te": "D11"}, inplace=True)
    
    if tech == "global":
        result = pd.merge(data, D11_df, left_index=True, right_index=True,
                          how="left")
        result["mlpi"] = ((1 + result.groupby(id_var)["D11"].shift(1)) / (1 +
            result["D11"]))
        result.drop(columns=["D11"], inplace=True)
        return result
    
    # D^t(x_{t+1}, y_{t+1}, b_{t+1})
    D12_list = []
    for t_val in t_values[1:]:
        eval_q = "{}=={}".format(t_var, t_val)
        
        if tech is None or tech == "com":
            ref_q = "{}=={}".format(t_var, t_values[t_values.index(t_val) - 1])
        elif tech == "seq":
            ref_q = "{}<={}".format(t_var, t_values[t_values.index(t_val) - 1])
        elif "window " in str(tech):
            h = int(tech.split(" ")[1])
            t_idx = t_values.index(t_val)
            t_min = t_values[max(0, t_idx - 1 - h)]
            t_max = t_values[min(len(t_values) - 1, t_idx - 1 + h)]
            ref_q = "{}<={}<={}".format(t_min, t_var, t_max)
        
        D12 = ddf(formula, data, eval_query=eval_q, ref_query=ref_q)[["te"]]
        D12_list.append(D12)
    
    D12_df = pd.concat(D12_list)
    D12_df.rename(columns={"te": "D12"}, inplace=True)
    
    # D^{t+1}(x_t, y_t, b_t)
    D21_list = []
    for t_val in t_values[1:]:
        eval_q = "{}=={}".format(t_var, t_values[t_values.index(t_val) - 1])
        
        if tech is None or tech == "com":
            ref_q = "{}=={}".format(t_var, t_val)
        elif tech == "seq":
            ref_q = "{}<={}".format(t_var, t_val)
        elif "window " in str(tech):
            h = int(tech.split(" ")[1])
            t_idx = t_values.index(t_val)
            t_min = t_values[max(0, t_idx - h)]
            t_max = t_values[min(len(t_values) - 1, t_idx + h)]
            ref_q = "{}<={}<={}".format(t_min, t_var, t_max)
        
        D21 = ddf(formula, data, eval_query=eval_q, ref_query=ref_q)[["te"]]
        D21_list.append(D21)
    
    D21_df = pd.concat(D21_list)
    D21_df.rename(columns={"te": "D21"}, inplace=True)
    
    # 合并并计算ML指数
    result_df = pd.concat([D11_df, D12_df, D21_df], axis=1)
    result = pd.merge(data, result_df, left_index=True, right_index=True,
                      how="left")
    
    # ML指数公式：sqrt[(1+D11.shift)/(1+D12) * (1+D21.shift)/(1+D11)]（shift 按个体分组）
    result["mlpi"] = np.sqrt((1 + result.groupby(id_var)["D11"].shift(1)) / (
        1 + result["D12"])
                             * (1 + result.groupby(id_var)["D21"].shift(1)) / (
                                 1 + result["D11"]))
    
    # Malmquist-Luenberger分解
    result["effch"] = (1 + result.groupby(id_var)["D11"].shift(1)) / (1 +
        result["D11"])
    result["techch"] = result["mlpi"] / result["effch"]
    
    result.drop(columns=["D11", "D12", "D21"], inplace=True)
    return result


# Markdown line 1757
import pandas as pd
from deabook.example.book_cases.loaders import china

ex4 = china(scaled=True)

# 计算ML指数
ml_result = mlpi("Y:CO2 = K L", ex4, "id", "t", tech="com")
print(ml_result)


# Markdown line 1788
def validate_efficiency_results(result_df, model_type="sbm", same_frontier=True):
    """
    验证DEA效率结果的合理性。
    
    Parameters
    ----------
    result_df : pd.DataFrame
        包含效率值的结果DataFrame
    model_type : str, default "sbm"
        模型类型："sbm"、"ddf"、"nddf"、"dea"
    
    Returns
    -------
    dict
        验证结果报告
    
    Raises
    ------
    AssertionError
        验证不通过时抛出警告
    """
    report = {"model_type": model_type, "checks_passed": True, "warnings": []}
    
    if model_type == "sbm":
        eff = result_df["obj"]
        # SBM效率应在 (0, 1] 区间
        assert (eff.dropna() > 0).all(), "SBM效率值应为正"
        assert (eff.dropna() <= 1.001).all(), "SBM效率值不应超过1（允许数值误差）"
        
        # 松弛变量应非负
        slack_cols = [c for c in result_df.columns if c.startswith("slack ")]
        for col in slack_cols:
            assert (result_df[col].dropna() >= -1e-6).all(), \
                f"松弛变量 {col} 出现负值"
    
    elif model_type in ["ddf", "nddf", "dea"]:
        eff = result_df["te"] if "te" in result_df.columns else result_df[
            "obj"]
        # DDF方向距离应 >= 0
        if same_frontier:
            assert (eff.dropna() >= -1e-6).all(), "同期方向距离不应为负"
    
    else:
        raise ValueError("未知模型类型")

    # 检查NaN比例
    nan_ratio = eff.isna().mean()
    if nan_ratio > 0.1:
        report["warnings"].append(
            f"警告: {nan_ratio*100:.1f}% 的DMU求解失败，请检查数据质量"
        )
        report["checks_passed"] = False
    
    return report


# Markdown line 1860
def cross_validate_with_deabook(data, xcol, ycol, model="SBM"):
    """
    使用deabook包交叉验证自编程结果。
    
    Parameters
    ----------
    data : pd.DataFrame
        投入产出数据
    xcol : list of str
        投入变量名
    ycol : list of str
        产出变量名
    model : str, default "SBM"
        模型类型
    
    Returns
    -------
    bool
        交叉验证是否通过
    """
    from deabook import DEA
    from deabook.constant import RTS_CRS, OPT_LOCAL
    if model.upper() != "SBM":
        raise ValueError("此函数只核验普通CRS SBM")
    formula = " ".join(ycol) + " = " + " ".join(xcol)
    benchmark = DEA.sbm(data, " ".join(xcol) + "=" + " ".join(ycol),
                        gx=[1]*len(xcol), gy=[1]*len(ycol),
                        rts=RTS_CRS).optimize(OPT_LOCAL, 'highs')
    actual = sbm(formula, data)["obj"]
    return bool(np.allclose(actual.to_numpy(), benchmark["rho"].to_numpy(),
                            atol=1e-6, rtol=1e-6))



# Markdown line 1906
def sensitivity_analysis_drop_dmu(formula, data, eval_query=None):
    """
    敏感性分析：逐一剔除每个DMU，观察其他DMU效率变化。
    
    若剔除某DMU后其他DMU效率大幅变化，则该DMU为"关键观测值"，
    对效率前沿有重要影响。
    
    Parameters
    ----------
    formula : str
        公式字符串
    data : pd.DataFrame
        投入产出数据
    eval_query : str or None, default None
        评价单元筛选条件
    
    Returns
    -------
    pd.DataFrame
        敏感性分析结果
    """
    # 基准效率
    base_result = sbm(formula, data, eval_query, eval_query)
    base_eff = base_result["obj"]
    
    # 准备数据
    if eval_query is not None:
        eval_data = data.query(eval_query)
    else:
        eval_data = data
    
    sensitivity = {}
    for dmu in eval_data.index:
        # 剔除该DMU
        data_drop = data.drop(index=dmu)
        result_drop = sbm(formula, data_drop, eval_query, eval_query)
        
        # 计算效率变化幅度（均方根变化）
        common_idx = base_eff.index.intersection(result_drop.index)
        if len(common_idx) > 0:
            rms_change = np.sqrt(
                np.mean((base_eff[common_idx] - result_drop.loc[common_idx,
                    "obj"])**2)
            )
            sensitivity[dmu] = rms_change
    
    sens_df = pd.DataFrame.from_dict(
        sensitivity, orient="index", columns=["rms_change"]
    )
    sens_df = sens_df.sort_values("rms_change", ascending=False)
    
    print("敏感性分析结果（影响最大的前5个DMU）：")
    print(sens_df.head())
    
    return sens_df


def sensitivity_analysis_perturbation(formula, data, eval_query=None,
                                      noise_level=0.01, n_replicates=10):
    """
    敏感性分析：对投入数据施加随机扰动，观察效率稳定性。
    
    Parameters
    ----------
    formula : str
        公式字符串
    data : pd.DataFrame
        投入产出数据
    eval_query : str or None, default None
        评价单元筛选条件
    noise_level : float, default 0.01
        扰动幅度（如0.01表示1%）
    n_replicates : int, default 10
        重复次数
    
    Returns
    -------
    pd.DataFrame
        每次重复的效率结果
    """
    _, _, xcol, ycol, _ = prepare_data(formula, data)
    
    all_results = []
    for rep in range(n_replicates):
        # 复制数据并施加随机扰动
        data_perturbed = data.copy()
        for col in xcol:
            noise = np.random.normal(0, noise_level, size=len(data_perturbed))
            data_perturbed[col] = data_perturbed[col] * (1 + noise)
        
        result = sbm(formula, data_perturbed, eval_query, eval_query)
        result["replicate"] = rep
        all_results.append(result[["obj", "replicate"]])
    
    all_df = pd.concat(all_results)
    
    # 计算变异系数（CV = std / mean）
    cv = all_df.groupby(all_df.index)["obj"].agg(["mean", "std"])
    cv["cv"] = cv["std"] / cv["mean"]
    
    print(f"效率值的平均变异系数: {cv['cv'].mean():.4f}")
    print(f"最大变异系数: {cv['cv'].max():.4f}")
    
    return all_df


# Markdown line 2077
"""DEA工具包——公共工具模块

本模块封装了DEA模型求解的通用功能，包括：
- 公式解析（parse_formula, normalize_formula）
- 数据验证（validate_data, prepare_data）
- 模型求解（solve_dea_model）
- 结果格式化（format_sbm_results, format_ddf_results）
"""

from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import pandas as pd
from deabook.example.book_cases.loaders import china
import numpy as np
import re


# ========== 公式解析 ==========

def normalize_formula(formula_str):
    """规范化公式字符串（详见12.3.3节）"""
    formula_str = formula_str.replace('：', ':').replace('＝', '=')
    formula_str = re.sub(r'\s+', ' ', formula_str).strip()
    return formula_str


def parse_formula(formula_str):
    """解析DEA公式字符串（详见12.3.3节）"""
    formula_str = normalize_formula(formula_str)
    parts = formula_str.split('=')
    if len(parts) != 2:
        raise ValueError(f"公式格式错误：'{formula_str}'")
    
    left_side = parts[0].strip()
    right_side = parts[1].strip()
    
    xcol = [v.strip() for v in right_side.split(' ') if v.strip()]
    
    if ':' in left_side:
        output_parts = left_side.split(':')
        ycol = [v.strip() for v in output_parts[0].split(' ') if v.strip()]
        bcol = [v.strip() for v in output_parts[1].split(' ') if v.strip()]
    else:
        ycol = [v.strip() for v in left_side.split(' ') if v.strip()]
        bcol = []
    
    return xcol, ycol, bcol


# ========== 数据验证与准备 ==========

def validate_data(df, xcol, ycol, bcol=None):
    """验证数据是否符合DEA模型要求（详见12.3.4节）"""
    cols = xcol + ycol
    if bcol:
        cols = cols + bcol
    
    assert df[cols].notna().all().all(), "数据包含缺失值"
    assert (df[xcol] > 0).all().all(), "投入数据必须严格为正"
    assert (df[ycol] > 0).all().all(), "期望产出数据必须严格为正"
    if bcol:
        assert (df[bcol] >= 0).all().all(), "非期望产出不能为负"


def prepare_data(formula, dataframe, eval_query=None, ref_query=None):
    """根据公式和查询条件准备数据（详见12.3.4节）"""
    xcol, ycol, bcol = parse_formula(formula)
    
    if eval_query is None:
        index_eval = dataframe.index
    else:
        index_eval = dataframe.query(eval_query).index
    
    if ref_query is None:
        index_ref = dataframe.index
    else:
        index_ref = dataframe.query(ref_query).index
    
    all_cols = xcol + ycol + (bcol if bcol else [])
    data = dataframe.loc[index_eval, all_cols]
    dataref = dataframe.loc[index_ref, all_cols]
    
    return data, dataref, xcol, ycol, bcol


# ========== 求解与结果提取 ==========

def solve_dea_model(model, solver_name="highs", timeout=60):
    """求解DEA模型并返回结果（详见12.3.5节）"""
    solver = SolverFactory(solver_name)
    
    if not solver.available(exception_flag=False):
        raise RuntimeError(f"求解器 '{solver_name}' 未安装或不可用。")
    
    results = solver.solve(model, timelimit=timeout, load_solutions=False)
    
    if results.solver.status != SolverStatus.ok:
        raise RuntimeError(f"求解器状态异常: {results.solver.status}")
    
    if results.solver.termination_condition != TerminationCondition.optimal:
        raise RuntimeError(
            f"求解未达到最优: {results.solver.termination_condition}"
        )
    
    model.solutions.load_from(results)
    return results


def format_sbm_results(obj_dict, slack_x_dict, slack_y_dict,
                       xcol, ycol, bcol=None, slack_b_dict=None):
    """格式化SBM结果为DataFrame（详见12.3.5节）"""
    obj_df = pd.DataFrame(obj_dict, index=["obj"]).T
    
    slack_x_df = pd.DataFrame(slack_x_dict, index=xcol).T
    slack_x_df.columns = slack_x_df.columns.map(lambda x: "slack " + str(x))
    
    slack_y_df = pd.DataFrame(slack_y_dict, index=ycol).T
    slack_y_df.columns = slack_y_df.columns.map(lambda y: "slack " + str(y))
    
    result = pd.concat([obj_df, slack_x_df, slack_y_df], axis=1)
    
    if bcol and slack_b_dict:
        slack_b_df = pd.DataFrame(slack_b_dict, index=bcol).T
        slack_b_df.columns = slack_b_df.columns.map(lambda b: "slack " + str(
            b))
        result = pd.concat([result, slack_b_df], axis=1)
    
    return result


def format_ddf_results(obj_dict, index_list):
    """格式化DDF结果为DataFrame（详见12.3.5节）"""
    return pd.DataFrame(obj_dict, index=["te"]).T


# Markdown line 2215
# 按本章顺序执行时，直接复用前文已经定义的模型函数
import pandas as pd
from deabook.example.book_cases.loaders import china

# 读取数据
data = china(scaled=True)

# SBM模型
result_sbm = sbm("Y = K L", data, "t==[1,2,3]", "t==[1,2,3]")

# DDF模型
result_ddf = ddf("Y:CO2 = K L", data)

# Malmquist指数
result_mpi = mpi("Y = K L", data, "id", "t", orient="oo", tech="com")


# Markdown line 2242
class DEAModel:
    """DEA模型的基类，展示面向对象设计的思路"""
    
    def __init__(self, formula, data, solver="highs"):
        self.formula = formula
        self.data = data
        self.solver = solver
        self.xcol, self.ycol, self.bcol = parse_formula(formula)
        validate_data(data, self.xcol, self.ycol, self.bcol)
    
    def fit(self, eval_query=None, ref_query=None):
        """拟合模型（子类必须实现）"""
        raise NotImplementedError
    
    def summary(self):
        """输出结果摘要（子类必须实现）"""
        raise NotImplementedError


class SBMModel(DEAModel):
    """SBM模型的面向对象实现"""
    
    def fit(self, eval_query=None, ref_query=None):
        # 复用函数式sbm的核心逻辑
        self.result = sbm(self.formula, self.data, eval_query, ref_query,
            self.solver)
        return self
    
    def summary(self):
        print(f"SBM模型结果 (N={len(self.result)}):")
        print(f"  平均效率: {self.result['obj'].mean():.4f}")
        print(f"  有效DMU数: {(self.result['obj'] >= 0.999).sum()}")


# Markdown line 2289
"""DEA模型单元测试"""
import pandas as pd
from deabook.example.book_cases.loaders import china
import numpy as np


class TestFormulaParser:
    """测试公式解析函数"""
    
    def test_basic_formula(self):
        xcol, ycol, bcol = parse_formula("Y = K L")
        assert xcol == ["K", "L"]
        assert ycol == ["Y"]
        assert bcol == []
    
    def test_undesirable_formula(self):
        xcol, ycol, bcol = parse_formula("Y:CO2 = K L")
        assert xcol == ["K", "L"]
        assert ycol == ["Y"]
        assert bcol == ["CO2"]
    
    def test_normalize_whitespace(self):
        xcol, ycol, bcol = parse_formula("  Y   =  K     L  ")
        assert xcol == ["K", "L"]
        assert ycol == ["Y"]


class TestSBMModel:
    """测试SBM模型"""
    
    def test_sbm_efficiency_range(self):
        """SBM效率值应在(0, 1]区间"""
        data = china(scaled=True).query("t==4").copy()
        result = sbm("Y = K L", data)
        assert (result["obj"] > 0).all()
        assert (result["obj"] <= 1 + 1e-7).all()
    
    def test_sbm_efficient_dmu(self):
        """效率前沿上的DMU效率应为1"""
        data = china(scaled=True).query("t==4").copy()
        result = sbm("Y = K L", data)
        assert np.isclose(result["obj"].max(), 1.0, atol=1e-7)

# 不依赖测试框架也可直接执行全部检查。
for test_class in (TestFormulaParser, TestSBMModel):
    instance = test_class()
    for name in dir(instance):
        if name.startswith("test_"):
            getattr(instance, name)()



# Markdown line 2366
from deabook.example.book_cases.run import FUN
FUN[12]()  # 重新计算本章，输出完整CSV、摘要表与PNG/PDF图

