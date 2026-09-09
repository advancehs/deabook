# 本文件从教材原文提取；在项目根目录运行。

# Markdown line 272
from deabook.example.book_cases.loaders import china
raw = china()
columns = ["K", "L", "E", "Y", "CO2"]
scale = raw[columns].mean()
data = raw.copy()
data[columns] = raw[columns] / scale


# Markdown line 285
"""
基于SBM对偶问题的影子价格计算模块。

本模块提供基于SBM模型对偶问题估计投入、产出和非期望产出
影子价格的函数。求解器默认采用 highs。
"""
from pyomo.environ import (
    ConcreteModel, Set, Var, Objective, Constraint, Reals,
    maximize, minimize, SolverFactory, value
)
import pandas as pd
from deabook.example.book_cases.loaders import china, china_lower
import numpy as np
import re
import warnings

warnings.filterwarnings("ignore")

# 读取数据
ex4 = china(scaled=True)
data = ex4[["K", "L", "Y", "CO2"]].copy()  # 提取投入、期望产出和非期望产出变量

# 设定变量个数
K = 2  # 投入变量个数
L = 1  # 期望产出变量个数
M = 1  # 非期望产出变量个数

# 准备参照集数据
xref = data.iloc[:, 0:K]
yref = data.iloc[:, K:K+L]
bref = data.iloc[:, K+L:K+L+M]

# 准备被评价DMU数据（以第1个DMU为例）
x = data.iloc[0, 0:K]
y = data.iloc[0, K:K+L]
b = data.iloc[0, K+L:K+L+M]

# 构建Pyomo模型
model = ConcreteModel()
model.I = Set(initialize=range(xref.shape[0]))  # 参照决策单元集合
model.K = Set(initialize=range(len(x)))          # 投入变量集合
model.L = Set(initialize=range(len(y)))          # 期望产出变量集合
model.M = Set(initialize=range(len(b)))          # 非期望产出变量集合

# 定义影子价格变量（施加非负约束以确保经济学含义）
model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
model.py = Var(model.L, bounds=(0, None), within=Reals, doc='期望产出的影子价格')
model.pb = Var(model.M, bounds=(0, None), within=Reals, doc='非期望产出的影子价格')
model.pomega = Var(bounds=(0, None), within=Reals, doc='SBM效率值（式(13.5)中的 xi）')

# 目标函数：最大化 pomega
def objective_rule(model):
    """返回目标函数表达式"""
    return model.pomega * 1

# 约束1：对偶约束（参照集加权影子价格之和非正）
def first_rule(model, i):
    """第i个参照DMU的对偶约束"""
    return (
        -sum(model.px[k] * xref.iloc[i, k] for k in model.K)
        + sum(model.py[l] * yref.iloc[i, l] for l in model.L)
        - sum(model.pb[m] * bref.iloc[i, m] for m in model.M)
        <= 0
    )

# 约束2：投入影子价格的下界约束
def second_rule(model, k):
    """第k个投入变量的影子价格下界"""
    return -model.px[k] <= -1 / (len(model.K) * x.iloc[k])

# 约束3：期望产出影子价格的上界约束
def third_rule(model, l):
    """第l个期望产出变量的影子价格上界"""
    return (
        model.pomega / (len(model.L) + len(model.M)) / y.iloc[l] - model.py[l]
        <= 0
    )

# 约束4：非期望产出影子价格的上界约束
def forth_rule(model, m):
    """第m个非期望产出变量的影子价格上界"""
    return (
        model.pomega / (len(model.L) + len(model.M)) / b.iloc[m] - model.pb[m]
        <= 0
    )

# 约束5：Charnes-Cooper变换的归一化约束
def fifth_rule(model):
    """归一化约束"""
    return (
        model.pomega
        + sum(model.px[k] * x.iloc[k] for k in model.K)
        - sum(model.py[l] * y.iloc[l] for l in model.L)
        + sum(model.pb[m] * b.iloc[m] for m in model.M)
        <= 1
    )

# 组装模型
model.obj = Objective(rule=objective_rule, sense=maximize, doc='目标函数')
model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
model.third = Constraint(model.L, rule=third_rule, doc='期望产出上界约束')
model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出上界约束')
model.fifth = Constraint(rule=fifth_rule, doc='归一化约束')

# 求解
opt = SolverFactory('highs')  # 指定 HiGHS 作为求解器
solution = opt.solve(model)   # 调用求解器求解

# 提取结果
px = np.asarray(list(model.px[:].value))
py = np.asarray(list(model.py[:].value))
pb = np.asarray(list(model.pb[:].value))
obj = value(model.obj)

print(f"投入影子价格 (px): {px}")
print(f"期望产出影子价格 (py): {py}")
print(f"非期望产出影子价格 (pb): {pb}")
print(f"最优目标值: {obj}")


# Markdown line 417
px, py, pb, obj = {}, {}, {}, {}

for j in range(data.shape[0]):
    x = data.iloc[j, 0:K]
    y = data.iloc[j, K:K+L]
    b = data.iloc[j, K+L:K+L+M]

    model = ConcreteModel()
    model.I = Set(initialize=range(xref.shape[0]))
    model.K = Set(initialize=range(len(x)))
    model.L = Set(initialize=range(len(y)))
    model.M = Set(initialize=range(len(b)))

    model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
    model.py = Var(model.L, bounds=(0, None), within=Reals, doc='期望产出的影子价格')
    model.pb = Var(model.M, bounds=(0, None), within=Reals, doc='非期望产出的影子价格')
    model.pomega = Var(bounds=(0, None), within=Reals,
        doc='SBM效率值（式(13.5)中的 xi）')

    def objective_rule(model):
        return model.pomega * 1

    def first_rule(model, i):
        return (
            -sum(model.px[k] * xref.iloc[i, k] for k in model.K)
            + sum(model.py[l] * yref.iloc[i, l] for l in model.L)
            - sum(model.pb[m] * bref.iloc[i, m] for m in model.M)
            <= 0
        )

    def second_rule(model, k):
        return -model.px[k] <= -1 / (len(model.K) * x.iloc[k])

    def third_rule(model, l):
        return (
            model.pomega / (len(model.L) + len(model.M)) / y.iloc[
                l] - model.py[l]
            <= 0
        )

    def forth_rule(model, m):
        return (
            model.pomega / (len(model.L) + len(model.M)) / b.iloc[
                m] - model.pb[m]
            <= 0
        )

    def fifth_rule(model):
        return (
            model.pomega
            + sum(model.px[k] * x.iloc[k] for k in model.K)
            - sum(model.py[l] * y.iloc[l] for l in model.L)
            + sum(model.pb[m] * b.iloc[m] for m in model.M)
            <= 1
        )

    model.obj = Objective(rule=objective_rule, sense=maximize, doc='目标函数')
    model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
    model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
    model.third = Constraint(model.L, rule=third_rule, doc='期望产出上界约束')
    model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出上界约束')
    model.fifth = Constraint(rule=fifth_rule, doc='归一化约束')

    opt = SolverFactory('highs')
    solution = opt.solve(model)

    # 终止条件判断：optimal（最优）、infeasible（不可行）、unbounded（无界）
    if solution.solver.termination_condition == "optimal":
        px[j] = np.asarray(list(model.px[:].value))
        py[j] = np.asarray(list(model.py[:].value))
        pb[j] = np.asarray(list(model.pb[:].value))
        obj[j] = value(model.obj)


# Markdown line 496
from pyomo.environ import (
    ConcreteModel, Set, Var, Objective, Constraint, Reals,
    maximize, SolverFactory, value
)
import pandas as pd
from deabook.example.book_cases.loaders import china, china_lower
import numpy as np
import re


def sbm_dual(dataframe, varname, num_inputs, num_outputs):
    """
    基于SBM模型对偶问题估计影子价格。

    求解SBM模型的对偶问题（乘数形式），获得各决策单元的投入、
    期望产出和非期望产出的影子价格。

    Parameters
    ----------
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框，变量按投入在前、
        期望产出在中、非期望产出在后的顺序排列。
    varname : list[str]
        投入变量、期望产出变量和非期望产出变量的名称列表，
        如 ["K", "L", "Y", "CO2"]。
    num_inputs : int
        投入变量的个数（P）。
    num_outputs : int
        期望产出变量的个数（Q）。

    Returns
    -------
    pd.DataFrame
        包含各DMU的影子价格估计结果的数据框。列包括：
        - 'obj': SBM效率值（式(13.5)中的 xi，取值于 (0, 1]）
        - 'shadow price <varname>': 各变量的影子价格

    Notes
    -----
    - 本函数对影子价格施加非负约束，以确保经济学含义的合理性。
    - 求解器默认使用 highs。
    - 关于SBM模型及其含非期望产出扩展的原问题形式，参见第8章（式(8.1)与式(8.3)）。

    Examples
    --------
    >>> ex4 = china(scaled=True)
    >>> result = sbm_dual(ex4, ["K", "L", "Y", "CO2"], 2, 1)
    >>> print(result.head())
    """
    data = dataframe.loc[:, varname]
    num_inputs_ = num_inputs
    num_outputs_ = num_outputs
    num_bad = data.shape[1] - num_inputs_ - num_outputs_

    xref = data.iloc[:, 0:num_inputs_]
    yref = data.iloc[:, num_inputs_:num_inputs_+num_outputs_]
    bref = data.iloc[:, num_inputs_+num_outputs_:num_inputs_+num_outputs_+num_bad]

    xcol = varname[0:num_inputs_]
    ycol = varname[num_inputs_:num_inputs_+num_outputs_]
    bcol = varname[num_inputs_+num_outputs_:num_inputs_+num_outputs_+num_bad]

    px, py, pb, obj = {}, {}, {}, {}

    for j in range(data.shape[0]):
        x = data.iloc[j, 0:num_inputs_]
        y = data.iloc[j, num_inputs_:num_inputs_+num_outputs_]
        b = data.iloc[j, num_inputs_+num_outputs_:num_inputs_+num_outputs_+num_bad]

        model = ConcreteModel()
        model.I = Set(initialize=xref.index)
        model.K = Set(initialize=range(len(x)))
        model.L = Set(initialize=range(len(y)))
        model.M = Set(initialize=range(len(b)))

        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
        model.py = Var(model.L, bounds=(0, None), within=Reals, doc='期望产出的影子价格')
        model.pb = Var(model.M, bounds=(0, None), within=Reals, doc='非期望产出的影子价格')
        model.pomega = Var(bounds=(0, None), within=Reals, doc='SBM效率值（式(13.5)中的 xi）')

        def objective_rule(model):
            return model.pomega * 1

        def first_rule(model, i):
            return (
                -sum(model.px[k] * xref.loc[i, xcol[k]] for k in model.K)
                + sum(model.py[l] * yref.loc[i, ycol[l]] for l in model.L)
                - sum(model.pb[m] * bref.loc[i, bcol[m]] for m in model.M)
                <= 0
            )

        def second_rule(model, k):
            return -model.px[k] <= -1 / (len(model.K) * x.iloc[k])

        def third_rule(model, l):
            return (
                model.pomega / (len(model.L) + len(model.M)) / y.iloc[l] - model.py[l]
                <= 0
            )

        def forth_rule(model, m):
            return (
                model.pomega / (len(model.L) + len(model.M)) / b.iloc[m] - model.pb[m]
                <= 0
            )

        def fifth_rule(model):
            return (
                model.pomega
                + sum(model.px[k] * x.iloc[k] for k in model.K)
                - sum(model.py[l] * y.iloc[l] for l in model.L)
                + sum(model.pb[m] * b.iloc[m] for m in model.M)
                <= 1
            )

        model.obj = Objective(rule=objective_rule, sense=maximize, doc='目标函数')
        model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
        model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
        model.third = Constraint(model.L, rule=third_rule, doc='期望产出上界约束')
        model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出上界约束')
        model.fifth = Constraint(rule=fifth_rule, doc='归一化约束')

        opt = SolverFactory('highs')
        solution = opt.solve(model)

        if solution.solver.termination_condition == "optimal":
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj)

    objdf = pd.DataFrame(obj, index=["obj"]).T
    pxdf = pd.DataFrame(px, index=xcol).T
    pxdf.columns = pxdf.columns.map(lambda x: "shadow price " + str(x))
    pydf = pd.DataFrame(py, index=ycol).T
    pydf.columns = pydf.columns.map(lambda y: "shadow price " + str(y))
    pbdf = pd.DataFrame(pb, index=bcol).T
    pbdf.columns = pbdf.columns.map(lambda b: "shadow price " + str(b))

    pdf = pd.concat([pxdf, pydf], axis=1)
    pdf = pd.concat([pdf, pbdf], axis=1)
    redata = pd.concat([objdf, pdf], axis=1)

    return redata


# Markdown line 645
import pandas as pd
from deabook.example.book_cases.loaders import china, china_lower

ex4 = china(scaled=True)
result = sbm_dual(ex4, ["K", "L", "Y", "CO2"], 2, 1)
print(result.head())


# Markdown line 658
def sbm_dual_select(dataframe, varname, num_inputs, num_outputs,
                    eval_query, ref_query):
    """
    基于SBM对偶问题估计影子价格（支持评价单元和参照集筛选）。

    Parameters
    ----------
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框。
    varname : list[str]
        投入、期望产出和非期望产出变量的名称列表，
        如 ["K", "L", "Y", "CO2"]。
    num_inputs : int
        投入变量的个数。
    num_outputs : int
        期望产出变量的个数。
    eval_query : str
        传入数据框.query()方法的参数，用于筛选待评价DMU，
        如 "dmu==1"、"dmu==[1,2,3]"。
    ref_query : str
        传入数据框.query()方法的参数，用于筛选参照集DMU，
        如 "dmu==1"、"dmu==[1,2,3]"。

    Returns
    -------
    pd.DataFrame
        包含各DMU影子价格估计结果的数据框。
    """
    index_list = dataframe.query(eval_query).index
    index_list_ref = dataframe.query(ref_query).index

    data = dataframe.loc[index_list, varname]
    dataref = dataframe.loc[index_list_ref, varname]

    num_inputs_ = num_inputs
    num_outputs_ = num_outputs
    num_bad = data.shape[1] - num_inputs_ - num_outputs_

    xref = dataref.iloc[:, 0:num_inputs_]
    yref = dataref.iloc[:, num_inputs_:num_inputs_+num_outputs_]
    bref = dataref.iloc[:, num_inputs_+num_outputs_:num_inputs_+num_outputs_+
        num_bad]

    xcol = varname[0:num_inputs_]
    ycol = varname[num_inputs_:num_inputs_+num_outputs_]
    bcol = varname[num_inputs_+num_outputs_:num_inputs_+num_outputs_+num_bad]

    px, py, pb, obj = {}, {}, {}, {}

    for j in data.index:
        x = data.loc[j, xcol]
        y = data.loc[j, ycol]
        b = data.loc[j, bcol]

        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)
        model.K = Set(initialize=range(len(x)))
        model.L = Set(initialize=range(len(y)))
        model.M = Set(initialize=range(len(b)))

        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
        model.py = Var(model.L, bounds=(0, None), within=Reals,
            doc='期望产出的影子价格')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,
            doc='非期望产出的影子价格')
        model.pomega = Var(bounds=(0, None), within=Reals,
            doc='SBM效率值（式(13.5)中的 xi）')

        def objective_rule(model):
            return model.pomega * 1

        def first_rule(model, i):
            return (
                -sum(model.px[k] * xref.loc[i, xcol[k]] for k in model.K)
                + sum(model.py[l] * yref.loc[i, ycol[l]] for l in model.L)
                - sum(model.pb[m] * bref.loc[i, bcol[m]] for m in model.M)
                <= 0
            )

        def second_rule(model, k):
            return -model.px[k] <= -1 / (len(model.K) * x.iloc[k])

        def third_rule(model, l):
            return (
                model.pomega / (len(model.L) + len(model.M)) / y.iloc[
                    l] - model.py[l]
                <= 0
            )

        def forth_rule(model, m):
            return (
                model.pomega / (len(model.L) + len(model.M)) / b.iloc[
                    m] - model.pb[m]
                <= 0
            )

        def fifth_rule(model):
            return (
                model.pomega
                + sum(model.px[k] * x.iloc[k] for k in model.K)
                - sum(model.py[l] * y.iloc[l] for l in model.L)
                + sum(model.pb[m] * b.iloc[m] for m in model.M)
                <= 1
            )

        model.obj = Objective(rule=objective_rule, sense=maximize, doc='目标函数')
        model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
        model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
        model.third = Constraint(model.L, rule=third_rule, doc='期望产出上界约束')
        model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出上界约束')
        model.fifth = Constraint(rule=fifth_rule, doc='归一化约束')

        opt = SolverFactory('highs')
        solution = opt.solve(model)

        if solution.solver.termination_condition == "optimal":
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj)

    objdf = pd.DataFrame(obj, index=["obj"]).T
    pxdf = pd.DataFrame(px, index=xcol).T
    pxdf.columns = pxdf.columns.map(lambda x: "shadow price " + str(x))
    pydf = pd.DataFrame(py, index=ycol).T
    pydf.columns = pydf.columns.map(lambda y: "shadow price " + str(y))
    pbdf = pd.DataFrame(pb, index=bcol).T
    pbdf.columns = pbdf.columns.map(lambda b: "shadow price " + str(b))

    pdf = pd.concat([pxdf, pydf], axis=1)
    pdf = pd.concat([pdf, pbdf], axis=1)
    redata = pd.concat([objdf, pdf], axis=1)

    return redata


# Markdown line 803
def sbm_dual_formula(formula, dataframe, eval_query=None, ref_query=None):
    """
    基于SBM对偶问题估计影子价格（支持公式字符串调用）。

    Parameters
    ----------
    formula : str
        公式字符串，格式为 "Y : CO2 = K L"。
        等号左边为产出变量（冒号前为期望产出，冒号后为非期望产出），
        等号右边为投入变量。
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框。
    eval_query : str or None, default None
        传入数据框.query()方法的参数，用于筛选待评价DMU。
        若为None，则评价所有DMU。
    ref_query : str or None, default None
        传入数据框.query()方法的参数，用于筛选参照集DMU。
        若为None，则使用全部DMU作为参照集。

    Returns
    -------
    pd.DataFrame
        包含各DMU影子价格估计结果的数据框。
    """
    px, py, pb, obj = {}, {}, {}, {}

    if eval_query is None:
        index_list = dataframe.index
    else:
        index_list = dataframe.query(eval_query).index

    if ref_query is None:
        index_list_ref = dataframe.index
    else:
        index_list_ref = dataframe.query(ref_query).index

    # 解析公式字符串
    input_vars = formula.split('=')[1].strip(' ')
    xcol = re.compile(' +').sub(' ', input_vars).split(' ')
    output_vars = formula.split('=')[0].split(':')[0].strip(' ')
    ycol = re.compile(' +').sub(' ', output_vars).split(' ')
    unoutput_vars = formula.split('=')[0].split(':')[1].strip(' ')
    bcol = re.compile(' +').sub(' ', unoutput_vars).split(' ')

    data = dataframe.loc[index_list, xcol + ycol + bcol]
    dataref = dataframe.loc[index_list_ref, xcol + ycol + bcol]

    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    bref = dataref.loc[:, bcol]

    for j in data.index:
        x = data.loc[j, xcol]
        y = data.loc[j, ycol]
        b = data.loc[j, bcol]

        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)
        model.K = Set(initialize=range(len(x)))
        model.L = Set(initialize=range(len(y)))
        model.M = Set(initialize=range(len(b)))

        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
        model.py = Var(model.L, bounds=(0, None), within=Reals,
            doc='期望产出的影子价格')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,
            doc='非期望产出的影子价格')
        model.pomega = Var(bounds=(0, None), within=Reals,
            doc='SBM效率值（式(13.5)中的 xi）')

        def objective_rule(model):
            return model.pomega * 1

        def first_rule(model, i):
            return (
                -sum(model.px[k] * xref.loc[i, xcol[k]] for k in model.K)
                + sum(model.py[l] * yref.loc[i, ycol[l]] for l in model.L)
                - sum(model.pb[m] * bref.loc[i, bcol[m]] for m in model.M)
                <= 0
            )

        def second_rule(model, k):
            return -model.px[k] <= -1 / (len(model.K) * x.iloc[k])

        def third_rule(model, l):
            return (
                model.pomega / (len(model.L) + len(model.M)) / y.iloc[
                    l] - model.py[l]
                <= 0
            )

        def forth_rule(model, m):
            return (
                model.pomega / (len(model.L) + len(model.M)) / b.iloc[
                    m] - model.pb[m]
                <= 0
            )

        def fifth_rule(model):
            return (
                model.pomega
                + sum(model.px[k] * x.iloc[k] for k in model.K)
                - sum(model.py[l] * y.iloc[l] for l in model.L)
                + sum(model.pb[m] * b.iloc[m] for m in model.M)
                <= 1
            )

        model.obj = Objective(rule=objective_rule, sense=maximize, doc='目标函数')
        model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
        model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
        model.third = Constraint(model.L, rule=third_rule, doc='期望产出上界约束')
        model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出上界约束')
        model.fifth = Constraint(rule=fifth_rule, doc='归一化约束')

        opt = SolverFactory('highs')
        solution = opt.solve(model)

        if solution.solver.termination_condition == "optimal":
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj)

    objdf = pd.DataFrame(obj, index=["obj"]).T
    pxdf = pd.DataFrame(px, index=xcol).T
    pxdf.columns = pxdf.columns.map(lambda x: "shadow price " + str(x))
    pydf = pd.DataFrame(py, index=ycol).T
    pydf.columns = pydf.columns.map(lambda y: "shadow price " + str(y))
    pbdf = pd.DataFrame(pb, index=bcol).T
    pbdf.columns = pbdf.columns.map(lambda b: "shadow price " + str(b))

    pdf = pd.concat([pxdf, pydf], axis=1)
    pdf = pd.concat([pdf, pbdf], axis=1)
    redata = pd.concat([objdf, pdf], axis=1)

    return redata


# Markdown line 944
import pandas as pd
from deabook.example.book_cases.loaders import china, china_lower

ex4 = china(scaled=True)
result = sbm_dual_formula("Y : CO2 = K L", ex4, "t==[1,2]", "t==[1,2,3]")
print(result)


# Markdown line 1081
from pyomo.environ import (
    ConcreteModel, Set, Var, Objective, Constraint, Reals,
    minimize, SolverFactory, value
)
import pandas as pd
from deabook.example.book_cases.loaders import china, china_lower
import numpy as np
import warnings

warnings.filterwarnings("ignore")

# 读取数据
ex5 = china_lower(scaled=True)
data = ex5[['labor', 'capital', 'energy', 'gdp', 'co2']]

# 设定变量个数
K = 3  # 投入变量个数（labor, capital, energy）
L = 1  # 期望产出变量个数（gdp）
M = 1  # 非期望产出变量个数（co2）

# 准备参照集数据
xref = data.iloc[:, 0:K]
yref = data.iloc[:, K:K+L]
bref = data.iloc[:, K+L:K+L+M]

# 准备被评价DMU数据（以第1个DMU为例）
x = data.iloc[0, 0:K]
y = data.iloc[0, K:K+L]
b = data.iloc[0, K+L:K+L+M]

# 变量名
xcol = xref.columns
ycol = yref.columns
bcol = bref.columns

# 设定方向向量（默认：投入减少、产出增加、非期望产出减少）
gx = [-1] * len(xcol)
gy = [1] * len(ycol)
gb = [-1] * len(bcol)

# 计算权重（默认等权重）
weight = []
denominator = int(gx[0] != 0) + int(gy[0] != 0) + int(gb[0] != 0)
for _ in range(len(xcol)):
    weight.append(1.0 / denominator / len(xcol))
for _ in range(len(ycol)):
    weight.append(1.0 / denominator / len(ycol))
for _ in range(len(bcol)):
    weight.append(1.0 / denominator / len(bcol))

iweight = weight[0:len(xcol)]
oweight = weight[len(xcol):len(xcol)+len(ycol)]
bweight = weight[len(xcol)+len(ycol):len(xcol)+len(ycol)+len(bcol)]

# 构建Pyomo模型
model = ConcreteModel()
model.I = Set(initialize=xref.index)             # 参照决策单元集合
model.K = Set(initialize=range(len(xcol)))       # 投入变量集合
model.L = Set(initialize=range(len(ycol)))       # 期望产出变量集合
model.M = Set(initialize=range(len(bcol)))       # 非期望产出变量集合

# 定义影子价格变量（施加非负约束）
model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
model.py = Var(model.L, bounds=(0, None), within=Reals, doc='期望产出的影子价格')
model.pb = Var(model.M, bounds=(0, None), within=Reals, doc='非期望产出的影子价格')
model.pomega = Var(bounds=(0, None), within=Reals, doc='方向距离函数值')


# 目标函数：最小化虚拟利润
def objective_rule(model):
    """返回目标函数表达式（虚拟利润）"""
    return (
        sum(model.px[k] * x.iloc[k] for k in model.K)
        - sum(model.py[l] * y.iloc[l] for l in model.L)
        + sum(model.pb[j] * b.iloc[j] for j in model.M)
        + model.pomega * 1
    )

# 约束1：对偶约束（参照集加权影子价格之和非负）
def first_rule(model, i):
    """第i个参照DMU的对偶约束"""
    return (
        sum(model.px[k] * xref.loc[i, xcol[k]] for k in model.K)
        - sum(model.py[l] * yref.loc[i, ycol[l]] for l in model.L)
        + sum(model.pb[m] * bref.loc[i, bcol[m]] for m in model.M)
        + model.pomega * 1
        >= 0
    )

# 约束2：投入影子价格的下界约束
def second_rule(model, k):
    """第k个投入变量的影子价格下界"""
    return -gx[k] * x.iloc[k] * model.px[k] >= iweight[k]

# 约束3：期望产出影子价格的下界约束
def third_rule(model, l):
    """第l个期望产出变量的影子价格下界"""
    return gy[l] * y.iloc[l] * model.py[l] >= oweight[l]

# 约束4：非期望产出影子价格的下界约束
def forth_rule(model, j):
    """第j个非期望产出变量的影子价格下界"""
    return -gb[j] * b.iloc[j] * model.pb[j] >= bweight[j]

# 组装模型
model.obj = Objective(rule=objective_rule, sense=minimize, doc='目标函数')
model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
model.third = Constraint(model.L, rule=third_rule, doc='期望产出下界约束')
model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出下界约束')

# 求解
opt = SolverFactory('highs')
solution = opt.solve(model)

# 提取结果
px = np.asarray(list(model.px[:].value))
py = np.asarray(list(model.py[:].value))
pb = np.asarray(list(model.pb[:].value))
obj = value(model.obj)

print(f"投入影子价格 (px): {px}")
print(f"期望产出影子价格 (py): {py}")
print(f"非期望产出影子价格 (pb): {pb}")
print(f"最优目标值: {obj}")


# Markdown line 1215
px, py, pb, obj = {}, {}, {}, {}

for j in range(data.shape[0]):
    x = data.iloc[j, 0:K]
    y = data.iloc[j, K:K+L]
    b = data.iloc[j, K+L:K+L+M]

    model = ConcreteModel()
    model.I = Set(initialize=xref.index)
    model.K = Set(initialize=range(len(xcol)))
    model.L = Set(initialize=range(len(ycol)))
    model.M = Set(initialize=range(len(bcol)))

    model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
    model.py = Var(model.L, bounds=(0, None), within=Reals, doc='期望产出的影子价格')
    model.pb = Var(model.M, bounds=(0, None), within=Reals, doc='非期望产出的影子价格')
    model.pomega = Var(bounds=(0, None), within=Reals, doc='方向距离函数值')

    def objective_rule(model):
        return (
            sum(model.px[k] * x.iloc[k] for k in model.K)
            - sum(model.py[l] * y.iloc[l] for l in model.L)
            + sum(model.pb[jj] * b.iloc[jj] for jj in model.M)
            + model.pomega * 1
        )

    def first_rule(model, i):
        return (
            sum(model.px[k] * xref.loc[i, xcol[k]] for k in model.K)
            - sum(model.py[l] * yref.loc[i, ycol[l]] for l in model.L)
            + sum(model.pb[m] * bref.loc[i, bcol[m]] for m in model.M)
            + model.pomega * 1
            >= 0
        )

    def second_rule(model, k):
        return -gx[k] * x.iloc[k] * model.px[k] >= iweight[k]

    def third_rule(model, l):
        return gy[l] * y.iloc[l] * model.py[l] >= oweight[l]

    def forth_rule(model, jj):
        return -gb[jj] * b.iloc[jj] * model.pb[jj] >= bweight[jj]

    model.obj = Objective(rule=objective_rule, sense=minimize, doc='目标函数')
    model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
    model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
    model.third = Constraint(model.L, rule=third_rule, doc='期望产出下界约束')
    model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出下界约束')

    opt = SolverFactory('highs')
    solution = opt.solve(model)

    if solution.solver.termination_condition == "optimal":
        px[j] = np.asarray(list(model.px[:].value))
        py[j] = np.asarray(list(model.py[:].value))
        pb[j] = np.asarray(list(model.pb[:].value))
        obj[j] = value(model.obj)


# Markdown line 1280
def nddf_dual(dataframe, varname, num_inputs, num_outputs):
    """
    基于NDDF模型对偶问题估计影子价格。

    求解NDDF模型的对偶问题（乘数形式），获得各决策单元的投入、
    期望产出和非期望产出的影子价格。

    Parameters
    ----------
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框，变量按投入在前、
        期望产出在中、非期望产出在后的顺序排列。
    varname : list[str]
        投入变量、期望产出变量和非期望产出变量的名称列表，
        如 ["labor", "capital", "energy", "gdp", "co2"]。
    num_inputs : int
        投入变量的个数。
    num_outputs : int
        期望产出变量的个数。

    Returns
    -------
    pd.DataFrame
        包含各DMU的影子价格估计结果的数据框。

    Notes
    -----
    - 默认方向向量：投入减少(-1)、期望产出增加(+1)、非期望产出减少(-1)。
    - 默认权重：等权重分配。
    - 如需自定义方向向量和权重，请使用 nddf_dual_advanced 函数。
    """
    data = dataframe.loc[:, varname]
    num_inputs_ = num_inputs
    num_outputs_ = num_outputs
    num_bad = data.shape[1] - num_inputs_ - num_outputs_

    xref = data.iloc[:, 0:num_inputs_]
    yref = data.iloc[:, num_inputs_:num_inputs_+num_outputs_]
    bref = data.iloc[:, num_inputs_+num_outputs_:num_inputs_+num_outputs_+
        num_bad]

    xcol = xref.columns
    ycol = yref.columns
    bcol = bref.columns

    # 默认方向向量
    gx = [-1] * len(xcol)
    gy = [1] * len(ycol)
    gb = [-1] * len(bcol)

    # 默认权重
    weight = []
    denominator = int(gx[0] != 0) + int(gy[0] != 0) + int(gb[0] != 0)
    for _ in range(len(xcol)):
        weight.append(1.0 / denominator / len(xcol))
    for _ in range(len(ycol)):
        weight.append(1.0 / denominator / len(ycol))
    for _ in range(len(bcol)):
        weight.append(1.0 / denominator / len(bcol))

    iweight = weight[0:len(xcol)]
    oweight = weight[len(xcol):len(xcol)+len(ycol)]
    bweight = weight[len(xcol)+len(ycol):len(xcol)+len(ycol)+len(bcol)]

    px, py, pb, obj = {}, {}, {}, {}

    for j in range(data.shape[0]):
        x = data.iloc[j, 0:num_inputs_]
        y = data.iloc[j, num_inputs_:num_inputs_+num_outputs_]
        b = data.iloc[j, num_inputs_+num_outputs_:num_inputs_+num_outputs_+
            num_bad]

        model = ConcreteModel()
        model.I = Set(initialize=xref.index)
        model.K = Set(initialize=range(len(xcol)))
        model.L = Set(initialize=range(len(ycol)))
        model.M = Set(initialize=range(len(bcol)))

        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
        model.py = Var(model.L, bounds=(0, None), within=Reals,
            doc='期望产出的影子价格')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,
            doc='非期望产出的影子价格')
        model.pomega = Var(bounds=(0, None), within=Reals, doc='方向距离函数值')

        def objective_rule(model):
            return (
                sum(model.px[k] * x.iloc[k] for k in model.K)
                - sum(model.py[l] * y.iloc[l] for l in model.L)
                + sum(model.pb[jj] * b.iloc[jj] for jj in model.M)
                + model.pomega * 1
            )

        def first_rule(model, i):
            return (
                sum(model.px[k] * xref.loc[i, xcol[k]] for k in model.K)
                - sum(model.py[l] * yref.loc[i, ycol[l]] for l in model.L)
                + sum(model.pb[m] * bref.loc[i, bcol[m]] for m in model.M)
                + model.pomega * 1
                >= 0
            )

        def second_rule(model, k):
            return -gx[k] * x.iloc[k] * model.px[k] >= iweight[k]


        def third_rule(model, l):
            return gy[l] * y.iloc[l] * model.py[l] >= oweight[l]

        def forth_rule(model, jj):
            return -gb[jj] * b.iloc[jj] * model.pb[jj] >= bweight[jj]

        model.obj = Objective(rule=objective_rule, sense=minimize, doc='目标函数')
        model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
        model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
        model.third = Constraint(model.L, rule=third_rule, doc='期望产出下界约束')
        model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出下界约束')

        opt = SolverFactory('highs')
        solution = opt.solve(model)

        if solution.solver.termination_condition == "optimal":
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj)

    objdf = pd.DataFrame(obj, index=["obj"]).T
    pxdf = pd.DataFrame(px, index=xcol).T
    pxdf.columns = pxdf.columns.map(lambda x: "shadow price " + str(x))
    pydf = pd.DataFrame(py, index=ycol).T
    pydf.columns = pydf.columns.map(lambda y: "shadow price " + str(y))
    pbdf = pd.DataFrame(pb, index=bcol).T
    pbdf.columns = pbdf.columns.map(lambda b: "shadow price " + str(b))

    pdf = pd.concat([pxdf, pydf], axis=1)
    pdf = pd.concat([pdf, pbdf], axis=1)
    redata = pd.concat([objdf, pdf], axis=1)

    return redata


# Markdown line 1427
def nddf_dual_advanced(formula, dataframe, gx=None, gy=None, gb=None,
                       weight=None, eval_query=None, ref_query=None):
    """
    基于NDDF对偶问题估计影子价格（高级版本）。

    支持自定义方向向量、权重，以及评价单元和参照集的筛选。

    Parameters
    ----------
    formula : str
        公式字符串，格式为 "期望产出:非期望产出=投入1 投入2 ..."，
        如 "gdp : co2 = labor capital energy"。
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框。
    gx : list or None, default None
        投入方向向量。默认为 [-1, -1, ...]，表示投入减少方向。
    gy : list or None, default None
        期望产出方向向量。默认为 [1, 1, ...]，表示期望产出增加方向。
    gb : list or None, default None
        非期望产出方向向量。默认为 [-1, -1, ...]，表示非期望产出减少方向。
    weight : list or None, default None
        各变量的权重向量。若为None，则使用等权重分配。
    eval_query : str or None, default None
        传入数据框.query()方法的参数，用于筛选待评价DMU。
        若为None，则评价所有DMU。
    ref_query : str or None, default None
        传入数据框.query()方法的参数，用于筛选参照集DMU。
        若为None，则使用全部DMU作为参照集。

    Returns
    -------
    pd.DataFrame
        包含各DMU影子价格估计结果的数据框。
    """
    px, py, pb, obj = {}, {}, {}, {}

    # 处理查询条件
    if eval_query is None:
        index_list = dataframe.index
    else:
        index_list = dataframe.query(eval_query).index

    if ref_query is None:
        index_list_ref = dataframe.index
    else:
        index_list_ref = dataframe.query(ref_query).index

    # 解析公式字符串
    input_vars = formula.split('=')[1].strip(' ')
    xcol = re.compile(' +').sub(' ', input_vars).split(' ')
    output_vars = formula.split('=')[0].split(':')[0].strip(' ')
    ycol = re.compile(' +').sub(' ', output_vars).split(' ')
    unoutput_vars = formula.split('=')[0].split(':')[1].strip(' ')
    bcol = re.compile(' +').sub(' ', unoutput_vars).split(' ')

    data = dataframe.loc[index_list, xcol + ycol + bcol]
    dataref = dataframe.loc[index_list_ref, xcol + ycol + bcol]

    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    bref = dataref.loc[:, bcol]

    # 处理方向向量（使用None判断默认值）
    if gx is None:
        gx = [-1] * len(xcol)
    if gy is None:
        gy = [1] * len(ycol)
    if gb is None:
        gb = [-1] * len(bcol)

    # 处理权重（使用None判断默认值）
    if weight is None:
        weight = []
        denominator = int(gx[0] != 0) + int(gy[0] != 0) + int(gb[0] != 0)
        for _ in range(len(xcol)):
            weight.append(1.0 / denominator / len(xcol))
        for _ in range(len(ycol)):
            weight.append(1.0 / denominator / len(ycol))
        for _ in range(len(bcol)):
            weight.append(1.0 / denominator / len(bcol))

    iweight = weight[0:len(xcol)]
    oweight = weight[len(xcol):len(xcol)+len(ycol)]
    bweight = weight[len(xcol)+len(ycol):len(xcol)+len(ycol)+len(bcol)]

    for j in data.index:
        x = data.loc[j, xcol]
        y = data.loc[j, ycol]
        b = data.loc[j, bcol]

        model = ConcreteModel()
        model.I = Set(initialize=xref.index)
        model.K = Set(initialize=range(len(xcol)))
        model.L = Set(initialize=range(len(ycol)))
        model.M = Set(initialize=range(len(bcol)))

        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='投入的影子价格')
        model.py = Var(model.L, bounds=(0, None), within=Reals,
            doc='期望产出的影子价格')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,
            doc='非期望产出的影子价格')
        model.pomega = Var(bounds=(0, None), within=Reals, doc='方向距离函数值')

        def objective_rule(model):
            return (
                sum(model.px[k] * x.iloc[k] for k in model.K)
                - sum(model.py[l] * y.iloc[l] for l in model.L)
                + sum(model.pb[jj] * b.iloc[jj] for jj in model.M)
                + model.pomega * 1
            )

        def first_rule(model, i):
            return (
                sum(model.px[k] * xref.loc[i, xcol[k]] for k in model.K)
                - sum(model.py[l] * yref.loc[i, ycol[l]] for l in model.L)
                + sum(model.pb[m] * bref.loc[i, bcol[m]] for m in model.M)
                + model.pomega * 1
                >= 0
            )

        def second_rule(model, k):
            return -gx[k] * x.iloc[k] * model.px[k] >= iweight[k]

        def third_rule(model, l):
            return gy[l] * y.iloc[l] * model.py[l] >= oweight[l]

        def forth_rule(model, jj):
            return -gb[jj] * b.iloc[jj] * model.pb[jj] >= bweight[jj]

        model.obj = Objective(rule=objective_rule, sense=minimize, doc='目标函数')
        model.first = Constraint(model.I, rule=first_rule, doc='对偶约束')
        model.second = Constraint(model.K, rule=second_rule, doc='投入下界约束')
        model.third = Constraint(model.L, rule=third_rule, doc='期望产出下界约束')
        model.forth = Constraint(model.M, rule=forth_rule, doc='非期望产出下界约束')

        opt = SolverFactory('highs')
        solution = opt.solve(model)

        if solution.solver.termination_condition == "optimal":
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj)

    objdf = pd.DataFrame(obj, index=["obj"]).T
    pxdf = pd.DataFrame(px, index=xcol).T
    pxdf.columns = pxdf.columns.map(lambda x: "shadow price " + str(x))
    pydf = pd.DataFrame(py, index=ycol).T
    pydf.columns = pydf.columns.map(lambda y: "shadow price " + str(y))
    pbdf = pd.DataFrame(pb, index=bcol).T
    pbdf.columns = pbdf.columns.map(lambda b: "shadow price " + str(b))

    pdf = pd.concat([pxdf, pydf], axis=1)
    pdf = pd.concat([pdf, pbdf], axis=1)
    redata = pd.concat([objdf, pdf], axis=1)

    return redata


# Markdown line 1591
import pandas as pd
from deabook.example.book_cases.loaders import china, china_lower

ex5 = china_lower(scaled=True)
result = nddf_dual_advanced(
    'gdp : co2 = labor capital energy',
    ex5,
    gx=[-1, -1, -1],
    gy=[1],
    gb=[-1],
    weight=None,
    eval_query=None,
    ref_query=None
)
print(result.head())


# Markdown line 1691
from deabook.example.book_cases.run import FUN
FUN[13]()  # 重新计算本章，输出完整CSV、摘要表与PNG/PDF图


# Markdown line 1729
import numpy as np
import pandas as pd

def resample_shadow_summary(dataframe, varname, num_inputs, num_outputs,
                            model_func, b=100, alpha=0.05, seed=42):
    if b < 2 or not 0 < alpha < 1 or dataframe.empty:
        raise ValueError("需至少两次重复、非空数据及合法alpha")
    rng = np.random.default_rng(seed)
    original = model_func(dataframe, varname, num_inputs, num_outputs)
    cols = [c for c in original if c.startswith('shadow price ')]
    estimates = []
    failures = 0
    for _ in range(b):
        sample = dataframe.iloc[rng.integers(len(dataframe), size=len(dataframe))]
        try:
            result = model_func(sample.reset_index(drop=True), varname,
                                num_inputs, num_outputs)
            if not np.isfinite(result[cols].to_numpy()).all():
                raise RuntimeError("重抽样结果含非有限值")
            estimates.append(result[cols].median())
        except RuntimeError:
            failures += 1
    if len(estimates) < 2:
        raise RuntimeError("有效重复不足，无法计算分位数")
    draws = pd.DataFrame(estimates)
    return pd.DataFrame({
        'sample_median': original[cols].median(),
        'resample_q_low': draws.quantile(alpha / 2),
        'resample_q_high': draws.quantile(1 - alpha / 2),
        'successful': len(estimates), 'failed': failures
    })

# 小规模教学运行；正式敏感性分析可提高重复次数。
resampling_result = resample_shadow_summary(
    china(scaled=True).query('t==4'), ['K', 'L', 'Y', 'CO2'],
    2, 1, sbm_dual, b=5)
print(resampling_result)


# Markdown line 1787
import numpy as np


def scale_for_shadow_price(data, input_cols, output_cols, bad_cols,
                           method="mean"):
    """
    按列缩放投入产出数据，以提高数值稳定性。

    只做缩放、不做中心化：SBM效率对列缩放严格不变，对平移则不然，
    因此不能使用z-score标准化（会把投入产出变成非正数）。

    Parameters
    ----------
    data : pd.DataFrame
        原始数据。
    input_cols, output_cols, bad_cols : list[str]
        投入、期望产出、非期望产出的列名。
    method : {"mean", "std", "max"}, default "mean"
        缩放基准：列均值、列标准差或列最大值，三者均为正常数。

    Returns
    -------
    scaled_data : pd.DataFrame
        缩放后的数据，各列量纲统一到1附近。
    scaler_params : dict
        缩放参数，键 'scale' 为各列除数（用于将影子价格反变换回
        原始尺度：price_original = price_scaled / scale）。
    """
    all_cols = input_cols + output_cols + bad_cols

    if method == "mean":
        scale = data[all_cols].mean()
    elif method == "std":
        scale = data[all_cols].std(ddof=0)
    elif method == "max":
        scale = data[all_cols].max()
    else:
        raise ValueError("method 只能取 'mean'、'std' 或 'max'")

    assert (scale > 0).all(), "缩放基准必须为正，请检查数据是否含非正列"

    scaled_data = data.copy()
    scaled_data[all_cols] = data[all_cols] / scale

    scaler_params = {
        'scale': scale,
        'input_cols': input_cols,
        'output_cols': output_cols,
        'bad_cols': bad_cols
    }

    return scaled_data, scaler_params


# Markdown line 1856
opt = SolverFactory('highs')

# 设置严格容差
opt.options['primal_feasibility_tolerance'] = 1e-9   # 原始可行性容差
opt.options['dual_feasibility_tolerance'] = 1e-9   # 对偶可行性容差
opt.options['ipm_optimality_tolerance'] = 1e-9  # 内点法最优性容差

solution = opt.solve(model)


# Markdown line 1871
def check_degeneracy(model, tol=1e-8):
    """
    检查DEA对偶问题是否存在退化。

    通过检查对偶约束的松弛变量是否接近零来判断退化。
    如果多个约束同时紧致，可能存在退化。

    Parameters
    ----------
    model : pyomo.environ.ConcreteModel
        已求解的Pyomo模型。
    tol : float, default 1e-8
        判断约束是否紧致的容差。

    Returns
    -------
    bool
        是否存在退化。
    """
    tight_count = 0
    for constr in model.first.values():
        slack = constr.uslack()  # 上界松弛
        if abs(slack) < tol:
            tight_count += 1

    # 如果超过维度数的约束同时紧致，可能存在退化
    is_degenerate = tight_count > len(model.K)
    return is_degenerate

