# 本文件从教材原文提取；在项目根目录运行。

# Markdown line 69
import numpy as np
from pyomo.environ import *
model = ConcreteModel()
model.I = Set(initialize = range(5)) # 采用range初始化技术参照决策单元个数的集合
model.K = Set(initialize = range(3)) # 采用range初始化投入变量个数的集合
model.L = Set(initialize = range(1)) # 采用range初始化产出变量个数的集合


# Markdown line 80
print(model.I.data())
print(model.K.data())
print(model.L.data())


# Markdown line 109
model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
model.intensity = Var(
    model.I, bounds=(0.0, None), doc='intensity variables'
)  # 强度变量（原lambda）


# Markdown line 122
def obj_rule(model):
    return model.theta
    # return model.theta*1 + sum(model.intensity[i] * 0 for i in model.I)
model.obj = Objective(rule=obj_rule, sense=maximize, doc='objective function')


# Markdown line 133
from deabook.example.book_cases.loaders import plants
data = plants().head(5)[["Labor", "Fuel", "Capital", "Output"]].to_numpy()
x = data[0, 0:3]      # 被评价电厂Coffeen的投入
y = data[0, 3:4]      # 被评价电厂Coffeen的产出
xref = data[:, 0:3]   # 所有参照决策单元的投入
yref = data[:, 3:4]   # 所有参照决策单元的产出

def input_rule(model, k):
    """投入约束：参照单元的加权投入不超过被评价单元的投入"""
    return (
        sum(model.intensity[i] * xref[i, k] for i in model.I)
        <= x[k]
    )

def output_rule(model, l):
    """产出约束：参照单元的加权产出不低于theta倍的被评价单元产出"""
    return (
        sum(model.intensity[i] * yref[i, l] for i in model.I)
        >= model.theta * y[l]
    )

model.input = Constraint(model.K, rule=input_rule, doc='input constraint')
model.output = Constraint(model.L, rule=output_rule, doc='output constraint')


# Markdown line 161
model2 = ConcreteModel()
model2.I = Set(initialize = range(4))
model2.J = Set(initialize = range(5))
model2.x = Var(model2.I, model2.J, within = Reals)
model2.c1 = Constraint(
    expr=2.0 <= model2.x[0, 0] + 4 * model2.x[0, 2]
)      # 不等式约束 2 <= x[0,0] + 4*x[0,2]
model2.c2 = Constraint(
    expr=3 * model2.x[2, 3] + 5 * model2.x[3, 4] == 10.0
)    # 等式约束 3*x[2,3] + 5*x[3,4] = 10


# Markdown line 176
def c3_rule(model2, i):
    if i == 1:
        return Constraint.Skip
    return (
        sum([model2.x[i, j] for j in model2.J])
        <= 2
    )
model2.c3 = Constraint(model2.I, rule=c3_rule)


# Markdown line 191
opt = SolverFactory('highs')       # 指定 HiGHS 作为求解器
solution = opt.solve(model)       # 调用求解器求解


# Markdown line 198
opt = SolverFactory('highs')
solution = opt.solve(model, report_timing=True)


# Markdown line 205
opt = SolverFactory('highs')
assert opt.available(exception_flag=False), "请安装 highspy"
solution = opt.solve(model)


# Markdown line 213
# 验证求解器是否可用
from pyomo.environ import SolverFactory
for solver in ['highs', 'glpk', 'mosek']:
    opt = SolverFactory(solver)
    print(f"{solver}: {'可用' if opt.available(exception_flag=False) else '不可用'}")


# Markdown line 223
theta = np.asarray(list(model.theta[:].value))          # 提取决策变量 theta（效率倒数）
intensity_vals = np.asarray(list(model.intensity[:].value))  # 提取强度变量 lambda
obj = value(model.obj)                                   # 提取目标函数值
print("optimum theta: \n {} ".format(theta))
print("optimum intensity: \n {} ".format(intensity_vals))
print("optimal objective: {}".format(obj))


# Markdown line 236
from pyomo.environ import *
import numpy as np

from deabook.example.book_cases.loaders import plants
data = plants().head(5)[["Labor", "Fuel", "Capital", "Output"]].to_numpy()

x = data[0, 0:3]       # 被评价单元A的投入
y = data[0, 3:4]       # 被评价单元A的产出
xref = data[:, 0:3]    # 参照单元投入
yref = data[:, 3:4]    # 参照单元产出

model = ConcreteModel()
model.I = Set(initialize=range(5))    # 技术参照决策单元集合
model.K = Set(initialize=range(3))    # 投入变量集合
model.L = Set(initialize=range(1))    # 产出变量集合

model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
model.intensity = Var(model.I, bounds=(0.0, None), doc='intensity variables')

def obj_rule(model):
    return model.theta
model.obj = Objective(rule=obj_rule, sense=maximize, doc='objective function')

def input_rule(model, k):
    return (
        sum(model.intensity[i] * xref[i, k] for i in model.I)
        <= x[k]
    )
def output_rule(model, l):
    return (
        sum(model.intensity[i] * yref[i, l] for i in model.I)
        >= model.theta * y[l]
    )

model.input = Constraint(model.K, rule=input_rule, doc='input constraint')
model.output = Constraint(model.L, rule=output_rule, doc='output constraint')

opt = SolverFactory('highs')
solution = opt.solve(model)

theta = np.asarray(list(model.theta[:].value))
intensity_vals = np.asarray(list(model.intensity[:].value))
obj = value(model.obj)
print("optimum theta: \n {} ".format(theta))
print("optimum intensity: \n {} ".format(intensity_vals))
print("optimal objective: {}".format(obj))


# Markdown line 289
from pyomo.environ import *
import numpy as np

from deabook.example.book_cases.loaders import plants
data = plants().head(5)[["Labor", "Fuel", "Capital", "Output"]].to_numpy()

theta_results = np.empty(5)    # 定义 theta_results 用于存储计算结果（原thetalt）

for j in range(5):             # 在 data 的行上循环
    x = data[j, 0:3]
    y = data[j, 3:4]
    xref = data[:, 0:3]
    yref = data[:, 3:4]
    
    model = ConcreteModel()
    model.I = Set(initialize=range(5))
    model.K = Set(initialize=range(3))
    model.L = Set(initialize=range(1))
    
    model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
    model.intensity = Var(model.I, bounds=(0.0, None),
        doc='intensity variables')
    
    def obj_rule(model):
        return model.theta
    model.obj = Objective(rule=obj_rule, sense=maximize,
        doc='objective function')
    
    def input_rule(model, k):
        return (
            sum(model.intensity[i] * xref[i, k] for i in model.I)
            <= x[k]
        )
    def output_rule(model, l):
        return (
            sum(model.intensity[i] * yref[i, l] for i in model.I)
            >= model.theta * y[l]
        )
    
    model.input = Constraint(model.K, rule=input_rule, doc='input constraint')
    model.output = Constraint(model.L, rule=output_rule,
        doc='output constraint')
    
    opt = SolverFactory('highs')
    solution = opt.solve(model)
    
    theta = np.asarray(list(model.theta[:].value))
    intensity_vals = np.asarray(list(model.intensity[:].value))
    obj = value(model.obj)
    
    theta_results[j] = theta[0]

te = 1 / theta_results
print(te)


# Markdown line 358
from pyomo.environ import *
import numpy as np

def dea1(data, dataref, numk):
    """
    计算DEA技术效率（基础版本）
    
    Parameters
    ----------
    data : np.ndarray
        待评价决策单元的投入产出数据
    dataref : np.ndarray
        参照决策单元的投入产出数据
    numk : int
        投入变量的个数
    
    Returns
    -------
    np.ndarray
        各决策单元的技术效率值（te）
    """
    n = data.shape[0]                    # 待评价决策单元个数
    theta_results = np.empty(n)          # 存储theta值的数组
    
    for j in range(n):                   # 在每个决策单元上循环
        x = data[j, 0:numk]
        y = data[j, numk:]
        xref = dataref[:, 0:numk]
        yref = dataref[:, numk:]
        
        model = ConcreteModel()
        model.I = Set(initialize=range(dataref.shape[0]))
        model.K = Set(initialize=range(numk))
        model.L = Set(initialize=range(data.shape[1] - numk))
        
        model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
        model.intensity = Var(model.I, bounds=(0.0, None),
            doc='intensity variables')
        
        def obj_rule(model):
            return model.theta
        model.obj = Objective(rule=obj_rule, sense=maximize,
            doc='objective function')
        
        def input_rule(model, k):
            return (
                sum(model.intensity[i] * xref[i, k] for i in model.I)
                <= x[k]
            )
        def output_rule(model, l):
            return (
                sum(model.intensity[i] * yref[i, l] for i in model.I)
                >= model.theta * y[l]
            )
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]
    
    te = 1.0 / theta_results
    return te


# Markdown line 439
from deabook.example.book_cases.loaders import plants
data = plants().head(5)[["Labor", "Fuel", "Capital", "Output"]].to_numpy()
dataref = data
numk = 3
te = dea1(data, dataref, numk)
print(te)


# Markdown line 456
from pyomo.environ import *
import pandas as pd
from deabook.example.book_cases.loaders import plants
import numpy as np

def dea2(dataframe, varname, numk):
    """
    计算DEA技术效率（支持DataFrame版本）
    
    Parameters
    ----------
    dataframe : pd.DataFrame
        待评价决策单元的投入产出数据框
    varname : list[str]
        投入和产出变量的列表，如["K","L","Y"]
    numk : int
        前numk个变量为投入数据
    
    Returns
    -------
    pd.DataFrame
        包含theta值和技术效率值的数据框
    """
    theta_results = {}           # 使用字典存储结果（保留原始索引）
    data = dataframe.loc[:, varname]
    dataref = dataframe.loc[:, varname]
    
    for j in range(data.shape[0]):
        x = data.iloc[j, 0:numk]
        y = data.iloc[j, numk:]
        xref = dataref.iloc[:, 0:numk]
        yref = dataref.iloc[:, numk:]
        xcol = xref.columns
        ycol = yref.columns
        
        model = ConcreteModel()
        model.I = Set(initialize=range(dataref.shape[0]))
        model.K = Set(initialize=range(numk))
        model.L = Set(initialize=range(len(ycol)))
        
        model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
        model.intensity = Var(model.I, bounds=(0.0, None),
            doc='intensity variables')
        
        def obj_rule(model):
            return model.theta
        model.obj = Objective(rule=obj_rule, sense=maximize,
            doc='objective function')
        
        def input_rule(model, k):
            return (
                sum(model.intensity[i] * xref.iloc[i, k] for i in model.I)
                <= x.iloc[k]
            )
        def output_rule(model, l):
            return (
                sum(model.intensity[i] * yref.iloc[i, l] for i in model.I)
                >= model.theta * y.iloc[l]
            )
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]
    
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    thetadf["te"] = 1 / thetadf["theta"]
    return thetadf


# Markdown line 541
import pandas as pd
from deabook.example.book_cases.loaders import plants
ex3 = plants()
te = dea2(ex3, ["Labor", "Fuel", "Capital", "Output"], 3)
print(te)


# Markdown line 555
def dea3(dataframe, varname, numk, eval_query):
    """
    计算DEA技术效率（支持查询条件指定评价单元）
    
    Parameters
    ----------
    eval_query : str or None
        传入数据框.query()方法的参数，如"dmu==1","dmu in [1,2,3]"
        若为None，则评价所有决策单元
    """
    theta_results = {}
    
    if eval_query is None:
        index_list = dataframe.index
    else:
        index_list = dataframe.query(eval_query).index
    
    data = dataframe.loc[index_list, varname]
    dataref = dataframe.loc[index_list, varname]
    
    for j in range(data.shape[0]):
        x = data.iloc[j, 0:numk]
        y = data.iloc[j, numk:]
        xref = dataref.iloc[:, 0:numk]
        yref = dataref.iloc[:, numk:]
        xcol = xref.columns
        ycol = yref.columns
        
        model = ConcreteModel()
        model.I = Set(initialize=range(dataref.shape[0]))
        model.K = Set(initialize=range(numk))
        model.L = Set(initialize=range(data.shape[1] - numk))
        
        model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
        model.intensity = Var(model.I, bounds=(0.0, None),
            doc='intensity variables')
        
        def obj_rule(model):
            return model.theta
        model.obj = Objective(rule=obj_rule, sense=maximize,
            doc='objective function')
        
        def input_rule(model, k):
            return (
                sum(model.intensity[i] * xref.iloc[i, k] for i in model.I)
                <= x.iloc[k]
            )
        def output_rule(model, l):
            return (
                sum(model.intensity[i] * yref.iloc[i, l] for i in model.I)
                >= model.theta * y.iloc[l]
            )
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]
    
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    thetadf["te"] = 1 / thetadf["theta"]
    return thetadf


# Markdown line 629
import pandas as pd
from deabook.example.book_cases.loaders import plants
ex3 = plants()
te = dea3(ex3, ["Labor", "Fuel", "Capital", "Output"], 3, "dmu in [1,2,3,4,5,6,7,8,9,10]")
print(te)


# Markdown line 639
te = dea3(ex3, ["Labor", "Fuel", "Capital", "Output"], 3, "dmu<=10")     # 评价前10个DMU
te = dea3(ex3, ["Labor", "Fuel", "Capital", "Output"], 3, None)           # 评价所有DMU


# Markdown line 648
def dea4(dataframe, varname, numk, eval_query, ref_query):
    """
    计算DEA技术效率（支持分别指定评价单元和参照单元）
    
    Parameters
    ----------
    eval_query : str or None
        筛选待评价决策单元的查询条件
    ref_query : str or None
        筛选参照决策单元的查询条件
    """
    theta_results = {}
    
    if eval_query is None:
        index_list = dataframe.index
    else:
        index_list = dataframe.query(eval_query).index
    
    if ref_query is None:
        index_list_ref = dataframe.index
    else:
        index_list_ref = dataframe.query(ref_query).index
    
    data = dataframe.loc[index_list, varname]
    dataref = dataframe.loc[index_list_ref, varname]
    
    for j in range(data.shape[0]):
        x = data.iloc[j, 0:numk]
        y = data.iloc[j, numk:]
        xref = dataref.iloc[:, 0:numk]
        yref = dataref.iloc[:, numk:]
        xcol = xref.columns
        ycol = yref.columns
        
        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)    # 使用参照单元的原始索引
        model.K = Set(initialize=range(numk))
        model.L = Set(initialize=range(data.shape[1] - numk))
        
        model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
        model.intensity = Var(model.I, bounds=(0.0, None),
            doc='intensity variables')
        
        def obj_rule(model):
            return model.theta
        model.obj = Objective(rule=obj_rule, sense=maximize,
            doc='objective function')
        
        def input_rule(model, k):
            return (
                sum(model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I)
                <= x.iloc[k]
            )
        def output_rule(model, l):
            return (
                sum(model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I)
                >= model.theta * y.iloc[l]
            )
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]
    
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    thetadf["te"] = 1 / thetadf["theta"]
    return thetadf


# Markdown line 726
import pandas as pd
from deabook.example.book_cases.loaders import plants
ex3 = plants()
te = dea4(ex3, ["Labor", "Fuel", "Capital", "Output"], 3, "dmu<=10", "dmu<=20")
print(te)


# Markdown line 738
te = dea4(ex3, ["Labor", "Fuel", "Capital", "Output"], 3, "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")


# Markdown line 746
def dea5(dataframe, varname, numk, eval_query, ref_query):
    """
    计算DEA技术效率（保留原始索引版本）
    """
    theta_results = {}
    
    if eval_query is None:
        index_list = dataframe.index
    else:
        index_list = dataframe.query(eval_query).index
    
    if ref_query is None:
        index_list_ref = dataframe.index
    else:
        index_list_ref = dataframe.query(ref_query).index
    
    xcol = varname[0:numk]
    ycol = varname[numk:]
    
    data = dataframe.loc[index_list, varname]
    dataref = dataframe.loc[index_list_ref, varname]
    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    
    for j in data.index:                    # 在data的原始索引上循环（关键改动）
        x = data.loc[j, xcol]               # 使用loc而非iloc
        y = data.loc[j, ycol]
        
        model = ConcreteModel()
        model.I = Set(initialize=dataref.index)
        model.K = Set(initialize=range(numk))
        model.L = Set(initialize=range(len(ycol)))
        
        model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
        model.intensity = Var(model.I, bounds=(0.0, None),
            doc='intensity variables')
        
        def obj_rule(model):
            return model.theta
        model.obj = Objective(rule=obj_rule, sense=maximize,
            doc='objective function')
        
        def input_rule(model, k):
            return (
                sum(model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I)
                <= x.loc[xcol[k]]
            )
        def output_rule(model, l):
            return (
                sum(model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I)
                >= model.theta * y.loc[ycol[l]]
            )
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]             # j为原始索引值
    
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    thetadf["te"] = 1 / thetadf["theta"]
    return thetadf


# Markdown line 824
import pandas as pd
from deabook.example.book_cases.loaders import plants
ex3 = plants()
te = dea5(ex3, ["Labor", "Fuel", "Capital", "Output"], 3, "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
print(te)


# Markdown line 840
import re

s = "Y = K     L  "
input_vars = s.split('=')[1].strip(' ')
input_vars = re.compile(' +').sub(' ', input_vars).split(' ')
output_vars = s.split('=')[0].strip(' ')
output_vars = re.compile(' +').sub(' ', output_vars).split(' ')


# Markdown line 852
from pyomo.environ import *
import pandas as pd
from deabook.example.book_cases.loaders import plants
import numpy as np
import re

def dea6(formula, dataframe, eval_query, ref_query):
    """
    计算DEA技术效率（支持公式字符串版本）
    
    Parameters
    ----------
    formula : str
        公式字符串，格式为 "产出 = 投入1 投入2 ..."
        例如："Output=Labor Fuel Capital"、"Y = K     L"
    """
    theta_results = {}
    
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
    output_vars = formula.split('=')[0].strip(' ')
    ycol = re.compile(' +').sub(' ', output_vars).split(' ')
    
    data = dataframe.loc[index_list, xcol + ycol]
    dataref = dataframe.loc[index_list_ref, xcol + ycol]
    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    
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
        
        def obj_rule(model):
            return model.theta
        model.obj = Objective(rule=obj_rule, sense=maximize,
            doc='objective function')
        
        def input_rule(model, k):
            return (
                sum(model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I)
                <= x.loc[xcol[k]]
            )
        def output_rule(model, l):
            return (
                sum(model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I)
                >= model.theta * y.loc[ycol[l]]
            )
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]
    
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    thetadf["te"] = 1 / thetadf["theta"]
    return thetadf


# Markdown line 939
import pandas as pd
from deabook.example.book_cases.loaders import plants
ex3 = plants()
te = dea6("Output=Labor Fuel Capital", ex3, "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
print(te)


# Markdown line 959
from pyomo.environ import *
import pandas as pd
from deabook.example.book_cases.loaders import plants
import numpy as np
import re

def dea7(formula, dataframe, rts, eval_query, ref_query):
    """
    计算DEA技术效率（支持CRS/VRS切换版本）
    
    Parameters
    ----------
    rts : str
        规模报酬假设："crs"（规模报酬不变）或"vrs"（规模报酬可变）
    """
    theta_results = {}
    
    if eval_query is None:
        index_list = dataframe.index
    else:
        index_list = dataframe.query(eval_query).index
    
    if ref_query is None:
        index_list_ref = dataframe.index
    else:
        index_list_ref = dataframe.query(ref_query).index
    
    # 解析公式
    input_vars = formula.split('=')[1].strip(' ')
    xcol = re.compile(' +').sub(' ', input_vars).split(' ')
    output_vars = formula.split('=')[0].strip(' ')
    ycol = re.compile(' +').sub(' ', output_vars).split(' ')
    
    data = dataframe.loc[index_list, xcol + ycol]
    dataref = dataframe.loc[index_list_ref, xcol + ycol]
    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    
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
        
        def obj_rule(model):
            return model.theta
        model.obj = Objective(rule=obj_rule, sense=maximize,
            doc='objective function')
        
        def input_rule(model, k):
            return (
                sum(model.intensity[i] * xref.loc[i, xcol[k]] for i in model.I)
                <= x.loc[xcol[k]]
            )
        def output_rule(model, l):
            return (
                sum(model.intensity[i] * yref.loc[i, ycol[l]] for i in model.I)
                >= model.theta * y.loc[ycol[l]]
            )
        def vrs_rule(model):
            """规模报酬可变约束：强度变量之和为1"""
            return sum(model.intensity[i] for i in model.I) == 1
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        # 根据rts参数条件添加VRS约束
        if rts == "vrs":
            model.vrs = Constraint(rule=vrs_rule, doc='VRS constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]
    
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    thetadf["te"] = 1 / thetadf["theta"]
    return thetadf


# Markdown line 1052
import pandas as pd
from deabook.example.book_cases.loaders import plants
ex3 = plants()
te = dea7("Output=Labor Fuel Capital", ex3, "vrs", "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
print(te)


# Markdown line 1075
from pyomo.environ import *
import pandas as pd
from deabook.example.book_cases.loaders import plants
import numpy as np
import re

def dea8(formula: str,
         dataframe: pd.DataFrame,
         rts: str = "crs",
         orient: str = "oo",
         eval_query: str | None = None,
         ref_query: str | None = None) -> pd.DataFrame:
    """
    计算DEA技术效率（完整版本：支持CRS/VRS、投入/产出导向）
    
    Parameters
    ----------
    formula : str
        公式字符串，格式为 "产出 = 投入1 投入2 ..."
    dataframe : pd.DataFrame
        包含投入产出数据的DataFrame
    rts : str, default "crs"
        规模报酬假设："crs"（规模报酬不变）或"vrs"（规模报酬可变）
    orient : str, default "oo"
        导向类型："oo"（产出导向，output-oriented）或"io"（投入导向，input-oriented）
    eval_query : str or None, default None
        筛选待评价决策单元的查询条件
    ref_query : str or None, default None
        筛选参照决策单元的查询条件
    
    Returns
    -------
    pd.DataFrame
        包含theta值和技术效率值的数据框
    """
    theta_results = {}
    
    # 处理查询条件（None表示使用全部数据）
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
    output_vars = formula.split('=')[0].strip(' ')
    ycol = re.compile(' +').sub(' ', output_vars).split(' ')
    
    data = dataframe.loc[index_list, xcol + ycol]
    dataref = dataframe.loc[index_list_ref, xcol + ycol]
    xref = dataref.loc[:, xcol]
    yref = dataref.loc[:, ycol]
    
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
            sense=(maximize if orient == "oo" else minimize),
            doc='objective function'
        )
        
        # 投入约束：根据导向变化
        def input_rule(model, k):
            if orient == "oo":           # 产出导向：投入不超过被评价单元
                return (
                    sum(model.intensity[i] * xref.loc[i, xcol[
                        k]] for i in model.I)
                    <= x.loc[xcol[k]]
                )
            elif orient == "io":         # 投入导向：投入不超过theta倍被评价单元
                return (
                    sum(model.intensity[i] * xref.loc[i, xcol[
                        k]] for i in model.I)
                    <= model.theta * x.loc[xcol[k]]
                )
        

        # 产出约束：根据导向变化
        def output_rule(model, l):
            if orient == "oo":           # 产出导向：产出不低于theta倍被评价单元
                return (
                    sum(model.intensity[i] * yref.loc[i, ycol[
                        l]] for i in model.I)
                    >= model.theta * y.loc[ycol[l]]
                )
            elif orient == "io":         # 投入导向：产出不低于被评价单元
                return (
                    sum(model.intensity[i] * yref.loc[i, ycol[
                        l]] for i in model.I)
                    >= y.loc[ycol[l]]
                )
        
        # VRS约束
        def vrs_rule(model):
            return sum(model.intensity[i] for i in model.I) == 1
        
        model.input = Constraint(model.K, rule=input_rule,
            doc='input constraint')
        model.output = Constraint(model.L, rule=output_rule,
            doc='output constraint')
        
        if rts == "vrs":
            model.vrs = Constraint(rule=vrs_rule, doc='VRS constraint')
        
        opt = SolverFactory('highs')
        solution = opt.solve(model)
        
        theta = np.asarray(list(model.theta[:].value))
        theta_results[j] = theta[0]
    
    thetadf = pd.DataFrame(theta_results, index=["theta"]).T
    # 技术效率计算：产出导向取倒数，投入导向直接取theta
    thetadf["te"] = (1 / thetadf["theta"] if orient == "oo" else thetadf[
        "theta"])
    return thetadf


# Markdown line 1223
import pandas as pd
from deabook.example.book_cases.loaders import plants
ex3 = plants()

# VRS投入导向
te = dea8("Output=Labor Fuel Capital", ex3, "vrs", "io", "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
print(te)


# Markdown line 1245
from pyomo.opt import SolverStatus, TerminationCondition

def solve_with_check(model, solver="highs"):
    """
    带状态检查的求解器调用
    
    检查求解器是否可用、求解是否成功，提供清晰的错误信息。
    """
    opt = SolverFactory(solver)
    
    # 检查求解器是否已安装
    if not opt.available(exception_flag=False):
        raise RuntimeError(
            f"求解器 '{solver}' 未安装或不可用。"
            f"请运行以下命令安装：\n"
            f"  conda install {solver}\n"
            f"或使用其他求解器：'highs'、'ipopt'"
        )
    
    # 调用求解器
    solution = opt.solve(model, load_solutions=False)
    
    # 检查求解状态
    if solution.solver.status != SolverStatus.ok:
        raise RuntimeError(f"求解器状态异常: {solution.solver.status}")
    
    # 检查求解终止条件
    if solution.solver.termination_condition != TerminationCondition.optimal:
        if (solution.solver.termination_condition ==
            TerminationCondition.infeasible):
            raise RuntimeError(
                "线性规划无可行解——请检查数据是否存在异常值（如0或负数），"
                "或投入产出变量是否设置正确。"
            )
        else:
            raise RuntimeError(
                f"求解未达到最优。终止条件: {solution.solver.termination_condition}"
            )
    
    model.solutions.load_from(solution)
    return solution


# Markdown line 1295
from pyomo.environ import *
from pyomo.opt import SolverStatus, TerminationCondition
import pandas as pd
from deabook.example.book_cases.loaders import plants
import numpy as np
import re


class DEAModel:
    """
    DEA效率评价模型（类封装版本）
    
    将DEA模型的构建、求解、结果提取封装为类，减少代码重复，
    提高代码的可维护性和可扩展性。
    
    Parameters
    ----------
    formula : str
        公式字符串，格式为 "产出 = 投入1 投入2 ..."
    dataframe : pd.DataFrame
        包含投入产出数据的DataFrame
    orient : str, default "oo"
        导向类型："oo"（产出导向）或"io"（投入导向）
    rts : str, default "crs"
        规模报酬假设："crs"或"vrs"
    """
    
    def __init__(self, formula, dataframe, orient="oo", rts="crs"):
        self.formula = formula
        self.dataframe = dataframe
        self.orient = orient
        self.rts = rts
        self._parse_formula()
    
    def _parse_formula(self):
        """解析公式字符串，提取投入变量名和产出变量名"""
        input_vars = self.formula.split('=')[1].strip(' ')
        self.xcol = re.compile(' +').sub(' ', input_vars).split(' ')
        output_vars = self.formula.split('=')[0].strip(' ')
        self.ycol = re.compile(' +').sub(' ', output_vars).split(' ')
    
    def _prepare_data(self, eval_query=None, ref_query=None):
        """
        根据查询条件准备评价数据和参考数据
        
        Returns
        -------
        data, dataref : pd.DataFrame
            评价数据和参照数据
        """
        if eval_query is None:
            index_list = self.dataframe.index
        else:
            index_list = self.dataframe.query(eval_query).index
        
        if ref_query is None:
            index_list_ref = self.dataframe.index
        else:
            index_list_ref = self.dataframe.query(ref_query).index
        
        cols = self.xcol + self.ycol
        data = self.dataframe.loc[index_list, cols]
        dataref = self.dataframe.loc[index_list_ref, cols]
        return data, dataref
    
    def _build_model(self, x, y, xref, yref):
        """
        构建Pyomo模型（通用部分）
        
        Parameters
        ----------
        x, y : pd.Series
            被评价单元的投入和产出
        xref, yref : pd.DataFrame
            参照单元的投入和产出
        
        Returns
        -------
        model : ConcreteModel
            构建好的Pyomo模型
        """
        model = ConcreteModel()
        model.I = Set(initialize=xref.index)
        model.K = Set(initialize=range(len(self.xcol)))
        model.L = Set(initialize=range(len(self.ycol)))
        
        model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
        model.intensity = Var(model.I, bounds=(0.0, None),
            doc='intensity variables')
        
        # 目标函数
        def obj_rule(model):
            return model.theta
        model.obj = Objective(
            rule=obj_rule,
            sense=(maximize if self.orient == "oo" else minimize)
        )
        
        # 投入约束
        def input_rule(model, k):
            if self.orient == "oo":
                return (
                    sum(model.intensity[i] * xref.loc[i, self.xcol[k]] 
                          for i in model.I)
                    <= x.loc[self.xcol[k]]
                )
            else:  # io
                return (
                    sum(model.intensity[i] * xref.loc[i, self.xcol[k]] 
                          for i in model.I)
                    <= model.theta * x.loc[self.xcol[k]]
                )
        
        # 产出约束
        def output_rule(model, l):

            if self.orient == "oo":
                return (
                    sum(model.intensity[i] * yref.loc[i, self.ycol[l]] 
                          for i in model.I)
                    >= model.theta * y.loc[self.ycol[l]]
                )
            else:  # io
                return (
                    sum(model.intensity[i] * yref.loc[i, self.ycol[l]] 
                          for i in model.I)
                    >= y.loc[self.ycol[l]]
                )
        
        model.input = Constraint(model.K, rule=input_rule)
        model.output = Constraint(model.L, rule=output_rule)
        
        # VRS约束
        if self.rts == "vrs":
            def vrs_rule(model):
                return sum(model.intensity[i] for i in model.I) == 1
            model.vrs = Constraint(rule=vrs_rule)
        
        return model
    
    def solve(self, eval_query=None, ref_query=None, solver="highs"):
        """
        求解DEA模型并返回效率结果
        
        Parameters
        ----------
        eval_query : str or None
            筛选待评价决策单元的查询条件
        ref_query : str or None
            筛选参照决策单元的查询条件
        solver : str, default "highs"
            求解器名称
        
        Returns
        -------
        pd.DataFrame
            包含theta值和技术效率值的数据框
        """
        data, dataref = self._prepare_data(eval_query, ref_query)
        xref = dataref.loc[:, self.xcol]
        yref = dataref.loc[:, self.ycol]
        
        theta_results = {}
        
        for j in data.index:
            x = data.loc[j, self.xcol]
            y = data.loc[j, self.ycol]
            
            model = self._build_model(x, y, xref, yref)
            
            # 调用带检查的求解
            opt = SolverFactory(solver)
            if not opt.available(exception_flag=False):
                raise RuntimeError(f"求解器 '{solver}' 未安装。请安装后重试。")
            
            solution = opt.solve(model)
            
            # 检查求解状态
            if (solution.solver.termination_condition 
                != TerminationCondition.optimal):
                print(f"警告: DMU {j} 的求解未达到最优，"
                      f"终止条件: {solution.solver.termination_condition}")
            
            theta = np.asarray(list(model.theta[:].value))
            theta_results[j] = theta[0]
        
        thetadf = pd.DataFrame(theta_results, index=["theta"]).T
        thetadf["te"] = (1 / thetadf["theta"] 
                         if self.orient == "oo" else thetadf["theta"])
        return thetadf


# Markdown line 1490
import pandas as pd
from deabook.example.book_cases.loaders import plants

# 读取数据
ex3 = plants()

# 创建模型实例（产出导向、VRS）
model = DEAModel("Output=Labor Fuel Capital", ex3, orient="oo", rts="vrs")

# 求解
te = model.solve(eval_query="dmu in [1,2,3,4,5]", ref_query="dmu<=20")
print(te)

# 切换到CRS投入导向，无需重新定义类
te_io = DEAModel("Output=Labor Fuel Capital", ex3, orient="io", rts="crs").solve()
print(te_io)


# Markdown line 1548
from deabook.example.book_cases.run import FUN
FUN[11]()  # 重新计算本章，输出完整CSV、摘要表与PNG/PDF图

