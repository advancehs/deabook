[TOC]

# 12. DEA程序编写的实战演练

在上一章中，我们详细介绍了如何在Python中一步一步的编写DEA程序。在本章中，我们将用6个例子进行程序编写的示范。

## 12.1 SBM模型的程序编写

Tone(2001)[^1]提出了基于松弛测度的效率模型（SBM）假设，决策单元以P种投入要素$\boldsymbol{X}$生产出Q种产出$\boldsymbol{Y}$ ，Tone(2001)定义了如下的技术效率测度模型：
$$
\begin{aligned} 
\min \quad & \rho = \frac{1 - \frac{1}{P}\sum_{p=1}^{P} s_{p}^{\boldsymbol{X}}/\boldsymbol{X}_{0p}}{1 + \frac{1}{Q}\sum_{j=1}^{Q} s_{j}^{\boldsymbol{Y}}/Y_{0j}} \\
\text{s.t.} \quad 
& \boldsymbol{X}_{0p} = \sum_{i=1}^{N} \lambda_{i}\boldsymbol{X}_{ip} + s_{p}^{\boldsymbol{X}}, & p = 1, \dots, P, \\ 
& Y_{0q} = \sum_{i=1}^{N} \lambda_{i}\boldsymbol{Y}_{iq} - s_{q}^{\boldsymbol{Y}}, & q = 1, \dots, Q, \\ 
& \lambda_{i}, s_{p}^{\boldsymbol{X}}, s_{q}^{\boldsymbol{Y}}, t \geq 0. 
\end{aligned}
$$
在上述表达式中，目标函数是非线性的形式，我们无法直接利用线性规划的技术进行求解，因此在求解之前我们需要对其进行如下线性转换:
$$
\begin{aligned} 
\min \quad & \tau = t - \frac{1}{P}\sum_{p=1}^{P} S_{p}^{\boldsymbol{X}}/\boldsymbol{X}_{0p} \\
\text{s.t.} \quad 
& 1 = t + \frac{1}{Q}\sum_{q=1}^{Q} S_{q}^{\boldsymbol{Y}}/\boldsymbol{Y}_{0q}, \\ 
& t\boldsymbol{X}_{0p} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{X}_{ip} + S_{p}^{\boldsymbol{X}}, & p = 1, \dots, P, \\ 
& t\boldsymbol{Y}_{0q} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{Y}_{iq} - S_{q}^{\boldsymbol{Y}}, & q = 1, \dots, Q, \\ 
& \Lambda_{i}, S_{p}^{\boldsymbol{X}}, S_{q}^{\boldsymbol{Y}} \geq 0,\ t > 0. 
\end{aligned}
$$
接下来我们开始程序编写的第一步工作，转化以上线性规划问题为Python可识别的模型。

在这个线性规划问题中，
(1)变量为{$t$, $S_{p}^{\boldsymbol{X}}$, $S_{q}^{\boldsymbol{Y}}$, $\Lambda_{i}$}，其中，$t$是一维的，$S_{p}^{\boldsymbol{X}}$的维度是$P$, $S_{q}^{\boldsymbol{Y}}$的维度是$Q$, $\Lambda_{i}$的维度是$N$；$\Lambda_{i}, S_{p}^{\boldsymbol{X}}, S_{q}^{\boldsymbol{Y}} \geq 0, t>0$；

(2)目标函数中变量$t$和$S_{p}^{\boldsymbol{X}}$对应的系数分别为1和$\frac{1}{P \cdot \boldsymbol{X}_{0p}}$；求最小值；

(3)约束条件有三个：
- 第一个约束对应的行数有$1$个（变量$t$的个数），其中$t$和$S_{q}^{\boldsymbol{Y}}$的系数分别为1和$\frac{1}{Q \cdot \boldsymbol{Y}_{0q}}$；
- 第二个约束对应的行数有$P$个（变量$S_{p}^{\boldsymbol{X}}$的个数），其中$t$, $\Lambda_{i}$, $S_{p}^{\boldsymbol{X}}$的系数分别为$\boldsymbol{X}_{0p}$, $\boldsymbol{X}_{ip}$和1；
- 第三个约束对应的行数有$Q$个（变量$S_{q}^{\boldsymbol{Y}}$的个数），其中$t$, $\Lambda_{i}$, $S_{q}^{\boldsymbol{Y}}$的系数分别为$\boldsymbol{Y}_{0q}$, $\boldsymbol{Y}_{iq}$和1；
三个约束的约束符号都为$=$。

下面给出求解SBM效率的Python代码。

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
import warnings
warnings.filterwarnings("ignore")
def sbm(formula, dataframe, evaquery, refquery):
    """
    formula: 产出变量=投入变量，如“ Y   = K     L ”
    dataframe: 待评价决策单元的投入产出数据框，按照投入和产出排列
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetax,thetay,obj = {},{}, {}        # 定义thetax,thetay,obj 用于存储计算结果，分别是slakx slacky，和obj
    indexlt = dataframe.query(evaquery).index
    indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0].strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars).split(' ')    
    data = dataframe.loc[indexlt,xcol+ycol]
    dataref = dataframe.loc[indexltref,xcol+ycol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.t = Var( domain =PositiveReals,  doc='CC trans')
        model.thetax = Var(model.K,bounds=(0.0, None), doc='slack x')
        model.thetay = Var(model.L,bounds=(0.0, None), doc='slack y')
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.t-sum(model.thetax[k]/x[k] for k in model.K) / len(model.K)
        def cctrans_rule(model):
            """Return the cctrans  constraint"""
            return 1==model.t+ (sum(model.thetay[l]/y[l] for l in model.L) )/(len(model.L))
        def input_rule(model, k):
            """Return the proper input constraint"""
            return sum(model.lamda[i] * xref.loc[i,xcol[k]] for i in model.I
                    ) == x[k]*model.t-model.thetax[k]
        def output_rule(model, l):
            """Return the proper output constraint"""
            return sum(model.lamda[i] * yref.loc[i,ycol[l]] for i in model.I
                    ) == y[l] * model.t+model.thetay[l]
        model.obj = Objective(rule=objective_rule, sense=minimize, doc='objective function')
        model.cctrans = Constraint(rule= cctrans_rule , doc='cctrans')
        model.input = Constraint(model.K, rule= input_rule , doc='input constraint')
        model.output = Constraint(model.L, rule= output_rule , doc='output constraint')
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            thetax[j] = np.asarray(list(model.thetax[:].value)) /np.asarray(list(model.t[:].value)) # 提取决策变量thetax
            thetay[j] = np.asarray(list(model.thetay[:].value)) /np.asarray(list(model.t[:].value)) # 提取决策变量thetay
            obj[j]= value(model.obj) # 提取目标函数 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        thetaxdf = pd.DataFrame(thetax,index=xcol).T
        thetaxdf.columns = thetaxdf.columns.map(lambda x : "slack "+ str(x) ) 
        thetaydf = pd.DataFrame(thetay,index=ycol).T
        thetaydf.columns = thetaydf.columns.map(lambda y : "slack "+ str(y) ) 
        thetadf=pd.concat([thetaxdf,thetaydf],axis=1)
        redata = pd.concat([objdf,thetadf],axis=1)
    return redata
```

以上SBM模型的Python函数代码基本上保留了第十二章的dea8函数的结构数据和变量传递方式，其主要改变的地方是对线性规划问题的目标函数和约束条件进行了重新定义。第46行声明线性规划是求最小值，第32行定义了线性规划问题的目标函数。第35-45行定义了线性规划问题的约束条件。第28-31行声明了变量及其范围，第28行限定了$t$为正实数，以满足其$>0$的范围。第52行对线性规划问题的求解状态是否收敛进行了判断，如果出现收敛，则取出松弛变量值和目标函数值。需要注意的是线性变换后松弛变量$s_{p}^{\boldsymbol{X}}=S_{p}^{\boldsymbol{X}}/t$，$s_{q}^{\boldsymbol{Y}}=S_{q}^{\boldsymbol{Y}}/t$，因此松弛变量值应除以t值。生成的松弛变量的值以slack的字符串为开头，以投入和产出的变量名为后缀。sbm命令的调用方式和运行结果如下：

```python
import pandas as pd
ex4 = pd.read_stata(r"../../data/Ex4.dta")
sbmte=sbm("Y =K L", ex4,"t==[1,2,3]","t==[1,2,3]" )
sbmte
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20260502163633126.png)
**图11.1：代码运行结果**

## 12.2 包含非合意产出的SBM模型的程序编写

Tone(2004)[^2]对原始SBM模型进行了扩展，使其可以应用到包含非合意产出情形下的技术效率测度，假设决策单元以$P$种投入要素$\boldsymbol X$生产出$Q$种产出$\boldsymbol Y$，并且伴随着产生了$R$种非合意产出$\boldsymbol B$，Tone(2004)[^2]定义了如下的技术效率测度模型：
$$
\begin{aligned} 
\min \quad & \rho=\frac{1-\frac{1}{P}\sum_{p=1}^{P}s_{p}^{\boldsymbol{X}}/\boldsymbol{X}_{0p}}
{1+\frac{1}{Q+R}\left(\sum_{j=1}^{Q}s_{j}^{\boldsymbol{Y}}/Y_{0j}+\sum_{j=1}^{R}s_{j}^{\boldsymbol{B}}/B_{0j}\right)} \\
\text{s.t.} \quad 
& \boldsymbol{X}_{0p}=\sum_{i=1}^{N}\lambda_{i}\boldsymbol{X}_{ip}+s_{p}^{\boldsymbol{X}}, \quad p=1,\dots,P, \\ 
& \boldsymbol{Y}_{0q}=\sum_{i=1}^{N}\lambda_{i}\boldsymbol{Y}_{iq}-s_{q}^{\boldsymbol{Y}}, \quad q=1,\dots,Q, \\ 
& \boldsymbol{B}_{0r}=\sum_{i=1}^{N}\lambda_{i}\boldsymbol{B}_{ir}+s_{r}^{\boldsymbol{B}}, \quad r=1,\dots,R, \\ 
& \lambda_{i},s_{p}^{\boldsymbol{X}},s_{q}^{\boldsymbol{Y}},s_{r}^{\boldsymbol{B}}\geq 0 . 
\end{aligned}
$$
利用线性转换方法，上述优化问题可以转换为如下线性规划问题：
$$
\begin{aligned} 
min \quad & \tau = t - \frac{1}{P}\sum_{p=1}^{P} S_{p}^{\boldsymbol{X}} / \boldsymbol{X}_{0p} \\
s.t. \quad 
& 1 = t + \frac{1}{Q+R} \left( \sum_{q=1}^{Q} S_{q}^{\boldsymbol{Y}} / \boldsymbol{Y}_{0q} + \sum_{r=1}^{R} S_{r}^{\boldsymbol{B}} / \boldsymbol{B}_{0r} \right), \\ 
& t\boldsymbol{X}_{0p} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{X}_{ip} + S_{p}^{\boldsymbol{X}}, \quad p = 1,...,P, \\ 
& t\boldsymbol{Y}_{0q} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{Y}_{iq} - S_{q}^{\boldsymbol{Y}}, \quad q = 1,...,Q, \\ 
& t\boldsymbol{B}_{0r} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{B}_{ir} + S_{r}^{\boldsymbol{B}}, \quad r = 1,...,R, \\ 
& \Lambda_{i}, S_{p}^{\boldsymbol{X}}, S_{q}^{\boldsymbol{Y}}, S_{r}^{\boldsymbol{B}} \geq 0, \quad t > 0. 
\end{aligned}
$$

在这个线性规划问题中，
(1)变量为{$t$, $S_{p}^{\boldsymbol{X}}$, $S_{q}^{\boldsymbol{Y}}$, $S_{r}^{\boldsymbol{B}}$,$\Lambda_{i}$}，其中，$t$是一维的，$S_{p}^{\boldsymbol{X}}$的维度是$P$, $S_{q}^{\boldsymbol{Y}}$的维度是$Q$, $S_{r}^{\boldsymbol{B}}$的维度是$R$，$\Lambda_{i}$的维度是$N$；$ t,\Lambda_{i},S_{p}^{\boldsymbol{X}},S_{q}^{\boldsymbol{Y}} \geq 0$；

(2)目标函数中变量$t$和$S_{p}^{\boldsymbol{X}}$对应的系数分别为1和$\frac{1}{P \cdot \boldsymbol{X}_{0p}}$；求最小值；

(3)约束条件有四个：
- 第一个约束对应的行数有$1$个（变量$t$的个数），其中$t$的系数为1；
- 第二个约束对应的行数有$P$个（变量$S_{p}^{\boldsymbol{X}}$的个数），其中$t$, $\Lambda_{i}$, $S_{p}^{\boldsymbol{X}}$的系数分别为$\boldsymbol{X}_{0p}$, $\boldsymbol{X}_{ip}$和1；
- 第三个约束对应的行数有$Q$个（变量$S_{q}^{\boldsymbol{Y}}$的个数），其中$t$, $\Lambda_{i}$, $S_{q}^{\boldsymbol{Y}}$的系数分别为$\boldsymbol{Y}_{0q}$, $\boldsymbol{Y}_{iq}$和1；
- 第四个约束对应的行数有$R$个（变量$S_{r}^{\boldsymbol{B}}$的个数），其中$t$, $\Lambda_{i}$, $S_{r}^{\boldsymbol{B}}$的系数分别为$\boldsymbol{B}_{0r}$, $\boldsymbol{B}_{ir}$和1；
四个约束的约束符号都为$=$。

由此可见，我们只需对上一节的SBM函数进行相应的修改便可。下面给出求解包含非合意产出的SBM效率的Python代码：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def sbm2( formula, dataframe, evaquery, refquery):
    """
    formula: 产出变量:非期望产出变量=投入变量，如“ Y :CO2  = K     L ”
    dataframe: 待评价决策单元的投入产出数据框，按照投入和产出排列
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetax,thetay,thetab,obj = {},{}, {} ,{} # 定义thetax,thetay,thetab,obj 用于存储计算结果，分别是slakx slacky slackb，和obj
    indexlt = dataframe.query(evaquery).index
    indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0] .split(':')[0] .strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars) .split(' ')
    unoutputvars = formula.split('=')[0] .split(':')[1] .strip(' ') 
    bcol=re.compile(' +').sub(' ',unoutputvars) .split(' ')   
    data = dataframe.loc[indexlt,xcol+ycol+bcol]
    dataref = dataframe.loc[indexltref,xcol+ycol+bcol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    bref=dataref.loc[:,bcol]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        b=data.loc[j,bcol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
        model.t = Var( domain =PositiveReals,  doc='CC trans')
        model.thetax = Var(model.K,bounds=(0.0, None), doc='slack x')
        model.thetay = Var(model.L,bounds=(0.0, None), doc='slack y')
        model.thetab = Var(model.M,bounds=(0.0, None), doc='slack b')
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.t-sum(model.thetax[k]/x[k] for k in model.K) / len(model.K)
        def cctrans_rule(model):
            """Return the cctrans  constraint"""
            return 1==model.t+ (sum(model.thetay[l]/y[l] for l in model.L
                    ) +sum(model.thetab[m]/b[m] for m in model.M))/(len(model.L)+len(model.M))        
        def input_rule(model, k):
            """Return the proper input constraint"""
            return sum(model.lamda[i] * xref.loc[i,xcol[k]] for i in model.I
                    ) == x[k]*model.t-model.thetax[k]
        def output_rule(model, l):
            """Return the proper output constraint"""
            return sum(model.lamda[i] * yref.loc[i,ycol[l]] for i in model.I
                    ) == y[l] * model.t+model.thetay[l]
        def undesirable_output_rule(model, m):
            """Return the proper undesirable output constraint"""
            return sum(model.lamda[i] * bref.loc[i,bcol[m]] for i in model.I
                    ) == b[m] * model.t-model.thetab[m]
        model.obj = Objective(rule=objective_rule, sense=minimize, doc='objective function')
        model.cctrans = Constraint(rule= cctrans_rule , doc='cctrans')
        model.input = Constraint(model.K, rule= input_rule , doc='input constraint')
        model.output = Constraint(model.L, rule= output_rule , doc='output constraint')
        model.undesirable_output = Constraint(model.M, rule= undesirable_output_rule , doc='undesirable output constraint')
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            thetax[j] = np.asarray(list(model.thetax[:].value)) /np.asarray(list(model.t[:].value)) # 提取决策变量thetax
            thetay[j] = np.asarray(list(model.thetay[:].value)) /np.asarray(list(model.t[:].value)) # 提取决策变量thetay
            thetab[j] = np.asarray(list(model.thetab[:].value)) /np.asarray(list(model.t[:].value)) # 提取决策变量thetab
            obj[j]= value(model.obj) # 提取目标函数 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        thetaxdf = pd.DataFrame(thetax,index=xcol).T
        thetaxdf.columns = thetaxdf.columns.map(lambda x : "slack "+ str(x) ) 
        thetaydf = pd.DataFrame(thetay,index=ycol).T
        thetaydf.columns = thetaydf.columns.map(lambda y : "slack "+ str(y) ) 
        thetabdf = pd.DataFrame(thetab,index=bcol).T
        thetabdf.columns = thetabdf.columns.map(lambda b : "slack "+ str(b) ) 
        thetadf=pd.concat([thetaxdf,thetaydf],axis=1)
        thetadf=pd.concat([thetadf,thetabdf],axis=1)
        redata = pd.concat([objdf,thetadf],axis=1)
    return redata
```

在sbm2命令中，我们需要分别传递投入变量、合意产出变量和非合意产出变量，因此我们采用如下表达式：投入变量=合意产出变量: 非合意产出变量（注意这里的冒号是英文的冒号），即：利用“=”和“:”对变量进行分割，从而实现准确传递。代码第13行使用“=”作为分隔符，对传入的sent进行分割，分割得到的第1部分是投入变量，使用strip方法和正则表达式去掉多余空格后，使用空格分割投入变量从而将其变为列表。第15行将使用“=”分割的第二部分再按照“:”分割，得到的第一部分是合意产出变量，第二部分是非合意产出变量。对这两部分同样使用strip方法和正则表达式去掉多余空格后，使用空格分割投入变量从而将其分别变为合意产出和非合意产出的列表。基于新定义的语法规则，我们便可以通过以下方式调用命令，且得到的结果如图所示：

```python
import pandas as pd
ex4 = pd.read_stata(r"../../data/Ex4.dta")
sbmte=sbm2("Y :CO2  = K     L", ex4,"t==[1,2,3]","t==[1,2,3]" )
sbmte
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20260502163726858.png)
**图11.2：代码运行结果**

## 12.3 DDF模型的程序编写

Chung et al. (1997) [^3] 提出的方向距离函数是测度技术效率的另一个常用工具。假设决策单元以$P$种投入要素$\boldsymbol X$ 生产出$Q$种产出$\boldsymbol Y$，并且伴随着产生了$R$种非合意产出$\boldsymbol B$，方向距离函数的定义如下：


$$
\begin{aligned}{D}(\boldsymbol{X}, \boldsymbol{Y}, \boldsymbol{B};\boldsymbol{g})=max{\beta:(\boldsymbol{X},\boldsymbol{Y},\boldsymbol{B})+ \beta \boldsymbol{g} \in T}\end{aligned}
$$
其中T是生产技术集，$g=(g_X,g_Y,g_B)'$是投入变量、合意产出变量和非合意产出变量调整的方向向量。利用数据包络分析的方法，方向距离函数可以由下面的线性规划问题求解：
$$
\begin{aligned}   & \mathop {D}(\boldsymbol{X}, \boldsymbol{Y}, \boldsymbol{B};\boldsymbol{g})=
\max\limits_{\beta, \lambda_1,...,\lambda_N}   \beta \\
\text { s.t. }    &\sum_{i=1}^{N} \lambda_{i} \boldsymbol{X_{i}} \leq \boldsymbol{X}_0+\beta \boldsymbol{g_X}  \\
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol{Y_{i}} \geq \boldsymbol{Y}_0 + \beta \boldsymbol{g_Y}  \\
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol{B_{i}}   = \boldsymbol{ B}_0+ \beta \boldsymbol{g_B}  \\
&\lambda_{i} \geq 0, i=1, \ldots, N  \end{aligned}
$$
在这个线性规划问题中，(1)变量为{$\beta$, $\lambda_{i}$}，其中，$\beta$是一维的， $\lambda_{i}$的维度是$N$；$ \lambda_{i}\geq 0$，$\beta$无约束；

(2)目标函数中变量$\beta$对应的系数为1；求最大值；

(3)约束条件有三个，第一个约束对应的行数有$P$个，$P$为变量$\boldsymbol{X}$的个数，第二个约束对应的行数有$Q$个，$Q$为变量$\boldsymbol{Y}$的个数，第二个约束对应的行数有$R$个，$R$为变量$\boldsymbol{B}$的个数；且第一个约束中变量$\lambda_{i}$和$\beta$对应的系数分别为$\boldsymbol{X_{i}}$和$\boldsymbol{g_X}$，第二个约束中变量$\lambda_{i}$和$\beta$对应的系数分别为$\boldsymbol{Y_{i}}$和$\boldsymbol{g_Y}$，第三个约束中变量$\lambda_{i}$和$\beta$对应的系数分别为$\boldsymbol{B_{i}}$和$\boldsymbol{g_B}$；第一个约束中约束符号为$\leq$，第二个约束中约束符号为$\geq$，第三个约束中约束符号为$=$。
下面给出求解DDF效率的Python代码：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def ddf(formula, dataframe, gx=None , gy=None , gb=None , evaquery=None, refquery=None ):
    """ddf: Directional distance function
    	formula: 产出变量:非期望产出变量=投入变量，如“ Y :CO2  = K     L ”
        dataframe: 待评价决策单元的投入产出数据框，按照投入和产出排列
        gx (list, optional): 投入方向向量. 默认为 [-1].
        gy (list, optional): 合意产出方向向量. 默认为 [1].
        gb (list, optional): 非合意产出方向向量. 默认为[-1].
        evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"。默认为全部
        refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]。默认为全部
    """
    obj = {}                # 定义obj 用于存储计算结果，是obj
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0] .split(':')[0] .strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars) .split(' ')
    unoutputvars = formula.split('=')[0] .split(':')[1] .strip(' ') 
    bcol=re.compile(' +').sub(' ',unoutputvars) .split(' ')   
    data = dataframe.loc[indexlt,xcol+ycol+bcol]
    dataref = dataframe.loc[indexltref,xcol+ycol+bcol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    bref=dataref.loc[:,bcol]
    if type(gx)==type(None):
        gx=[-1]*len(xcol)
    if type(gy)==type(None):
        gy=[1]*len(ycol)
    if type(gb)==type(None):
        gb=[-1]*len(bcol)
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        b=data.loc[j,bcol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index)    # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(bounds=(None, None), within=Reals,doc='directional distance')
        model.lamda = Var(model.I , bounds=(0.0, None),within=Reals, doc='intensity variables')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.theta *1  
        def input_rule(model, k):
            """Return the proper input constraint"""
            return sum(model.lamda[i] * xref.loc[i,xcol[k]] for i in model.I
                    ) - model.theta*gx[k]*x[k] <= x[k]
        def output_rule(model, l):
            """Return the proper output constraint"""
            return -sum(model.lamda[i] * yref.loc[i,ycol[l]] for i in model.I
                    ) + model.theta*gy[l] *y[l]<= -y[l]
        def undesirable_output_rule(model, m):
            """Return the proper undesirable output constraint"""
            return sum(model.lamda[i] * bref.loc[i,bcol[m]] for i in model.I
                    ) -model.theta*gb[m]*b[m]==  b[m]
        model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
        model.input = Constraint(model.K,  rule=input_rule, doc='input constraint')
        model.output = Constraint(model.L,  rule=output_rule, doc='output constraint')
        model.undesirable_output = Constraint(model.M, rule=undesirable_output_rule, doc='undesirable output constraint')
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            obj[j]= value(model.obj) # 提取目标函数 
        objdf= pd.DataFrame(obj,index=["te"]).T
    return objdf
```

上述ddf函数代码比sbm2函数多了三个输入项，用以传递方向向量的变量，$\boldsymbol{g_X}$， $\boldsymbol{g_Y}$，$\boldsymbol{g_B}$。第48-49行对决策变量进行了定义。第50-52行对线性规划问题的目标函数进行了声明，第53-64行对线性规划问题的投入变量、合意产出变量和非合意产出变量对应的不等式约束进行了声明。在命令语法上，ddf比sbm2函数多了gx，gy，gb三个选项用于传递方向向量的变量。gx，gy，gb的默认值为None，表示他们为可选项。当它们省缺时，第33=38行进行了省缺值时的相关赋值。当投入变量的方向向量缺省时，我们令$\boldsymbol{g_X}=-\boldsymbol X$；当合意产出变量的方向向量缺省时，我们令$\boldsymbol{g_Y}=\boldsymbol Y$；当非合意产出变量的方向向量缺省时，我们令$\boldsymbol{g_B}=-\boldsymbol B$。以下代码说明了命令ddf的使用方法：

```python
import pandas as pd
ex4 = pd.read_stata(r"../../data/Ex4.dta")
ddfte=ddf("Y:    CO2 =K     L  ",ex4,    )
ddfte
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20260502163845791.png)
**图11.3：代码运行结果**


## 12.4 NDDF模型的程序编写

Zhou et al. (2012)[^4]提出的非径向方向距离函数(NDDF)允许投入和产出变量进行不同比例的调整。假设决策单元以$P$种投入要素$\boldsymbol X$生产出$Q$种产出$\boldsymbol Y$，并且伴随着产生了$R$种非合意产出$\boldsymbol B$，NDDF的定义如下：
$$
\vec{D}(\boldsymbol X, \boldsymbol Y, \boldsymbol B;\boldsymbol g)=max\{\boldsymbol{w'\beta}:(\boldsymbol X, \boldsymbol Y, \boldsymbol B)+\boldsymbol{\beta'}diag(\boldsymbol g) \in T\}
$$
其中$T$是生产技术集；$\boldsymbol g=(\boldsymbol{g_X,g_Y,g_B})'$是投入变量，合意产出变量和非合意产出变量调整的方向向量，$\boldsymbol{\beta}=({\beta_X,\beta_Y,\beta_B})'$是投入和产出变量对应的调整因子向量；$\boldsymbol w=(\omega_X,\omega_Y,\omega_B)'$是各个调整因子对应的权重向量。利用数据包络分析的方法，非径向方向距离函数可以由以下线性规划问题求解：
$$
\begin{array}{l}
 & \vec{D}(\boldsymbol x, \boldsymbol y, \boldsymbol b;\boldsymbol g) = \max   \boldsymbol{w'\beta} \\
\text { s.t. } \quad 
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol X_{i} \leq \boldsymbol x+diag(\boldsymbol g_{X})\boldsymbol {\beta_{X}} \\
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol Y_{i} \geq \boldsymbol y+diag(\boldsymbol g_{Y})\boldsymbol {\beta_{Y}} \\
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol  B_{i}=\boldsymbol b+diag(\boldsymbol g_{B})\boldsymbol {\beta_{B}} \\
&\lambda_{n} \geq 0, n=1, \ldots, N ,\boldsymbol{\beta} \geq0
\end{array}
$$
在这个线性规划问题中，(1)变量为{$\beta$, $\lambda_{i}$}，其中，$\beta$的维度是一， $\lambda_{i}$的维度是$N$；$ \lambda_{i}\geq 0$，$\beta\geq0$；

(2)目标函数中变量$\beta$对应的系数为$w$；求最大值；

(3)约束条件有三个，第一个约束对应的行数有$P$个，$P$为变量$\boldsymbol{X}$的个数，第二个约束对应的行数有$Q$个，$Q$为变量$\boldsymbol{Y}$的个数，第二个约束对应的行数有$R$个，$R$为变量$\boldsymbol{B}$的个数；且第一个约束中变量$\lambda_{i}$和$\beta$对应的系数分别为$\boldsymbol{X_{i}}$和$\boldsymbol{g_X}$，第二个约束中变量$\lambda_{i}$和$\beta$对应的系数分别为$\boldsymbol{Y_{i}}$和$\boldsymbol{g_Y}$，第三个约束中变量$\lambda_{i}$和$\beta$对应的系数分别为$\boldsymbol{B_{i}}$和$\boldsymbol{g_B}$；第一个约束中约束符号为$\leq$，第二个约束中约束符号为$\geq$，第三个约束中约束符号为$=$。
下面给出求解NDDF效率的Python代码：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def nddf(formula, dataframe, gx=None , gy=None, gb=None, weight =None, evaquery=None, refquery=None ):
    """nddf: Non-radial Directional distance function
    	formula: 产出变量:非期望产出变量=投入变量，如“ Y :CO2  = K     L ”
        dataframe: 待评价决策单元的投入产出数据框，按照投入和产出排列
        gx (list, optional): 投入方向向量. 默认为 [-1].
        gy (list, optional): 合意产出方向向量. 默认为 [1].
        gb (list, optional): 非合意产出方向向量. 默认为[-1].
        weight(list, optional): weght matrix
        evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"。默认为全部
        refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]。默认为全部
    """
    thetax,thetay,thetab,obj = {},{},{},{}      # 定义thetax,thetay,thetab,obj 用于存储计算结果， 分别是betax,betay,betab,obj
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0] .split(':')[0] .strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars) .split(' ')
    unoutputvars = formula.split('=')[0] .split(':')[1] .strip(' ') 
    bcol=re.compile(' +').sub(' ',unoutputvars) .split(' ')   
    data = dataframe.loc[indexlt,xcol+ycol+bcol]
    dataref = dataframe.loc[indexltref,xcol+ycol+bcol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    bref=dataref.loc[:,bcol]
    if type(gx)==type(None):
        gx=[-1]*len(xcol)
    if type(gy)==type(None):
        gy=[1]*len(ycol)
    if type(gb)==type(None):
        gb=[-1]*len(bcol)
    if type(weight) == type(None):
        weight=[]
        fenmu = 1*int(gx[0]!=0) + 1*int(gy[0]!=0) + 1*int(gb[0]!=0)
        for _ in range(len(xcol)):
            weight.append(1/fenmu/len(xcol))
        for _ in range(len(ycol)):
            weight.append(1/fenmu/len(ycol))
        for _ in range(len(bcol)):
            weight.append(1/fenmu/len(bcol))  
    iweight = weight[0:len(xcol)]
    oweight = weight[len(xcol):len(xcol)+len(ycol)]
    bweight = weight[len(xcol)+len(ycol):len(xcol)+len(ycol)+len(bcol)]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        b=data.loc[j,bcol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
        model.thetax = Var(model.K,bounds=(0.0, None),within=Reals, doc='scale factor x')
        model.thetay = Var(model.L,bounds=(0.0, None),within=Reals, doc='scale factor y')
        model.thetab = Var(model.M,bounds=(0.0, None),within=Reals, doc='scale factor b')
        model.lamda = Var(model.I , bounds=(0.0, None),within=Reals, doc='intensity variables')
        def objective_rule(model):
            """Return the proper objective function"""
            return -sum( iweight[k]*gx[k]* model.thetax[k] for k in model.K
                ) + sum(oweight[l]*gy[l]* model.thetay[l] for l in model.L
                ) - sum(bweight[m]*gb[m]* model.thetab[m] for m in model.M)
        def input_rule(model, k):
            """Return the proper input constraint"""
            return sum(model.lamda[i] * xref.loc[i,xcol[k]] for i in model.I
                ) - gx[k]*x[k]*model.thetax[k]<=  x[k]
        def output_rule(model, l):
            """Return the proper output constraint"""
            return -1* sum(model.lamda[i] * yref.loc[i,ycol[l]] for i in model.I
                ) + gy[l]*y[l]*model.thetay[l] <= -1*y[l]
        def undesirable_output_rule(model, m):
            """Return the proper undesirable output constraint"""
            return sum(model.lamda[i] * bref.loc[i,bcol[m]] for i in model.I
                ) - gb[m]*b[m]*model.thetab[m] == b[m]
        model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
        model.input = Constraint(model.K,  rule=input_rule, doc='input constraint')
        model.output = Constraint(model.L,  rule=output_rule, doc='output constraint')
        model.undesirable_output = Constraint(model.M, rule=undesirable_output_rule, doc='undesirable output constraint')
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            thetax[j] = np.asarray(list(model.thetax[:].value))  # 提取决策变量thetax
            thetay[j] = np.asarray(list(model.thetay[:].value))  # 提取决策变量thetay
            thetab[j] = np.asarray(list(model.thetab[:].value))   # 提取决策变量thetab
            obj[j]= value(model.obj) # 提取目标函数 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        thetaxdf = pd.DataFrame(thetax,index=xcol).T
        thetaxdf.columns = thetaxdf.columns.map(lambda x : "scale "+ str(x) ) 
        thetaydf = pd.DataFrame(thetay,index=ycol).T
        thetaydf.columns = thetaydf.columns.map(lambda y : "scale "+ str(y) ) 
        thetabdf = pd.DataFrame(thetab,index=bcol).T
        thetabdf.columns = thetabdf.columns.map(lambda b : "scale "+ str(b) ) 
        thetadf=pd.concat([thetaxdf,thetaydf],axis=1)
        thetadf=pd.concat([thetadf,thetabdf],axis=1)
        redata = pd.concat([objdf,thetadf],axis=1)
    return redata
```

上述nddf函数代码比ddf函数代码多了一个输入项weight，用以传递权重向量$\boldsymbol w$。第61-64行对决策变量进行了定义。第65-69行对线性规划问题的目标函数进行了声明，第70-81行对线性规划问题的投入变量、合意产出变量和非合意产出变量对应的不等式约束进行了声明。如果线性规划问题求解收敛，则在第88-102行将求得的目标函数值值以及$\beta$数据框redata中，并返回。weight的默认值为None，表示他们为可选项。当它省缺时，第40=51行进行了省缺值时的相关赋值。以下代码说明了命令nddf的使用方法：

```python
import pandas as pd
ex4 = pd.read_stata(r"Ex4.dta")
nddfte=nddf(" Y:    CO2   =K     L  ",  ex4,  )
nddfte
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20260502163918820.png)
**图11.4：代码运行结果**


## 12.5 Malmquist生产率指数的程序编写
### 12.5.1 产出方向的Malmquist生产率指数
产出方向的Malmquist生产率指数定义如下：
$$
M(\boldsymbol X_{t+1}, \boldsymbol Y_{t+1}, \boldsymbol X_{t}, \boldsymbol Y_{t})=
\left[\frac{D^{t}(\boldsymbol X_{t+1}, \boldsymbol Y_{t+1})}{D^{t}(\boldsymbol X_{t}, \boldsymbol Y_{t})} 
\frac{D^{t+1}(\boldsymbol X_{t+1}, \boldsymbol Y_{t+1})}{D^{t+1}(\boldsymbol X_{t}, \boldsymbol Y_{t})}\right]^{1 / 2}
$$
从这个定义中，我们可以看到计算Malmquist生产率指数，需要估计4个谢泼德产出距离函数：${D^{t}(\boldsymbol X_{t+1}, \boldsymbol Y_{t+1})}$、${D^{t}(\boldsymbol X_{t}, \boldsymbol Y_{t})}$ 、${D^{t+1}(\boldsymbol X_{t+1}, \boldsymbol Y_{t+1})}$和${D^{t+1}(\boldsymbol X_{t}, \boldsymbol Y_{t})}$。对此，我们可以使用dea8命令进行相应的估计。

下面我们以ex4.dta数据集为例进行Malmquist生产率指数的计算。

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def dea8(formula, dataframe, rts, orient, evaquery=None, refquery=None):
    """
    formula: 产出变量 = 投入变量，如“ Y    = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    rts: 传入表示可变规模报酬或不变规模报酬的字符串，"crs","vrs" 二选一
    orient: 传入表示产出方向或投入方向的字符串，"oo","io" 二选一
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0].strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars).split(' ')    
    data = dataframe.loc[indexlt,xcol+ycol]
    dataref = dataframe.loc[indexltref,xcol+ycol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense= (maximize if orient == "oo" else minimize),doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
            elif orient == "io":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= model.theta*x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
            elif orient == "io":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=y[l]
        def vrs_rule(model):
            """Return the various return to scale constraint"""
            return sum(model.lamda[i] for i in model.I) == 1
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
        if rts == "vrs":
            model.vrs = Constraint(rule=vrs_rule, doc='various return to scale rule') # 定义与可变规模报酬有关的约束条件
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
        lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
        obj = value(model.obj) # 提取目标函数 
        #print("optimum theta: \n {} ".format(theta))
        #print("optimum lamda: \n {} ".format(lamda))
        #print("optimal objective: {}".format(obj))
        thetalt[j]=theta[0]  
        thetadf =pd.DataFrame(thetalt,index=["theta"]).T
        thetadf["te"] = (1/thetadf["theta"] if orient == "oo" else  thetadf["theta"] )
    return thetadf

"""D^{t}(X_{t},Y_{t})""" 
D11_list = []
for t in range(1,4):
    D11 = dea8("Y=K L ", ex4, "crs", "oo", "t=={}".format(t), "t=={}".format(t))[["te"]]
    D11_list.append(D11)
D11df = pd.concat(D11_list)
D11df.rename(columns={"te":"D11"}, inplace=True)
"""D^{t+1}(X_{t+1},Y_{t+1})""" 
D22_list = []
for t in range(2,4):
    D22 = dea8("Y=K L ", ex4, "crs", "oo", "t=={}".format(t), "t=={}".format(t))[["te"]]
    D22_list.append(D22)
D22df = pd.concat(D22_list)
D22df.rename(columns={"te":"D22"}, inplace=True)
"""D^{t}(X_{t+1},Y_{t+1})""" 
D12_list = []
for t in range(2,4):
    D12 = dea8("Y=K L ", ex4, "crs", "oo", "t=={}".format(t), "t=={}".format(t-1))[["te"]]
    D12_list.append(D12)
D12df = pd.concat(D12_list)
D12df.rename(columns={"te":"D12"}, inplace=True)
"""D^{t+1}(X_{t},Y_{t})""" 
D21_list = []
for t in range(2,4):
    D21 = dea8("Y=K L ", ex4, "crs", "oo", "t=={}".format(t-1), "t=={}".format(t))[["te"]]
    D21_list.append(D21)
D21df = pd.concat(D21_list)
D21df.rename(columns={"te":"D21"}, inplace=True)

df = pd.concat([D11df, D22df, D12df, D21df], axis=1)
ex42 = pd.merge(ex4, df, left_index=True, right_index=True, how="left")
ex42["mpi"] = ex42["D12"] / ex42["D11"].shift(1) * ex42["D11"] / ex42["D21"].shift(1)
```

以上代码呈现了Malmquist生产率指数的具体计算过程。接下来，我们将此过程封装为python命令。

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def dea8(formula, dataframe, rts, orient, evaquery=None, refquery=None):
    """
    formula: 产出变量 = 投入变量，如“ Y    = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    rts: 传入表示可变规模报酬或不变规模报酬的字符串，"crs","vrs" 二选一
    orient: 传入表示产出方向或投入方向的字符串，"oo","io" 二选一
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0].strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars).split(' ')    
    data = dataframe.loc[indexlt,xcol+ycol]
    dataref = dataframe.loc[indexltref,xcol+ycol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense= (maximize if orient == "oo" else minimize),doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
            elif orient == "io":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= model.theta*x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
            elif orient == "io":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=y[l]
        def vrs_rule(model):
            """Return the various return to scale constraint"""
            return sum(model.lamda[i] for i in model.I) == 1
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
        if rts == "vrs":
            model.vrs = Constraint(rule=vrs_rule, doc='various return to scale rule') # 定义与可变规模报酬有关的约束条件
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
        lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
        obj = value(model.obj) # 提取目标函数 
        #print("optimum theta: \n {} ".format(theta))
        #print("optimum lamda: \n {} ".format(lamda))
        #print("optimal objective: {}".format(obj))
        thetalt[j]=theta[0]  
        thetadf =pd.DataFrame(thetalt,index=["theta"]).T
        thetadf["te"] = (1/thetadf["theta"] if orient == "oo" else  thetadf["theta"] )
    return thetadf
def mpi( formula ,data, id, t):
    """
    formula: 产出变量=投入变量 ，如“Y=K L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    id: 个体变量
    t: 时间变量
    """
    tlt = pd.Series(data[t]).drop_duplicates().sort_values() 
    D11_list = []
    for tindex in tlt.index:
        print(t, tlt.iloc[tindex])
        D11 = dea8(
            formula, data, "crs", "oo",
            "{}=={}".format(t, tlt.iloc[tindex]),
            "{}=={}".format(t, tlt.iloc[tindex])
        )[["te"]]
        D11_list.append(D11)
    D11df = pd.concat(D11_list)
    D11df.rename(columns={"te": "D11"}, inplace=True)
    """D^{t}(X_{t+1},Y_{t+1})""" 
    D12_list = []
    for tindex in tlt.index[1:]:
        D12 = dea8(
            formula, data, "crs", "oo",  # 关键修复：ex4 -> data
            "{}=={}".format(t, tlt.iloc[tindex]),
            "{}=={}".format(t, tlt.iloc[tindex-1])
        )[["te"]]
        D12_list.append(D12)
    D12df = pd.concat(D12_list)
    D12df.rename(columns={"te": "D12"}, inplace=True)
    """D^{t+1}(X_{t},Y_{t})""" 
    D21_list = []
    for tindex in tlt.index[1:]:
        D21 = dea8(
            formula, data, "crs", "oo",  # 关键修复：ex4 -> data
            "{}=={}".format(t, tlt.iloc[tindex-1]),
            "{}=={}".format(t, tlt.iloc[tindex])
        )[["te"]]
        D21_list.append(D21)
    D21df = pd.concat(D21_list)
    D21df.rename(columns={"te": "D21"}, inplace=True)
    df = pd.concat([D11df,D12df],axis=1)
    df = pd.concat([df,D21df],axis=1)
    data2 = pd.merge(data,df,left_index=True,right_index=True,how="left")
    data2["mpi"] = data2["D12"] / data2["D11"].shift(1) * data2["D11"]/ data2["D21"].shift(1)
    data2.drop(columns = ["D11","D12","D21"],inplace =True)
    return data2
```

mpi命令的id和t分别传递了决策单元的名称变量和时间变量；我们用mpi命令重新计算决策单元的Malmquist生产率指数。具体代码如下：

```python
import pandas as pd
ex4 = pd.read_stata(r"Ex4.dta")
data2=mpi( "Y=K L " , ex4, "id", "t" )
data2
```

在mpi命令中，我们采用了同期生产技术的假定，其生产前沿的构造只使用一个时期决策单元的投入和产出数据。接下来，我们在mpi命令的基础上进行修改，以便考虑各类不同的生产技术设定。首先，我们考虑时序生产技术下的Malmquist生产率指数。时序生产技术假定过去的所有生产技术都是可以被当期使用的，因此，其生产前沿的构造使用了决策单元当期及之前的所有投入和产出数据。对于计算时序生产技术下的Malmquist生产率指数，我们只需将mpi命令代码中的第85行、第91行和第97行中的"=="于修改为"<="。其次，我们考虑视窗（window）生产技术下的Malmquist生产率指数。视窗生产技术以一定时间窗口内的观测值来构造生产前沿。例如设定窗宽为$h$，在时期$t$的生产技术采用了${t-h,...t+h}$时期的决策单元投入和产出变量的观测值作为技术参照集。因此我们需要将mpi命令中的第85行、第91行和第97行代码分别修改为：

```python
D11 = dea8(formula, data, "crs", "oo","{}=={}".format(t,tlt.iloc[tindex]) ,"{}<={}<={}".format(tlt.iloc[0 if tindex-h<0 else tindex-h],t,tlt.iloc[tlt.index.max() if tindex+h>tlt.index.max() else tindex+h]) )[["te"]]
D12 = dea8(formula, ex4, "crs", "oo","{}=={}".format(t,tlt.iloc[tindex]) ,"{}<={}<={}".format(tlt.iloc[0 if tindex-1-h<0 else tindex-1-h],t,tlt.iloc[tlt.index.max() if tindex-1+h>tlt.index.max() else tindex-1+h]) )[["te"]] 
D21 = dea8(formula, ex4, "crs", "oo","{}=={}".format(t,tlt.iloc[tindex-1]) ,"{}<={}<={}".format(tlt.iloc[0 if tindex-h<0 else tindex-h],t,tlt.iloc[tlt.index.max() if tindex+h>tlt.index.max() else tindex+h])  ) [["te"]]
```

现在我们将时序生产技术和视窗生产技术作为两个选项，添加到mpi命令中，具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def dea8(formula, dataframe, rts, orient, evaquery=None, refquery=None):
    """
    formula: 产出变量 = 投入变量，如“ Y    = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    rts: 传入表示可变规模报酬或不变规模报酬的字符串，"crs","vrs" 二选一
    orient: 传入表示产出方向或投入方向的字符串，"oo","io" 二选一
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0].strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars).split(' ')    
    data = dataframe.loc[indexlt,xcol+ycol]
    dataref = dataframe.loc[indexltref,xcol+ycol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense= (maximize if orient == "oo" else minimize),doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
            elif orient == "io":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= model.theta*x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
            elif orient == "io":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=y[l]
        def vrs_rule(model):
            """Return the various return to scale constraint"""
            return sum(model.lamda[i] for i in model.I) == 1
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
        if rts == "vrs":
            model.vrs = Constraint(rule=vrs_rule, doc='various return to scale rule') # 定义与可变规模报酬有关的约束条件
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
        lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
        obj = value(model.obj) # 提取目标函数 
        #print("optimum theta: \n {} ".format(theta))
        #print("optimum lamda: \n {} ".format(lamda))
        #print("optimal objective: {}".format(obj))
        thetalt[j]=theta[0]  
        thetadf =pd.DataFrame(thetalt,index=["theta"]).T
        thetadf["te"] = (1/thetadf["theta"] if orient == "oo" else  thetadf["theta"] )
    return thetadf
def mpi2(formula, data, id, t, tech=None ):
    """
    formula: 投入变量=产出变量，如“ Y=K     L  ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    id: 个体变量
    t: 时间变量
    tecn: 从当期生产技术（None/"com"）、时序生产技术（"seq"）、视窗生产技术中选择（"window 4"）"window"后面添加视窗的大小。
    """
    tlt = pd.Series(data[t]).drop_duplicates().sort_values() 
    """D^{t}(X_{t},Y_{t})""" 
    D11df = pd.DataFrame()
    for tindex in tlt.index:
        evaquery11="{}=={}".format(t,tlt.iloc[tindex])
        if (type(tech)==type(None)) or (tech=="com"):
            refquery11="{}=={}".format(t,tlt.iloc[tindex])
        elif tech=="seq":
            refquery11="{}<={}".format(t,tlt.iloc[tindex])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery11="{}<={}<={}".format(
                    tlt.iloc[0 if tindex-h<0 else tindex-h],t,tlt.iloc[tlt.index.max() if tindex+h>tlt.index.max() else tindex+h])
        D11 = dea8(formula,data, "crs", "oo",evaquery11 ,refquery11 )[["te"]]
        D11df = D11df.append(D11 )
    D11df.rename(columns = {"te":"D11"},inplace=True)
    """D^{t}(X_{t+1},Y_{t+1})""" 
    D12df = pd.DataFrame()
    for tindex in tlt.index[1:]:
        evaquery12="{}=={}".format(t,tlt.iloc[tindex])
        if (type(tech)==type(None)) or tech=="com":
            refquery12="{}=={}".format(t,tlt.iloc[tindex-1])
        elif tech=="seq":
            refquery12="{}<={}".format(t,tlt.iloc[tindex-1])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery12="{}<={}<={}".format(
                tlt.iloc[0 if tindex-1-h<0 else tindex-1-h],t,tlt.iloc[tlt.index.max() if tindex-1+h>tlt.index.max() else tindex-1+h])
        D12 = dea8(formula,data, "crs", "oo",evaquery12 ,refquery12 )[["te"]] 
        D12df = D12df.append(D12 )
    D12df.rename(columns = {"te":"D12"},inplace=True)
    """D^{t+1}(X_{t},Y_{t})""" 
    D21df = pd.DataFrame()
    for tindex in tlt.index[1:]:
        evaquery21="{}=={}".format(t,tlt.iloc[tindex-1])
        if (type(tech)==type(None)) or tech=="com":
            refquery21="{}=={}".format(t,tlt.iloc[tindex])
        elif tech=="seq":
            refquery21="{}<={}".format(t,tlt.iloc[tindex])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery21="{}<={}<={}".format(
                    tlt.iloc[0 if tindex-h<0 else tindex-h],t,tlt.iloc[tlt.index.max() if tindex+h>tlt.index.max() else tindex+h])
        D21 = dea8(formula,data, "crs", "oo",evaquery21 ,refquery21  ) [["te"]]
        D21df = D21df.append(D21)
    D21df.rename(columns = {"te":"D21"},inplace=True)
    df = pd.concat([D11df,D12df],axis=1)
    df = pd.concat([df,D21df],axis=1)
    data2 = pd.merge(data,df,left_index=True,right_index=True,how="left")
    data2["mpi"] = data2["D12"] / data2["D11"].shift(1) * data2["D11"]/ data2["D21"].shift(1)
    data2.drop(columns = ["D11","D12","D21"],inplace =True)
    return data2
```


### 12.5.2 全局生产技术下的Malmquist生产率指数
接下来，我们考虑全局生产技术下的Malmquist生产率指数。全局生产技术使用所有时期决策单元投入和产出变量的观测值来构造技术参照集。具体而言，全局生产技术下的Malmquist生产率指数定义为：
$$
M(\boldsymbol X_{t+1}, \boldsymbol Y_{t+1}, \boldsymbol X_{t}, \boldsymbol Y_{t})=
\frac{D^{G}(\boldsymbol X_{t+1}, \boldsymbol Y_{t+1})}{D^{G}(\boldsymbol X_{t}, \boldsymbol Y_{t})}
$$
其中$D^{G}(·)$表示全局生产技术下的产出距离函数，可以通过以下线性规划问题求解：
$$
\begin{aligned} & D^{G}(X^*,Y^*)^{-1} = \max\limits_{\theta,\lambda_{11},...,\lambda_{NT1}} \theta
\\ 
\text {s.t.} & \sum_ {t=1}^{T} \sum_{i=1}^{N} \lambda_{it} \boldsymbol {X}_{kit} \leq \boldsymbol { X}_k^* , k=1,2,...,K \\
& \sum_ {t=1}^{T} \sum_{i=1}^{N} \lambda_{it}\boldsymbol {Y} _ {lit} \geq \boldsymbol {Y}_l^* , l=1,2,...,L\\ 
& { \lambda_{it} } \geq 0 ,i=1,2,...,N;t=1,2,...,T;-\infty \leq \theta \leq +\infty \end{aligned}
$$
基于以上表达式，我们可以使用dea8命令进行计算。具体代码如下：

```python
import pandas as pd
ex4 = pd.read_stata(r"Ex4.dta")
DG_list = [] 
for t in range(1, 4):
    DG = dea8(
        "Y=K L ", 
        ex4, 
        "crs", 
        "oo",
        evaquery="t=={}".format(t), 
        refquery=None               
    )[["te"]]
    DG_list.append(DG)
DGdf = pd.concat(DG_list)  
ex42 = pd.concat([ex4,DGdf],axis=1)
ex42["mpi"] = ex42["te"] / ex42["te"].shift(1)   
ex42.drop(columns = ["te"],inplace = True)
ex42
```

我们通过global选项，把以上全局技术情形下的Malmquist生产率指数添加到原mpi2命令中，具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def dea8(formula, dataframe, rts, orient, evaquery=None, refquery=None):
    """
    formula: 产出变量 = 投入变量，如“ Y    = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    rts: 传入表示可变规模报酬或不变规模报酬的字符串，"crs","vrs" 二选一
    orient: 传入表示产出方向或投入方向的字符串，"oo","io" 二选一
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0].strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars).split(' ')    
    data = dataframe.loc[indexlt,xcol+ycol]
    dataref = dataframe.loc[indexltref,xcol+ycol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense= (maximize if orient == "oo" else minimize),doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
            elif orient == "io":
                return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= model.theta*x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            if orient == "oo":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
            elif orient == "io":
                return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=y[l]
        def vrs_rule(model):
            """Return the various return to scale constraint"""
            return sum(model.lamda[i] for i in model.I) == 1
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
        if rts == "vrs":
            model.vrs = Constraint(rule=vrs_rule, doc='various return to scale rule') # 定义与可变规模报酬有关的约束条件
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
        lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
        obj = value(model.obj) # 提取目标函数 
        #print("optimum theta: \n {} ".format(theta))
        #print("optimum lamda: \n {} ".format(lamda))
        #print("optimal objective: {}".format(obj))
        thetalt[j]=theta[0]  
        thetadf =pd.DataFrame(thetalt,index=["theta"]).T
        thetadf["te"] = (1/thetadf["theta"] if orient == "oo" else  thetadf["theta"] )
    return thetadf
def mpi3( formula, data, id, t, tech=None ):
    """
    formula: 投入变量=产出变量，如“ Y=K     L  ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    id: 个体变量
    t: 时间变量
    tech: 从当期生产技术（None/"com"）、时序生产技术（"seq"）、视窗生产技术中选择（"window 4"）"window"后面添加视窗的大小和全局生产技术（"global"）。    """
    tlt = pd.Series(data[t]).drop_duplicates().sort_values() 
    """D^{t}(X_{t},Y_{t})"""
    D11_list = []
    for tindex in tlt.index:
        evaquery11 = "{}=={}".format(t, tlt.iloc[tindex])
        if (tech is None) or (tech == "com"):
            refquery11 = "{}=={}".format(t, tlt.iloc[tindex])
        elif tech == "seq":
            refquery11 = "{}<={}".format(t, tlt.iloc[tindex])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery11 = "{}<={}<={}".format(
                tlt.iloc[0 if tindex-h < 0 else tindex-h], t,
                tlt.iloc[tlt.index.max() if tindex+h > tlt.index.max() else tindex+h]
            )
        elif tech == "global":
            refquery11 = None
        D11 = dea8(formula, data, "crs", "oo", evaquery11, refquery11)[["te"]]
        D11_list.append(D11)
    D11df = pd.concat(D11_list)
    D11df.rename(columns={"te": "D11"}, inplace=True)

    """D^{t}(X_{t+1},Y_{t+1})"""
    D12_list = []
    for tindex in tlt.index[1:]:
        evaquery12 = "{}=={}".format(t, tlt.iloc[tindex])
        if (tech is None) or (tech == "com"):
            refquery12 = "{}=={}".format(t, tlt.iloc[tindex-1])
        elif tech == "seq":
            refquery12 = "{}<={}".format(t, tlt.iloc[tindex-1])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery12 = "{}<={}<={}".format(
                tlt.iloc[0 if tindex-1-h < 0 else tindex-1-h], t,
                tlt.iloc[tlt.index.max() if tindex-1+h > tlt.index.max() else tindex-1+h]
            )
        elif tech == "global":
            refquery12 = None
        D12 = dea8(formula, data, "crs", "oo", evaquery12, refquery12)[["te"]]  # 关键修复：ex4 -> data
        D12_list.append(D12)
    D12df = pd.concat(D12_list)
    D12df.rename(columns={"te": "D12"}, inplace=True)

    """D^{t+1}(X_{t},Y_{t})"""
    D21_list = []
    for tindex in tlt.index[1:]:
        evaquery21 = "{}=={}".format(t, tlt.iloc[tindex-1])
        if (tech is None) or (tech == "com"):
            refquery21 = "{}=={}".format(t, tlt.iloc[tindex])
        elif tech == "seq":
            refquery21 = "{}<={}".format(t, tlt.iloc[tindex])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery21 = "{}<={}<={}".format(
                tlt.iloc[0 if tindex-h < 0 else tindex-h], t,
                tlt.iloc[tlt.index.max() if tindex+h > tlt.index.max() else tindex+h]
            )
        elif tech == "global":
            refquery21 = None
        D21 = dea8(formula, data, "crs", "oo", evaquery21, refquery21)[["te"]]  # 关键修复：ex4 -> data
        D21_list.append(D21)
    D21df = pd.concat(D21_list)
    D21df.rename(columns={"te": "D21"}, inplace=True)
    df = pd.concat([D11df,D12df],axis=1)
    df = pd.concat([df,D21df],axis=1)
    data2 = pd.merge(data,df,left_index=True,right_index=True,how="left")
    data2["mpi"] = data2["D12"] / data2["D11"].shift(1) * data2["D11"]/ data2["D21"].shift(1)
    data2.drop(columns = ["D11","D12","D21"],inplace =True)
    return data2
```

至此我们便完成了一个估计Malmquist生产率指数的相对完整命令。这一命令通过选项的方式将多种生产技术情形综合起来，体现了Python命令在使用上的灵活性和扩展性。下面我们给出几个调用mpi3命令的例子：

```python
import pandas as pd
ex4 = pd.read_stata(r"Ex4.dta")
data2=mpi3( "Y=K L ", ex4, "id", "t",  tech="com")
data3=mpi3( "Y=K L ",ex4, "id", "t",  tech="seq")
data4=mpi3( "Y=K L ", ex4, "id", "t", tech="window 2")
data5=mpi3( "Y=K L ",ex4, "id", "t",   tech="global")
```

## 12.6 Malmquist-Luenberger生产率指数的程序编写

产出方向的Malmquist-Luenberger生产率指数的定义如下：

$$
M L_{t}^{t+1}=\left[\frac{1+D^{t}\left(\boldsymbol{x}_{t}, \boldsymbol{y}_{t}, \boldsymbol{b}_{t} ; \boldsymbol{g}\right)}{1+D^{t}\left(\boldsymbol{x}_{t+1}, \boldsymbol{y}_{t+1}, \boldsymbol{b}_{t+1} ; \boldsymbol{g}\right)} \times \frac{1+D^{t+1}\left(\boldsymbol{x}_{t}, \boldsymbol{y}_{t}, \boldsymbol{b}_{t} ; \boldsymbol{g}\right)} {1+D^{t+1}\left(\boldsymbol{x}_{t+1}, \boldsymbol{y}_{t+1}, \boldsymbol{b}_{t+1} ; \boldsymbol{g}\right)}  \right]^{1 / 2}
$$
从这个定义中，我们可以看到Malmquist-Luenberger生产率指数的计算需要估计4个方向距离函数：$D^{t}\left(\boldsymbol{x}_{t}, \boldsymbol{y}_{t}, \boldsymbol{b}_{t} ; \boldsymbol{g}\right)$、$D^{t+1}\left(\boldsymbol{x}_{t+1}, \boldsymbol{y}_{t+1}, \boldsymbol{b}_{t+1} ; \boldsymbol{g}\right)$、$D^{t+1}\left(\boldsymbol{x}_{t}, \boldsymbol{y}_{t}, \boldsymbol{b}_{t} ; \boldsymbol{g}\right)$和$D^{t}\left(\boldsymbol{x}_{t+1}, \boldsymbol{y}_{t+1}, \boldsymbol{b}_{t+1} ; \boldsymbol{g}\right)$。对此我们可以使用ddf命令进行相应的估计。
下面我们以Ex4.dta数据集为例来展示如何计算Malmquist-Luenberger生产率指数。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def ddf(formula, dataframe, gx=None , gy=None , gb=None , evaquery=None, refquery=None ):
    """ddf: Directional distance function
    	formula: 产出变量:非期望产出变量=投入变量，如“ Y :CO2  = K     L ”
        dataframe: 待评价决策单元的投入产出数据框，按照投入和产出排列
        gx (list, optional): 投入方向向量. 默认为 [-1].
        gy (list, optional): 合意产出方向向量. 默认为 [1].
        gb (list, optional): 非合意产出方向向量. 默认为[-1].
        evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"。默认为全部
        refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]。默认为全部
    """
    obj = {}                # 定义obj 用于存储计算结果，是obj
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0] .split(':')[0] .strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars) .split(' ')
    unoutputvars = formula.split('=')[0] .split(':')[1] .strip(' ') 
    bcol=re.compile(' +').sub(' ',unoutputvars) .split(' ')   
    data = dataframe.loc[indexlt,xcol+ycol+bcol]
    dataref = dataframe.loc[indexltref,xcol+ycol+bcol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    bref=dataref.loc[:,bcol]
    if type(gx)==type(None):
        gx=[-1]*len(xcol)
    if type(gy)==type(None):
        gy=[1]*len(ycol)
    if type(gb)==type(None):
        gb=[-1]*len(bcol)
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        b=data.loc[j,bcol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index)    # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(bounds=(None, None), within=Reals,doc='directional distance')
        model.lamda = Var(model.I , bounds=(0.0, None),within=Reals, doc='intensity variables')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.theta *1  
        def input_rule(model, k):
            """Return the proper input constraint"""
            return sum(model.lamda[i] * xref.loc[i,xcol[k]] for i in model.I
                    ) - model.theta*gx[k]*x[k] <= x[k]
        def output_rule(model, l):
            """Return the proper output constraint"""
            return -sum(model.lamda[i] * yref.loc[i,ycol[l]] for i in model.I
                    ) + model.theta*gy[l] *y[l]<= -y[l]
        def undesirable_output_rule(model, m):
            """Return the proper undesirable output constraint"""
            return sum(model.lamda[i] * bref.loc[i,bcol[m]] for i in model.I
                    ) -model.theta*gb[m]*b[m]==  b[m]
        model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
        model.input = Constraint(model.K,  rule=input_rule, doc='input constraint')
        model.output = Constraint(model.L,  rule=output_rule, doc='output constraint')
        model.undesirable_output = Constraint(model.M, rule=undesirable_output_rule, doc='undesirable output constraint')
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            obj[j]= value(model.obj) # 提取目标函数 
        objdf= pd.DataFrame(obj,index=["te"]).T
    return objdf

"""D^{t}(X_{t},Y_{t})""" 
D11_list = []
for t in range(1, 4):
    D11 = ddf("Y:CO2=K L", ex4, evaquery="t=={}".format(t), refquery="t=={}".format(t))
    D11_list.append(D11)
D11df = pd.concat(D11_list)
D11df.rename(columns={"te": "D11"}, inplace=True)

"""D^{t+1}(X_{t+1},Y_{t+1})""" 
D22_list = []
for t in range(2, 4):
    D22 = ddf("Y:CO2=K L", ex4, evaquery="t=={}".format(t), refquery="t=={}".format(t))
    D22_list.append(D22)
D22df = pd.concat(D22_list)
D22df.rename(columns={"te": "D22"}, inplace=True)

"""D^{t}(X_{t+1},Y_{t+1})""" 
D12_list = []
for t in range(2, 4):
    D12 = ddf("Y:CO2=K L", ex4, evaquery="t=={}".format(t), refquery="t=={}".format(t-1))
    D12_list.append(D12)
D12df = pd.concat(D12_list)
D12df.rename(columns={"te": "D12"}, inplace=True)

"""D^{t+1}(X_{t},Y_{t})""" 
D21_list = []
for t in range(2, 4):
    D21 = ddf("Y:CO2=K L", ex4, evaquery="t=={}".format(t-1), refquery="t=={}".format(t))
    D21_list.append(D21)
D21df = pd.concat(D21_list)
D21df.rename(columns={"te": "D21"}, inplace=True)
df = pd.concat([D11df,D22df],axis=1)
df = pd.concat([df,D12df],axis=1)
df = pd.concat([df,D21df],axis=1)
ex42 = pd.merge(ex4,df,left_index=True,right_index=True,how="left")
ex42["mpi"] = (1+ex42["D11"].shift(1))/(1+ex42["D12"])   * (1+ ex42["D21"].shift(1))/ (1+ex42["D22"])
```

从以上代码中我们不难发现，整个计算过程和Malmquist生产率指数的计算过程相似。正如我们此前所讨论的，DEA编程的基本框架包含两部分：一是DEA模型对应的线性规划问题求解，二是不同技术参照集的数据选择。参考mpi3命令我们接下来编写一个计算Malmquist-Luenberger生产率指数的mlpi命令。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def ddf(formula, dataframe, gx=None , gy=None , gb=None , evaquery=None, refquery=None ):
    """ddf: Directional distance function
    	formula: 产出变量:非期望产出变量=投入变量，如“ Y :CO2  = K     L ”
        dataframe: 待评价决策单元的投入产出数据框，按照投入和产出排列
        gx (list, optional): 投入方向向量. 默认为 [-1].
        gy (list, optional): 合意产出方向向量. 默认为 [1].
        gb (list, optional): 非合意产出方向向量. 默认为[-1].
        evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"。默认为全部
        refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]。默认为全部
    """
    obj = {}                # 定义obj 用于存储计算结果，是obj
    if type(evaquery)==type(None):
        indexlt = dataframe.index
    else:
        indexlt = dataframe.query(evaquery).index
    if type(refquery)==type(None):
        indexltref = dataframe. index
    else:
        indexltref = dataframe.query(refquery).index
    inputvars = formula.split('=')[1].strip(' ') 
    xcol = re.compile(' +').sub(' ',inputvars).split(' ')
    outputvars = formula.split('=')[0] .split(':')[0] .strip(' ') 
    ycol = re.compile(' +').sub(' ',outputvars) .split(' ')
    unoutputvars = formula.split('=')[0] .split(':')[1] .strip(' ') 
    bcol=re.compile(' +').sub(' ',unoutputvars) .split(' ')   
    data = dataframe.loc[indexlt,xcol+ycol+bcol]
    dataref = dataframe.loc[indexltref,xcol+ycol+bcol]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    bref=dataref.loc[:,bcol]
    if type(gx)==type(None):
        gx=[-1]*len(xcol)
    if type(gy)==type(None):
        gy=[1]*len(ycol)
    if type(gb)==type(None):
        gb=[-1]*len(bcol)
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        b=data.loc[j,bcol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index)    # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(bounds=(None, None), within=Reals,doc='directional distance')
        model.lamda = Var(model.I , bounds=(0.0, None),within=Reals, doc='intensity variables')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.theta *1  
        def input_rule(model, k):
            """Return the proper input constraint"""
            return sum(model.lamda[i] * xref.loc[i,xcol[k]] for i in model.I
                    ) - model.theta*gx[k]*x[k] <= x[k]
        def output_rule(model, l):
            """Return the proper output constraint"""
            return -sum(model.lamda[i] * yref.loc[i,ycol[l]] for i in model.I
                    ) + model.theta*gy[l] *y[l]<= -y[l]
        def undesirable_output_rule(model, m):
            """Return the proper undesirable output constraint"""
            return sum(model.lamda[i] * bref.loc[i,bcol[m]] for i in model.I
                    ) -model.theta*gb[m]*b[m]==  b[m]
        model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
        model.input = Constraint(model.K,  rule=input_rule, doc='input constraint')
        model.output = Constraint(model.L,  rule=output_rule, doc='output constraint')
        model.undesirable_output = Constraint(model.M, rule=undesirable_output_rule, doc='undesirable output constraint')
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            obj[j]= value(model.obj) # 提取目标函数 
    objdf= pd.DataFrame(obj,index=["te"]).T
    return objdf


def mlpi( formula, data, id, t, tech=None ):
    """
    formula: 产出变量:非期望产出变量=投入变量，如“ Y :CO2  = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    id: 个体变量
    t: 时间变量
    tech: 从当期生产技术（None/"com"）、时序生产技术（"seq"）、视窗生产技术中选择（"window 4"）"window"后面添加视窗的大小和全局生产技术（"global"）。    """
    tlt = pd.Series(data[t]).drop_duplicates().sort_values() 
    """D^{t}(X_{t},Y_{t})""" 
    D11df = pd.DataFrame()
    for tindex in tlt.index:
        evaquery11="{}=={}".format(t,tlt.iloc[tindex])
        if (type(tech)==type(None)) or (tech=="com"):
            refquery11="{}=={}".format(t,tlt.iloc[tindex])
        elif tech=="seq":
            refquery11="{}<={}".format(t,tlt.iloc[tindex])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery11="{}<={}<={}".format(
                    tlt.iloc[0 if tindex-h<0 else tindex-h],t,tlt.iloc[tlt.index.max() if tindex+h>tlt.index.max() else tindex+h])
        elif tech =="global":
            refquery11=None
        D11 = ddf(formula, data, evaquery=evaquery11 ,refquery=refquery11 )[["te"]]
        D11df = D11df.append(D11 )
    D11df.rename(columns = {"te":"D11"},inplace=True)
    """D^{t}(X_{t+1},Y_{t+1})""" 
    D12df = pd.DataFrame()
    for tindex in tlt.index[1:]:
        evaquery12="{}=={}".format(t,tlt.iloc[tindex])
        if (type(tech)==type(None)) or tech=="com":
            refquery12="{}=={}".format(t,tlt.iloc[tindex-1])
        elif tech=="seq":
            refquery12="{}<={}".format(t,tlt.iloc[tindex-1])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery12="{}<={}<={}".format(
                tlt.iloc[0 if tindex-1-h<0 else tindex-1-h],t,tlt.iloc[tlt.index.max() if tindex-1+h>tlt.index.max() else tindex-1+h])
        elif tech =="global":
            break
        D12 = ddf(formula, data, evaquery=evaquery12 ,refquery=refquery12 )[["te"]] 
        D12df = D12df.append(D12 )
    try:
        D12df.rename(columns = {"te":"D12"},inplace=True)
    except:
        pass
    """D^{t+1}(X_{t},Y_{t})""" 
    D21df = pd.DataFrame()
    for tindex in tlt.index[1:]:
        evaquery21="{}=={}".format(t,tlt.iloc[tindex-1])
        if (type(tech)==type(None)) or tech=="com":
            refquery21="{}=={}".format(t,tlt.iloc[tindex])
        elif tech=="seq":
            refquery21="{}<={}".format(t,tlt.iloc[tindex])
        elif "window " in tech:
            h = int(tech.split(" ")[1])
            refquery21="{}<={}<={}".format(
                    tlt.iloc[0 if tindex-h<0 else tindex-h],t,tlt.iloc[tlt.index.max() if tindex+h>tlt.index.max() else tindex+h])
        elif tech =="global":
            break
        D21 = ddf(formula, data, evaquery= evaquery21 ,refquery=refquery21  ) [["te"]]
        D21df = D21df.append(D21)
    try:
        D21df.rename(columns = {"te":"D21"},inplace=True)
    except:
        pass
    if tech !="global":
        df = pd.concat([D11df,D12df],axis=1)
        df = pd.concat([df,D21df],axis=1)
        data2 = pd.merge(data,df,left_index=True,right_index=True,how="left")
        data2["mlpi"] = (1+data2["D11"].shift(1))/(1+data2["D12"])   * (1+ data2["D21"].shift(1))/ (1+data2["D11"])
        data2.drop(columns = ["D11","D12","D21"],inplace =True)
    else:
        data2 = pd.merge(data,D11df,left_index=True,right_index=True,how="left")
        data2["mlpi"] = (1+data2["D11"].shift(1) )/(1+ data2["D11"]   ) 
    return data2
```

 mlpi命令保留了mpi3命令的代码结构，主要是将谢泼德产出距离函数替换为方向距离函数，即将mpi3命令调用的dea8命令更改为ddf命令。



## 参考文献
[^1]:Tone, K. (2001). A slacks-based measure of efficiency in data envelopment analysis. European Journal of Operational Research, 130(3), 498–509. https://doi.org/10.1016/s0377-2217(99)00407-5
[^2]:Tone, K. (2004). Dealing with undesirable outputs in DEA: a slacks-based measure (SBM) approach. 日本オペレーションズ・リサーチ学会春季研究発表会アブストラクト集, 2004, 44–45.
[^3]:Chung, Y. H., Färe, R., & Grosskopf, S. (1997). Productivity and Undesirable Outputs: A Directional Distance Function Approach. Journal of Environmental Management, 51(3), 229–240. https://doi.org/10.1006/jema.1997.0146
[^4]:Zhou, P., Ang, B. W., & Wang H.. (2012). Energy and CO2 emission performance in electricity generation: A non-radial directional distance function approach. European Journal of Operational Research, 90(1), 625–635. https://doi.org/10.1016/j.ejor.2012.04.022









