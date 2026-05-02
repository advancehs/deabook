[TOC]

# 13. 基于DEA对偶问题的影子价格计算

在前面的章节中，我们通过线性规划对DEA模型进行了求解。每一个线性规划问题都存在一个对偶问题。通过求解对偶问题，我们可以进一步计算生产资源（投入要素和非期望产出）的影子价格。当我们在DEA模型中讨论影子价格，我们在讨论的是：在保持当前生产效率不变的前提下，每增加或减少一单位投入要素（如劳动力、资本）或非期望产出（如污染、废弃物）时，对效率值的边际影响。其经济学含义是资源在最优配置下的隐含价值。本文将介绍如何应用Python编写基于DEA模型的影子价格计算程序。

## 13.1 线性规划的对偶问题

考虑如下问题，企业A有3中投入要素，每种投入要素的最大可以使用量分别为8、16 和12。企业计划生产2种商品，第一种商品每单位需要投入的生产要素为分别为1、4和0，第二种商品每单位需要投入的生产要素为分别为2、0和4。两商品的市场价格分别为2和3。企业该如何安排生产才能使其利润实现最大化？

上述问题可以总结为表11.1：

|       | 商品1 | 商品2 |      |
| ----- | ----- | ----- | ---- |
| 要素1 | 1     | 2     | 8    |
| 要素2 | 4     | 0     | 16   |
| 要素3 | 0     | 4     | 12   |
|       | 2     | 3     |      |

设企业生产的商品分别为$x_1$和$x_2$，则企业的决策问题可以表示为：
$$
\begin{align}
& \max z=2x_{1}+3 x_{2}  \\
\text { s.t. } 
& \ x_{1}+2 x_{2} \leq 8     \\
& \ 4 x_{1}        \leq 16     \\
& \ 4 x_{2} \leq 12   \\
& \ x_{1} \geq 0, x_{2} \geq 0  
\end{align}
$$
现在我们从另一个角度进行分析，假设企业B想将企业A的所有生产资源收买过来，其需要支付的最小代价是多少呢？换言之，企业A放弃其拥有的生产资源需要得到多少的补偿呢？虽然这,就要求企业A出让生产一单位某种商品的资源所获得的回报大于等于该商品的价格，即$\ y_1+4y_2 \geq 2, \ 2y_1+4y_3      \geq 4  $（$y_1,y_2,y_3$为一单位生产资源的价格）。进而，我们可以得到以下线性规划问题：
$$
\begin{align}
& \min w=8y_1+16y_2+12y_3    \\
\text { s.t. } 
& \ y_1+4y_2 \geq 2      \\
& \ 2y_1+4y_3      \geq 4     \\
& \ y_{1} \geq 0, y_{2} \geq 0   , y_{3} \geq 0   
\end{align}
$$
以上便是线性规划的对偶问题。

接下来，我们从数学表达式的角度分析原问题的对偶问题。我们知道上述实际问题的原问题的表达式如下的数学问题：
$$
\begin{align}
& \max z=2x_{1}+3 x_{2} \tag{1}  \\
\text { s.t. } 
& \ x_{1}+2 x_{2} \leq 8    \tag{2} \\
& \ 4 x_{1}        \leq 16    \tag{3} \\
& \ 4 x_{2} \leq 12 \tag{4} \\
& \ x_{1} \geq 0, x_{2} \geq 0   \tag{5}
\end{align}
$$
我们的目的是获得目标函数式(1)式的最大值，可以对约束条件式(2)、式(3)和式(4)的组合对目标函数进行放缩，从而获得目标函数式(1)的上限。如$(2)+(3) \times1/4$得$2x_{1}+3 x_{2}=z \leq15$，可知在这种组合情况下，目标函数式(1)的上限为15。那么怎样找到一个最佳的组合方式？我们可以设定约束条件式(2)、式(3)和式(4)分别乘上$y_1$、$y_2$和$y_3$，且$y_1\geq0,y_2\geq0,y_3\geq0$，即考虑$(2)\times y_1+(3) \times y_2+(4)\times y_3$，可得式(6)：
$$
(y_1+4y_2)\times x_1+(2y_1+4y_3)\times x_2\leq 8y_1+16y_2+12y_3  \tag{6}
$$
要想的到式(1)，我们需要限定式(6)中$x_1$的系数不小于式(1)中$x_1$的系数，式(6)中$x_2$的系数不小于式(1)中$x_2$的系数，可得：
$$
z=2x_1+3x_2\leq(y_1+4y_2)\times x_1+(2y_1+4y_3)\times x_2\leq 8y_1+16y_2+12y_3 \tag{7}
$$
如果知道了$8y_1+16y_2+12y_3$的最小值，则$z$的最大值就可以求出来。此时我们得到了另一个线性规划问题：
$$
\begin{align}
& \min w=8y_1+16y_2+12y_3    \\
\text { s.t. } 
& \ y_1+4y_2 \geq 2      \\
& \ 2y_1+4y_3      \geq 4     \\
& \ y_{1} \geq 0, y_{2} \geq 0   , y_{3} \geq 0   
\end{align}
$$
此时，得到了一个线性规划原问题的对偶问题。

基于此，我们总结得到以下要点：
（1）原问题有几个约束，对偶问题就有几个决策变量（用于在约束两边同时乘该决策变量）；原问题有几个决策变量，对偶问题就有几个约束；

（2）原问题的决策变量的系数对应对偶问题的右端项；原问题的右端项对应对偶问题的决策变量的系数（原问题的右端项乘了对偶问题的决策变量）；

（3）原问题的系数矩阵的转置对应对偶问题的系数矩阵（原问题的左端想乘了对偶问题的决策变量，合并同类项后便可得到对偶问题的系数矩阵）；

（4）原问题是最大化问题，对偶问题就是最小化问题；原问题是最小化问题，对偶问题就是最大化问题；

（5）原问题与对偶问题之间约束条件符号方向和决策变量符号范围的对应关系。其探讨的核心逻辑围绕“正常情况”展开，并通过约束条件与变量符号的映射关系建立原问题与对偶问题的联系。

首先，线性规划问题的“正常情况”是其最基础形式，为对偶转换提供基准，主要包括：
- **最大化问题**的“正常情况”要求所有决策变量非负（\(x_j \geq 0\)），且约束条件为小于等于型（\(\leq\)）。这一设计的逻辑在于：最大化目标需上界约束以防止解趋向无穷大，而非负变量符合实际资源的物理意义（如投入不可为负）。  
- **最小化问题**的“正常情况”则要求变量非负（\(x_j \geq 0\)），且约束为大于等于型（\(\geq\)），以确保资源投入满足最低需求的下界限制。

其次，原问题的约束类型按照以下规则，直接决定对偶变量的符号范围：
- **原问题约束为\(\leq\)（最大化）或\(\geq\)（最小化）**：对偶变量\(y_i \geq 0\)。此时约束代表资源的有效限制，影子价格反映其边际贡献或成本节省，非负性符合经济学直觉。  
- **原问题约束为等式（\(=\)）**：对偶变量\(y_i\)无约束（可正可负）。等式约束意味着资源被严格固定，其影子价格可能因冗余或短缺呈现不同符号。

此时，我们可以进一步讨论原问题变量的符号范围如何决定对偶约束的方向：  
（1）**原问题变量非负（\(x_j \geq 0\)）**：  
   - 若原问题为**最大化**，对偶约束为\(\leq c_j\)（资源消耗的边际收益上限）；  
   - 若原问题为**最小化**，对偶约束为\(\geq c_j\)（资源投入的边际成本下限）。  

（2）**原问题变量非正（\(x_j \leq 0\)）**：需通过变量替换（如\(x_j' = -x_j \geq 0\)）转化为非负形式。替换后原问题的目标系数与约束矩阵相应列符号反转，对偶约束方向随之反向。  
（3）**原问题变量无约束**：对偶约束必须为等式（\(=\)）。无约束变量可自由取值，其对偶约束需严格平衡目标函数与资源消耗的关系，避免解无界。  


除上述情况外，我们还可能遇到原问题包含非标准约束或变量符号的情况。此时需通过数学变形转化为正常情况：  
（1）**约束条件符号不符合正常情况**：  
   - 最大化问题中出现\(\geq\)约束时，等式两边乘以\(-1\)转化为\(\leq\)；  
   - 最小化问题中出现\(\leq\)约束时，同样乘以\(-1\)转化为\(\geq\)。  
   - 转化后，对偶变量仍按正常情况的符号规则处理。 

例如，原问题为最大化但约束为 \(\geq\) 时，乘以 \(-1\) 后约束变为 \(\leq\)，对应对偶变量 \(y_i \geq 0\)。但由于原约束方向非正常，实际对偶变量应为 \(y_i \leq 0\)。

（2）**变量符号非正常**：  
   - 非正变量（\(x_j \leq 0\)）通过替换\(x_j = -x_j'\)（\(x_j' \geq 0\)）转化为非负形式，此时原问题的目标系数与约束矩阵对应列符号反转。  
   - 变量无约束时直接对应对偶约束为等式，无需额外调整。  

上述规则的本质是通过数学对称性揭示资源的隐含价值（影子价格）。我们可将原问题与对偶问题的符号对应归纳为以下原则。
（1）**原问题约束类型**：  
   - \(\leq\)（最大化）或\(\geq\)（最小化）→ 对偶变量非负（\(y_i \geq 0\)）；  
   - \(=\) → 对偶变量无约束。  

（2）**原问题变量符号**：  
   - \(x_j \geq 0\) → 对偶约束为\(\leq\)（最大化）或\(\geq\)（最小化）；  
   - \(x_j \leq 0\) → 对偶约束方向反向；  
   - \(x_j\)无约束 → 对偶约束为等式。  


## 13.2 SBM模型的对偶问题

了解了如何从原问题得出对偶问题后，我们再来看一下Tone(2004)[^1]的SBM模型的线性规划问题：
$$
\begin{aligned} 
\min \quad & \tau = t - \frac{1}{P} \sum_{p=1}^{P} S_{p}^{\boldsymbol{X}} / \boldsymbol{X}_{0p} \\
\text{s.t.} \quad & 1 = t + \frac{1}{Q+R} \left( \sum_{q=1}^{Q} S_{q}^{\boldsymbol{Y}} / \boldsymbol{Y}_{0q} + \sum_{r=1}^{R} S_{r}^{\boldsymbol{B}} / \boldsymbol{B}_{0r} \right) \\ 
& t\boldsymbol{X}_{0p} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{X}_{ip} + S_{p}^{\boldsymbol{X}}, \quad p=1,\dots,P \\ 
& t\boldsymbol{Y}_{0q} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{Y}_{iq} - S_{q}^{\boldsymbol{Y}}, \quad q=1,\dots,Q \\ 
& t\boldsymbol{B}_{0r} = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{B}_{ir} + S_{r}^{\boldsymbol{B}}, \quad r=1,\dots,R \\ 
& \Lambda_{i}, S_{p}^{\boldsymbol{X}}, S_{q}^{\boldsymbol{Y}} \geq 0, \quad t > 0 
\end{aligned}
$$
首先，由于$t>0$暗含在第一个约束条件中，所以我们可以去掉原线性规划中$t>0$的约束。此外，为了方便从原问题写对偶问题，我们可以将原问题的目标函数和约束条件中的决策变量按照固定的顺序排列，系数为0的可以补0。从而得到如下线性规划问题：
$$
\begin{aligned} 
\min \quad & \tau = \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{0} - \frac{1}{P} \sum_{p=1}^{P} \frac{S_{p}^{\boldsymbol{X}}}{\boldsymbol{X}_{0p}} + S_{q}^{\boldsymbol{Y}}\boldsymbol{0} + S_{r}^{\boldsymbol{B}}\boldsymbol{0} + t \\
\text{s.t.} \quad & 
\sum_{i=1}^{N} \Lambda_{i}\boldsymbol{0} + S_{p}^{\boldsymbol{X}}\boldsymbol{0} + \frac{1}{Q+R} \sum_{q=1}^{Q} \frac{S_{q}^{\boldsymbol{Y}}}{\boldsymbol{Y}_{0q}} + \frac{1}{Q+R} \sum_{r=1}^{R} \frac{S_{r}^{\boldsymbol{B}}}{\boldsymbol{B}_{0r}} + t = 1 \\ 
& \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{X}_{ip} + S_{p}^{\boldsymbol{X}} + S_{q}^{\boldsymbol{Y}}\boldsymbol{0} + S_{r}^{\boldsymbol{B}}\boldsymbol{0} - t\boldsymbol{X}_{0p} = 0, \quad p=1,\dots,P \\ 
& \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{Y}_{iq} + S_{p}^{\boldsymbol{X}}\boldsymbol{0} - S_{q}^{\boldsymbol{Y}} + S_{r}^{\boldsymbol{B}}\boldsymbol{0} - t\boldsymbol{Y}_{0q} = 0, \quad q=1,\dots,Q \\ 
& \sum_{i=1}^{N} \Lambda_{i}\boldsymbol{B}_{ir} + S_{p}^{\boldsymbol{X}}\boldsymbol{0} + S_{q}^{\boldsymbol{Y}}\boldsymbol{0} + S_{r}^{\boldsymbol{B}} - t\boldsymbol{B}_{0r} = 0, \quad r=1,\dots,R \\ 
& \Lambda_{i}, S_{p}^{\boldsymbol{X}}, S_{q}^{\boldsymbol{Y}} \geq 0 
\end{aligned}
$$
下面我们根据原问题写出对偶问题：
$$
\begin{aligned} \quad & \xi^\ast= max \quad   \xi
\\\mathrm{\ s.t.\ }
& v\boldsymbol X+py\boldsymbol Y_{iq}+pb\boldsymbol B_{ir}\leq 0 \\
& v_p\leq -\frac{1}{P}\left[1/x_{op}\right] , p=1,...,P \\
& \frac{\xi}{R+Q}\left[1/y_{oq}^g\right]  - u^g_q \leq 0 , q=1,...,Q \\
& \frac{\xi}{R+Q}\left[1/y_{or}^b\right]  + u^b_q \leq 0,r=1,...,R \\
& \xi+\boldsymbol{vx_o}-\boldsymbol {u^g y_o^g}+\boldsymbol{u^b y_o^b}=1 \\
& \xi, v, u^g,u^b 无约束
\end{aligned}.
$$

## 13.3 基于SBM模型对偶问题的影子价格计算程序

接下来，我们以Ex4.dta数据为例，写出求解以上线性规划问题的Python代码。具体代码如下：

```python
# 导入优化建模和数据处理库
from pyomo.environ import *  # 导入Pyomo优化建模环境
import pandas as pd ; import numpy as np ; import re  # 导入数据处理、数值计算和正则表达式库
import warnings  # 导入警告控制模块
warnings.filterwarnings("ignore")  # 忽略所有警告信息
ex4 = pd.read_stata(r"../../data/Ex4.dta")  # 读取Stata数据文件
data  = ex4[["K","L","Y","CO2"]]  # 提取投入、产出和非期望产出变量
K=2
L=1
M=1
xref=data.iloc[:,0:K]
yref=data.iloc[:,K:K+L]
bref=data.iloc[:,K+L:K+L+M]
x=data.iloc[0,0:K]
y=data.iloc[0,K:K+L]
b=data.iloc[0,K+L:K+L+M]
model = ConcreteModel()
model.I = Set(initialize = range(xref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
model.K = Set(initialize = range(len(x))) # 采用列表初始化投入变量个数的集合
model.L = Set(initialize = range(len(y))) # 采用列表初始化产出变量个数的集合
model.M = Set(initialize = range(len(b))) # 采用列表初始化产出变量个数的集合
model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
def objective_rule(model):
    """Return the proper objective function"""
    return model.pomega*1
def first_rule(model, i):
    """Return the proper input constraint"""
    return -sum(model.px[k] * xref.iloc[i,k] for k in model.K
            ) + sum(model.py[l] * yref.iloc[i,l] for l in model.L
            ) - sum(model.pb[m] * bref.iloc[i,m] for m in model.M)    <= 0
def second_rule(model,k):
    """Return the proper second constraint"""
    return -model.px[k] <= -1 / (len(model.K)*x[k] )
def third_rule(model,l):
    """Return the proper third constraint"""
    return model.pomega / (len(model.L)+len(model.M)) /y[l]-  model.py[l]<=0
def forth_rule(model,m):
    """Return the proper forth constraint"""
    return model.pomega / (len(model.L)+len(model.M)) /b[m] - model.pb[m]<=0
def fifth_rule(model):
    """Return the proper fifth constraint"""
    return model.pomega + sum(model.px[k] * x[k] for k in model.K
            ) - sum(model.py[l] * y[l] for l in model.L
            ) + sum(model.pb[m] * b[m] for m in model.M)<=1                           
model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')
model.fifth = Constraint(         rule= fifth_rule , doc='fifth constraint')             
opt = SolverFactory('mosek') # 指定 mosek 作为求解器
solution = opt.solve(model) # 调用求解器求解
px = np.asarray(list(model.px[:].value))
py = np.asarray(list(model.py[:].value))
pb = np.asarray(list(model.pb[:].value))
obj = value(model.obj) 
print("optimum shadow price x: \n {} ".format(px))
print("optimum shadow price y: \n {} ".format(py))
print("optimum shadow price b: \n {} ".format(pb))
print("optimal object: {}".format(obj))
```

需要特别指出的是，在以上代码的第19-22行，我们对决策变量施加了下界为0的约束。而在原来的对偶问题中，决策变量是没有上下界限制的。这主要是因为从经济学角度来讲，影子价格不太可能小于零，为了计算得到的影子价格具有良好的经济学含义，我们进一步施加了这个约束。以上代码与求解DEA模型的代码没有太大的差异，这里就不再进行详细的解释了。运行上述代码后，我们便可以得到各个变量对应的影子价格。Python运行结果会显示第1个决策单元投入和产出变量$(\boldsymbol {X,Y,B})$对应的影子价格。

接下来我们对所有的样本数据进行迭代求解，获取所有决策单元生产要素的影子价格。在以下代码中，我们使用for语句进行循环。第12=15行，在每次循环中取出一个决策单元的投入、期望产出和非期望产出数据。第55行对线性规划的最优解是否存在进行判断，如果存在，则提出最优解的参数值。运行这部分代码，所有决策单元生产要素的影子价格会存储在变量$px$、$py$和$pb$中。

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
ex4 = pd.read_stata(r"../../data/Ex4.dta")
data  = ex4[["K","L","Y","CO2"]]
K=2
L=1
M=1
xref=data.iloc[:,0:K]
yref=data.iloc[:,K:K+L]
bref=data.iloc[:,K+L:K+L+M]
px,py,pb,obj = {},{},{},{}             
for j in range(data.shape[0]):
    x=data.iloc[j,0:K]
    y=data.iloc[j,K:K+L]
    b=data.iloc[j,K+L:K+L+M]
    model = ConcreteModel()
    model.I = Set(initialize = range(xref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
    model.K = Set(initialize = range(len(x))) # 采用列表初始化投入变量个数的集合
    model.L = Set(initialize = range(len(y))) # 采用列表初始化产出变量个数的集合
    model.M = Set(initialize = range(len(b))) # 采用列表初始化产出变量个数的集合
    model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
    model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
    model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
    model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
    def objective_rule(model):
        """Return the proper objective function"""
        return model.pomega*1
    def first_rule(model, i):
        """Return the proper input constraint"""
        return -sum(model.px[k] * xref.iloc[i,k] for k in model.K
                ) + sum(model.py[l] * yref.iloc[i,l] for l in model.L
                ) - sum(model.pb[m] * bref.iloc[i,m] for m in model.M)    <= 0
    def second_rule(model,k):
        """Return the proper second constraint"""
        return -model.px[k] <= -1 / (len(model.K)*x[k] )
    def third_rule(model,l):
        """Return the proper third constraint"""
        return model.pomega / (len(model.L)+len(model.M)) /y[l]-  model.py[l]<=0
    def forth_rule(model,m):
        """Return the proper forth constraint"""
        return model.pomega / (len(model.L)+len(model.M)) /b[m] - model.pb[m]<=0
    def fifth_rule(model):
        """Return the proper fifth constraint"""
        return model.pomega + sum(model.px[k] * x[k] for k in model.K
                ) - sum(model.py[l] * y[l] for l in model.L
                ) + sum(model.pb[m] * b[m] for m in model.M)<=1                           
    model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
    model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
    model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
    model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
    model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')
    model.fifth = Constraint(         rule= fifth_rule , doc='fifth constraint')             
    opt = SolverFactory('mosek') # 指定 mosek 作为求解器
    solution = opt.solve(model) # 调用求解器求解
    if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
        px[j] = np.asarray(list(model.px[:].value))
        py[j] = np.asarray(list(model.py[:].value))
        pb[j] = np.asarray(list(model.pb[:].value))
        obj[j] = value(model.obj) 
```

通过以上两个步骤，我们便较为完整的实现了SBM模型对偶问题的求解。为了代码使用的方便，我们接着进行Python函数编写。在以上代码的基础上，我们需要进一步考虑Python函数的变量和参数的传递问题。在变量和参数的传递上，我们需要传递投入和产出数据集、投入变量的个数和期望产出变量的个数。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def sbmdual( dataframe, varname ,P,Q):
    """
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    varname: 投入变量、期望产出变量和非期望产出变量的列表，如["K","L","Y","CO2"]。
    P: 投入变量的个数。
    Q:期望产出变量的个数。
    """
    data  = dataframe.loc[:,varname]
    K=P
    L=Q
    M=data.shape[1]-P-Q
    xref=data.iloc[:,0:K]
    yref=data.iloc[:,K:K+L]
    bref=data.iloc[:,K+L:K+L+M]
    xcol=varname[0:K]
    ycol=varname[K:K+L]
    bcol=varname[K+L:K+L+M]
    px,py,pb,obj = {},{},{},{}             
    for j in range(data.shape[0]):
        x=data.iloc[j,0:K]
        y=data.iloc[j,K:K+L]
        b=data.iloc[j,K+L:K+L+M]
        model = ConcreteModel()
        model.I = Set(initialize = range(xref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(x))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(y))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(b))) # 采用列表初始化产出变量个数的集合
        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
        model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
        model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.pomega*1
        def first_rule(model, i):
            """Return the proper input constraint"""
            return -sum(model.px[k] * xref.iloc[i,k] for k in model.K
                    ) + sum(model.py[l] * yref.iloc[i,l] for l in model.L
                    ) - sum(model.pb[m] * bref.iloc[i,m] for m in model.M)    <= 0
        def second_rule(model,k):
            """Return the proper second constraint"""
            return -model.px[k] <= -1 / (len(model.K)*x[k] )
        def third_rule(model,l):
            """Return the proper third constraint"""
            return model.pomega / (len(model.L)+len(model.M)) /y[l]-  model.py[l]<=0
        def forth_rule(model,m):
            """Return the proper forth constraint"""
            return model.pomega / (len(model.L)+len(model.M)) /b[m] - model.pb[m]<=0
        def fifth_rule(model):
            """Return the proper fifth constraint"""
            return model.pomega + sum(model.px[k] * x[k] for k in model.K
                    ) - sum(model.py[l] * y[l] for l in model.L
                    ) + sum(model.pb[m] * b[m] for m in model.M)<=1                           
        model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
        model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
        model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
        model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
        model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')
        model.fifth = Constraint(         rule= fifth_rule , doc='fifth constraint')             
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj) 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        pxdf = pd.DataFrame(px,index=xcol).T
        pxdf.columns = pxdf.columns.map(lambda x : "shadow price "+ str(x) ) 
        pydf = pd.DataFrame(py,index=ycol).T
        pydf.columns = pydf.columns.map(lambda y : "shadow price "+ str(y) ) 
        pbdf = pd.DataFrame(pb,index=bcol).T
        pbdf.columns = pbdf.columns.map(lambda b : "shadow price "+ str(b) ) 
        pdf=pd.concat([pxdf,pydf],axis=1)
        pdf=pd.concat([pdf,pbdf],axis=1)
        redata = pd.concat([objdf,pdf],axis=1)
    return redata
```

下面我们给出一个应用的例子：

```python
import pandas as pd
ex4 = pd.read_stata(r"../../data/Ex4.dta")
data = sbmdual( ex4, ["K","L","Y","CO2"] ,2,1)
data
```

现在我们可以对sbmdual函数增加样本和技术参照集的限定选项。为此，我们需要额外增加两个参数：evaquery和refquery，分别用来筛选待评价单元和筛选技术参照集。以下是具体的代码：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def sbmdual2( dataframe, varname ,P,Q,evaquery,refquery):
    """
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    varname: 投入变量、期望产出变量和非期望产出变量的列表，如["K","L","Y","CO2"]。
    P: 投入变量的个数。
    Q:期望产出变量的个数。
    evaquery: 传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery: 传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    indexlt = dataframe.query(evaquery).index
    indexltref = dataframe.query(refquery).index
    data = dataframe.loc[indexlt,varname]
    dataref = dataframe.loc[indexltref,varname]
    K=P
    L=Q
    M=data.shape[1]-P-Q
    xref=dataref.iloc[:,0:K]
    yref=dataref.iloc[:,K:K+L]
    bref=dataref.iloc[:,K+L:K+L+M]
    xcol=varname[0:K]
    ycol=varname[K:K+L]
    bcol=varname[K+L:K+L+M]
    px,py,pb,obj = {},{},{},{}             
    for j in data.index:
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        b=data.loc[j,bcol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(x))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(y))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(b))) # 采用列表初始化产出变量个数的集合
        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
        model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
        model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.pomega*1
        def first_rule(model, i):
            """Return the proper input constraint"""
            return -sum(model.px[k] * xref.iloc[i,k] for k in model.K
                    ) + sum(model.py[l] * yref.iloc[i,l] for l in model.L
                    ) - sum(model.pb[m] * bref.iloc[i,m] for m in model.M)    <= 0
        def second_rule(model,k):
            """Return the proper second constraint"""
            return -model.px[k] <= -1 / (len(model.K)*x[k] )
        def third_rule(model,l):
            """Return the proper third constraint"""
            return model.pomega / (len(model.L)+len(model.M)) /y[l]-  model.py[l]<=0
        def forth_rule(model,m):
            """Return the proper forth constraint"""
            return model.pomega / (len(model.L)+len(model.M)) /b[m] - model.pb[m]<=0
        def fifth_rule(model):
            """Return the proper fifth constraint"""
            return model.pomega + sum(model.px[k] * x[k] for k in model.K
                    ) - sum(model.py[l] * y[l] for l in model.L
                    ) + sum(model.pb[m] * b[m] for m in model.M)<=1                           
        model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
        model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
        model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
        model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
        model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')
        model.fifth = Constraint(         rule= fifth_rule , doc='fifth constraint')             
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj) 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        pxdf = pd.DataFrame(px,index=xcol).T
        pxdf.columns = pxdf.columns.map(lambda x : "shadow price "+ str(x) ) 
        pydf = pd.DataFrame(py,index=ycol).T
        pydf.columns = pydf.columns.map(lambda y : "shadow price "+ str(y) ) 
        pbdf = pd.DataFrame(pb,index=bcol).T
        pbdf.columns = pbdf.columns.map(lambda b : "shadow price "+ str(b) ) 
        pdf=pd.concat([pxdf,pydf],axis=1)
        pdf=pd.concat([pdf,pbdf],axis=1)
        redata = pd.concat([objdf,pdf],axis=1)
    return redata
```

最后，我们使用另一种调用命令的形式，即使用“期望产出变量 : 非期望产出变量=投入变量”的语法调用sbmdual函数。关于这部分，我们可以参考上一章中sbm命令的代码。

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def sbmdual3( formula, dataframe, evaquery, refquery):
    """
    formula: 产出变量:非期望产出变量=投入变量，如“ Y :CO2  = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    evaquery: 传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery: 传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    px,py,pb,obj = {},{},{},{}  
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
    for j in data.index:
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        b=data.loc[j,bcol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(x))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(y))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(b))) # 采用列表初始化产出变量个数的集合
        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
        model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
        model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
        def objective_rule(model):
            """Return the proper objective function"""
            return model.pomega*1
        def first_rule(model, i):
            """Return the proper input constraint"""
            return -sum(model.px[k] * xref.iloc[i,k] for k in model.K
                    ) + sum(model.py[l] * yref.iloc[i,l] for l in model.L
                    ) - sum(model.pb[m] * bref.iloc[i,m] for m in model.M)    <= 0
        def second_rule(model,k):
            """Return the proper second constraint"""
            return -model.px[k] <= -1 / (len(model.K)*x[k] )
        def third_rule(model,l):
            """Return the proper third constraint"""
            return model.pomega / (len(model.L)+len(model.M)) /y[l]-  model.py[l]<=0
        def forth_rule(model,m):
            """Return the proper forth constraint"""
            return model.pomega / (len(model.L)+len(model.M)) /b[m] - model.pb[m]<=0
        def fifth_rule(model):
            """Return the proper fifth constraint"""
            return model.pomega + sum(model.px[k] * x[k] for k in model.K
                    ) - sum(model.py[l] * y[l] for l in model.L
                    ) + sum(model.pb[m] * b[m] for m in model.M)<=1                           
        model.obj = Objective(rule=objective_rule, sense=maximize, doc='objective function')
        model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
        model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
        model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
        model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')
        model.fifth = Constraint(         rule= fifth_rule , doc='fifth constraint')             
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj) 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        pxdf = pd.DataFrame(px,index=xcol).T
        pxdf.columns = pxdf.columns.map(lambda x : "shadow price "+ str(x) ) 
        pydf = pd.DataFrame(py,index=ycol).T
        pydf.columns = pydf.columns.map(lambda y : "shadow price "+ str(y) ) 
        pbdf = pd.DataFrame(pb,index=bcol).T
        pbdf.columns = pbdf.columns.map(lambda b : "shadow price "+ str(b) ) 
        pdf=pd.concat([pxdf,pydf],axis=1)
        pdf=pd.concat([pdf,pbdf],axis=1)
        redata = pd.concat([objdf,pdf],axis=1)
    return redata
```

下面是调用sbmdual3命令的例子：

```python
import pandas as pd
ex4 = pd.read_stata(r"../../data/Ex4.dta")
data = sbmdual3( " Y:    CO2 = K     L ", ex4, "t==[1,2]","t==[1,2,3]")
data
```

## 13.4 NDDF模型的对偶问题

首先，给出一个用以估计NDDF模型的线性规划问题如下：
$$
\begin{array}{l}
\vec{D}(\boldsymbol{x}, \boldsymbol{y}, \boldsymbol{b};\boldsymbol{g})=\max  &\boldsymbol{\omega_{x}' \beta_{X}}+\boldsymbol{\omega_{y}' \beta_{Y}}+\boldsymbol{\omega_{b}' \beta_{B}} \\
\text { s.t. } \quad  &\sum_{i=1}^{N} \lambda_{i} \boldsymbol{X_{i}} \leq \boldsymbol{x}+diag(\boldsymbol{g_{X}})\boldsymbol{\beta_{X}} \\
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol{Y_{i}} \geq \boldsymbol{y}+diag(\boldsymbol{g_{Y}})\boldsymbol{\beta_{Y}} \\
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol{B_{i}}=\boldsymbol{b}+diag(\boldsymbol{g_{B}})\boldsymbol{\beta_{B}} \\
&\lambda_{i} \geq 0, \beta\geq0 ,i=1, \ldots, N 
\end{array}
$$
为了便于确定对偶问题，我们把在目标函数中的决策变量$\lambda_{i}$的系数补上$\boldsymbol0$，将目标函数和约束条件中的决策变量按顺序排列，并将约束条件的符号变为$\leq$。因此，上述模型可以变换为：
$$
\begin{array}{l}
\vec{D}(\boldsymbol{x}, \boldsymbol{y}, \boldsymbol{b};\boldsymbol{g})=\max  &\sum_{i=1}^{N} \lambda_{i}\boldsymbol0+\boldsymbol{\omega_{x}' \beta_{X}}+\boldsymbol{\omega_{y}' \beta_{Y}}+\boldsymbol{\omega_{b}' \beta_{B}} \\
\text { s.t. } \quad 
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol{X_{i}}-diag(\boldsymbol{g_{X}})\boldsymbol{\beta_{X}}+\boldsymbol{0}\boldsymbol{\beta_{Y}} +\boldsymbol{0}\boldsymbol{\beta_{B}}  \leq \boldsymbol{x} \\
&-\sum_{i=1}^{N} \lambda_{i} \boldsymbol{Y_{i}}+\boldsymbol{0}\boldsymbol{\beta_{X}} +diag(\boldsymbol{g_{Y}})\boldsymbol{\beta_{Y}} +\boldsymbol{0}\boldsymbol{\beta_{B}}  \leq -\boldsymbol{y}\\
&\sum_{i=1}^{N} \lambda_{i} \boldsymbol{B_{i}}+\boldsymbol{0}\boldsymbol{\beta_{X}} +\boldsymbol{0}\boldsymbol{\beta_{Y}} -diag(\boldsymbol{g_{B}})\boldsymbol{\beta_{B}}=\boldsymbol{b} \\
&\lambda_{i} \geq 0, \beta\geq0 ,i=1, \ldots, N 
\end{array}
$$
按照线性规划的原问题与对偶问题的关系，我们便可以写出对偶问题的具体表达式，即：
$$
\begin{array}{l}
\min & \boldsymbol{x'v}_x -  \boldsymbol{y'u}_y +  \boldsymbol{b'u}_b  \\
\text { s.t. } \quad  
&\boldsymbol{Xv_X}- \boldsymbol{Yu_Y} + \boldsymbol{Bu_B} \geq \boldsymbol{0} \\
&-diag(\boldsymbol{g_{X}})\boldsymbol{v_{X}} \geq \boldsymbol{\omega_X} \\
&diag(\boldsymbol{g_{Y}})\boldsymbol{u_{Y}} \geq \boldsymbol{\omega_Y} \\
&-diag(\boldsymbol{g_{B}})\boldsymbol{u_{B}} \geq \boldsymbol{\omega_B} \\
&\boldsymbol{v}_x ,\boldsymbol{u}_y,\boldsymbol{u}_b\geq 0
\end{array}
$$

## 13.5 基于NDDF模型对偶问题的影子价格计算程序

我们以Ex5.dta数据为例，编写求解以上线性规划问题的Python代码。具体代码如下：

```python
ex5 = pd.read_stata(r"../../data/Ex5.dta")
data  = ex5[['labor','capital','energy','gdp','co2']]
K=3
L=1
M=1
xref=data.iloc[:,0:K]
yref=data.iloc[:,K:K+L]
bref=data.iloc[:,K+L:K+L+M]
x=data.iloc[0,0:K]
y=data.iloc[0,K:K+L]
b=data.iloc[0,K+L:K+L+M]
xcol=xref.columns
ycol=yref.columns
bcol=bref.columns
gx=[-1]*len(xcol)
gy=[ 1]*len(ycol)
gb=[-1]*len(bcol)
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
model = ConcreteModel()
model.I = Set(initialize = range(xref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
def objective_rule(model):
    """Return the proper objective function"""
    return sum(model.px[k]*x[k] for k in model.K  
        ) - sum(model.py[l]*y[l] for l in model.L 
        ) + sum(model.pb[j]*b[j] for j in model.M 
        ) + model.pomega*1 
def first_rule(model, i):
    """Return the proper first constraint"""
    return sum(model.px[k] * xref.loc[i,xcol[k]] for k in model.K
        ) - sum(model.py[l] * yref.loc[i,ycol[l]] for l in model.L
        ) + sum(model.pb[m] * bref.loc[i,bcol[m]] for m in model.M
        ) + model.pomega*1    >=0
def second_rule(model, k):
    """Return the proper second constraint"""
    return  -gx[k]*x[k] * model.px[k] >= iweight[k] 
def third_rule(model, l):
    """Return the proper third constraint"""
    return  gy[l]*y[l] * model.py[l] >= oweight[l]  
def forth_rule(model, j):
    """Return the proper forth constraint"""
    return  -gb[j]*b[j] * model.pb[j] >= bweight[j]                           
model.obj = Objective(rule=objective_rule, sense=minimize, doc='objective function')
model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')         
opt = SolverFactory('mosek') # 指定 mosek 作为求解器
solution = opt.solve(model) # 调用求解器求解
px = np.asarray(list(model.px[:].value))
py = np.asarray(list(model.py[:].value))
pb = np.asarray(list(model.pb[:].value))
obj = value(model.obj) 
print("optimum shadow price x: \n {} ".format(px))
print("optimum shadow price y: \n {} ".format(py))
print("optimum shadow price b: \n {} ".format(pb))
print("optimal object: {}".format(obj))
```

在以上代码中第38=58行是按照该线性规划问题的目标函数和不等式约束逐步设定的，这里就不再展开解释了。接下来我们在上述代码中加上循环，对决策单元的数据进行迭迭代，用以求解所有决策单元投入和产出变量的影子价格。

```python
ex5 = pd.read_stata(r"../../data/Ex5.dta")
data  = ex5[['labor','capital','energy','gdp','co2']]
K=3
L=1
M=1
xref=data.iloc[:,0:K]
yref=data.iloc[:,K:K+L]
bref=data.iloc[:,K+L:K+L+M]
xcol=xref.columns
ycol=yref.columns
bcol=bref.columns
gx=[-1]*len(xcol)
gy=[ 1]*len(ycol)
gb=[-1]*len(bcol)
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
px,py,pb,obj={},{},{},{}
for j in range(data.shape[0]):
    x=data.iloc[j,0:K]
    y=data.iloc[j,K:K+L]
    b=data.iloc[j,K+L:K+L+M]
    model = ConcreteModel()
    model.I = Set(initialize = range(xref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
    model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
    model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
    model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
    model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
    model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
    model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
    model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
    def objective_rule(model):
        """Return the proper objective function"""
        return sum(model.px[k]*x[k] for k in model.K  
            ) - sum(model.py[l]*y[l] for l in model.L 
            ) + sum(model.pb[j]*b[j] for j in model.M 
            ) + model.pomega*1 
    def first_rule(model, i):
        """Return the proper first constraint"""
        return sum(model.px[k] * xref.loc[i,xcol[k]] for k in model.K
            ) - sum(model.py[l] * yref.loc[i,ycol[l]] for l in model.L
            ) + sum(model.pb[m] * bref.loc[i,bcol[m]] for m in model.M
            ) + model.pomega*1    >=0
    def second_rule(model, k):
        """Return the proper second constraint"""
        return  -gx[k]*x[k] * model.px[k] >= iweight[k] 
    def third_rule(model, l):
        """Return the proper third constraint"""
        return  gy[l]*y[l] * model.py[l] >= oweight[l]  
    def forth_rule(model, j):
        """Return the proper forth constraint"""
        return  -gb[j]*b[j] * model.pb[j] >= bweight[j]                           
    model.obj = Objective(rule=objective_rule, sense=minimize, doc='objective function')
    model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
    model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
    model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
    model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')         
    opt = SolverFactory('mosek') # 指定 mosek 作为求解器
    solution = opt.solve(model) # 调用求解器求解
    if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
        px[j] = np.asarray(list(model.px[:].value))
        py[j] = np.asarray(list(model.py[:].value))
        pb[j] = np.asarray(list(model.pb[:].value))
        obj[j] = value(model.obj) 
```

以上代码中，第27=30行，在每次循环中取出一个决策单元的投入、期望产出和非期望产出数据。第68行对线性规划的最优解是否存在进行判断，如果存在，则提出最优解的参数值。运行这部分代码，所有决策单元生产要素的影子价格会存储在变量$px$、$py$和$pb$中。通过以上两个步骤，我们便较为完整的实现了NDDF模型对偶问题的求解。为了代码使用的方便，我们接着进行Python函数编写。在以上代码的基础上，我们便可封装成Python函数。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def nddfdual( dataframe, varname ,P,Q):
    """
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    varname: 投入变量、期望产出变量和非期望产出变量的列表，如["K","L","Y","CO2"]。
    P: 投入变量的个数。
    Q:期望产出变量的个数。
    """
    data  = dataframe.loc[:,varname]
    K=P
    L=Q
    M=data.shape[1]-P-Q
    xref=data.iloc[:,0:K]
    yref=data.iloc[:,K:K+L]
    bref=data.iloc[:,K+L:K+L+M]
    xcol=xref.columns
    ycol=yref.columns
    bcol=bref.columns
    gx=[-1]*len(xcol)
    gy=[ 1]*len(ycol)
    gb=[-1]*len(bcol)
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
    px,py,pb,obj={},{},{},{}
    for j in range(data.shape[0]):
        x=data.iloc[j,0:K]
        y=data.iloc[j,K:K+L]
        b=data.iloc[j,K+L:K+L+M]
        model = ConcreteModel()
        model.I = Set(initialize = range(xref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
        model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
        model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
        def objective_rule(model):
            """Return the proper objective function"""
            return sum(model.px[k]*x[k] for k in model.K  
                ) - sum(model.py[l]*y[l] for l in model.L 
                ) + sum(model.pb[j]*b[j] for j in model.M 
                ) + model.pomega*1 
        def first_rule(model, i):
            """Return the proper first constraint"""
            return sum(model.px[k] * xref.loc[i,xcol[k]] for k in model.K
                ) - sum(model.py[l] * yref.loc[i,ycol[l]] for l in model.L
                ) + sum(model.pb[m] * bref.loc[i,bcol[m]] for m in model.M
                ) + model.pomega*1    >=0
        def second_rule(model, k):
            """Return the proper second constraint"""
            return  -gx[k]*x[k] * model.px[k] >= iweight[k] 
        def third_rule(model, l):
            """Return the proper third constraint"""
            return  gy[l]*y[l] * model.py[l] >= oweight[l]  
        def forth_rule(model, j):
            """Return the proper forth constraint"""
            return  -gb[j]*b[j] * model.pb[j] >= bweight[j]                           
        model.obj = Objective(rule=objective_rule, sense=minimize, doc='objective function')
        model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
        model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
        model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
        model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')         
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj) 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        pxdf = pd.DataFrame(px,index=xcol).T
        pxdf.columns = pxdf.columns.map(lambda x : "shadow price "+ str(x) ) 
        pydf = pd.DataFrame(py,index=ycol).T
        pydf.columns = pydf.columns.map(lambda y : "shadow price "+ str(y) ) 
        pbdf = pd.DataFrame(pb,index=bcol).T
        pbdf.columns = pbdf.columns.map(lambda b : "shadow price "+ str(b) ) 
        pdf=pd.concat([pxdf,pydf],axis=1)
        pdf=pd.concat([pdf,pbdf],axis=1)
        redata = pd.concat([objdf,pdf],axis=1)
    return redata
```

现在我们可以对nddfdual函数增加样本和技术参照集的限定选项。为此，我们需要额外增加两个参数：evaquery和refquery，分别用来筛选待评价单元和筛选技术参照集。另外，我们需要添加对方向向量和权重自定义的选项。为此，我们需要额外增加参数gx、gy、gb和weight。另外，我们使用“投入变量=期望产出变量：非期望产出变量”的语法传递参数。以下是具体的代码：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def nddfdual2(formula, dataframe, gx=None, gy=None, gb=None,weight =None, evaquery=None, refquery=None ):
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
    px,py,pb,obj={},{},{},{}
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
        model.I = Set(initialize = range(xref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(len(xcol))) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.M = Set(initialize = range(len(bcol))) # 采用列表初始化产出变量个数的集合
        model.px = Var(model.K, bounds=(0, None), within=Reals, doc='shadow price of x')
        model.py = Var(model.L, bounds=(0, None), within=Reals, doc='shadow price of y')
        model.pb = Var(model.M, bounds=(0, None), within=Reals,doc='shadow price of b')
        model.pomega = Var( bounds=(0, None), within=Reals,doc='shadow price of 1')
        def objective_rule(model):
            """Return the proper objective function"""
            return sum(model.px[k]*x[k] for k in model.K  
                ) - sum(model.py[l]*y[l] for l in model.L 
                ) + sum(model.pb[j]*b[j] for j in model.M 
                ) + model.pomega*1 
        def first_rule(model, i):
            """Return the proper first constraint"""
            return sum(model.px[k] * xref.loc[i,xcol[k]] for k in model.K
                ) - sum(model.py[l] * yref.loc[i,ycol[l]] for l in model.L
                ) + sum(model.pb[m] * bref.loc[i,bcol[m]] for m in model.M
                ) + model.pomega*1    >=0
        def second_rule(model, k):
            """Return the proper second constraint"""
            return  -gx[k]*x[k] * model.px[k] >= iweight[k] 
        def third_rule(model, l):
            """Return the proper third constraint"""
            return  gy[l]*y[l] * model.py[l] >= oweight[l]  
        def forth_rule(model, j):
            """Return the proper forth constraint"""
            return  -gb[j]*b[j] * model.pb[j] >= bweight[j]                           
        model.obj = Objective(rule=objective_rule, sense=minimize, doc='objective function')
        model.first = Constraint(model.I, rule= first_rule, doc='first constraint')
        model.second =Constraint(model.K, rule= second_rule , doc='second constraint')
        model.third = Constraint(model.L, rule= third_rule , doc='third constraint')
        model.forth = Constraint(model.M, rule= forth_rule , doc='forth constraint')         
        opt = SolverFactory('mosek') # 指定 mosek 作为求解器
        solution = opt.solve(model) # 调用求解器求解
        if solution.solver.termination_condition == "optimal":     # 终止条件 一般包括三种 optimal, feasible, infeasible
            px[j] = np.asarray(list(model.px[:].value))
            py[j] = np.asarray(list(model.py[:].value))
            pb[j] = np.asarray(list(model.pb[:].value))
            obj[j] = value(model.obj) 
        objdf= pd.DataFrame(obj,index=["obj"]).T
        pxdf = pd.DataFrame(px,index=xcol).T
        pxdf.columns = pxdf.columns.map(lambda x : "shadow price "+ str(x) ) 
        pydf = pd.DataFrame(py,index=ycol).T
        pydf.columns = pydf.columns.map(lambda y : "shadow price "+ str(y) ) 
        pbdf = pd.DataFrame(pb,index=bcol).T
        pbdf.columns = pbdf.columns.map(lambda b : "shadow price "+ str(b) ) 
        pdf=pd.concat([pxdf,pydf],axis=1)
        pdf=pd.concat([pdf,pbdf],axis=1)
        redata = pd.concat([objdf,pdf],axis=1)
    return redata
```

我们可以看到以上代码与nddf命令基本一致，因为nddf命令和nddfdual命令需要传递的变量和参数是一样的。不同的地方在于求解的线性规划问题不同，主要体现在第65-85行定义了不同的目标函数和约束条件。下面是命令调用的例子：

```python
import pandas as pd
ex5 = pd.read_stata(r"../../data/Ex5.dta")
data = nddfdual2(' gdp : co2=labor capital energy ', ex5, gx=[-1,-1,-1],gy=[1],gb=[-1],weight=None, )
data
```


## 参考文献
[^1]:Tone, K. (2004). A conventional scheme for coping with negative output data in DEA : a slacks-based measure (SBM) approach(DEA(1)). 日本オペレーションズ・リサーチ学会秋季研究発表会アブストラクト集, 2004, 298–299.



