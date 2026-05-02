[TOC]

# 11. DEA模型编程的基本思路

在前面的章节，我们介绍了多个效率评价的DEA模型。虽然这些DEA模型形式各不相同，但他们的基本结构是一样的，都由以下两部分组成：技术集和效率评价准则。在数学表达式上，他们共同形成了一个线性规划问题。效率评价准则主要体现在线性规划问题的目标函数上，而技术集主要体现在线性规划问题的约束条件上。从本章开始，我们将逐步介绍编写一个DEA模型的Python函数。

## 11.1 DEA模型编写的基本步骤

根据DEA模型的结构和Python软件的特点，我们可以将编写DEA程序的过程划分为以下三个步骤：
步骤1：将一个决策单元的DEA模型转化为软件的线性规划求解器能识别的线性规划模型。
步骤2：将将步骤1在所有待评价决策单元的数据集上进行循环迭代。
步骤3：将步骤2封装成Python函数。

## 11.2 步骤1：转化一个决策单元的线性规划问题

在DEA模型中，通常将效率评价的测度转化为一个线性规划问题。我们可以使用Python的pyomo包来写一个线性规划问题，并求解。读者可以访问[pyomo的文档](https://pyomo.readthedocs.io/en/stable/installation.html)来查看说明文档和技术细节。我们下面以求解基于谢泼德产出导向距离函数的技术效率为例，进行详细说明。假定给定的数据集如表12.1所示。

表11.1 给定数据集

| 决策单元 | 投入  |       |       | 产出  |       |
| -------- | ----- | ----- | ----- | ----- | ----- |
|          | $x_1$ | $x_2$ | $x_3$ | $y_1$ | $y_2$ |
| A        | 5     | 6     | 7     | 23    | 67    |
| B        | 5     | 9     | 7     | 14    | 34    |
| C        | 8     | 6     | 87    | 45    | 12    |
| D        | 34    | 67    | 32    | 11    | 37    |
| E        | 9     | 6     | 19    | 31    | 29    |

我们将不变规模报酬技术下的谢泼德产出导向距离函数表示为如下线性规划问题：
$$
\begin{aligned} D(\mathbf{X}^*,\mathbf{Y}^*)^{-1}&=\max_{\theta,\lambda_1,\lambda_2,...,\lambda_5}  \quad  \theta  \\s.t.& \sum_{i=1}^{5} \lambda_{i}X_{ki} \leq  X_k^*,  k=1,2,3 \\  & \sum_{i=1}^{5} \lambda_{i}Y_{i}  \geq  \theta Y_k^* ,l=1,2  \\& \lambda_{i} \geq 0, n=1, 2,\ldots, 5 ;-\infty \leq \theta\leq +\infty \end{aligned}
$$

其中，$X^* \in R^3_+$和$Y^* \in R^2_+$分别为评价决策单元的投入和产出向量。

具体而言，当被评价决策单元为A，技术参照决策单元为A、B、C、D、E全体时，我们可以通过以下流程对DEA模型的线性规划问题进行转换：

(1)确定决策变量，包括确定决策变量的维度，确定维度是一维还是多维的，以及决策变量的上下界；

(2)确定目标函数，包括求最大值还是最小值，以及变量对应的系数；

(3)确定约束条件，包括约束条件的维度，确定维度是一维的还是多维的，以及变量对应的系数和约束符号（$\geq$、$\leq$和$=$）。

在上述谢泼德产出导向距离函数中：

(1)变量为$\theta$和$\lambda_i$，其中，$\theta$是一维的，$\lambda_i$是多维的，$\lambda_i$的维度为$1×5$；$\theta$无约束，且$\lambda_i\geq0$；

(2)目标函数中变量$\theta$和$\lambda_i$对应的系数分别为1和$\mathbf{0}_{1\times 5}$；求最小值；

(3)约束条件有两个，第一个约束对应的行数有$3$个，$3$为变量$X$的个数，第二个约束对应的行数有$2$个，$2$为变量$Y$的个数；且第一个约束中变量$\lambda_i$对应的系数为全部决策单元的投入变量，即$x_1、x_2、x_3$的第1-5行，$\theta$对应的系数为被评价决策单元A的投入变量，即$x_1、x_2、x_3$的第1行，第二个约束中变量$\lambda_i$对应的系数为全部决策单元的产出变量，即$y_1、y_2$的第1-5行；第一个约束中约束符号为$\leq$，第二个约束中约束符号为$\geq$。

在确定上述决策变量、约束和目标函数后，我们就可以用代码实现优化模型，我们需要将变量、约束和目标函数转换为pyomo定义的语言。在转换之前，需要考虑另外一个不可或缺的、在pyomo模块中的要素：索引集合。索引集合通常被用作在定义变量和约束时，声明变量和约束的索引，从而可以通过索引调用变量和约束的具体内容。索引集合可以采用如下方法调用：

```python
import numpy as np
from pyomo.environ import *
model = ConcreteModel()
model.I = Set(initialize = range(5)) # 采用列表初始化技术参照决策单元个数的集合
model.K = Set(initialize = range(3)) # 采用列表初始化投入变量个数的集合
model.L = Set(initialize = range(2)) # 采用列表初始化产出变量个数的集合
```

常用参数initialize表示初始化参数，可以采用列表，元组和集合来给一个Set初始化。上面是使用列表集合进行初始化。可以通过data属性查看索引集合的内容。

```
print(model.I.data())
print(model.K.data())
print(model.L.data())
```

(1)决策变量

在定义完索引后，使用Var方法定义决策变量。常用的参数有四个：index_set表示决策变量的索引，这个索引可以是上面我们定义的索引集合；bounds表示决策变量上下界；within表示决策变量类型；doc为该决策变量取名。bounds传入一个元组，元组的第一个值表示变量的下界，第二个值表示变量的上界，且都包含等号，$\infty$对应传入None；within 默认为实数，可以选择的类型有：

|       传入值        |     含义     |
| :-----------------: | :----------: |
|        Reals        |     实数     |
|    PositiveReals    |    正实数    |
|  NonPositiveReals   |   非正实数   |
|    NegativeReals    |    负实数    |
|  NonNegativeReals   |   非负实数   |
|   PercentFraction   | 0到1之间实数 |
|    UnitInterval     | 0到1之间实数 |
|      Integers       |     整数     |
|  PositiveIntegers   |    正整数    |
| NonPositiveIntegers |   非正整数   |
|  NegativeIntegers   |    负整数    |
| NonNegativeIntegers |   非负整数   |
|       Boolean       |   布尔变量   |
|       Binary        |   0,1变量    |



由于DEA对应的规划问题为线性规划问题，且变量的范围为正实数或非负实数，所以within一般不需要传入参数，只需要在bounds中传入变量的上下界即可。决策变量是一维的还是多维的对应Var函数调用的两种方法：当决策变量是一维变量时，需要只传入一个索引变量；而当决策变量是多维变量时，需要根据决策变量的维度传入多个索引变量。比如，下面梁行代码分别定义了一维变量和二维变量：

```python
model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency')
model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables')
```

(2)目标函数

使用Objective方法定义决策变量。常用参数有三个：obj_rule表示目标函数定义规则；sense表示最大化问题还是最小化问题，maximize代表是极大化问题，minimize代表极小化问题；doc为该目标函数取名。下面的代码定义了上述DEA模型的目标函数。由于$\lambda_i$的系数为0，所以可以省略。

```python
def obj_rule(model):
    return model.theta
    #return model.theta*1  + sum(model.lamda[i] *0 for i in model.I)
model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function')
```

(3)约束条件

使用Constraint方法定义决策变量。参数有三个：index_set表示约束索引集合，可以是上面定义的索引集合；rule表示约束定义规则，可以传入一个函数表达式，一般用来定义多维约束；expr表示约束表达式，一般用来定义一维约束。下面使用rule定义上述DEA模型的约束，这是一个多维约束。此时需要用到技术参照决策单元和被评价决策单元的投入产出变量作为决策变量的系数。

```python
data=np.array([[ 5,  6,  7, 23, 67],
               [ 5,  9,  7, 14, 34],
               [ 8,  6, 87, 45, 12],
               [34, 67, 32, 11, 37],
               [ 9,  6, 19, 31, 29]])
x=data[0,0:3]
y=data[0,3:5]
xref=data[:,0:3]
yref=data[:,3:5]
def input_rule(model,k):
	"""Return the proper input constraint"""
	return sum(model.lamda[i]*xref[i,k] for i in model.I) <= x[k]
def output_rule(model,l):
    """Return the proper output constraint"""
    return sum(model.lamda[i]*yref[i,l] for i in model.I) >=model.theta*y[l]
model.input = Constraint(model.K, rule=input_rule, doc='input constraint')
model.output = Constraint(model.L, rule=output_rule, doc='output constraint')
```

下面的例子定义模型model2中的一维约束。模型model2中变量的维度为$4\times 5$，第一个一维约束为$x_{11}+4x_{13} \geq 2$，第二个一维约束为$3x_{34}+4x_{45} = 10$。需要注意python中的下标从0开始。

```python
model2 = ConcreteModel()
model2.I = Set(initialize = range(4))
model2.J = Set(initialize = range(5))
model2.x = Var(model2.I, model2.J, within = Reals)
model2.c1 = Constraint(expr=2.0 <= model2.x[0,0]+4*model2.x[0,2]) # 不等式约束 2 <= x[0,0] + 4*x[0,2]
model2.c2 = Constraint(expr=3*model2.x[2,3]+5*model2.x[3,4]==10.0) # 等式约束 3*x[2,3] + 5*x[3,4] = 10
```

可以使用Constraint.Skip方法跳过使用某一个索引值对应的变量来组成约束。如下面定义的多维约束为$\sum_{j,j \neq 1 } x_{ij} \leq 2, i=1,2,3,4 $：

```python
def c3_rule(model2, i):
    if j==1:
        return Objective.Skip
    return sum([model2.p[i,j]*model2.x[i,j] for j in model2.J]) <= 2
model2.c3= Constraint(model2.I, rule = c3_rule)
```

(4)求解

Pyomo 只能完成建模的任务，求解的需要将 Pyomo 的模型导入到专业的优化求解器来实现，例如glpk、 ipopt、Cplex、SCIP、mosek等。如下代码就是调用glpk求解器来实现求解：

```python
opt = SolverFactory('glpk') # 指定 glpk 作为求解器
solution = opt.solve(model) # 调用求解器求解
```

另外，可以在solve方法中传入“report_timing=True”来查看模型运行的时间：

```python
opt = SolverFactory('glpk') # 指定 glpk 作为求解器
solution = opt.solve(model,report_timing=True) # 调用求解器求解
```

如下表格展示了各种常见求解器求解一个被评价决策单元的效率所需的时间。可以发现mosek求解器消耗的时间分别是mosek、ipopt、scip求解器的1/6、1/8和1/10，因此mosek求解器被推荐使用。

```python
opt = SolverFactory('mosek') # 指定 mosek 作为求解器
solution = opt.solve(model,report_timing=True) # 调用求解器求解

opt = SolverFactory('ipopt') # 指定 ipopt 作为求解器
solution = opt.solve(model,report_timing=True) # 调用求解器求解

opt = SolverFactory('scip') # 指定 scip 作为求解器
solution = opt.solve(model,report_timing=True) # 调用求解器求解

opt = SolverFactory('mosek') # 指定 mosek 作为求解器
solution = opt.solve(model,report_timing=True) # 调用求解器求解
```

| 求解器 | 时间                                                         |      |
| ------ | ------------------------------------------------------------ | ---- |
| glpk   | ![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220919231604.png) |      |
| ipopt  | ![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220919231509.png) |      |
| scip   | ![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220919231635.png) |      |
| mosek  | ![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220920103822.png) |      |



在第9章中，我们有时需要使用非线性求解器，虽然免费的ipopt求解器可以求解非线性规划问题，但是其求解速度非常慢，因此，需要使用其他求解器来进行操作。AMPL（A Mathematical Programming Language）可以提供30天的学术许可证供我们使用，过期后可以使用学生身份申请增加期限。AMPL是一种专为数学规划问题设计的高级代数建模语言，旨在通过抽象的代数表达式简化复杂优化模型的构建。AMPL本身不直接求解问题，而是将模型转换为标准格式（如`.nl`文件）后调用外部求解器。支持的求解器包括Knitro、Gurobi、MINOS、IPOPT等，涵盖线性、非线性、整数规划等多种问题类型。我们使用AMPL提供的Knitro求解器来求解非线性规划问题。

AMPL的下载步骤如下：

第一步，注册并登录。

进入链接https://portal.ampl.com/account/ampl/，结果如下图所示，使用教育邮箱注册并登录。

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20250521182014492.png)

第二步，登录完成后，可以在链接https://portal.ampl.com/user/ampl/license/list中找到ampl在python中的安装步骤，如下图所示。

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20250521184644150.png)

```shell
# Install Python API for AMPL:
$ python -m pip install amplpy --upgrade

# Install solver modules:
$ python -m amplpy.modules install highs gurobi knitro

# Activate your license (e.g., free ampl.com/ce or ampl.com/courses licenses):
$ python -m amplpy.modules activate <your-license-uuid>

```

可以发现，改网站提供了上述安装步骤，只需在anaconda自带的`Anaconda Powershell Prompt`中输入上述命令，即可安装好非线性规划求解器knitro。

另外，还需要手动安装glpk以及ipopt求解器。

使用如下命令安装glpk求解器：

```python
conda install glpk
```

在Windows.上，安装ipopt很容易。只需要从ipopt的官方网站上(https::/github.com/coin-or/Ipopt)下载预编译好的 ipopt二进制文件，然后将它添加到环境变量中即可。在Linux和Mac OS X.上，可以从源代码编译Ipopt。编译方法请参考官方提供的文档。

如果没有安装好ipopt，则直接不实用ipopt求解器也是可以的，因为这个求解器对于非线性规划问题求解起来非常慢，而对于线性规划问题，我们可以直接使用mosek求解器。

对于scip求解器，可以参照https://blog.csdn.net/zx_glave/article/details/128397314进行安装。如果没有安装好scip求解器，也可以直接跳过这一步骤，因为我们可以使用mosek替代scip求解器。


(5)优化结果查看

```python
theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
obj = value(model.obj) # 提取目标函数 
print("optimum theta: \n {} ".format(theta))
print("optimum lamda: \n {} ".format(lamda))
print("optimal objective: {}".format(obj))
```

![图11.1 代码运行结果](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914160108.png)
**图11.1：代码运行结果**

一个决策单元的DEA线性规划问题的代码如下：

```python
from pyomo.environ import *
import numpy as np
data=np.array([[ 5,  6,  7, 23, 67],
               [ 5,  9,  7, 14, 34],
               [ 8,  6, 87, 45, 12],
               [34, 67, 32, 11, 37],
               [ 9,  6, 19, 31, 29]])
x=data[0,0:3]
y=data[0,3:5]
xref=data[:,0:3]
yref=data[:,3:5]
model = ConcreteModel()
model.I = Set(initialize = range(5)) # 采用列表初始化技术参照决策单元个数的集合
model.K = Set(initialize = range(3)) # 采用列表初始化投入变量个数的集合
model.L = Set(initialize = range(2)) # 采用列表初始化产出变量个数的集合
model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
def obj_rule(model):
    return model.theta 
    #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
def input_rule(model,k):
    """Return the proper input constraint"""
    return sum(model.lamda[i]*xref[i,k] for i in model.I) <= x[k]
def output_rule(model,l):
    """Return the proper output constraint"""
    return sum(model.lamda[i]*yref[i,l] for i in model.I) >=model.theta*y[l]
model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
opt = SolverFactory('mosek') # 指定 mosek 作为求解器
solution = opt.solve(model) # 调用求解器求解
theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
obj = value(model.obj) # 提取目标函数 
print("optimum theta: \n {} ".format(theta))
print("optimum lamda: \n {} ".format(lamda))
print("optimal objective: {}".format(obj))
```

## 11.3 步骤2：循环迭代

在前面的步骤中，我们已经完成了决策单元A的技术效率求解。接下来我们将求解决策单元B的技术效率。决策单元B的投入和产出数据在nparray的第二行。事实上我们只要改变被评价决策单元的投入和产出数据（X和Y所在的行），即将步骤1代码中的第8=9行修改为 “x=data[1,0:3]，y=data[1,3:5]”，求解得到的就是第2个决策单元的技术效率。求解所有决策单元的技术效率的基本思路，就是让步骤1的线性规划求解在所有决策单元上进行迭代。我们可以用for循环来实现这一操作。具体代码如下：

```python
from pyomo.environ import *
import numpy as np
data=np.array([[ 5,  6,  7, 23, 67],
               [ 5,  9,  7, 14, 34],
               [ 8,  6, 87, 45, 12],
               [34, 67, 32, 11, 37],
               [ 9,  6, 19, 31, 29]])
thetalt = np.empty(5)                # 定义thetalt 用于存储计算结果
for j in range(5):              # 在data的行上循环
    x=data[j,0:3]
    y=data[j,3:5]
    xref=data[:,0:3]
    yref=data[:,3:5]
    model = ConcreteModel()
    model.I = Set(initialize = range(5)) # 采用列表初始化技术参照决策单元个数的集合
    model.K = Set(initialize = range(3)) # 采用列表初始化投入变量个数的集合
    model.L = Set(initialize = range(2)) # 采用列表初始化产出变量个数的集合
    model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
    model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
    def obj_rule(model):
        return model.theta 
        #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
    model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
    def input_rule(model,k):
        """Return the proper input constraint"""
        return sum(model.lamda[i]*xref[i,k] for i in model.I) <= x[k]
    def output_rule(model,l):
        """Return the proper output constraint"""
        return sum(model.lamda[i]*yref[i,l] for i in model.I) >=model.theta*y[l]
    model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
    model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
    opt = SolverFactory('mosek') # 指定 mosek 作为求解器
    solution = opt.solve(model) # 调用求解器求解
    theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
    lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
    obj = value(model.obj) # 提取目标函数 
    #print("optimum theta: \n {} ".format(theta))
    #print("optimum lamda: \n {} ".format(lamda))
    #print("optimal objective: {}".format(obj))
    thetalt[j]=theta[0]  
te = 1/thetalt
print(te)
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914161931.png)
**图11.2：代码运行结果**

在上述代码中，第8行提前定义了一个变量thetalt（其元素值初始化为缺省值），用以存储下面循环计算得到的结果。我们可以看到这部分代码与步骤1中的核心计算代码完全相同，只是增加了一些for循环的语句。选中以上代码运行之后，te变量就存储了5个决策单元的技术效率值。运行后的结果如图11.2所示。

## 11.4 步骤3：封装成Python函数

我们可以将以上代码写成一个函数，以提高代码应用的灵活性，实现代码的可重复使用。具体代码如下：

```python
from pyomo.environ import *
import numpy as np
data=np.array([[ 5,  6,  7, 23, 67],
               [ 5,  9,  7, 14, 34],
               [ 8,  6, 87, 45, 12],
               [34, 67, 32, 11, 37],
               [ 9,  6, 19, 31, 29]])
thetalt = np.empty(5)                # 定义thetalt 用于存储计算结果
for j in range(5):              # 在data的行上循环
    x=data[j,0:3]
    y=data[j,3:5]
    xref=data[:,0:3]
    yref=data[:,3:5]
    model = ConcreteModel()
    model.I = Set(initialize = range(5)) # 采用列表初始化技术参照决策单元个数的集合
    model.K = Set(initialize = range(3)) # 采用列表初始化投入变量个数的集合
    model.L = Set(initialize = range(2)) # 采用列表初始化产出变量个数的集合
    model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
    model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
    def obj_rule(model):
        return model.theta 
        #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
    model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
    def input_rule(model,k):
        """Return the proper input constraint"""
        return sum(model.lamda[i]*xref[i,k] for i in model.I) <= x[k]
    def output_rule(model,l):
        """Return the proper output constraint"""
        return sum(model.lamda[i]*yref[i,l] for i in model.I) >=model.theta*y[l]
    model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
    model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
    opt = SolverFactory('mosek') # 指定 mosek 作为求解器
    solution = opt.solve(model) # 调用求解器求解
    theta = np.asarray(list(model.theta[:].value)) # 提取决策变量theta
    lamda = np.asarray(list(model.lamda[:].value)) # 提取决策变量lamda
    obj = value(model.obj) # 提取目标函数 
    print("optimum theta: \n {} ".format(theta))
    #print("optimum lamda: \n {} ".format(lamda))
    #print("optimal objective: {}".format(obj))
    thetalt[j]=theta[0]  
te = 1./thetalt
print(te)
```

现在我们对以上代码进行简要的解释。第3行def dea1()定义了一个函数，其名称为dea1。def dea1后面括号中有三项，分别是data、dataref和numk，分别向函数传入评价单元的投入产出数据，参照单元的投入产出数据和投入变量的个数。第9行的np.empty定义一个空的数组。np.empty括号中需要传入一个表示数组长度的数字，使用data.shape[0]确定。shape属性可以获得data的长度（行数）和宽度（列数），[0]表示选择长度。第10行对评价单元的长度进行循环。第11-14行使用k对评价单元和参照单元的投入产出数据进行索引。第16行使用shape方法判断参照单元的个数，第17行用numk定义投入变量个数，第18行用评价单元的总的宽度减去numk来定义产出变量个数。第43行将生成的效率变量返回，即调用dea1的函数后，返回一个储存计数效率值的变量。

在定义完dea1函数后，我们便可以调用dea1函数来实现先前例子中的计数效率计算，具体代码如下:

```python
data=np.array([[ 5,  6,  7, 23, 67],
               [ 5,  9,  7, 14, 34],
               [ 8,  6, 87, 45, 12],
               [34, 67, 32, 11, 37],
               [ 9,  6, 19, 31, 29]])
dataref=data
numk=3
te=dea1(data,dataref,numk)
print(te)
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914171928.png)
**图11.3：代码运行结果**

以上dea1函数可以很好地实现我们第一个DEA模型的计算，但我们可以进一步进行优化。在dea1函数调用之前中，我们使用数组定义被评价决策单元和参照决策单元的投入产出数据。但是数组无法将变量名称包含在内。所以，我们接下来考虑使用pandas的Datarame（数据框）传入投入产出数据，并用数据框的变量名确定投入产出变量。下面我们以stata数据集EX3.dta（见图11.4）作为例子进行说明。该数据集中包含了20个决策单元，在投入和产出变量上有两种投入要素K和L，一种产出Y。

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914174006.png)
**图11.4：数据集EX3.dta**

我们定义另一个函数dea2，用来传入使用数据框定义的投入产出数据。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np
def dea2(dataframe,varname,numk):
    """
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    varname: 字符串列表，表示投入和产出变量的列表，如["K","L","Y"]
    k：表示前k个变量为投入数据
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    data = dataframe.loc[:,varname]
    dataref = dataframe.loc[:,varname]
    for j in range(data.shape[0]):              # 在data的行上循环
        x=data.iloc[j,0:numk]
        y=data.iloc[j,numk:]
        xref=dataref.iloc[:,0:numk]
        yref=dataref.iloc[:,numk:]
        xcol=xref.columns
        ycol=yref.columns
        model = ConcreteModel()
        model.I = Set(initialize = range(dataref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(numk)) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(len(ycol))) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
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
    thetadf["te"] = 1/thetadf["theta"]
    return thetadf
```

以上代码中dea2()函数与dea1()函数主要有以下不同。dea2()函数中的第10-11、13-16行替代了dea1()中的第10-11、13-14行。第10行使用数据框的loc方法取出EX3数据中的投入产出数据。在确定评价单元和参照单元时，使用数据框的iloc方法取出想要的行和列。第31行和34行代码是在约束条件中对xref和yref的调用。在dea1中xref、yref是数组，使用[]可以直接调用，而dea2中xref、yref是数据框，需要使用变量名称和索引值调用。因为k和l分别表示输入变量和输出变量的第k个和第l个表里，所以需要在第17、18行用columns()方法取到变量名的列表，进而使用k和l选取变量名。dea2的第9行和45-47行将dea1中原本使用数组存储优化结果的方法改变为使用字典存储优化结果，并进一步将字典变成数据框。

我们可以调用dea2()函数来实现先前例子中的技术效率的计算。具体代码如下：

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea2(ex3,["K","L","Y"],2)
te
```

计算结果展示在图11.5中。

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914182143.png)
**图11.5：代码运行结果**

接下来，我们对以上函数进行如下扩展：通过数据框中常用的切片方法（如$index \geq 3$、$year==2010$）来指定待评价的决策单元。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np
def dea3(dataframe, varname, numk, evaquery):
    """
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    varname: 字符串列表，表示投入和产出变量的列表，如["K","L","Y"]
    k：表示前k个变量为投入数据
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    indexlt = dataframe.query(evaquery).index
    data = dataframe.loc[indexlt,varname]
    dataref = dataframe.loc[indexlt,varname]
    for j in range(data.shape[0]):              # 在data的行上循环
        x=data.iloc[j,0:numk]
        y=data.iloc[j,numk:]
        xref=dataref.iloc[:,0:numk]
        yref=dataref.iloc[:,numk:]
        xcol=xref.columns
        ycol=yref.columns
        model = ConcreteModel()
        model.I = Set(initialize = range(dataref.shape[0])) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(numk)) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(data.shape[1]-numk)) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
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
    thetadf["te"] = 1/thetadf["theta"]
    return thetadf
```

我们可以调用dea3()函数来实现先前例子中的技术效率的计算。具体代码如下：

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea3(ex3,["K","L","Y"],2,"dmu==[1,2,3,4,5,6,7,8,9,10]")
te
```

计算结果展示在图11.6中。

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220916225103.png)
**图11.6：代码运行结果**

evaquery的其他调用方法：

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea3(ex3,["K","L","Y"],2,"dmu<=10")
te
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220916225201.png)
**图11.7：代码运行结果**

接下来，我们继续对以上函数进行如下扩展：通过数据框中常用的切片方法（如$index \geq 3$、$year==2010$）来指定参照决策单元。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np
def dea4(dataframe, varname, numk, evaquery, refquery):
    """
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    varname: 字符串列表，表示投入和产出变量的列表，如["K","L","Y"]
    k：表示前k个变量为投入数据
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    indexlt = dataframe.query(evaquery).index
    indexltref = dataframe.query(refquery).index
    data = dataframe.loc[indexlt,varname]
    dataref = dataframe.loc[indexltref,varname]
    for j in range(data.shape[0]):              # 在data的行上循环
        x=data.iloc[j,0:numk]
        y=data.iloc[j,numk:]
        xref=dataref.iloc[:,0:numk]
        yref=dataref.iloc[:,numk:]
        xcol=xref.columns
        ycol=yref.columns
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(numk)) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(data.shape[1]-numk)) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
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
    thetadf["te"] = 1/thetadf["theta"]
    return thetadf
```

我们可以调用dea4()函数来实现先前例子中的技术效率的计算。具体代码如下：

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea4(ex3,["K","L","Y"],2,"dmu<=10","dmu<=20")
te
```

执行结果如图11.8所示。

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914202111.png)
**图11.8：代码运行结果**

evaquery和refquery还拥有其他调用方法。具体代码如下：

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea4(ex3,["K","L","Y"],2,"dmu==[10,11,2,3,4,5,6,18]","dmu<=20")
te
```

代码执行结果在11.9中。

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914202458.png)
**图11.9：代码运行结果**

从上图我们可以发现，我们选择的是第10、11、2 、3、4、5 、6 、18个dmu作为待评价单元，但是得到的te的索引却是0-7。我们期望得到和Ex3数据框中第10、11、2 、3、4、5 、6 、18个dmu一样的索引。这样做的好处是使数据框的索引和原始数据的索引相同，从而可以将优化结果与原始数据合并。这需要对dea4进行修改。具体代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np
def dea5(dataframe, varname, numk, evaquery, refquery):
    """
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    varname: 字符串列表，表示投入和产出变量的列表，如["K","L","Y"]
    k：表示前k个变量为投入数据
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
    indexlt = dataframe.query(evaquery).index
    indexltref = dataframe.query(refquery).index
    data = dataframe.loc[indexlt,varname]
    dataref = dataframe.loc[indexltref,varname]
    xcol=varname[0:numk]
    ycol=varname[numk:]
    xref=dataref.loc[:,xcol]
    yref=dataref.loc[:,ycol]
    for j in data.index:              # 在data的索引上循环
        x=data.loc[j,xcol]
        y=data.loc[j,ycol]
        model = ConcreteModel()
        model.I = Set(initialize = dataref.index) # 采用列表初始化技术参照决策单元个数的集合
        model.K = Set(initialize = range(numk)) # 采用列表初始化投入变量个数的集合
        model.L = Set(initialize = range(data.shape[1]-numk)) # 采用列表初始化产出变量个数的集合
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
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
    thetadf["te"] = 1/thetadf["theta"]
    return thetadf
```

对比dea5和dea4，可以发现只有第16-22行的代码发生了变化，变化的原因由dea5的第20行的改变导致。由于我们需要data的索引，所以需要直接在data的索引上进行循环，这样得到的$j$的值就是data的索引。这导致x的切片的方法从iloc变成了loc。iloc被用来索引相对位置，不管索引是什么，使用iloc后都会从0开始数。而loc被用来索引数据框原来的索引值。因此，iloc的第二个元素需要传入变量的位置，而loc的第二个元素需要传入变量名，iloc[j,0:numk]被改成了loc[j,xcol]。dea5的第16-19行对应的变量原来在循环中，现在被放在循环外。且xcol和ycol使用varname的切片来获得。

我们可以通过以下方式进行调用，分析上面相同dmu的技术效率。执行结果展示在图11.10。可以发现，te的索引变成了在Ex3中相同的索引。我们可以通过索引将te和Ex3合并了。

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea5(ex3,["K","L","Y"],2,"dmu==[10,11,2,3,4,5,6,18]","dmu<=20")
te
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220914204918.png)
**图11.10：代码运行结果**

dea5可以很好的实现我们的程序应用目的。在以上函数中我们通过选项varname和numk传递投入变量、产出变量以及投入变量的个数。接下来我们对dea5函数的输入项进行优化，实现以一个R语言中使用的formula表达式格式来传递投入和产出变量。例如我们以表达式"Y=K L"来实现相关变量的传递。在这个表达式中，我们以"="作为分隔符等式，等式左边是产出变量，等式右边是投入变量。对此我们可以使用Python中字符串的split方法来实现投入和产出变量的获取。对于字符串“Y = K L”，我们可以用下面的代码进行分割。

```python
import re
s = "Y = K     L  "
inputvars = s.split('=')[1].strip(' ') 
inputvars = re.compile(' +').sub(' ',inputvars).split(' ')
outputvars = s.split('=')[0]    .strip(' ') 
outputvars = re.compile(' +').sub(' ',outputvars).split(' ')
```

第3行代码中的split("=")将s切分成一个列表，元素分别为”Y“和” K L “。[1]表示选择第2个元素”K     L “。strip(' ') 将”K     L “两端的空格去掉，生成”K     L“。接下来第4行使用正则表达式的sub方法，将1个或多个空格替换为1个空格。.split(' ')表示根据1个空格切分字符串，得到列表['K', 'L']。第5-6行代码同理，得到列表['Y']。这两个列表正好可以用来对数据框的列进行索引。我们可以利用上述代码对dea5进行升级为dea6。下面我们给出dea6的代码：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def dea6(formula, dataframe, evaquery, refquery):
    """
    formula: 产出变量=投入变量，如“ Y   = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
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
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
        model.input = Constraint(model.K, rule=input_rule, doc='input constraint') # 定义与投入变量有关的约束条件
        model.output =Constraint(model.L, rule=output_rule, doc='output constraint') # 定义与产出变量有关的约束条件
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
    thetadf["te"] = 1/thetadf["theta"]
    return thetadf
```

dea6函数的第13-16行使用上述split和正则表达式进行切分，得到xcol和ycol。另外，第26-27行的投入变量和产出变量个数的集合的代码也直接使用xcol和ycol的长度进行定义。其他代码没有发生改变。

dea6函数调用方式和代码运行结果如下：

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea6("Y=K L",ex3,"dmu==[10,11,2,3,4,5,6,18]","dmu<=20")
te
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220922173138.png)
**图11.11：代码运行结果**

上述dea6可以完成对产出方向且规模报酬不变的谢泼德距离函数的技术效率计算，如果想要将技术设定为规模报酬可变，则需要多添加一个约束条件：$\sum_{i=1}^{5} \lambda_i=1$。公式如下：
$$
\begin{aligned} D(\mathbf{X}^*,\mathbf{Y}^*)^{-1}&=\max_{\theta,\lambda_1,\lambda_2,...,\lambda_5}  \quad  \theta  \\s.t.& \sum_{i=1}^{5} \lambda_{i}X_{ki} \leq  X_k^*,  k=1,2,3 \\  & \sum_{i=1}^{5} \lambda_{i}Y_{i}  \geq  \theta Y_k^* ,l=1,2  \\ &\sum_{i=1}^{5} \lambda_i=1 \\& \lambda_{i} \geq 0, n=1, 2,\ldots, 5 ;-\infty \leq \theta\leq +\infty  \end{aligned}
$$
可以将dea6进一步完善为dea7。dea7的代码如下：

```python
from pyomo.environ import *
import pandas as pd ; import numpy as np ; import re
def dea7(formula, dataframe, rts,evaquery, refquery):
    """
    formula: 产出变量=投入变量，如“ Y   = K     L ”
    data: 待评价决策单元的投入产出数据框，按照投入和产出排列
    rts: 传入表示可变规模报酬或不变规模报酬的字符串，如"crs"，"vrs"
    evaquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    refquery:传入数据框.query()方法中的参数，如"dmu==1","dmu==[1,2,3]"
    """
    thetalt = {}          # 定义thetalt 用于存储计算结果
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
        model.theta = Var(within=Reals,bounds=(None, None), doc='efficiency') # 定义决策变量theta
        model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables') # 定义决策变量lamda
        def obj_rule(model):
            return model.theta 
            #return model.theta[0]*1  + sum(model.lamda[i] *0 for i in model.I)
        model.obj = Objective(rule = obj_rule, sense=maximize,doc='objective function') # 定义目标函数
        def input_rule(model,k):
            """Return the proper input constraint"""
            return sum(model.lamda[i]*xref.loc[i,xcol[k]] for i in model.I) <= x[k]
        def output_rule(model,l):
            """Return the proper output constraint"""
            return sum(model.lamda[i]*yref.loc[i,ycol[l]] for i in model.I) >=model.theta*y[l]
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
    thetadf["te"] = 1/thetadf["theta"]
    return thetadf
```

dea7与dea6的不同点在于第41-43行和第46-47行，并新加入了rts参数。我们加入了判断，如果rts参数传入“vrs”，那么在模型中加入vrs_rule。dea7的执行方式和代码运行结果如下：

```python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea7("Y=K L",ex3,"vrs","dmu==[10,11,2,3,4,5,6,18]","dmu<=20")
te
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220922173209.png)
**图11.12：代码运行结果**

接下来，我们以投入导向进行研究。给出不变规模报酬的投入方向谢泼德距离函数的公式为：
$$
\begin{aligned} D(\mathbf{X}^*,\mathbf{Y}^*)&=\min_{\theta,\lambda_1,\lambda_2,...,\lambda_5}  \quad  \theta  \\s.t.& \sum_{i=1}^{5} \lambda_{i}X_{ki} \leq  \theta X_k^*,  k=1,2,3 \\  & \sum_{i=1}^{5} \lambda_{i}Y_{i}  \geq  Y_k^* ,l=1,2  \\& \lambda_{i} \geq 0, n=1, 2,\ldots, 5 ;-\infty \leq \theta\leq +\infty \end{aligned}
$$
给出可变规模报酬的投入方向谢泼德距离函数的公式为：
$$
\begin{aligned} D(\mathbf{X}^*,\mathbf{Y}^*)&=\min_{\theta,\lambda_1,\lambda_2,...,\lambda_5}  \quad  \theta  \\s.t.& \sum_{i=1}^{5} \lambda_{i}X_{ki} \leq  \theta X_k^*,  k=1,2,3 \\  & \sum_{i=1}^{5} \lambda_{i}Y_{i}  \geq  Y_k^* ,l=1,2  \\ &\sum_{i=1}^{5} \lambda_i=1 \\& \lambda_{i} \geq 0, n=1, 2,\ldots, 5 ;-\infty \leq \theta\leq +\infty \end{aligned}
$$

我们可以改进dea7得到dea8，从而加入计算投入方向的谢泼德距离函数。具体代码如下：

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
```

dea8与dea7的区别在于参数添加了orient以定义产出方向或投入方向。第一，第41行的model.obj 的sense加入了判断，从而使最优化方向随着投入或产出的不同方向改变。第二，第42-53行的input_rule和output_rule中加入了判断，从而使不同的方向的约束条件不同。第三，第71行，te的计算方法也加入了判断。第四，第13-20行判断evaquery和refquery是否使用默认值None。当取默认值None时，被评价单元和参考决策单元使用所有的数据。dea8的执行方式和代码运行结果如下：

```Python
import pandas as pd
ex3 = pd.read_stata(r"Ex3.dta")
te=dea8("Y=K L",ex3,"vrs","io","dmu==[10,11,2,3,4,5,6,18]","dmu<=20")
te
```

![](https://advance-markdown.oss-cn-shenzhen.aliyuncs.com/img/20220922173228.png)
**图11.13：代码运行结果**