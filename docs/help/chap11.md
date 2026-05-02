# 第十一章：Pyomo 建模基础与 DEA 从头实现

**主题**：使用 Pyomo 优化建模语言从零构建 DEA 模型，从基础语法到完整函数封装。

**数据文件**：
- `Ex3.dta`（20个 DMU，3个投入(K,L)，1个产出(Y)）

**变量含义**：
- `K`：资本, `L`：劳动, `Y`：产出

---

## Pyomo 介绍与安装

Pyomo 是一个基于 Python 的开源优化建模语言，支持线性规划、非线性规划等多种优化问题。

### 安装

```bash
pip install pyomo
```

还需要安装至少一个求解器：
- **mosek**（推荐，商业）：`pip install mosek`
- **glpk**（免费）：`conda install -c conda-forge glpk`

### 本章需要的导入

```python
from pyomo.environ import *
import pandas as pd
import numpy as np
import re
```

### Pyomo 核心概念

| 概念 | Pyomo 类 | 说明 |
|------|---------|------|
| 集合 | `Set` | 索引集合，如 DMU 集合、投入/产出变量集合 |
| 决策变量 | `Var` | 优化变量，如 theta、lamda |
| 目标函数 | `Objective` | 最大化或最小化目标 |
| 约束条件 | `Constraint` | 线性/非线性约束 |
| 求解器 | `SolverFactory` | 调用求解器求解模型 |

### 建模步骤

```python
# 1. 创建模型
model = ConcreteModel()

# 2. 定义集合
model.I = Set(initialize=range(5))  # DMU 集合
model.K = Set(initialize=range(3))  # 投入变量集合

# 3. 定义决策变量
model.theta = Var(within=Reals, bounds=(None, None), doc='efficiency')
model.lamda = Var(model.I, bounds=(0.0, None), doc='intensity variables')

# 4. 定义目标函数
def obj_rule(model):
    return model.theta
model.obj = Objective(rule=obj_rule, sense=maximize)

# 5. 定义约束条件
def input_rule(model, k):
    return sum(model.lamda[i] * xref[i, k] for i in model.I) <= x[k]
model.input = Constraint(model.K, rule=input_rule)

# 6. 求解
opt = SolverFactory('mosek')
solution = opt.solve(model)

# 7. 提取结果
theta = np.asarray(list(model.theta[:].value))
```

---

### 参数详解（dea8 最终版）

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `formula` | str | `"Y=K L"` | `=` 右边为投入，左边为产出 |
| `dataframe` | DataFrame | — | DMU 数据 |
| `rts` | str | `"crs"` | 不变规模报酬（无额外约束） |
| | str | `"vrs"` | 可变规模报酬（`sum(lamda) == 1`） |
| `orient` | str | `"oo"` | 产出导向：maximize theta，`te = 1/theta` |
| | str | `"io"` | 投入导向：minimize theta，`te = theta` |
| `evaquery` | str | `None` | 评价 DMU 筛选，默认全部 |
| | str | `"dmu<=10"` | 评价 DMU 1-10 |
| | str | `"t==2"` | 评价第 2 期 |
| `refquery` | str | `None` | 参考集 DMU 筛选，默认全部 |

> `evaquery` / `refquery` 使用 `DataFrame.query()` 语法。

### 约束规则

**产出导向 (`"oo"`)**：
- 投入约束：`sum(lamda * xref) <= x`（投入不增）
- 产出约束：`sum(lamda * yref) >= theta * y`（产出最大比例扩张）

**投入导向 (`"io"`)**：
- 投入约束：`sum(lamda * xref) <= theta * x`（投入最大比例缩减）
- 产出约束：`sum(lamda * yref) >= y`（产出不减）

---

## 函数演进路径

| 函数 | 新增特性 | 输入格式 |
|------|---------|---------|
| `dea1` | 基础 DEA 函数 | numpy array |
| `dea2` | 支持 DataFrame + 列名 | DataFrame + varname 列表 |
| `dea3` | 支持 evaquery 筛选 | + evaquery 参数 |
| `dea4` | 支持 refquery 参考集 | + refquery 参数 |
| `dea5` | 支持索引对齐（保留原始 index） | 同 dea4 |
| `dea6` | 支持 formula 公式语法 | formula 字符串 |
| `dea7` | 支持 rts 规模报酬 | + rts 参数 |
| `dea8` | 支持 orient 导向（完整版） | + orient 参数 |

---

## dea1：基础 DEA 函数

```python
def dea1(data, dataref, numk):
    """
    data: 待评价 DMU 的投入产出数据（numpy array）
    dataref: 参考技术 DMU 的投入产出数据
    numk: 前k个变量为投入
    """
    thetalt = np.empty(data.shape[0])
    for j in range(data.shape[0]):
        x = data[j, 0:numk]
        y = data[j, numk:]
        xref = dataref[:, 0:numk]
        yref = dataref[:, numk:]
        # ... Pyomo 建模 + 求解 ...
        thetalt[j] = theta[0]
    te = 1 / thetalt
    return te
```

```python
data = np.array([[5,6,7,23,67], [5,9,7,14,34], ...])
te = dea1(data, data, 3)
```

---

## dea2：支持 DataFrame

```python
def dea2(dataframe, varname, numk):
    """
    dataframe: 数据框
    varname: 变量名列表，如 ["K","L","Y"]
    numk: 前k个变量为投入
    """
```

```python
ex3 = pd.read_stata("Ex3.dta")
te = dea2(ex3, ["K","L","Y"], 2)
```

> 输入从 numpy array 改为 DataFrame，返回含 `theta` 和 `te` 列的 DataFrame。

---

## dea3：支持评价 DMU 筛选

```python
def dea3(dataframe, varname, numk, evaquery):
    """
    evaquery: DataFrame.query() 参数，如 "dmu==[1,2,3]"
    """
```

```python
te = dea3(ex3, ["K","L","Y"], 2, "dmu<=10")
```

---

## dea4：支持参考集筛选

```python
def dea4(dataframe, varname, numk, evaquery, refquery):
    """
    refquery: 参考集 DMU 筛选，如 "dmu<=20"
    """
```

```python
te = dea4(ex3, ["K","L","Y"], 2, "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
```

---

## dea5：索引对齐

```python
def dea5(dataframe, varname, numk, evaquery, refquery):
```

> 使用 DataFrame 原始索引（`data.index`）而非 `range()`，返回结果保留原始 index。

---

## dea6：formula 公式语法

```python
def dea6(formula, dataframe, evaquery, refquery):
    """
    formula: "Y=K L"（产出=投入）
    """
```

```python
te = dea6("Y=K L", ex3, "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
```

> 支持空格不规范的写法（内部用 `re.compile(' +')` 处理）。

---

## dea7：支持规模报酬

```python
def dea7(formula, dataframe, rts, evaquery, refquery):
    """
    rts: "crs" 或 "vrs"
    """
```

```python
te = dea7("Y=K L", ex3, "vrs", "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
```

---

## dea8：支持导向选择（完整版）

```python
def dea8(formula, dataframe, rts, orient, evaquery=None, refquery=None):
    """
    rts: "crs" 或 "vrs"
    orient: "oo"（产出导向）或 "io"（投入导向）
    """
```

```python
# 产出导向 + VRS
te = dea8("Y=K L", ex3, "vrs", "oo", "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")

# 投入导向 + VRS
te = dea8("Y=K L", ex3, "vrs", "io", "dmu==[10,11,2,3,4,5,6,18]", "dmu<=20")
```

---

## 调试方法

### model.pprint()

打印模型所有组件（集合、变量、目标函数、约束）：

```python
model.pprint()
```

### 单步调试

在循环中加入 `break` 可以只构建第一个 DMU 的模型进行检查：

```python
for j in data.index:
    # ... 建模 ...
    break  # 只看第一个 DMU 的模型

model.pprint()  # 查看模型结构
```

### 求解计时

```python
solution = opt.solve(model, report_timing=True)
```

### 求解器选择

| 求解器 | 类型 | 适用场景 | 安装方式 |
|--------|------|---------|---------|
| `"mosek"` | 商业 | 线性/非线性规划（推荐） | `pip install mosek` |
| `"glpk"` | 免费 | 线性规划 | `conda install glpk` |
| `"ipopt"` | 免费 | 非线性规划 | `conda install ipopt` |
| `"scip"` | 免费 | 混合整数规划 | `conda install scip` |
