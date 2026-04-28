# deabook 使用帮助文档

本文档按模块详细说明 `deabook` 包中所有类、方法、参数和常量定义。

---

## 1. 常量模块 (constant)

```python
from deabook.constant import *
```

### 1.1 复合误差项类型 (CET)

| 常量 | 值 | 含义 |
|------|----|------|
| `CET_ADDI` | `"addi"` | 加性复合误差项 |
| `CET_MULT` | `"mult"` | 乘性复合误差项 |

### 1.2 前沿类型 (FUN)

| 常量 | 值 | 含义 |
|------|----|------|
| `FUN_PROD` | `"prod"` | 生产前沿 |
| `FUN_COST` | `"cost"` | 成本前沿 |

### 1.3 规模报酬类型 (RTS)

| 常量 | 值 | 含义 |
|------|----|------|
| `RTS_VRS1` | `"vrs1"` | 可变规模报酬（统一缩减因子 theta，theta bounds=[0,1]） |
| `RTS_VRS2` | `"vrs2"` | 可变规模报酬（差异化缩减因子 mu，mu bounds=[0, +inf]，按参考单元索引） |
| `RTS_CRS` | `"crs"` | 不变规模报酬（无 VRS 约束） |

### 1.4 残差分解方法 (RED)

| 常量 | 值 | 含义 |
|------|----|------|
| `RED_MOM` | `"MOM"` | 矩估计法 (Method of Moments) |
| `RED_QLE` | `"QLE"` | 拟似然估计法 (Quasi-Likelihood Estimation) |
| `RED_KDE` | `"KDE"` | 核密度估计法 (Kernel Density Estimation)，使用 Richardson-Lucy 盲反卷积 |

### 1.5 求解器选项 (OPT)

| 常量 | 值 | 含义 |
|------|----|------|
| `OPT_LOCAL` | `"local"` | 本地求解 |
| `OPT_DEFAULT` | `None` | 使用默认求解器 |

### 1.6 技术类型 (TECH)

| 常量 | 值 | 含义 |
|------|----|------|
| `TOTAL` | `"Global production technology"` | 全局生产技术（所有期参考集合并） |
| `CONTEMPORARY` | `"Contemporary production technolog"` | 当期生产技术（仅当期参考集合） |

### 1.7 指数类型

| 常量 | 值 | 含义 |
|------|----|------|
| `MAL` | `" malquist prodcutivity index or malquist-luenberger prodcutivity index"` | Malmquist 或 Malmquist-Luenberger 指数 |
| `LUE` | `"luenberger prodcutivity index"` | Luenberger 指数 |

---

## 2. DEA 模块 (DEA)

```python
from deabook.DEA import DEA, DEA2, DDF, DDF2, DEADUAL, DDFDUAL
```

### 2.1 DEA（批量包络形式）

经典数据包络分析（DEA）的径向效率模型，一次性求解所有决策单元。

#### 构造函数

```python
DEA(data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 包含投入和产出变量的数据框 |
| `sent` | str | 必填 | 变量声明字符串，格式为 `"inputvars=outputvars"`，例如 `"K L = Y"` |
| `gy` | list | `[1]` | 产出距离向量。`gy[k]=1` 表示第 k 个产出参与径向缩放，`gy[k]=0` 表示不参与。通过 `sum(gy)>=1` 且 `sum(gx)==0` 判断为产出导向；`sum(gx)>=1` 且 `sum(gy)==0` 判断为投入导向 |
| `gx` | list | `[0]` | 投入距离向量。`gx[j]=1` 表示第 j 个投入参与径向缩放，`gx[j]=0` 表示不参与 |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_VRS1`（可变规模报酬）或 `RTS_CRS`（不变规模报酬） |
| `baseindex` | str | `None` | 被评价决策单元的筛选条件，例如 `"Year=[2009,2010]"` |
| `refindex` | str | `None` | 参考决策单元的筛选条件，例如 `"Year=[2010]"` |

#### 方向与优化目标

| 方向 | 条件 | 目标 | rho 含义 |
|------|------|------|----------|
| 投入导向 | `sum(gx)>=1` 且 `sum(gy)==0` | 最小化 rho | 径向缩减因子，te = rho，值域 [0,1] |
| 产出导向 | `sum(gy)>=1` 且 `sum(gx)==0` | 最大化 rho | 径向扩张因子，te = 1/rho，值域 [1, +inf) |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | `email`: 远程优化邮箱地址；`solver`: 求解器名称 | 无 | 求解模型 |
| `display_status()` | 无 | 无 | 打印优化状态 |
| `display_rho()` | 无 | 无 | 打印 rho 值 |
| `display_lamda()` | 无 | 无 | 打印 lamda（强度变量）值 |
| `get_status()` | 无 | int | 返回优化状态码 |
| `get_rho()` | 无 | numpy.ndarray | 返回 rho 值数组 |
| `get_lamda()` | 无 | numpy.ndarray | 返回 lamda 值矩阵（DMU x 参考单元） |

---

### 2.2 DEA2（逐 DMU 包络形式）

经典 DEA 径向效率模型，逐个求解每个决策单元。

#### 构造函数

```python
DEA2(data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None)
```

参数与 DEA 完全相同，但内部为每个 DMU 构建独立的优化模型。额外支持超导向（hyper）模式：`sum(gx)>=1` 且 `sum(gy)>=1`。

#### 方向与优化目标

| 方向 | 条件 | 目标 | rho 约束 | 返回列 |
|------|------|------|----------|--------|
| 投入导向 | `sum(gx)>=1, sum(gy)==0` | 最小化 rho | [0, 1] | `optimization_status`, `rho`, `te`（te=rho） |
| 产出导向 | `sum(gy)>=1, sum(gx)==0` | 最大化 rho | [1, +inf) | `optimization_status`, `rho`, `te`（te=1/rho） |
| 超导向 + CRS | `sum(gx)>=1, sum(gy)>=1, rts=RTS_CRS` | 最小化 rho | [0, 1] | `optimization_status`, `rho`, `te`（te=sqrt(rho)） |
| 超导向 + VRS1 | `sum(gx)>=1, sum(gy)>=1, rts=RTS_VRS1` | 最大化 rho | [0, +inf) | `optimization_status`, `rho`, `tei`（1-rho）, `teo`（1/(1+rho)） |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | 同 DEA | pandas.DataFrame | 返回包含优化状态和效率值的 DataFrame |
| `display_status()` | 无 | 无 | 打印每个 DMU 的优化状态 |
| `display_rho()` | 无 | 无 | 打印每个 DMU 的 rho 值 |
| `display_lamda()` | 无 | 无 | 打印每个 DMU 的 lamda 值 |
| `get_status()` | 无 | dict | 返回优化状态字典 |
| `get_rho()` | 无 | pandas.Series | 返回 rho 值序列 |
| `get_lamda()` | 无 | pandas.DataFrame | 返回 lamda 值数据框 |
| `info(dmu="all")` | `dmu`: DMU 索引标签或列表 | 无 | 打印指定 DMU 的完整模型信息 |

---

### 2.3 DDF（批量方向距离函数）

方向距离函数（DDF）模型，一次性求解所有决策单元。继承自 DEA。

#### 构造函数

```python
DDF(data, sent, gy=[1], gx=[1], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 数据框 |
| `sent` | str | 必填 | 变量声明字符串，格式同 DEA |
| `gy` | list | `[1]` | 产出方向向量。每个元素为对应产出的方向权重 |
| `gx` | list | **`[1]`** | 投入方向向量。**注意：默认值与 DEA 不同，DDF 默认 gx=[1]** |
| `rts` | str | `RTS_VRS1` | 规模报酬 |
| `baseindex` | str | `None` | 被评价单元筛选条件 |
| `refindex` | str | `None` | 参考单元筛选条件 |

**关键区别**：DDF 的目标始终是最大化 rho（方向距离），投入约束为 `x - rho*gx*x`，产出约束为 `y + rho*gy*y`。

#### 方法

继承 DEA 的所有方法（`optimize`, `display_status`, `display_rho`, `display_lamda`, `get_status`, `get_rho`, `get_lamda`）。优化调用 `optimize()` 后通过 `get_rho()` 获取方向距离值。

---

### 2.4 DDF2（逐 DMU 方向距离函数）

方向距离函数模型，逐个求解每个决策单元。继承自 DEA2。

#### 构造函数

```python
DDF2(data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `gx` | list | `[0]` | **注意：DDF2 默认 gx=[0]，与 DDF 的 gx=[1] 不同** |
| 其余参数 | | | 同 DDF |

#### 返回列

| 方向 | 条件 | 返回列 |
|------|------|--------|
| 投入方向 | `sum(gx)>=1, sum(gy)==0` | `optimization_status`, `rho`, `tei`（1-rho） |
| 产出方向 | `sum(gy)>=1, sum(gx)==0` | `optimization_status`, `rho`, `teo`（1/(1+rho)） |
| 超方向 | `sum(gx)>=1, sum(gy)>=1` | `optimization_status`, `rho`, `tei`（1-rho）, `teo`（1/(1+rho)） |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | 同上 | pandas.DataFrame | 返回包含优化状态和效率值的 DataFrame |
| `display_status()` | 无 | 无 | 打印每个 DMU 的优化状态 |
| `display_rho()` | 无 | 无 | 打印每个 DMU 的 rho 值 |
| `display_lamda()` | 无 | 无 | 打印每个 DMU 的 lamda 值 |
| `get_status()` | 无 | dict | 返回优化状态字典 |
| `get_rho()` | 无 | pandas.Series | 返回 rho 值序列 |
| `get_lamda()` | 无 | pandas.DataFrame | 返回 lamda 值数据框 |
| `info(dmu="all")` | `dmu`: DMU 索引 | 无 | 打印指定 DMU 的完整模型信息 |

---

### 2.5 DEADUAL（乘子形式 / 对偶形式）

DEA 的乘子（对偶）形式，用于估计影子价格。继承自 DEA。

#### 构造函数

```python
DEADUAL(data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None)
```

参数与 DEA 完全相同。

#### 额外方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | 同 DEA | 无 | 求解乘子模型 |
| `display_delta()` | 无 | 无 | 打印 delta（投入乘子/影子价格） |
| `display_gamma()` | 无 | 无 | 打印 gamma（产出乘子/影子价格） |
| `display_alpha()` | 无 | 无 | 打印 alpha（VRS 截距项），仅 RTS_VRS1 |
| `get_delta()` | 无 | numpy.ndarray | 返回 delta 值矩阵（DMU x 投入） |
| `get_gamma()` | 无 | numpy.ndarray | 返回 gamma 值矩阵（DMU x 产出） |
| `get_alpha()` | 无 | numpy.ndarray | 返回 alpha 值数组，仅 RTS_VRS1 |
| `get_efficiency()` | 无 | numpy.ndarray | 返回效率值数组。投入导向 CRS: sum(delta*y)；投入导向 VRS1: sum(delta*y)+alpha；产出导向 CRS: sum(gamma*x)；产出导向 VRS1: sum(gamma*x)+alpha |

---

### 2.6 DDFDUAL（方向距离函数乘子形式）

DDF 的乘子（对偶）形式。继承自 DEADUAL。

#### 构造函数

```python
DDFDUAL(data, sent, gy=[1], gx=[1], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `gx` | **`[1]`** | 与 DDF 一致，默认 gx=[1] |
| 其余参数 | | 同 DEADUAL |

DDF 乘子形式的标准化约束为 `sum(delta*gx*x) + sum(gamma*gy*y) == 1`。

#### 方法

继承 DEADUAL 的所有方法（`optimize`, `display_delta`, `display_gamma`, `display_alpha`, `get_delta`, `get_gamma`, `get_alpha`, `get_efficiency`）。

---

## 3. DEAweak 模块 (DEAweak)

```python
from deabook.DEAweak import DEAweak, DEAweak2, DDFweak, DDFweak2, NDDFweak2
```

弱处置性（Weak Disposability）DEA 模型，处理包含非期望产出（如 CO2 排放）的效率评价。

### 3.1 DEAweak（批量弱处置性 DEA）

#### 构造函数

```python
DEAweak(data, sent, gy, gx, gb, rts, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 数据框 |
| `sent` | str | 必填 | 变量声明字符串，格式为 `"inputvars=outputvars:unoutputvars"`，例如 `"K L = Y : CO2"` |
| `gy` | list | **必填，无默认值** | 产出距离向量 |
| `gx` | list | **必填，无默认值** | 投入距离向量 |
| `gb` | list | **必填，无默认值** | 非期望产出方向向量 |
| `rts` | str | **必填，无默认值** | 规模报酬：`RTS_VRS1`, `RTS_VRS2`, `RTS_CRS` |
| `baseindex` | str | `None` | 被评价单元筛选条件 |
| `refindex` | str | `None` | 参考单元筛选条件 |

**注意**：`gy`, `gx`, `gb`, `rts` 均为必填位置参数，无默认值。

#### RTS 行为

| RTS | 额外变量 | VRS 约束 |
|-----|---------|---------|
| `RTS_VRS1` | `theta`（bounds=[0,1]，统一缩减因子） | `sum(lamda) == theta` |
| `RTS_VRS2` | `mu`（bounds=[0,+inf]，按参考单元索引） | `sum(lamda + mu) == 1` |
| `RTS_CRS` | 无 | 无 VRS 约束 |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | 同 DEA | 无 | 求解模型 |
| `display_status()` | 无 | 无 | 打印优化状态 |
| `display_theta()` | 无 | 无 | 打印 theta 值（仅 RTS_VRS1） |
| `display_rho()` | 无 | 无 | 打印 rho 值 |
| `display_lamda()` | 无 | 无 | 打印 lamda 值 |
| `get_status()` | 无 | int | 返回优化状态码 |
| `get_theta()` | 无 | numpy.ndarray | 返回 theta 值数组（仅 RTS_VRS1） |
| `get_rho()` | 无 | numpy.ndarray | 返回 rho 值数组 |
| `get_lamda()` | 无 | numpy.ndarray | 返回 lamda 值矩阵 |

---

### 3.2 DEAweak2（逐 DMU 弱处置性 DEA）

#### 构造函数

```python
DEAweak2(data, sent, gy=[1], gx=[0], gb=[1], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 数据框 |
| `sent` | str | 必填 | 变量声明字符串，格式为 `"inputvars=outputvars:unoutputvars"` |
| `gy` | list | `[1]` | 产出距离向量 |
| `gx` | list | `[0]` | 投入距离向量 |
| `gb` | list | `[1]` | 非期望产出方向向量 |
| `rts` | str | `RTS_VRS1` | 规模报酬 |
| `baseindex` | str | `None` | 被评价单元筛选条件 |
| `refindex` | str | `None` | 参考单元筛选条件 |

#### 方向定义

| 方向 | 条件 |
|------|------|
| 投入导向 (input_oriented) | `sum(gx)>=1, sum(gy)==0, sum(gb)==0` |
| 产出导向 (output_oriented) | `sum(gy)>=1, sum(gx)==0, sum(gb)==0` |
| 非期望产出导向 (unoutput_oriented) | `sum(gb)>=1, sum(gx)==0, sum(gy)==0` |
| 超导向 yx (hyper_orientedyx) | `sum(gx)>=1, sum(gy)>=1, sum(gb)==0` |
| 超导向 yb (hyper_orientedyb) | `sum(gb)>=1, sum(gy)>=1, sum(gx)==0` |

#### RTS + 方向组合的 ValueError 限制

| 方向 | RTS | 是否支持 |
|------|-----|---------|
| 投入导向 | `RTS_VRS1` | **不支持**，抛出 ValueError |
| 超导向 yx | `RTS_VRS1` | **不支持**，抛出 ValueError |
| 产出导向 | 所有 RTS | 支持 |
| 非期望产出导向 | 所有 RTS | 支持 |
| 超导向 yb | 所有 RTS | 支持 |

#### optimize() 返回列

| 方向 | RTS | 返回列 |
|------|-----|--------|
| 投入导向 / 非期望产出导向 | CRS / VRS2 | `optimization_status`, `rho`, `te`（te=rho） |
| 产出导向 | 所有 RTS | `optimization_status`, `rho`, `te`（te=1/rho） |
| 超导向 yx | CRS | `optimization_status`, `rho`, `te`（te=sqrt(rho)） |
| 超导向 yx | VRS2 | `optimization_status`, `rho`, `tei`（1-rho）, `teo`（1/(1+rho)） |
| 超导向 yb | CRS | `optimization_status`, `rho`, `te`（te=sqrt(rho)） |
| 超导向 yb | VRS1 | `optimization_status`, `rho`, `teuo`（1-rho）, `teo`（1/(1+rho)） |
| 超导向 yb | VRS2 | `optimization_status`, `rho`, `teuo`（1-rho）, `teo`（1/(1+rho)） |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | 同 DEA | pandas.DataFrame | 返回效率结果 |
| `display_status()` | 无 | 无 | 打印每个 DMU 的优化状态 |
| `display_rho()` | 无 | 无 | 打印每个 DMU 的 rho 值 |
| `display_theta()` | 无 | 无 | 打印每个 DMU 的 theta 值（仅 RTS_VRS1） |
| `display_lamda()` | 无 | 无 | 打印每个 DMU 的 lamda 值 |
| `display_mu()` | 无 | 无 | 打印每个 DMU 的 mu 值（仅 RTS_VRS2） |
| `get_status()` | 无 | dict | 返回优化状态字典 |
| `get_rho()` | 无 | pandas.Series | 返回 rho 值序列 |
| `get_theta()` | 无 | pandas.Series | 返回 theta 值序列（仅 RTS_VRS1） |
| `get_lamda()` | 无 | pandas.DataFrame | 返回 lamda 值数据框 |
| `get_mu()` | 无 | pandas.DataFrame | 返回 mu 值数据框（仅 RTS_VRS2） |
| `info(dmu="all")` | `dmu`: DMU 索引 | 无 | 打印指定 DMU 的完整模型信息 |

---

### 3.3 DDFweak（批量弱处置性方向距离函数）

#### 构造函数

```python
DDFweak(data, sent, gy=[1], gx=[1], gb=[1], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `gx` | **`[1]`** | 与 DEAweak 不同，DDFweak 默认 gx=[1] |
| `gb` | **`[1]`** | 默认 gb=[1] |
| 其余参数 | | 同 DEAweak |

DDFweak 始终最大化 rho。投入约束为 `x - rho*gx*x`，产出约束为 `y + rho*gy*y`，非期望产出约束为 `b - rho*gb*b`。

#### 方法

继承 DEAweak 的所有方法。

---

### 3.4 DDFweak2（逐 DMU 弱处置性方向距离函数）

#### 构造函数

```python
DDFweak2(data, sent, gy=[1], gx=[0], gb=[1], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `gx` | `[0]` | 与 DDFweak 不同，DDFweak2 默认 gx=[0] |
| 其余参数 | | 同 DDFweak |

#### 方向定义（7 种）

| 方向 | 条件 |
|------|------|
| 投入导向 (input_oriented) | `sum(gx)>=1, sum(gy)==0, sum(gb)==0` |
| 产出导向 (output_oriented) | `sum(gy)>=1, sum(gx)==0, sum(gb)==0` |
| 非期望产出导向 (unoutput_oriented) | `sum(gb)>=1, sum(gx)==0, sum(gy)==0` |
| 超导向 yx (hyper_orientedyx) | `sum(gx)>=1, sum(gy)>=1, sum(gb)==0` |
| 超导向 yb (hyper_orientedyb) | `sum(gb)>=1, sum(gy)>=1, sum(gx)==0` |
| 超导向 xb (hyper_orientedxb) | `sum(gb)>=1, sum(gx)>=1, sum(gy)==0` |
| 超导向 yxb (hyper_orientedyxb) | `sum(gb)>=1, sum(gx)>=1, sum(gy)>=1` |

#### RTS + 方向组合的 ValueError 限制

| 方向 | RTS | 是否支持 |
|------|-----|---------|
| 投入导向 | `RTS_VRS1` | **不支持** |
| 超导向 yx | `RTS_VRS1` | **不支持** |
| 超导向 xb | `RTS_VRS1` | **不支持** |
| 超导向 yxb | `RTS_VRS1` | **不支持** |
| 其余组合 | | 支持 |

#### optimize() 返回列

| 方向 | RTS | 返回列 |
|------|-----|--------|
| 投入导向 | CRS / VRS2 | `optimization_status`, `rho`, `tei`（1-rho） |
| 产出导向 | 所有 RTS | `optimization_status`, `rho`, `teo`（1/(1+rho)） |
| 非期望产出导向 | 所有 RTS | `optimization_status`, `rho`, `tei`（1-rho） |
| 超导向 yx | CRS / VRS2 | `optimization_status`, `rho`, `tei`（1-rho）, `teo`（1/(1+rho)） |
| 超导向 yb | 所有 RTS | `optimization_status`, `rho`, `tei`（1-rho）, `teo`（1/(1+rho)） |
| 超导向 xb | CRS / VRS2 | `optimization_status`, `rho`, `tei`（1-rho）, `teo`（1/(1+rho)） |
| 超导向 yxb | CRS / VRS2 | `optimization_status`, `rho`, `tei`（1-rho）, `teo`（1/(1+rho)） |

#### 方法

继承 DEAweak2 的所有方法（`optimize`, `display_status`, `display_rho`, `display_theta`, `display_lamda`, `display_mu`, `get_status`, `get_rho`, `get_theta`, `get_lamda`, `get_mu`, `info`）。

---

### 3.5 NDDFweak2（非径向方向距离函数）

非径向 DDF（Non-radial DDF），为每个变量独立计算效率比率，而非统一 rho。

**注意**：NDDFweak2 类在源码中被 MQDEAweak 模块引用（`from .DEAweak import ... NDDFweak2`），具体参数和返回值与 DDFweak2 类似但每个变量使用独立的松弛变量代替统一 rho。

---

## 4. MQDEA 模块 (MQDEA)

```python
from deabook.MQDEA import MQDEA, MQDDF
```

Malmquist 生产率指数模块，用于动态效率变化分解。

### 4.1 MQDEA（Malmquist 指数）

#### 构造函数

```python
MQDEA(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 面板数据框，必须包含 `id` 和 `year` 指定的列 |
| `id` | str | 必填 | 决策单元标识列名，例如 `"Firm"` |
| `year` | str | 必填 | 时间列名，例如 `"Year"` |
| `sent` | str | `"inputvar=outputvar"` | 变量声明字符串，格式同 DEA |
| `gy` | list | `[1]` | 产出距离向量 |
| `gx` | list | `[0]` | 投入距离向量 |
| `rts` | str | `RTS_VRS1` | 规模报酬 |
| `tech` | str | `TOTAL` | 技术类型：`TOTAL`（全局生产技术）或 `CONTEMPORARY`（当期生产技术） |
| `email` | str | `OPT_LOCAL` | 远程优化邮箱 |
| `solver` | str | `OPT_DEFAULT` | 求解器名称 |

**重要**：计算在 `__init__` 中自动完成，`optimize()` 仅返回已计算的结果。

#### 内部计算逻辑

- **TOTAL 模式**：对每个年份，评价期 DMU 以所有年份 DMU 为参考集（D11）。计算 `mqpi = D11(t) / D11(t-1)`。
- **CONTEMPORARY 模式**：计算三个距离函数：
  - D11：t 期 DMU 以 t 期参考集（同期）
  - D12：t 期 DMU 以 t-1 期参考集
  - D21：t-1 期 DMU 以 t 期参考集
  - 计算 `MQ = sqrt(D12/D11_prev * D11/D21_prev)`
  - 计算 `MEFFCH = D11 / D11_prev`（技术效率变化）
  - 计算 `MTECHCH = sqrt((D12/D11) * (D11_prev/D21_prev))`（技术变化）

#### optimize() 返回 DataFrame 列

**TOTAL 模式**：

| 方向/RTS | 返回列 |
|----------|--------|
| 投入导向 / 产出导向 / 超导向+CRS | `D11`, `mqpi` |
| 超导向+VRS1 | `D11_tei`, `D11_teo`, `mqpi_tei`, `mqpi_teo` |

**CONTEMPORARY 模式**：

| 方向/RTS | 返回列 |
|----------|--------|
| 投入导向 / 产出导向 / 超导向+CRS | `D11`, `D12`, `D21`, `MQ`, `MEFFCH`, `MTECHCH` |
| 超导向+VRS1 | `D11_tei`, `D11_teo`, `D12_tei`, `D12_teo`, `D21_tei`, `D21_teo`, `MQ_tei`, `MQ_teo`, `MEFFCH_tei`, `MEFFCH_teo`, `MTECHCH_tei`, `MTECHCH_teo` |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize()` | 无 | pandas.DataFrame | 返回计算结果（计算已在初始化时完成） |

---

### 4.2 MQDDF（Malmquist-Luenberger 指数）

使用方向距离函数（DDF2）计算 Malmquist-Luenberger 指数。

#### 构造函数

```python
MQDDF(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)
```

参数与 MQDEA 完全相同。内部使用 DDF2 而非 DEA2 进行效率计算。

#### 与 MQDEA 的关键差异

- 内部使用 `DDF2` 而非 `DEA2`
- 从 DDF2 结果中提取 `tei`/`teo` 列（而非 `te` 列）
- 距离函数值的含义为方向距离而非径向效率

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize()` | 无 | pandas.DataFrame | 返回计算结果 |

---

## 5. MQDEAweak 模块 (MQDEAweak)

```python
from deabook.MQDEAweak import MQDEAweak, MQDDFweak, MQNDDFweak
```

弱处置性 Malmquist 生产率指数模块，处理包含非期望产出的动态效率分解。

### 5.1 MQDEAweak

#### 构造函数

```python
MQDEAweak(data, id, year, sent="inputvar=outputvar:unoutputvar", gy=[1], gx=[0], gb=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 面板数据框 |
| `id` | str | 必填 | 决策单元标识列名 |
| `year` | str | 必填 | 时间列名 |
| `sent` | str | `"inputvar=outputvar:unoutputvar"` | 变量声明字符串，格式为 `"K L = Y : CO2"` |
| `gy` | list | `[1]` | 产出距离向量 |
| `gx` | list | `[0]` | 投入距离向量 |
| `gb` | list | `[0]` | 非期望产出方向向量 |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_VRS1`, `RTS_VRS2`, `RTS_CRS` |
| `tech` | str | `TOTAL` | 技术类型：`TOTAL` 或 `CONTEMPORARY` |
| `email` | str | `OPT_LOCAL` | 远程优化邮箱 |
| `solver` | str | `OPT_DEFAULT` | 求解器名称 |

#### 方向定义（5 种）

| 方向 | 条件 |
|------|------|
| 投入导向 (input_oriented) | `sum(gx)>=1, sum(gy)==0, sum(gb)==0` |
| 产出导向 (output_oriented) | `sum(gy)>=1, sum(gx)==0, sum(gb)==0` |
| 非期望产出导向 (undesirable_oriented) | `sum(gb)>=1, sum(gx)==0, sum(gy)==0` |
| 超导向 yx (hyper_orientedyx) | `sum(gx)>=1, sum(gy)>=1, sum(gb)==0` |
| 超导向 yb (hyper_orientedyb) | `sum(gb)>=1, sum(gy)>=1, sum(gx)==0` |

#### optimize() 返回 DataFrame 列

**TOTAL 模式**：

| 方向/RTS | 返回列 |
|----------|--------|
| 投入/产出/非期望产出导向 / 超导向+CRS | `D11`, `mqpi` |
| 超导向 yb + VRS1/VRS2 | `D11_teuo`, `D11_teo`, `mqpi_teuo`, `mqpi_teo` |
| 超导向 yx + VRS2 | `D11_tei`, `D11_teo`, `mqpi_tei`, `mqpi_teo` |

**CONTEMPORARY 模式**：

| 方向/RTS | 返回列 |
|----------|--------|
| 标准方向 / 超导向+CRS | `D11`, `D12`, `D21`, `MQ`, `MEFFCH`, `MTECHCH` |
| 超导向 yx + VRS2 | `D11_tei`, `D11_teo`, `D12_tei`, `D12_teo`, `D21_tei`, `D21_teo`, `MQ_tei`, `MQ_teo`, `MEFFCH_tei`, `MEFFCH_teo`, `MTECHCH_tei`, `MTECHCH_teo` |
| 超导向 yb + VRS1/VRS2 | 同上，使用 `teuo`/`teo` 后缀 |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize()` | 无 | pandas.DataFrame | 返回计算结果（计算已在初始化时完成） |

---

### 5.2 MQDDFweak

弱处置性 Malmquist-Luenberger 指数，内部使用 DDFweak2。

### 5.3 MQNDDFweak

非径向 DDF 的 Malmquist 指数变体，内部使用 NDDFweak2。

---

## 6. MB 模块 (MB)

```python
from deabook.MB import MB, MBx, MBx2, MBxy, MBxy2, MB2
```

物质平衡（Material Balance）分解模块，将运营效率（OE）分解为技术效率（TE）、潜在投入分配效率（PIAE）、非潜在投入分配效率（NPIAE）和非潜在产出分配效率（NPOAE）等组分。

### 6.1 MB（基础物质平衡分解）

#### 构造函数

```python
MB(data, sent="inputvar_np + inputvar_p = outputvar_np + outputvar_p : unoutputvar",
   sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], level=5,
   rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 数据框 |
| `sent` | str | `"inputvar_np + inputvar_p =outputvar_np + outputvar_p:unoutputvar"` | 变量声明字符串。使用 `+` 区分污染和非污染变量。格式：`"非污染投入 + 污染投入 = 非污染产出 + 污染产出 : 非期望产出"`。例如 `"K L + E = Y : CO2"` |
| `sx` | list | `[[1,1,1],[1,1,1]]` | 投入包含污染物质系数矩阵。外层维度对应非期望产出个数，内层维度对应投入个数 |
| `sy` | list | `[[1],[1]]` | 期望产出包含污染物质系数矩阵。外层维度对应非期望产出个数，内层维度对应产出个数 |
| `level` | int | `5` | 分解层级：1=仅变量化 b；2=再变量化 x_p；3=再变量化 x_np；4=再变量化 y_p（无则变量化 y_np）；5=再变量化 y_np |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_VRS1` 或 `RTS_CRS` |
| `baseindex` | str | `None` | 被评价单元筛选条件 |
| `refindex` | str | `None` | 参考单元筛选条件 |

#### sent 格式说明

sent 字符串支持四种组合（取决于数据中是否存在非污染产出/污染产出）：

| 格式 | 示例 | 含义 |
|------|------|------|
| 有非污染产出 + 污染产出 | `"K L + E = Y + Y2 : CO2"` | 使用 MB1111 内部类 |
| 仅有非污染产出 | `"K L + E = Y : CO2"` | 使用 MB1110 内部类 |
| 仅有污染产出 | `"K L + E = + Y2 : CO2"` | 使用 MB1101 内部类 |
| 无产出 | `"K L + E = : CO2"` | 使用 MB1100 内部类 |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(solver=OPT_DEFAULT, dmu="1")` | `solver`: 求解器；`dmu`: 显示信息的 DMU 编号 | `(data3, info)` | 返回元组：`data3` 为结果 DataFrame（包含 obj, var Undesirable, Input slack, Output slack, Undesirable Output slack 等列）；`info` 为模型详细信息 |
| `info(dmu="all")` | `dmu`: DMU 编号 | 无 | 打印模型信息 |

---

### 6.2 MBx（投入端物质平衡分解 + 面板）

逐 DMU 求解，包含年份维度，仅对投入端进行分解。

#### 构造函数

```python
MBx(data, year, sent="inputvar=outputvar", sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 面板数据框 |
| `year` | str 或 array-like | 必填 | 时间列名或时间序列 |
| 其余参数 | | | 同 MB（但不含 `level` 参数） |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(solver=OPT_DEFAULT)` | `solver`: 求解器 | pandas.DataFrame | 返回包含 obj, var Undesirable, 各变量 slack 的 DataFrame |
| `info(dmu="all")` | `dmu`: DMU 编号 | 无 | 打印模型信息 |

---

### 6.3 MBx2（投入端物质平衡分解 + 面板变体）

与 MBx 功能相同，为变体实现。参数和方法完全一致。

```python
MBx2(data, year, sent="inputvar=outputvar", sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None, refindex=None)
```

---

### 6.4 MBxy（投入+产出端物质平衡分解 + 面板）

同时包含投入和产出的物质平衡分解，增加了 `objy` 变量（目标产出）。

#### 构造函数

```python
MBxy(data, year, sent="inputvar=outputvar", sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None, refindex=None)
```

参数与 MBx 相同。

**与 MBx 的关键差异**：
- 增加了 `objy` 变量（目标产出）
- 产出约束：`sum(lamda*yref) - thetay == objy`（而非等于当前 DMU 的产出）
- 物质平衡约束中包含产出松弛项

#### 方法

与 MBx 相同（`optimize`, `info`）。

---

### 6.5 MBxy2（投入+产出端物质平衡分解变体）

与 MBxy 类似，但：
- 目标函数改为**最大化** `thetab`（非期望产出松弛）
- 非期望产出约束：`sum(lamda*bref) + thetab == b_current`（不使用 theta 变量）
- 无 `var Undesirable` 列

```python
MBxy2(data, year, sent="inputvar=outputvar", sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None, refindex=None)
```

---

### 6.6 MB2（简化物质平衡分解）

无 `objx`/`objy` 目标变量，直接使用当前 DMU 的投入/产出数据。

#### 构造函数

```python
MB2(data, year, sent="inputvar=outputvar", sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None, refindex=None)
```

**关键差异**：
- 无 `objx`/`objy` 变量
- 投入约束：`sum(lamda*xref) + thetax == x_current`
- 产出约束：`sum(lamda*yref) - thetay == y_current`
- 非期望产出约束：`sum(lamda*bref) + thetab == b_current`
- 目标函数为**最大化** `sum(thetab)`
- 物质平衡约束：`sum(sx*thetax) + sum(sy*thetay) == thetab`

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(solver=OPT_DEFAULT)` | `solver`: 求解器 | pandas.DataFrame | 返回结果 |
| `info(dmu="all")` | `dmu`: DMU 编号 | 无 | 打印模型信息 |

---

## 7. CNLSSDFDDFweak 模块 (CNLSSDFDDFweak)

```python
from deabook.CNLSSDFDDFweak import CNLSSDweak, CNLSDDFweak
```

凸非参数最小二乘（CNLS）随机前沿估计模块，支持 Shephard 距离函数和方向距离函数。

### 7.1 CNLSSDweak（Shephard 距离函数 CNLS）

#### 构造函数

```python
CNLSSDweak(data, sent, z=None, gy=[1], gx=[0], gb=[0], cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 数据框 |
| `sent` | str | 必填 | 变量声明字符串，格式为 `"inputvars=outputvars:unoutputvars"` |
| `z` | list 或 None | `None` | 环境变量（Contextual Variables），用于控制异质性 |
| `gy` | list | `[1]` | 产出距离向量。`sum(gy)>=1` 且 `sum(gx)==0` 且 `sum(gb)==0` 为产出方向 |
| `gx` | list | `[0]` | 投入距离向量。`sum(gx)>=1` 且 `sum(gy)==0` 且 `sum(gb)==0` 为投入方向 |
| `gb` | list | `[0]` | 非期望产出方向向量。`sum(gb)>=1` 且 `sum(gx)==0` 且 `sum(gy)==0` 为非期望产出方向 |
| `cet` | str | `CET_MULT` | **必须为 `CET_MULT`**。若传入 `CET_ADDI` 将抛出 ValueError |
| `fun` | str | `FUN_PROD` | 前沿类型：`FUN_PROD`（生产前沿）或 `FUN_COST`（成本前沿） |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_VRS1`, `RTS_VRS2`, `RTS_CRS` |

#### 方向与残差符号

| 方向 | 条件 | 回归方程残差符号 | `get_residual()` 符号处理 |
|------|------|----------------|--------------------------|
| 投入方向 | `sum(gx)>=1` | `+epsilon` | `+1*v` |
| 非期望产出方向 | `sum(gb)>=1` | `+epsilon` | `+1*v` |
| 产出方向 | `sum(gy)>=1` | `-epsilon` | `-1*v` |

#### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | `email`: 邮箱；`solver`: 求解器 | 无 | 求解模型 |
| `display_status()` | 无 | 无 | 打印优化状态 |
| `display_alpha()` | 无 | 无 | 打印 alpha 值（仅 VRS） |
| `display_delta()` | 无 | 无 | 打印 delta（投入系数） |
| `display_gamma()` | 无 | 无 | 打印 gamma（产出系数） |
| `display_kappa()` | 无 | 无 | 打印 kappa（非期望产出系数） |
| `display_omega()` | 无 | 无 | 打印 omega（环境变量系数），仅当 `z` 非 None |
| `display_residual()` | 无 | 无 | 打印残差 |
| `get_status()` | 无 | int | 返回优化状态码 |
| `get_alpha()` | 无 | numpy.ndarray | 返回 alpha 值数组（仅 VRS） |
| `get_delta()` | 无 | numpy.ndarray | 返回 delta 值矩阵（DMU x 投入） |
| `get_gamma()` | 无 | numpy.ndarray | 返回 gamma 值矩阵（DMU x 产出） |
| `get_kappa()` | 无 | numpy.ndarray | 返回 kappa 值矩阵（DMU x 非期望产出） |
| `get_omega()` | 无 | numpy.ndarray | 返回 omega 值数组，仅当 `z` 非 None |
| `get_residual()` | 无 | numpy.ndarray | 返回残差数组（已按方向调整符号） |
| `get_adjusted_residual()` | 无 | numpy.ndarray | 返回平移后的残差（`residual - max(residual)`） |
| `get_adjusted_alpha()` | 无 | numpy.ndarray | 返回平移后的 alpha（`alpha + max(residual)`） |

---

### 7.2 CNLSDDFweak（方向距离函数 CNLS）

#### 构造函数

```python
CNLSDDFweak(data, sent, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS1,
            deltap=None, gammap=None, kappap=None, epsilonp=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | pandas.DataFrame | 必填 | 数据框 |
| `sent` | str | 必填 | 变量声明字符串 |
| `z` | list 或 None | `None` | 环境变量 |
| `gy` | list | `[1]` | 产出方向向量 |
| `gx` | list | **`[1]`** | 投入方向向量（默认 [1]，与 CNLSSDweak 的 [0] 不同） |
| `gb` | list | **`[1]`** | 非期望产出方向向量（默认 [1]） |
| `fun` | str | `FUN_PROD` | 前沿类型 |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_VRS1`, `RTS_VRS2`, `RTS_CRS` |
| `deltap` | list 或 None | `None` | delta 变量边界，格式为 `[lower, upper]`。若为 None 则 bounds=[0, +inf) |
| `gammap` | list 或 None | `None` | gamma 变量边界，格式为 `[lower, upper]`。若为 None 则 bounds=[0, +inf) |
| `kappap` | list 或 None | `None` | kappa 变量边界，格式为 `[lower, upper]`。若为 None 则 bounds=[0, +inf) |
| `epsilonp` | list 或 None | `None` | epsilon 变量边界（当前未生效，bounds 固定为 None） |

**与 CNLSSDweak 的关键差异**：
- 无 `cet` 参数（始终为加性误差项 CET_ADDI）
- 默认 `gx=[1]`, `gb=[1]`（CNLSSDweak 默认 `gx=[0]`, `gb=[0]`）
- 支持 `deltap`/`gammap`/`kappap` 参数自定义变量边界

#### 方向与残差符号

| 方向 | 条件 | `get_residual()` 符号处理 |
|------|------|--------------------------|
| 投入方向 | `sum(gx)>=1` | `-1*v` |
| 非期望产出方向 | `sum(gb)>=1` | `-1*v` |
| 产出方向 | `sum(gy)>=1` | `+1*v` |

**注意**：残差符号处理与 CNLSSDweak **相反**。

#### 方法

继承 CNLSSDweak 的所有方法，但 `get_residual()` 的符号处理逻辑不同（见上表）。`get_adjusted_residual()` 和 `get_adjusted_alpha()` 同样可用。

---

## 8. StoNED 模块 (StoNED)

```python
from deabook.StoNED import StoNED
```

随机非参数包络数据（StoNED）模块，用于对 CNLS 模型的残差进行分解，估计技术效率和噪声。

### 8.1 StoNED

#### 构造函数

```python
StoNED(model)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `model` | CNLSSDweak / CNLSDDFweak / CNLSSD / CNLSDDF / CNLSSDweakmeta / CNLSDDFweakmeta | 必填 | 已优化的 CNLS 模型对象 |

#### 构造函数行为

根据 `model.__class__.__name__` 初始化不同的属性：

| 模型类名 | 额外初始化属性 |
|----------|---------------|
| `CNLSSDweak`, `CNLSDDFweak` | `basexy`, `gb`（来自 model） |
| `CNLSSD`, `CNLSDDF` | `basexy`（来自 model） |
| `CNLSSDweakmeta`, `CNLSDDFweakmeta` | `basexy`, `gb`, `gddf_er`, `gddf`, `gresidual`, `basexy_old` |

通用属性：`y`, `gy`, `gx`, `epsilonhat`（残差）。

#### 方法

##### get_mean_of_inefficiency

```python
get_mean_of_inefficiency(method=RED_MOM)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `method` | str | `RED_MOM` | 残差分解方法：`RED_MOM`（矩估计）、`RED_QLE`（拟似然估计）、`RED_KDE`（核密度估计） |

**返回值**：`self.mu`（float）-- 无条件期望非效率值。

该方法同时设置以下属性：
- `self.mu`：无条件期望非效率
- `self.sigma_u`：非效率项标准差
- `self.sigma_v`：噪声项标准差
- `self.residual_minus`：调整后残差（`residual - mu` 或 `residual + mu`，取决于 fun）

##### get_technical_inefficiency

```python
get_technical_inefficiency(method=RED_MOM)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `method` | str | `RED_MOM` | 仅支持 `RED_QLE` 和 `RED_KDE` |

**返回值**：numpy.ndarray -- 每个 DMU 的技术非效率估计值。

| 方法 | 计算方式 |
|------|---------|
| `RED_QLE` | JLMS 估计：`sigmas * norm.pdf(mus/sigmas) / norm.cdf(mus/sigmas) + mus`。对 CET_ADDI 和 CET_MULT 均使用相同公式 |
| `RED_KDE` | Richardson-Lucy 盲反卷积算法，返回 `u_final`（非负非效率估计） |

**注意**：若传入 `RED_MOM`，将抛出 ValueError（MOM 方法无法估计个体非效率）。

##### get_technical_efficiency

```python
get_technical_efficiency(method=RED_MOM)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `method` | str | `RED_MOM` | 仅支持 `RED_QLE` 和 `RED_KDE` |

**返回值**：numpy.ndarray -- 每个 DMU 的技术效率值，计算方式为 `exp(-技术非效率)`。

对 CET_ADDI 和 CET_MULT 均返回 `exp(-get_technical_inefficiency(method))`。

##### get_technical_efficiency_ratio

```python
get_technical_efficiency_ratio(method=RED_MOM)
```

仅适用于方向距离函数模型（CNLSDDFweak, CNLSDDF, CNLSSDweakmeta, CNLSDDFweakmeta）。

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `method` | str | `RED_MOM` | 仅支持 `RED_QLE` 和 `RED_KDE` |

**返回值**：

对于 CNLSDDFweak / CNLSDDF：

| 方向 | TE 计算公式 |
|------|-----------|
| `sum(gy)>=1` | `TE = basexy / (basexy + ddfhat)` |
| `sum(gx)>=1` | `TE = (-basexy + epsilonhat) / (-basexy + epsilonhat + ddfhat)` |
| `sum(gb)>=1` | `TE = (-basexy + epsilonhat) / (-basexy + epsilonhat + ddfhat)` |

返回 numpy.ndarray。

对于 CNLSSDweakmeta / CNLSDDFweakmeta：

返回 pandas.DataFrame，包含以下列：

| 列名 | 含义 |
|------|------|
| `TGR` | 技术差距比（Technology Gap Ratio） |
| `GTE` | 组技术效率（Group Technical Efficiency） |
| `MTE` | 元前沿技术效率（Metafrontier Technical Efficiency）= TGR * GTE |

##### richardson_lucy_blind_corrected

```python
richardson_lucy_blind_corrected(method)
```

内部方法，由 `get_technical_inefficiency(method=RED_KDE)` 调用。使用修正的 Richardson-Lucy 盲反卷积算法从残差中分离非效率 u 和噪声 v。

| 参数 | 默认值 | 含义 |
|------|--------|------|
| `kernel_size` | 5 | 噪声核大小 |
| `max_m` | 500 | 最大盲迭代次数 |
| `max_j` | 500 | 每次盲迭代中的 RL 迭代次数 |
| `tol` | 1e-16 | 收敛容差 |

---

## 9. 使用示例

### 9.1 基本 DEA 分析

```python
import pandas as pd
from deabook.DEA import DEA, DEA2

# 准备数据
data = pd.read_csv("data.csv")

# 投入导向 DEA（批量）
model = DEA(data, sent="K L = Y", gy=[0], gx=[1, 1], rts="vrs1")
model.optimize(solver="mosek")
rho = model.get_rho()

# 投入导向 DEA（逐 DMU）
model2 = DEA2(data, sent="K L = Y", gy=[0], gx=[1, 1], rts="vrs1")
results = model2.optimize(solver="mosek")
print(results[['optimization_status', 'rho', 'te']])
```

### 9.2 方向距离函数

```python
from deabook.DEA import DDF, DDFDUAL

# DDF 包络形式
model = DDF(data, sent="K L = Y", gy=[1], gx=[1], rts="vrs1")
model.optimize(solver="mosek")
ddf_value = model.get_rho()

# DDF 乘子形式（影子价格）
dual_model = DDFDUAL(data, sent="K L = Y", gy=[1], gx=[1], rts="vrs1")
dual_model.optimize(solver="mosek")
delta = dual_model.get_delta()  # 投入影子价格
gamma = dual_model.get_gamma()  # 产出影子价格
```

### 9.3 弱处置性 DEA（非期望产出）

```python
from deabook.DEAweak import DEAweak2, DDFweak2

# 投入导向弱处置性 DEA（逐 DMU）
model = DEAweak2(data, sent="K L = Y : CO2",
                 gy=[0], gx=[1, 1], gb=[0],
                 rts="crs")
results = model.optimize(solver="mosek")
print(results[['optimization_status', 'rho', 'te']])

# 方向距离函数弱处置性
model2 = DDFweak2(data, sent="K L = Y : CO2",
                  gy=[1], gx=[1], gb=[1],
                  rts="vrs2")
results2 = model2.optimize(solver="mosek")
print(results2[['optimization_status', 'rho', 'tei', 'teo']])
```

### 9.4 Malmquist 生产率指数

```python
from deabook.MQDEA import MQDEA

# 全局技术 Malmquist 指数
mq = MQDEA(data, id="Firm", year="Year",
           sent="K L = Y", gy=[0], gx=[1, 1],
           rts="vrs1", tech="Global production technology")
result = mq.optimize()
print(result[['D11', 'mqpi']])

# 当期技术 Malmquist 指数（完整分解）
mq2 = MQDEA(data, id="Firm", year="Year",
            sent="K L = Y", gy=[0], gx=[1, 1],
            rts="vrs1", tech="Contemporary production technolog")
result2 = mq2.optimize()
print(result2[['D11', 'D12', 'D21', 'MQ', 'MEFFCH', 'MTECHCH']])
```

### 9.5 弱处置性 Malmquist 指数

```python
from deabook.MQDEAweak import MQDEAweak

mq = MQDEAweak(data, id="Firm", year="Year",
               sent="K L = Y : CO2", gy=[0], gx=[1, 1], gb=[0],
               rts="crs", tech="Global production technology")
result = mq.optimize()
print(result[['D11', 'mqpi']])
```

### 9.6 物质平衡分解

```python
from deabook.MB import MB, MBx, MBxy

# 基础 MB 分解
model = MB(data, sent="K L + E = Y : CO2",
           sx=[[1,1,1],[1,1,1]], sy=[[1],[1]],
           level=5, rts="vrs1")
data3, info = model.optimize(solver="mosek")

# 面板数据 MB 分解（投入端）
model_x = MBx(data, year="Year", sent="K L + E = Y : CO2",
              sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts="vrs1")
result_x = model_x.optimize(solver="mosek")

# 面板数据 MB 分解（投入+产出端）
model_xy = MBxy(data, year="Year", sent="K L + E = Y : CO2",
                sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts="vrs1")
result_xy = model_xy.optimize(solver="mosek")
```

### 9.7 CNLS 随机前沿 + StoNED 分解

```python
from deabook.CNLSSDFDDFweak import CNLSSDweak, CNLSDDFweak
from deabook.StoNED import StoNED

# Shephard 距离函数 CNLS
cnls_model = CNLSSDweak(data, sent="K L = Y : CO2",
                        z=None, gy=[0], gx=[1, 1], gb=[0],
                        cet="mult", fun="prod", rts="vrs1")
cnls_model.optimize(solver="mosek")

# StoNED 残差分解（拟似然估计）
stoned = StoNED(cnls_model)
mu = stoned.get_mean_of_inefficiency(method="QLE")
efficiency = stoned.get_technical_efficiency(method="QLE")

# 方向距离函数 CNLS
ddf_model = CNLSDDFweak(data, sent="K L = Y : CO2",
                        z=None, gy=[1], gx=[1], gb=[1],
                        fun="prod", rts="vrs1")
ddf_model.optimize(solver="mosek")

# StoNED DDF 效率比率
stoned_ddf = StoNED(ddf_model)
te_ratio = stoned_ddf.get_technical_efficiency_ratio(method="QLE")

# StoNED KDE 方法（Richardson-Lucy 盲反卷积）
kde_inefficiency = stoned.get_technical_inefficiency(method="KDE")
kde_efficiency = stoned.get_technical_efficiency(method="KDE")
```
