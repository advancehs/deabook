# 第十三章：SBM、DDF、NDDF、MPI 从零实现

**主题**：使用 Pyomo 从零实现 SBM、DDF、NDDF、Malmquist 指数和 Malmquist-Luenberger 指数。

**数据文件**：
- `Ex4.dta`（K, L, Y, CO2）
- `Ex5.dta`（labor, capital, energy, gdp, co2）

**说明**：本章不使用 deabook 包，所有函数为 notebook 内自定义实现。

---

## 导入

```python
import pandas as pd
import numpy as np
from pyomo.environ import ConcreteModel, Set, Var, Objective, Constraint, Reals, minimize, maximize
```

---

## 公式语法

与 dea8 相同，使用 `"Y~K L"` 格式，扩展支持非期望产出：

```
"Y~K L"           # 无非期望产出
"Y:CO2~K L"       # 含非期望产出（: 分隔）
```

- `~` 左边 `:` 前：期望产出（Y）
- `~` 左边 `:` 后：非期望产出（CO2）
- `~` 右边：投入（K, L，空格分隔）

---

## 函数 1：sbm — SBM 模型（无非期望产出）

```python
def sbm(formula, dataframe, evaquery=None, refquery=None):
    ...
```

### 示例

```python
res = sbm("Y~K L", data)
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y~K L"` |
| `dataframe` | DataFrame |
| `evaquery` | 评价 DMU 筛选 |
| `refquery` | 参考 DMU 筛选 |

### 返回值

DataFrame，包含效率值。

---

## 函数 2：sbm2 — SBM 模型（含非期望产出）

```python
def sbm2(formula, dataframe, evaquery=None, refquery=None):
    ...
```

### 示例

```python
res = sbm2("Y:CO2~K L", data)
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y:CO2~K L"`（含非期望产出） |
| `dataframe` | DataFrame |
| `evaquery` | 评价 DMU 筛选 |
| `refquery` | 参考 DMU 筛选 |

### 返回值

DataFrame，包含效率值。

### sbm vs sbm2

| 特征 | sbm | sbm2 |
|------|-----|------|
| 非期望产出 | 不支持 | 支持 |
| formula | `"Y~K L"` | `"Y:CO2~K L"` |
| 松弛变量 | 仅投入/产出 | 额外包含非期望产出松弛 |

---

## 函数 3：ddf — 方向距离函数

```python
def ddf(formula, dataframe, gx=None, gy=None, gb=None, evaquery=None, refquery=None):
    ...
```

### 示例

```python
res = ddf("Y:CO2~K L", data, gx=[1,1], gy=[1], gb=[1])
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y:CO2~K L"` |
| `gx` | 投入方向向量（如 `[1,1]`） |
| `gy` | 期望产出方向向量（如 `[1]`） |
| `gb` | 非期望产出方向向量（如 `[1]`） |
| `evaquery` | 评价 DMU 筛选 |
| `refquery` | 参考 DMU 筛选 |

### 返回值

DataFrame，包含方向距离函数值。

---

## 函数 4：nddf — 非径向方向距离函数

```python
def nddf(formula, dataframe, gx=None, gy=None, gb=None, weight=None, evaquery=None, refquery=None):
    ...
```

### 示例

```python
res = nddf("Y:CO2~K L", data, gx=[1,1], gy=[1], gb=[1], weight=[1/5,1/5,1/5,1/5,1/5])
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y:CO2~K L"` |
| `gx` | 投入方向向量 |
| `gy` | 期望产出方向向量 |
| `gb` | 非期望产出方向向量 |
| `weight` | 各变量权重向量（长度 = 投入数 + 产出数 + 非期望产出数） |
| `evaquery` | 评价 DMU 筛选 |
| `refquery` | 参考 DMU 筛选 |

### nddf vs ddf 区别

| 特征 | ddf | nddf |
|------|-----|------|
| 类型 | 径向 DDF | 非径向 DDF |
| 效率度量 | 单一 theta | 逐变量独立 rho |
| 权重 | 无 | `weight` 参数 |
| 返回粒度 | 整体效率 | 每个变量的独立效率 |

---

## 函数 5：dea8 — DEA 完整版

同第十二章的 `dea8`，此处作为 MPI 计算的底层模型。

```python
res = dea8("Y~K L", data, rts="crs", orient="oo")
```

---

## 函数 6：mpi — Malmquist 指数

```python
def mpi(formula, data, id, t):
    ...
```

### 示例

```python
res = mpi("Y~K L", data, id="province", t="year")
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y~K L"` |
| `data` | 面板数据 DataFrame |
| `id` | DMU 标识列名 |
| `t` | 时间列名 |

### 返回值

DataFrame，包含：
- `D11`：当期效率
- `MQ`：Malmquist 指数
- `MEFFCH`：效率变化
- `MTECHCH`：技术变化

---

## 函数 7：mpi2 — 增强版 Malmquist（多种技术选项）

```python
def mpi2(formula, data, id, t, tech=None):
    ...
```

### tech 参数

| 值 | 含义 |
|----|------|
| `None` | 默认技术 |
| `"com"` | 当期技术（Contemporary） |
| `"seq"` | 顺序技术（Sequential） |
| `"window N"` | 窗口技术，N 为窗口大小，如 `"window 3"` |

---

## 函数 8：mpi3 — 增强版 Malmquist（增加 global 技术）

```python
def mpi3(formula, data, id, t, tech=None):
    ...
```

### tech 参数

在 mpi2 基础上增加：

| 值 | 含义 |
|----|------|
| `"global"` | 全局技术（所有期 DMU 构成统一参考集） |
| 其余同 mpi2 | `"com"`, `"seq"`, `"window N"` |

---

## 函数 9：mlpi — Malmquist-Luenberger 指数

```python
def mlpi(formula, data, id, t, tech=None):
    ...
```

### 示例

```python
res = mlpi("Y:CO2~K L", data, id="province", t="year")
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y:CO2~K L"`（含非期望产出） |
| `data` | 面板数据 |
| `id` | DMU 标识列名 |
| `t` | 时间列名 |
| `tech` | 技术选项（同 mpi2） |

### mlpi vs mpi 区别

| 特征 | mpi | mlpi |
|------|-----|------|
| 底层模型 | dea8（DEA） | ddf（方向距离函数） |
| 非期望产出 | 不支持 | 支持 |
| formula | `"Y~K L"` | `"Y:CO2~K L"` |
| 效率距离 | DEA 径向 | DDF 方向距离 |

### 返回值

DataFrame，包含：
- `D11`：当期距离函数值
- `MLPI`：Malmquist-Luenberger 指数
- `MEFFCH`：效率变化
- `MTECHCH`：技术变化
