# 第十四章：影子价格 / 对偶模型从零实现

**主题**：使用 Pyomo 构建 SBM 和 NDDF 的对偶模型，估计投入/产出的影子价格。

**数据文件**：
- `Ex4.dta`（K, L, Y, CO2）
- `Ex5.dta`（labor, capital, energy, gdp, co2）

**说明**：本章不使用 deabook 包，所有函数为 notebook 内自定义实现。

---

## 导入

```python
import pandas as pd
import numpy as np
from pyomo.environ import ConcreteModel, Set, Var, Objective, Constraint, Reals
```

---

## 函数 1：sbmdual — SBM 对偶模型（基础版）

```python
def sbmdual(dataframe, varname, P, Q):
    ...
```

### 示例

```python
res = sbmdual(data, "DMU_name", P=2, Q=1)
```

### 参数

| 参数 | 类型 | 含义 |
|------|------|------|
| `dataframe` | DataFrame | 数据 |
| `varname` | str | DMU 名称列名 |
| `P` | int | 投入变量个数 |
| `Q` | int | 期望产出变量个数 |

### 返回值

DataFrame，包含列：
- `obj`：目标函数值
- `shadow price K`：资本影子价格
- `shadow price L`：劳动影子价格
- `shadow price Y`：产出影子价格
- `shadow price CO2`：CO2 影子价格（如含非期望产出）

---

## 函数 2：sbmdual2 — SBM 对偶模型（含评价/参考分离）

```python
def sbmdual2(dataframe, varname, P, Q, evaquery, refquery):
    ...
```

### 示例

```python
res = sbmdual2(data, "DMU_name", P=2, Q=1, evaquery="year==2017", refquery="year<=2017")
```

### 参数

| 参数 | 含义 |
|------|------|
| `dataframe` | DataFrame |
| `varname` | DMU 名称列 |
| `P` | 投入个数 |
| `Q` | 期望产出个数 |
| `evaquery` | 评价 DMU 筛选 |
| `refquery` | 参考 DMU 筛选 |

### sbmdual vs sbmdual2

| 特征 | sbmdual | sbmdual2 |
|------|---------|---------|
| 评价/参考分离 | 不支持 | 支持 evaquery/refquery |
| 灵活性 | 所有 DMU 同时评价和参考 | 可分别指定 |

---

## 函数 3：sbmdual3 — SBM 对偶模型（公式版，含非期望产出）

```python
def sbmdual3(formula, dataframe, evaquery=None, refquery=None):
    ...
```

### 示例

```python
res = sbmdual3("Y:CO2~K L", data)
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y:CO2~K L"` |
| `dataframe` | DataFrame |
| `evaquery` | 评价 DMU 筛选 |
| `refquery` | 参考 DMU 筛选 |

### sbmdual3 vs sbmdual/sbmdual2

| 特征 | sbmdual/sbmdual2 | sbmdual3 |
|------|-----------------|---------|
| 变量指定 | 显式 P, Q | 公式语法 |
| 非期望产出 | 隐含 | 显式（`:` 语法） |
| 列名 | 固定格式 `shadow price X` | 根据变量名动态命名 |

### 返回值

DataFrame，列名根据公式中的变量名动态生成：
- `obj`：目标函数值
- `shadow price K`、`shadow price L`：投入影子价格
- `shadow price Y`：期望产出影子价格
- `shadow price CO2`：非期望产出影子价格

---

## 函数 4：nddfdual — NDDF 对偶模型（基础版）

```python
def nddfdual(dataframe, varname, P, Q):
    ...
```

### 示例

```python
res = nddfdual(data, "DMU_name", P=2, Q=1)
```

### 参数

| 参数 | 含义 |
|------|------|
| `dataframe` | DataFrame |
| `varname` | DMU 名称列 |
| `P` | 投入个数 |
| `Q` | 期望产出个数 |

### 返回值

与 sbmdual 相同格式的 DataFrame。

### nddfdual vs sbmdual

| 特征 | sbmdual | nddfdual |
|------|---------|---------|
| 原始模型 | SBM（松弛模型） | NDDF（非径向方向距离函数） |
| 对偶结构 | 基于 SBM LP | 基于 NDDF LP |
| 影子价格含义 | SBM 松弛的边际值 | NDDF 方向距离的边际值 |

---

## 函数 5：nddfdual2 — NDDF 对偶模型（公式版）

```python
def nddfdual2(formula, dataframe, gx=None, gy=None, gb=None, weight=None, evaquery=None, refquery=None):
    ...
```

### 示例

```python
res = nddfdual2("Y:CO2~K L", data, gx=[1,1], gy=[1], gb=[1], weight=[1/4,1/4,1/4,1/4])
```

### 参数

| 参数 | 含义 |
|------|------|
| `formula` | `"Y:CO2~K L"` |
| `gx` | 投入方向向量 |
| `gy` | 期望产出方向向量 |
| `gb` | 非期望产出方向向量 |
| `weight` | 各变量权重 |
| `evaquery` | 评价 DMU 筛选 |
| `refquery` | 参考 DMU 筛选 |

### 返回值

DataFrame，列名根据变量名动态生成：
- `obj`：目标函数值
- 逐变量影子价格列

---

## 影子价格的经济含义

影子价格表示在最优解处，每增加一个单位的约束资源对目标函数的边际贡献：

| 变量类型 | 影子价格含义 |
|---------|------------|
| 投入 | 增加一单位投入的边际成本 |
| 期望产出 | 增加一单位产出的边际收益 |
| 非期望产出 | 减少一单位排放的边际减排成本 |

---

## 公式语法汇总

所有公式版函数使用统一语法：

| 语法 | 含义 |
|------|------|
| `"Y~K L"` | 产出 Y，投入 K 和 L |
| `"Y:CO2~K L"` | 产出 Y，非期望产出 CO2，投入 K 和 L |
| `"Y1 Y2~K L E"` | 多产出，多投入 |
