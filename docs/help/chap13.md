# 第十三章：影子价格与对偶模型

**主题**：使用 Pyomo 构建 SBM 和 NDDF 的对偶模型，估计投入、产出和非期望产出的影子价格。

**数据文件**：
- `Ex4.dta`（30个中国省份，3期，含 K,L,Y,CO2）
- `Ex5.dta`（含 labor,capital,energy,gdp,co2）

**变量含义**：
- `K`：资本, `L`：劳动, `Y`：产出(GDP), `CO2`：碳排放（非期望产出）
- `labor`：劳动, `capital`：资本, `energy`：能源, `gdp`：GDP, `co2`：CO2

---

## 安装与导入

本章所有模型使用 Pyomo 从零构建对偶模型，不依赖 deabook 包。

```bash
pip install pyomo mosek
```

```python
from pyomo.environ import *
import pandas as pd
import numpy as np
import re
import warnings
warnings.filterwarnings("ignore")
```

---

## 对偶模型概念

对偶模型（影子价格模型）是 SBM/NDDF 原始模型的对偶问题，通过求解对偶问题获得各变量的边际价值（影子价格）：

| 模型 | 原始问题 | 对偶问题 |
|------|---------|---------|
| SBM | 最小化效率损失 | 最大化 pomega |
| NDDF | 最大化方向距离 | 最小化总成本 |

---

### 参数详解

#### SBM 对偶系列

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `dataframe` | DataFrame | — | DMU 数据（sbmdual / sbmdual2） |
| `varname` | list | `["K","L","Y","CO2"]` | 变量名列表 |
| `P` | int | `2` | 投入变量个数 |
| `Q` | int | `1` | 期望产出变量个数 |
| `formula` | str | `"Y:CO2=K L"` | 公式语法（仅 sbmdual3） |
| `evaquery` | str | `"t==[1,2]"` | 评价 DMU 筛选 |
| `refquery` | str | `"t==[1,2,3]"` | 参考集 DMU 筛选 |

#### NDDF 对偶系列

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `dataframe` | DataFrame | — | DMU 数据（nddfdual） |
| `varname` | list | `["labor","capital","energy","gdp","co2"]` | 变量名列表 |
| `P` | int | `3` | 投入变量个数 |
| `Q` | int | `1` | 期望产出变量个数 |
| `formula` | str | `"gdp:co2=labor capital energy"` | 公式语法（仅 nddfdual2） |
| `gx` | list | `[-1,-1,-1]` | 投入方向向量 |
| `gy` | list | `[1]` | 期望产出方向向量 |
| `gb` | list | `[-1]` | 非期望产出方向向量 |
| `weight` | list | `None` | 权重矩阵 |
| `evaquery` | str | `None` | 评价 DMU 筛选 |
| `refquery` | str | `None` | 参考集 DMU 筛选 |

---

### 返回值

#### SBM 对偶

| 返回列 | 含义 |
|--------|------|
| `obj` | 目标函数值 |
| `shadow price K`, `shadow price L` 等 | 投入影子价格 |
| `shadow price Y` 等 | 期望产出影子价格 |
| `shadow price CO2` 等 | 非期望产出影子价格 |

#### NDDF 对偶

| 返回列 | 含义 |
|--------|------|
| `obj` | 目标函数值 |
| `shadow price labor`, `shadow price capital` 等 | 投入影子价格 |
| `shadow price gdp` 等 | 期望产出影子价格 |
| `shadow price co2` 等 | 非期望产出影子价格 |

### 影子价格解释

| 影子价格 | 含义 |
|---------|------|
| `px`（投入） | 增加一单位投入的边际成本 |
| `py`（期望产出） | 增加一单位期望产出的边际收益 |
| `pb`（非期望产出） | 减少一单位非期望产出的边际成本（即排放成本） |

---

## SBM 对偶 vs NDDF 对偶

| 特征 | SBM 对偶 | NDDF 对偶 |
|------|---------|-----------|
| 目标函数 | maximize pomega | minimize sum(px*x) - sum(py*y) + sum(pb*b) + pomega |
| 权重 | 无 | 使用 weight 参数 |
| 方向向量 | 无 | 使用 gx/gy/gb |
| 约束数 | 5 类 | 4 类 |

---

## 模型 1：sbmdual — 基于列名

```python
ex4 = pd.read_stata("Ex4.dta")
data = sbmdual(ex4, ["K","L","Y","CO2"], 2, 1)
```

---

## 模型 2：sbmdual2 — 支持筛选

```python
data = sbmdual2(ex4, ["K","L","Y","CO2"], 2, 1, "t==[1,2]", "t==[1,2,3]")
```

---

## 模型 3：sbmdual3 — formula 语法

```python
data = sbmdual3("Y:CO2=K L", ex4, "t==[1,2]", "t==[1,2,3]")
```

---

## 模型 4：nddfdual — 基于列名

```python
ex5 = pd.read_stata("Ex5.dta")
data = nddfdual(ex5, ["labor","capital","energy","gdp","co2"], 3, 1)
```

---

## 模型 5：nddfdual2 — formula + 方向向量

```python
data = nddfdual2("gdp:co2=labor capital energy", ex5, gx=[-1,-1,-1], gy=[1], gb=[-1])
```

---

## 函数速查表

| 函数 | 输入格式 | 特点 |
|------|---------|------|
| `sbmdual` | varname 列表 + P, Q | SBM 对偶，全样本 |
| `sbmdual2` | varname 列表 + P, Q + evaquery/refquery | SBM 对偶，支持筛选 |
| `sbmdual3` | formula 字符串 | SBM 对偶，公式语法 |
| `nddfdual` | varname 列表 + P, Q | NDDF 对偶，全样本 |
| `nddfdual2` | formula + gx/gy/gb/weight | NDDF 对偶，完整参数 |
