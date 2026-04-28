# 第十二章：Pyomo DEA 从零实现（dea1-dea8 渐进式开发）

**主题**：使用 Pyomo 从零开始渐进式实现 DEA 模型，从最简版本到完整的公式化 DEA。

**数据文件**：
- `Ex3.dta`（Stata 格式数据）

**说明**：本章不使用 deabook 包，所有函数为 notebook 内自定义实现。

---

## 导入

```python
import pandas as pd
import numpy as np
from pyomo.environ import ConcreteModel, Set, Var, Objective, Constraint, Reals
```

---

## 函数演进路线

本章按 8 个版本渐进实现 DEA：

| 版本 | 函数名 | 新增功能 |
|------|--------|---------|
| 1 | `dea1` | numpy 基础版本，固定变量 |
| 2 | `dea2` | 使用 DataFrame 输入/输出 |
| 3 | `dea3` | 增加 `evaquery`（评价 DMU 筛选） |
| 4 | `dea4` | 增加 `refquery`（参考 DMU 筛选） |
| 5 | `dea5` | 索引化重构（内部实现变化） |
| 6 | `dea6` | 公式语法 `"Y~K L"` 代替显式变量列表 |
| 7 | `dea7` | 增加 `rts` 参数（`"crs"`/`"vrs"`） |
| 8 | `dea8` | 增加 `orient` 参数（`"oo"`/`"io"`） |

---

## dea1：numpy 基础版本

```python
def dea1(data, xcol, ycol):
    ...
```

- 输入：DataFrame + 显式列名列表
- 最简单的 DEA LP 构建

---

## dea2：DataFrame 输入输出

```python
def dea2(data, xcol, ycol):
    ...
```

- 返回带 DataFrame 索引的结果

---

## dea3：增加 evaquery

```python
def dea3(data, xcol, ycol, evaquery=None):
    ...
```

### 参数

| 参数 | 含义 |
|------|------|
| `evaquery` | 评价 DMU 筛选条件，使用 DataFrame.query() 语法，如 `"year==2017"` |

---

## dea4：增加 refquery

```python
def dea4(data, xcol, ycol, evaquery=None, refquery=None):
    ...
```

### 参数

| 参数 | 含义 |
|------|------|
| `evaquery` | 评价 DMU 筛选条件 |
| `refquery` | 参考 DMU 筛选条件，如 `"year<=2017"` |

---

## dea5：索引化重构

```python
def dea5(data, xcol, ycol, evaquery=None, refquery=None):
    ...
```

内部实现改为基于索引的方式，API 不变。

---

## dea6：公式语法

```python
def dea6(formula, data, evaquery=None, refquery=None):
    ...
```

### 参数

| 参数 | 格式 | 含义 |
|------|------|------|
| `formula` | `"Y~K L"` | `~` 左边为产出，右边为投入（空格分隔） |
| `data` | DataFrame | 数据 |
| `evaquery` | str/None | 评价 DMU 筛选 |
| `refquery` | str/None | 参考 DMU 筛选 |

### formula 语法

```
"Y~K L"
```

- `~` 左边：产出变量（Y）
- `~` 右边：投入变量（K, L，空格分隔）

---

## dea7：增加 rts 参数

```python
def dea7(formula, data, evaquery=None, refquery=None, rts="crs"):
    ...
```

### 参数

| 参数 | 值 | 含义 |
|------|-----|------|
| `rts` | `"crs"` | 不变规模报酬（默认） |
| `rts` | `"vrs"` | 可变规模报酬（添加 `sum(lamda) = 1` 约束） |

---

## dea8：完整版（增加 orient 参数）

```python
def dea8(formula, data, evaquery=None, refquery=None, rts="crs", orient="oo"):
    ...
```

### 参数

| 参数 | 值 | 含义 |
|------|-----|------|
| `orient` | `"oo"` | 产出导向（Output Oriented）— 最大化产出 |
| `orient` | `"io"` | 投入导向（Input Oriented）— 最小化投入 |
| `rts` | `"crs"`/`"vrs"` | 规模报酬 |
| `formula` | `"Y~K L"` | 公式 |
| `evaquery` | str/None | 评价 DMU |
| `refquery` | str/None | 参考 DMU |

### 返回值

DataFrame，包含列：
- `optimization_status`：求解状态
- `efficiency`：效率值
- （产出导向：效率 ≥ 1，越大越有效）
- （投入导向：效率 ≤ 1，越小越有效）

### dea8 完整 LP 结构

**产出导向**：
- 目标：最大化 `theta`（产出缩放因子）
- 约束：`theta * y_i ≤ sum(lamda_j * y_j)`
- 投入约束：`x_i ≥ sum(lamda_j * x_j)`

**投入导向**：
- 目标：最小化 `theta`（投入缩放因子）
- 约束：`theta * x_i ≥ sum(lamda_j * x_j)`
- 产出约束：`y_i ≤ sum(lamda_j * y_j)`

**VRS 约束**：
- `sum(lamda_j) = 1`

---

## evaquery / refquery 语法

使用 pandas DataFrame.query() 兼容语法：

```python
evaquery = "year==2017"
refquery = "year<=2017"
evaquery = "province in ['北京','上海']"
```
