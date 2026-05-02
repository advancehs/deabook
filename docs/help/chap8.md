# 第八章：基于松弛测度的效率模型（SBM）

**主题**：使用 SBM（Slacks-Based Measure）模型测量非径向效率，包含投入导向、产出导向、非期望产出导向和综合非导向 SBM。

**数据文件**：
- `oecd data.xlsx`（示例数据，含 K, L, E, Y, CO2, year 列）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出（期望产出）, `CO2`：碳排放（非期望产出）

---

## deabook 包介绍与安装

deabook 是一个用于数据包络分析（DEA）的 Python 包，提供径向 DEA、方向距离函数（DDF）、Malmquist 指数、弱可处置性模型、物质平衡分解、CNLS 随机前沿、SBM 等功能。

### 安装

```bash
pip install deabook
```

### 本章需要的导入

```python
from deabook.SBM import SBM
from deabook.constant import RTS_CRS, RTS_VRS1
import pandas as pd
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表。



---

### 参数详解

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | 包含 sent 中指定列的数据 |
| `sent` | str | `"K L=Y"` | 不含非期望产出：`=` 左边为投入，`=` 右边为期望产出 |
| | str | `"K L E=Y:CO2"` | 含非期望产出：`:` 后为非期望产出 |
| | | | |
| `gy` | list[int] | `[0]` | 不调整期望产出 |
| | list[int] | `[1]` | 沿期望产出方向扩大（松弛进入目标函数） |
| | | | |
| `gx` | list[int] | `[0,0]` | 不调整投入 |
| | list[int] | `[1,0]` | 仅第一个投入允许松弛 |
| | list[int] | `[1,1]` | 所有投入允许松弛 |
| | | | |
| `gb` | list[int] | `[0]` | 不调整非期望产出 |
| | list[int] | `[1]` | 沿非期望产出方向缩减（松弛进入目标函数） |
| | | | |
| **gy/gx/gb 组合** | **方向模式** | **条件** | **说明** |
| | `input_oriented` | gx≥1, gy\==0, gb==0 | 投入导向 SBM（仅投入松弛进入目标函数） |
| | `output_oriented` | gy≥1, gx\==0, gb==0 | 期望产出导向 SBM |
| | `unoutput_oriented` | gb≥1, gx\==0, gy==0 | 非期望产出导向 SBM |
| | `hyper_orientedyx` | gx≥1, gy≥1, gb==0 | 投入+期望产出综合 SBM |
| | `hyper_orientedyb` | gb≥1, gy≥1, gx==0 | 期望产出+非期望产出综合 SBM |
| | `hyper_orientedxb` | gx≥1, gb≥1, gy==0 | 投入+非期望产出综合 SBM |
| | `hyper_orientedyxb` | gx≥1, gy≥1, gb≥1 | 全方向综合 SBM |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬 |
| | | | |
| `baseindex` | str | `None` | 被评价 DMU 筛选条件，如 `"t=[1,2,3]"` |
| `refindex` | str | `None` | 参考 DMU 筛选条件，如 `"t=[1,2,3]"` |
| | | | |
| `solver` | str | `"mosek"` | 商业求解器 |
| | str | `"glpk"` | 开源求解器 |

### 返回值

`optimize()` 返回 DataFrame，列结构取决于 **方向模式** 和 **模型类型**：

| 模型类型 | 触发条件 | 返回列 | 效率公式 |
|----------|----------|--------|---------|
| `ratio` | gx≥1 且 (gy≥1 或 gb≥1) | `optimization_status`, `direction`, `objective_value`, `rho`, `t` | `rho = τ*`（Charnes-Cooper 变换后目标值） |
| `input_only` | gx≥1, gy\==0, gb==0 | `optimization_status`, `direction`, `objective_value`, `rho` | `rho = 1 - avg(s⁻/x₀)` |
| `output_bad_only` | gx==0, (gy≥1 或 gb≥1) | `optimization_status`, `direction`, `objective_value`, `rho`, `phi` | `phi = 1 + avg(s⁺/y₀)`, `rho = 1/phi` |

所有方向均追加松弛列：

| 列 | 含义 |
|----|------|
| `slack_x_变量名` | 投入冗余（如 `slack_x_K`, `slack_x_L`） |
| `slack_y_变量名` | 期望产出不足（如 `slack_y_Y`） |
| `slack_b_变量名` | 非期望产出过剩（如 `slack_b_CO2`） |

### 辅助方法

| 方法 | 返回值 |
|------|--------|
| `get_rho()` | SBM 效率值 Series |
| `get_slacks()` | 所有松弛变量 DataFrame |
| `get_lamda()` | 强度变量（λ）DataFrame |
| `get_status()` | 求解状态 Series |
| `info(dmu="all")` | 打印指定 DMU 的 Pyomo 模型详情 |

### SBM 与径向 DEA 的区别

| 特征 | 径向 DEA / DDF | SBM |
|------|----------------|-----|
| 调整方式 | 按统一比例因子同时收缩/扩张 | 每个变量独立松弛 |
| 松弛处理 | 松弛不影响效率值 | 松弛直接进入效率值 |
| 效率值范围 | 投入导向 0-1 | 0-1（更严格，通常 ≤ 径向效率） |
| 适用场景 | 比例性改进 | 存在非比例性冗余/不足时 |

---

## 模型 1：投入导向 SBM — CRS

```python
model = SBM(data, sent="K L=Y", gy=[0], gx=[1,1], rts=RTS_CRS,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

---

## 模型 2：投入导向 SBM — VRS

```python
model = SBM(data, sent="K L=Y", gy=[0], gx=[1,1], rts=RTS_VRS1,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

与模型 1 唯一区别：`rts=RTS_VRS1`（可变规模报酬，增加凸性约束 $\sum_n \lambda_n = 1$）。

---

## 模型 3：期望产出导向 SBM — CRS

```python
model = SBM(data, sent="K L=Y", gy=[1], gx=[0,0], rts=RTS_CRS,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

---

## 模型 4：期望产出导向 SBM — VRS

```python
model = SBM(data, sent="K L=Y", gy=[1], gx=[0,0], rts=RTS_VRS1,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

---

## 模型 5：投入与期望产出综合 SBM — CRS

```python
model = SBM(data, sent="K L=Y", gy=[1], gx=[1,1], rts=RTS_CRS,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

此方向对应 `hyper_orientedyx`，同时考虑投入冗余和产出不足，返回 `rho`（分式 SBM 效率）和 `t`（Charnes-Cooper 变换变量）。

---

## 模型 6：投入与期望产出综合 SBM — VRS

```python
model = SBM(data, sent="K L=Y", gy=[1], gx=[1,1], rts=RTS_VRS1,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

---

## 模型 7：含非期望产出的非期望产出导向 SBM — CRS

```python
model = SBM(data, sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1], rts=RTS_CRS,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

仅非期望产出松弛进入目标函数，衡量污染物可削减空间。

---

## 模型 8：含非期望产出的非期望产出导向 SBM — VRS

```python
model = SBM(data, sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1], rts=RTS_VRS1,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

---

## 模型 9：含非期望产出的全方向综合 SBM — CRS

```python
model = SBM(data, sent="K L E=Y:CO2", gy=[1], gx=[1,1,1], gb=[1], rts=RTS_CRS,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

对应 `hyper_orientedyxb`，投入、期望产出和非期望产出的松弛全部进入目标函数。

---

## 模型 10：含非期望产出的全方向综合 SBM — VRS

```python
model = SBM(data, sent="K L E=Y:CO2", gy=[1], gx=[1,1,1], gb=[1], rts=RTS_VRS1,
            baseindex=None, refindex=None)
res = model.optimize(solver="mosek")
```

---

## 数据要求

SBM 的目标函数中存在 $s/x_0$、$s/y_0$、$s/b_0$ 等相对松弛比率，因此所有投入、期望产出和非期望产出变量**必须为严格正数**（> 0）。若数据中存在 0 或负数，需先进行平移或重新定义变量。

---

**参考文献**：
- Tone, K. (2001). A slacks-based measure of efficiency in data envelopment analysis. *European Journal of Operational Research*, 130(3), 498-509.
- Tone, K. (2004). Dealing with undesirable outputs in DEA: A slacks-based measure (SBM) approach. *GRIPS Research Report Series*.

---

## 本章说明

本章介绍了 SBM 模型的数学原理、参数设定和 Python 实现。通过 `deabook.SBM` 模块，可以灵活设置投入导向、产出导向、非期望产出导向和综合非导向 SBM 模型。
