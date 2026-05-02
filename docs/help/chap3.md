# 第三章：DEA 静态效率分析（资本效率 / 能源效率）

**主题**：使用径向 DEA、超效率 DEA 和 DDF 测量静态效率。

**数据文件**：
- `oecd data.xlsx`（33个OECD国家，2016-2019）
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本 (Capital)
- `L`：劳动 (Labor)
- `E`：能源 (Energy)
- `Y`：产出 (GDP)
- `CO2`：二氧化碳排放（非期望产出，本章不使用）

---

## deabook 包介绍与安装

deabook 是一个用于数据包络分析（DEA）的 Python 包，提供径向 DEA、方向距离函数（DDF）、Malmquist 指数、弱可处置性模型、物质平衡分解、CNLS 随机前沿等功能。

### 安装

```bash
pip install deabook
```

### 本章需要的导入

```python
from deabook import DEA
from deabook.constant import RTS_VRS1, RTS_CRS, OPT_LOCAL
import pandas as pd
import numpy as np
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表。

---

### 参数详解

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | 包含 sent 中指定列的数据，如K, L, E, Y  |
| `sent` | str | `"投入1 投入2 ...=产出1 产出2 ..."` | `=` 左边为投入变量，`=`右边为产出变量 |
| | | | |
| `gy` | list[int] | `[0]` | 不沿该产出方向调整 |
| | list[int] | `[1]` | 沿该产出方向扩大产出；与 `gx`≥1 同时使用时（Hyper 模式），表示沿产出方向同步扩大 |
| | | | |
| `gx` | list[int] | `[0]` | 不沿该投入方向调整 |
| | list[int] | `[1]` | 沿该投入方向缩减 |
| | | | |
| **gy/gx 组合** | **判断规则** | **`sum(gx)>=1`** | 至少缩减一个投入 |
| | | **`sum(gy)>=1`** | 至少扩大一个产出 |
| | | | |
| | **组合** | ✗ ✓ → 产出导向 | 不缩减投入，只扩大产出 |
| | | ✓ ✗ → 投入导向 | 缩减投入中 gx=1 对应的变量 |
| | | ✓ ✓ → 超效率(Hyper) | 同时缩减投入、扩大产出 |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬 |
| | | | |
| `baseindex` | str | `None` | 被评价 DMU 筛选条件，如 `"t==2019"`，`None` 表示全部评价 |
| | str | `"t=[2017,2018,2019]"` | 仅评价指定时期 |
| | | | |
| `refindex` | str | `None` | 参考技术筛选条件，如 `"t==2019"`，`None` 表示使用全部 DMU 构建前沿 |
| | str | `"t=[2017,2018,2019]"` | 仅用指定时期构建参考前沿 |
| | | | |
| `email` | str | `OPT_LOCAL` | 本地求解 |
| | str | 邮箱地址 | 远程求解通知 |
| | | | |
| `solver` | str | `"mosek"` | 商业求解器，速度快、精度高 |
| | str | `"glpk"` | 开源求解器 |
| | str | `"cbc"` | 开源求解器 |

### 返回值

`optimize()` 返回 DataFrame，列结构取决于 **gy/gx 导向** 和 **rts**：

#### DEA2 返回值

| 情况 | 返回列 | 效率公式 |
|------|--------|---------|
| 投入导向（gx≥1, gy==0） | `optimization_status`, `rho`, `te` | `te = rho`（范围 0-1，1 表示有效） |
| 产出导向（gy≥1, gx==0） | `optimization_status`, `rho`, `te` | `te = 1/rho` |
| 超效率 + CRS（gx≥1, gy≥1） | `optimization_status`, `rho`, `te` | `te = sqrt(rho)` |
| 超效率 + VRS（gx≥1, gy≥1） | `optimization_status`, `rho`, `tei`, `teo` | `tei = 1 - rho`（投入端），`teo = 1/(1+rho)`（产出端） |

> **Hyper 不是只能搭配 VRS**：CRS 和 VRS 都支持 Hyper。区别在于 CRS 返回单一 `te = sqrt(rho)`，VRS 返回 `tei` + `teo` 双列。

#### DDF2 返回值

| 情况 | 返回列 | 效率公式 |
|------|--------|---------|
| 投入导向 | `optimization_status`, `rho`, `tei` | `tei = 1 - rho` |
| 产出导向 | `optimization_status`, `rho`, `teo` | `teo = 1/(1+rho)` |
| 超效率 | `optimization_status`, `rho`, `tei`, `teo` | 同上 |

#### DEA2 vs DDF2 对比

| 特征 | DEA2 | DDF2 |
|------|------|------|
| 类型 | 径向 DEA | 方向距离函数 |
| 投入导向效率 | `te = rho` | `tei = 1 - rho` |
| 产出导向效率 | `te = 1/rho` | `teo = 1/(1+rho)` |

---

## 模型 1：投入导向 DEA — CRS（纯资本效率 PEO）

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

---

## 模型 2：投入导向 DEA — VRS（纯资本效率 PEO）

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

与模型 1 唯一区别：`rts=RTS_VRS1`（可变规模报酬）。

---

## 模型 3：超效率 DEA (Hyperbolic) — CRS

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[1], gx=[1,0,0], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

---

## 模型 4：超效率 DEA (Hyperbolic) — VRS

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[1], gx=[1,0,0], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

---

## 模型 5：DDF 投入导向 — CRS

```python
model = DEA.DDF2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

---

## 模型 6：DDF 投入导向 — VRS

```python
model = DEA.DDF2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

与模型 5 唯一区别：`rts=RTS_VRS1`（可变规模报酬）。

---
