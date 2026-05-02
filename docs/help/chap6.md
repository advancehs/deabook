# 第六章：环境绩效测度方法（非期望产出）

**主题**：使用弱可处置性 DEA 模型测量环境效率，包含 DEAweak、DDFweak、NDDFweak。

**数据文件**：
- `oecd data.xlsx`（33个OECD国家，2016-2019）
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP), `CO2`：碳排放（非期望产出）

---

## deabook 包介绍与安装

deabook 是一个用于数据包络分析（DEA）的 Python 包，提供径向 DEA、方向距离函数（DDF）、Malmquist 指数、弱可处置性模型、物质平衡分解、CNLS 随机前沿等功能。

### 安装

```bash
pip install deabook
```

### 本章需要的导入

```python
from deabook import DEAweak
from deabook.constant import RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_LOCAL
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表。



---

### 参数详解

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | 包含 sent 中指定列的数据 |
| `sent` | str | `"K L E=Y:CO2"` | `=` 左边：投入（K, L, E）<br>`=` 右边 `:` 前：期望产出（Y）<br>`:` 后：非期望产出（CO2） |
| | | | |
| `gy` | list[int] | `[0]` | 不调整期望产出 |
| | list[int] | `[1]` | 沿期望产出方向扩大 |
| | | | |
| `gx` | list[int] | `[0,0,0]` | 不调整投入 |
| | list[int] | `[1,0,0]` | 沿第一个投入缩减 |
| | list[int] | `[1,1,1]` | 沿所有投入缩减 |
| | | | |
| `gb` | list[int] | `[0]` | 不调整非期望产出 |
| | list[int] | `[1]` | 沿非期望产出方向缩减（污染导向） |
| | | | |
| **gy/gx/gb 组合** | **方向模式** | **条件** | **说明** |
| | `unoutput_oriented` | gb≥1, gx\==0, gy==0 | 非期望产出导向 |
| | `input_oriented` | gx≥1, gy\==0, gb==0 | 投入导向 |
| | `output_oriented` | gy≥1, gx\==0, gb==0 | 产出导向 |
| | `hyper_orientedyx` | gx≥1, gy≥1, gb==0 | 投入+产出超效率 |
| | `hyper_orientedyb` | gb≥1, gy≥1, gx==0 | 产出+非期望超效率 |
| | `hyper_orientedxb` | gb≥1, gx≥1, gy==0 | 投入+非期望超效率 |
| | `hyper_orientedyxb` | gb≥1, gx≥1, gy≥1 | 全方向超效率 |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬（相同缩减因子） |
| | str | `RTS_VRS2` | 可变规模报酬（不同缩减因子，仅弱可处置性模型支持） |
| | | | |
| `email` | str | `OPT_LOCAL` | 本地求解 |
| | str | 邮箱地址 | 远程求解通知 |
| | | | |
| `solver` | str | `"mosek"` | 商业求解器 |
| | str | `"glpk"` | 开源求解器 |
| | | | |
| `baseindex` | str | `None` | 被评价 DMU 筛选条件，`None` 表示全部 |
| `refindex` | str | `None` | 参考技术筛选条件，`None` 表示全部 |

### 返回值

#### DEAweak2 返回值

所有方向均返回基础列 `optimization_status`, `rho`, `objective_value`，再根据方向追加效率列：

| 方向模式 | 追加列 | 效率公式 |
|----------|--------|---------|
| `input_oriented` | `tei` | `tei = 1 - rho` |
| `output_oriented` | `teo` | `teo = 1/(1+rho)` |
| `unoutput_oriented` | `te` | `te = rho` |
| `hyper_orientedyx` | `tei`, `teo` | 同上 |
| `hyper_orientedyb` | `teuo`, `teo` | 同上 |
| `hyper_orientedxb` | `tei`, `teuo` | 同上 |
| `hyper_orientedyxb` | `tei`, `teuo`, `teo` | 同上 |

#### DDFweak2 返回值

| 方向模式 | 返回列 | 效率公式 |
|----------|--------|---------|
| `input_oriented` | `optimization_status`, `rho`, `tei` | `tei = 1 - rho` |
| `output_oriented` | `optimization_status`, `rho`, `teo` | `teo = 1/(1+rho)` |
| `unoutput_oriented` | `optimization_status`, `rho`, `te` | `te = rho` |
| `hyper_orientedyx` | `optimization_status`, `rho`, `tei`, `teo` | 同上 |
| `hyper_orientedyb` | `optimization_status`, `rho`, `teuo`, `teo` | 同上 |
| `hyper_orientedxb` | `optimization_status`, `rho`, `tei`, `teuo` | 同上 |
| `hyper_orientedyxb` | `optimization_status`, `rho`, `tei`, `teuo`, `teo` | 同上 |

#### NDDFweak2 返回值

NDDF（非径向方向距离函数）为每个变量生成**独立的效率比率**：

| 列 | 含义 |
|----|------|
| `rhoY`, `rhoCO2` | 各变量缩减/扩张比率 |
| `teY`, `teCO2` | 各变量效率 |
| `weight_teY`, `weight_teCO2` | 加权效率 |
| `teo` | 产出端聚合效率 |
| `teuo` | 非期望产出端聚合效率 |
| `teuo2o` | 复合效率 `= teuo / teo` |
| `objective_value` | 目标函数值 |

#### RTS 差异

| RTS | 额外变量 | 约束 | 含义 |
|-----|---------|------|------|
| `RTS_CRS` | 无 | 无 VRS 约束 | 所有 DMU 同一规模 |
| `RTS_VRS1` | `theta`（标量） | `sum(lamda) = theta` | 相同缩减因子 |
| `RTS_VRS2` | `mu`（向量） | `sum(lamda+mu) = 1` | 不同缩减因子 |

---

## 模型 1：DEAweak2 非期望产出导向 — CRS/VRS1/VRS2

### CRS

```python
model = DEAweak.DEAweak2(data, sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

### VRS1

```python
model = DEAweak.DEAweak2(data, sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

### VRS2

```python
model = DEAweak.DEAweak2(data, sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1], rts=RTS_VRS2)
res = model.optimize(solver="mosek")
```

---

## 模型 2：DEAweak2 超效率(Hyper) — CRS/VRS1/VRS2

### CRS

```python
model = DEAweak.DEAweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

### VRS1

```python
model = DEAweak.DEAweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

### VRS2

```python
model = DEAweak.DEAweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS2)
res = model.optimize(solver="mosek")
```

与模型 1 唯一区别：`gy=[1]`，触发 Hyper 模式。

---

## 模型 3：DDFweak2 — CRS/VRS1/VRS2

```python
# CRS
model = DEAweak.DDFweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_CRS)
res = model.optimize(solver="mosek")

# VRS1
model = DEAweak.DDFweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS1)
res = model.optimize(solver="mosek")

# VRS2
model = DEAweak.DDFweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS2)
res = model.optimize(solver="mosek")
```

---

## 模型 4：NDDFweak2 — CRS/VRS1/VRS2

```python
# CRS
model = DEAweak.NDDFweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_CRS)
res = model.optimize(solver="mosek")

# VRS1
model = DEAweak.NDDFweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS1)
res = model.optimize(solver="mosek")

# VRS2
model = DEAweak.NDDFweak2(data, sent="K L E=Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS2)
res = model.optimize(solver="mosek")
```
