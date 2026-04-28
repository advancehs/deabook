# 第五章：环境效率（含非期望产出的静态分析）

**主题**：使用弱可处置性 DEA 模型测量环境效率，包含 DEAweak、DDFweak、NDDFweak。

**数据文件**：
- `oecd data.xlsx`（33个OECD国家，2016-2019）
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP), `CO2`：碳排放（非期望产出）

---

## 导入

```python
from deabook import DEAweak
from deabook.constant import RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_LOCAL
```

**导入说明**：
- `DEAweak`：弱可处置性 DEA 模块，包含 `DEAweak2`、`DDFweak2`、`NDDFweak2`
- `RTS_VRS2`：可变规模报酬（不同缩减因子）— 仅弱可处置性模型支持

---

## sent 公式格式

弱可处置性模型使用 `:` 分隔非期望产出：

```
"K L E=Y:CO2"
```

- `=` 左边：投入（K, L, E）
- `=` 右边 `:` 前：期望产出（Y）
- `=` 右边 `:` 后：非期望产出（CO2）

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

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `sent` | `"K L E=Y:CO2"` | 含非期望产出 |
| `gy` | `[0]` | 不调整期望产出 |
| `gx` | `[0,0,0]` | 不调整投入 |
| `gb` | `[1]` | 沿非期望产出方向调整（污染导向） |
| `rts` | `RTS_CRS`/`RTS_VRS1`/`RTS_VRS2` | 规模报酬 |

### 返回值

非期望产出导向（`sum(gb)>=1`, `sum(gx)==0`, `sum(gy)==0`）：
- `optimization_status`：求解状态
- `rho`：效率原始值
- `te`：技术效率（`te = rho`）

### RTS 规模报酬差异

| RTS | 额外变量 | 约束 | 含义 |
|-----|---------|------|------|
| `RTS_CRS` | 无 | 无 VRS 约束 | 所有 DMU 同一规模 |
| `RTS_VRS1` | `theta`（标量，bounds 0-1） | `sum(lamda) = theta` | 相同缩减因子 |
| `RTS_VRS2` | `mu`（向量，per ref DMU） | `sum(lamda+mu) = 1` | 不同缩减因子 |

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

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `gy` | `[1]` | 同时沿产出方向调整 |
| `gb` | `[1]` | 同时沿非期望产出方向调整 |
| 方向判定 | `hyper_orientedyb` | `sum(gb)>=1` 且 `sum(gy)>=1` 且 `sum(gx)==0` |

### 返回值

`hyper_orientedyb` 模式：
- CRS: `optimization_status`, `rho`, `te`（`te = sqrt(rho)`）
- VRS1: `optimization_status`, `rho`, `teuo`（`teuo = 1-rho`）, `teo`（`teo = 1/(1+rho)`）
- VRS2: `optimization_status`, `rho`, `teuo`（`teuo = 1-rho`）, `teo`（`teo = 1/(1+rho)`）

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

### DDFweak2 vs DEAweak2 区别

| 特征 | DEAweak2 | DDFweak2 |
|------|----------|----------|
| 类型 | 径向 DEA | 方向距离函数 |
| 效率计算 | `te = rho` 或 `te = 1/rho` | `tei/teo/teuo = f(rho)` |
| 方向模式数 | 5 种 | 7 种 |
| 额外方向 | — | `hyper_orientedxb`, `hyper_orientedyxb` |

### 返回值

`hyper_orientedyb` 模式（与 DEAweak2 相同参数）：
- `optimization_status`, `rho`, `objective_value`
- `teuo`（`teuo = 1-rho`）
- `teo`（`teo = 1/(1+rho)`）

### DDFweak2 完整方向模式

| 方向 | 条件 | 返回列 |
|------|------|--------|
| `input_oriented` | gx≥1, gy==0, gb==0 | `tei` |
| `output_oriented` | gy≥1, gx==0, gb==0 | `teo` |
| `unoutput_oriented` | gb≥1, gx==0, gy==0 | `teuo` |
| `hyper_orientedyx` | gx≥1, gy≥1, gb==0 | `tei`, `teo` |
| `hyper_orientedyb` | gb≥1, gy≥1, gx==0 | `teuo`, `teo` |
| `hyper_orientedxb` | gb≥1, gx≥1, gy==0 | `tei`, `teuo` |
| `hyper_orientedyxb` | gb≥1, gx≥1, gy≥1 | `tei`, `teuo`, `teo` |

---

## 模型 4：NDDFweak2 — CRS/VRS1/VRS2

```python
# CRS
model = DEAweak.NDDFweak2(data, sent="K L E = Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_CRS)
res = model.optimize(solver="mosek")

# VRS1
model = DEAweak.NDDFweak2(data, sent="K L E = Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS1)
res = model.optimize(solver="mosek")

# VRS2
model = DEAweak.NDDFweak2(data, sent="K L E = Y:CO2", gy=[1], gx=[0,0,0], gb=[1], rts=RTS_VRS2)
res = model.optimize(solver="mosek")
```

### NDDFweak2 特点

NDDF（非径向方向距离函数）为每个变量生成**独立的效率比率**，而非单一 rho 值。

### 返回值

`hyper_orientedyb` 模式返回 DataFrame，包含：
- 产出端逐变量列：`rhoY`, `teY`, `weight_teY`, 聚合为 `teo`
- 非期望产出端逐变量列：`rhoCO2`, `teCO2`, `weight_teCO2`, 聚合为 `teuo`
- 复合效率列：`teuo2o = teuo / teo`
- `objective_value`：目标函数值

### 逐变量效率列说明

对于每个变量（以产出 Y 为例）：
- `rhoY`：该变量的缩减/扩张比率
- `teY`：该变量效率，`teo`方向：`te = 1/(1+rho)`；`tei`方向：`te = 1-rho`
- `weight_teY`：`teY × 权重`
- `teo`：所有产出端 `weight_te` 之和

### NDDFweak2 完整方法列表

| 方法 | 返回值 | 说明 |
|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | DataFrame | 求解并返回逐变量效率 |
| `display_status()` | 打印 | 求解状态 |
| `display_objective_value()` | 打印 | 目标函数值 |
| `display_theta()` | 打印 | theta 值（VRS1） |
| `display_lamda()` | 打印 | 强度变量 |
| `display_mu()` | 打印 | mu 值（VRS2） |
| `display_rhox()` | 打印 | 投入端逐变量 rho |
| `display_rhoy()` | 打印 | 产出端逐变量 rho |
| `display_rhob()` | 打印 | 非期望产出端逐变量 rho |
| `get_objective_value()` | Series | 目标函数值 |
| `get_theta()` | Series | theta 值 |
| `get_lamda()` | DataFrame | 强度变量 |
| `get_mu()` | DataFrame | mu 值 |
| `get_rhox()` | DataFrame | 投入端 rho |
| `get_rhoy()` | DataFrame | 产出端 rho |
| `get_rhob()` | DataFrame | 非期望产出端 rho |
| `get_results_df()` | DataFrame | 完整结果 DataFrame |
