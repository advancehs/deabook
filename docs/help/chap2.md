# 第二章：DEA 静态效率分析（资本效率 / 能源效率）

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

## 导入

```python
from deabook import DEA
from deabook.constant import RTS_VRS1, RTS_CRS, OPT_LOCAL
import pandas as pd
import numpy as np
```

**导入说明**：
- `DEA`：核心 DEA 模块，包含 `DEA2`（径向 DEA）、`DDF2`（方向距离函数）
- `RTS_VRS1`：可变规模报酬（相同缩减因子）
- `RTS_CRS`：不变规模报酬
- `OPT_LOCAL`：本地求解（默认值）

---

## 模型 1：投入导向 DEA — CRS（纯资本效率 PEO）

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `data` | DataFrame | 包含 K、L、E、Y 列的数据 |
| `sent` | `"K L E=Y"` | `=` 左边为投入（K, L, E），右边为产出（Y） |
| `gy` | `[0]` | 产出方向向量。`[0]` 表示不沿产出方向调整（投入导向） |
| `gx` | `[1,0,0]` | 投入方向向量。`[1,0,0]` 表示只沿第一个投入（K，资本）缩减 |
| `rts` | `RTS_CRS` | 不变规模报酬 |

### 返回值

`res` 为 DataFrame，包含列：
- `optimization_status`：求解状态（`"ok"` 表示成功）
- `rho`：效率原始值
- `te`：技术效率（投入导向：`te = rho`，范围 0-1，1 表示有效）

### gy/gx 方向向量说明

`gy` 和 `gx` 是**布尔型方向向量**，长度分别等于产出和投入的个数：
- `1` 表示沿该变量方向调整
- `0` 表示不调整

通过设置不同的方向向量组合实现不同导向：

| 导向 | gy | gx | 说明 |
|------|----|----|------|
| 投入导向 | `[0]` | `[1,0,0]` | 只缩减资本 K，L 和 E 不变 |
| 产出导向 | `[1]` | `[1,0,0]` | 同时缩减 K、扩大 Y |
| 超效率(Hyper) | `[1]` | `[1,0,0]` | 投入和产出同时调整 |

---

## 模型 2：投入导向 DEA — VRS（纯资本效率 PEO）

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

与模型 1 唯一区别：`rts=RTS_VRS1`（可变规模报酬）。
返回列相同：`optimization_status`, `rho`, `te`。

---

## 模型 3：超效率 DEA (Hyperbolic) — CRS

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[1], gx=[1,0,0], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `gy` | `[1]` | 沿产出方向调整（`sum(gy)>=1`） |
| `gx` | `[1,0,0]` | 沿资本方向调整（`sum(gx)>=1`） |
| `rts` | `RTS_CRS` | 不变规模报酬 |

### 返回值

`sum(gx)>=1` 且 `sum(gy)>=1` 触发**超效率(Hyperbolic)**模式：
- CRS 时返回列：`optimization_status`, `rho`, `te`（`te = sqrt(rho)`）

---

## 模型 4：超效率 DEA (Hyperbolic) — VRS

```python
model = DEA.DEA2(data, sent="K L E=Y", gy=[1], gx=[1,0,0], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

### 返回值

VRS + Hyperbolic 模式返回**两个效率指标**：
- `optimization_status`：求解状态
- `rho`：距离函数原始值
- `tei`：投入端技术效率（`tei = 1 - rho`）
- `teo`：产出端技术效率（`teo = 1 / (1 + rho)`）

---

## 模型 5：DDF 投入导向 — CRS

```python
model = DEA.DDF2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_CRS)
res = model.optimize(solver="mosek")
```

### DEA2 vs DDF2 区别

| 特征 | DEA2 | DDF2 |
|------|------|------|
| 类型 | 径向 DEA | 方向距离函数 |
| 效率计算 | `te = rho` 或 `te = 1/rho` | `tei = 1 - rho` 或 `teo = 1/(1+rho)` |
| 投入导向返回 | `te` | `tei` |
| 产出导向返回 | `te` | `teo` |
| 超效率返回 | `tei`, `teo` | `tei`, `teo` |

### 返回值

投入导向（`sum(gx)>=1`, `sum(gy)==0`）：
- `optimization_status`, `rho`, `tei`（`tei = 1 - rho`）

---

## 模型 6：DDF 投入导向 — VRS

```python
model = DEA.DDF2(data, sent="K L E=Y", gy=[0], gx=[1,0,0], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

返回列与模型 5 相同：`optimization_status`, `rho`, `tei`。

---

## DEA2 类完整方法列表

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)` | email: 求解邮箱; solver: 求解器名 | DataFrame | 求解并返回效率结果 |
| `display_status()` | 无 | 打印 | 显示求解状态 |
| `display_rho()` | 无 | 打印 | 显示 rho 值 |
| `display_lamda()` | 无 | 打印 | 显示强度变量 |
| `get_status()` | 无 | 求解状态 | 获取优化状态 |
| `get_rho()` | 无 | numpy array | 获取效率原始值 |
| `get_lamda()` | 无 | DataFrame | 获取强度变量 |
| `info(dmu)` | dmu: DMU 索引 | 打印 | 查看指定 DMU 的 LP 模型 |

## DDF2 类完整方法列表

与 DEA2 相同：optimize, display_status, display_rho, display_lamda, get_status, get_rho, get_lamda, info(dmu)。

---

## 求解器说明

| 求解器 | 适用场景 | 说明 |
|--------|---------|------|
| `"mosek"` | DEA、DDF 线性规划 | 商业求解器，速度快、精度高 |
| `"glpk"` | 线性规划（教学） | 开源求解器，速度较慢 |
| `"knitro"` | 非线性规划 | 用于 CNLS/StoNED 模型 |
