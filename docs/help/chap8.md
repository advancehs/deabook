# 第八章：CNLS/StoNED 随机前沿分析（碳排放效率）

**主题**：使用 CNLS（凸非参数最小二乘）和 StoNED（随机半非参数估计）测量碳排放效率。

**数据文件**：
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP), `CO2`：碳排放（非期望产出）

---

## 导入

```python
from deabook import CNLSSDFDDFweak, StoNED
from deabook.constant import FUN_PROD, OPT_LOCAL, RTS_VRS1, RTS_VRS2, RTS_CRS, CET_ADDI, CET_MULT, RED_QLE, RED_MOM, RED_KDE
import pandas as pd
```

**导入说明**：
- `CNLSSDFDDFweak`：CNLS 模块，包含 `CNLSSDweak`（Shephard 距离函数）和 `CNLSDDFweak`（方向距离函数）
- `StoNED`：随机半非参数效率估计模块
- `CET_MULT`：乘法复合误差项（CNLSSDweak 必须使用）
- `FUN_PROD`：生产前沿
- `RED_QLE`：准似然估计残差分解方法
- `RED_MOM`：矩方法残差分解
- `RED_KDE`：核密度估计残差分解

---

## sent 公式格式

与弱可处置性模型相同：

```
"K L E=Y:CO2"
```

---

## 第一部分：CNLSSDweak — Shephard 距离函数

CNLSSDweak 使用乘法复合误差项（`cet=CET_MULT`），求解器为 **knitro**（非线性规划）。

### 模型 1：CRS

```python
model = CNLSSDFDDFweak.CNLSSDweak(
    data, sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    cet=CET_MULT, rts=RTS_CRS, fun=FUN_PROD
)
model.optimize(solver="knitro")

rd = StoNED.StoNED(model)
res = rd.get_technical_efficiency(RED_QLE)
```

### 模型 2：VRS1

```python
model = CNLSSDFDDFweak.CNLSSDweak(
    data, sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    cet=CET_MULT, rts=RTS_VRS1, fun=FUN_PROD
)
model.optimize(solver="knitro")

rd = StoNED.StoNED(model)
res = rd.get_technical_efficiency(RED_QLE)
```

### 模型 3：VRS2

```python
model = CNLSSDFDDFweak.CNLSSDweak(
    data, sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    cet=CET_MULT, rts=RTS_VRS2, fun=FUN_PROD
)
model.optimize(solver="knitro")

rd = StoNED.StoNED(model)
res = rd.get_technical_efficiency(RED_QLE)
```

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `sent` | `"K L E=Y:CO2"` | 含非期望产出的公式 |
| `gy` | `[0]` | 不调整期望产出 |
| `gx` | `[0,0,0]` | 不调整投入 |
| `gb` | `[1]` | 沿非期望产出方向调整（碳排放导向） |
| `cet` | `CET_MULT` | 乘法复合误差项（**必须**，设为其他值会报错） |
| `rts` | `RTS_CRS`/`RTS_VRS1`/`RTS_VRS2` | 规模报酬 |
| `fun` | `FUN_PROD` | 生产前沿 |
| `solver` | `"knitro"` | 非线性求解器（CNLS 必须用 knitro） |

### 返回值

- `get_technical_efficiency(RED_QLE)` → numpy array：每个 DMU 的技术效率 TE（0-1，越大越好）

### TE 计算方式

```
TE = exp(-technical_inefficiency)
```

通过 StoNED 残差分解（QLE/矩方法/核密度）估计无效率项 u，再取指数变换得到 TE。

---

## 第二部分：CNLSDDFweak — 方向距离函数

CNLSDDFweak 不需要指定 `cet` 参数（内部使用加法模型），求解器为 **mosek**（线性规划）。

### 模型 4：CRS

```python
model = CNLSSDFDDFweak.CNLSDDFweak(
    data, sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    fun=FUN_PROD, rts=RTS_CRS
)
model.optimize(solver="mosek")

rd = StoNED.StoNED(model)
res = rd.get_technical_efficiency_ratio(RED_QLE)
```

### 模型 5：VRS1

```python
model = CNLSSDFDDFweak.CNLSDDFweak(
    data, sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    fun=FUN_PROD, rts=RTS_VRS1
)
model.optimize(solver="mosek")

rd = StoNED.StoNED(model)
res = rd.get_technical_efficiency_ratio(RED_QLE)
```

### 模型 6：VRS2

```python
model = CNLSSDFDDFweak.CNLSDDFweak(
    data, sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    fun=FUN_PROD, rts=RTS_VRS2
)
model.optimize(solver="mosek")

rd = StoNED.StoNED(model)
res = rd.get_technical_efficiency_ratio(RED_QLE)
```

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `gy` | `[0]` | 不调整期望产出 |
| `gx` | `[0,0,0]` | 不调整投入 |
| `gb` | `[1]` | 沿非期望产出方向调整 |
| `fun` | `FUN_PROD` | 生产前沿 |
| `rts` | `RTS_CRS`/`RTS_VRS1`/`RTS_VRS2` | 规模报酬 |
| `solver` | `"mosek"` | 线性求解器 |

### 返回值

- `get_technical_efficiency_ratio(RED_QLE)` → numpy array：每个 DMU 的效率比率 TE

### TE ratio 计算方式（按导向）

```
gb≥1（非期望产出导向）：
    TE = (-basexy + epsilon) / (-basexy + epsilon + ddfhat)

gy≥1（产出导向）：
    TE = basexy / (basexy + ddfhat)

gx≥1（投入导向）：
    TE = (-basexy + epsilon) / (-basexy + epsilon + ddfhat)
```

---

## CNLSSDweak vs CNLSDDFweak 区别

| 特征 | CNLSSDweak | CNLSDDFweak |
|------|-----------|-------------|
| 类型 | Shephard 距离函数 | 方向距离函数 |
| cet 参数 | **必须** `CET_MULT` | 无（内部加法） |
| 求解器 | `"knitro"`（非线性） | `"mosek"`（线性） |
| 目标函数 | 最小化残差平方和（log 域） | 最小化残差平方和 |
| StoNED 方法 | `get_technical_efficiency()` | `get_technical_efficiency_ratio()` |
| 返回值 | TE = exp(-u) | TE = ratio 公式 |

---

## StoNED 完整方法列表

```python
rd = StoNED.StoNED(model)
```

| 方法 | 返回值 | 说明 |
|------|--------|------|
| `get_technical_efficiency(method)` | numpy array | 技术效率（CNLSSDweak 用） |
| `get_technical_efficiency_ratio(method)` | numpy array | 效率比率（CNLSDDFweak 用） |
| `get_technical_inefficiency(method)` | numpy array | 技术无效率值 |
| `get_mean_of_inefficiency(method)` | float | 无效率均值 μ |

### method 参数

| 值 | 全称 | 含义 |
|----|------|------|
| `RED_MOM` | Method of Moments | 矩方法（默认） |
| `RED_QLE` | Quasi-Likelihood Estimation | 准似然估计 |
| `RED_KDE` | Kernel Density Estimation | 核密度估计（含 Richardson-Lucy 盲反卷积） |

---

## CNLSSDweak 完整参数列表

```python
CNLSSDFDDFweak.CNLSSDweak(data, sent, z=None, gy=[1], gx=[0], gb=[0], cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | 必填 | 数据 |
| `sent` | str | 必填 | `"K L E=Y:CO2"` |
| `z` | list | `None` | 上下文变量 |
| `gy` | list | `[1]` | 期望产出方向向量 |
| `gx` | list | `[0]` | 投入方向向量 |
| `gb` | list | `[0]` | 非期望产出方向向量 |
| `cet` | str | `CET_MULT` | 误差类型（**必须 CET_MULT**） |
| `fun` | str | `FUN_PROD` | 前沿类型：`FUN_PROD`/`FUN_COST` |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_CRS`/`RTS_VRS1`/`RTS_VRS2` |

## CNLSDDFweak 完整参数列表

```python
CNLSSDFDDFweak.CNLSDDFweak(data, sent, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS1, deltap=None, gammap=None, kappap=None, epsilonp=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | 必填 | 数据 |
| `sent` | str | 必填 | `"K L E=Y:CO2"` |
| `z` | list | `None` | 上下文变量 |
| `gy` | list | `[1]` | 期望产出方向向量 |
| `gx` | list | `[1]` | 投入方向向量 |
| `gb` | list | `[1]` | 非期望产出方向向量 |
| `fun` | str | `FUN_PROD` | 前沿类型 |
| `rts` | str | `RTS_VRS1` | 规模报酬 |
| `deltap` | tuple | `None` | delta 变量 bounds `(min, max)` |
| `gammap` | tuple | `None` | gamma 变量 bounds `(min, max)` |
| `kappap` | tuple | `None` | kappa 变量 bounds `(min, max)` |
| `epsilonp` | tuple | `None` | epsilon 变量 bounds `(min, max)` |

---

## CNLSSDweak 方法列表

| 方法 | 返回值 | 说明 |
|------|--------|------|
| `optimize(email, solver)` | — | 求解模型 |
| `display_status()` | 打印 | 求解状态 |
| `display_alpha()` | 打印 | alpha 值（VRS） |
| `display_delta()` | 打印 | delta 值（投入系数） |
| `display_gamma()` | 打印 | gamma 值（产出系数） |
| `display_kappa()` | 打印 | kappa 值（非期望产出系数） |
| `display_residual()` | 打印 | 残差值 |
| `get_status()` | int | 求解状态 |
| `get_alpha()` | array | alpha 值 |
| `get_delta()` | DataFrame | delta 值 |
| `get_gamma()` | DataFrame | gamma 值 |
| `get_kappa()` | DataFrame | kappa 值 |
| `get_residual()` | array | 残差值 |
| `get_adjusted_residual()` | array | 调整后残差 |

## CNLSDDFweak 方法列表

与 CNLSSDweak 类似，额外支持 `deltap`/`gammap`/`kappap` bounds 参数。

---

## 求解器说明

| 求解器 | 适用场景 | 说明 |
|--------|---------|------|
| `"knitro"` | CNLSSDweak（非线性） | 商业非线性求解器，必须用于 CET_MULT 模型 |
| `"mosek"` | CNLSDDFweak（线性） | 商业线性求解器，速度快 |
