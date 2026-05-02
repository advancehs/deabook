# 第十章：CNLS/StoNED 随机前沿分析（碳排放效率）

**主题**：使用 CNLS（凸非参数最小二乘）和 StoNED（随机半非参数估计）测量碳排放效率。

**数据文件**：
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
from deabook import CNLSSDFDDFweak, StoNED
from deabook.constant import FUN_PROD, OPT_LOCAL, RTS_VRS1, RTS_VRS2, RTS_CRS, CET_ADDI, CET_MULT, RED_QLE, RED_MOM, RED_KDE
import pandas as pd
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表。



### 参数详解

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | DMU 数据，需包含 sent 中指定的列 |
| `sent` | str | `"K L E=Y:CO2"` | `=` 左边：投入（K, L, E）<br>`=` 右边 `:` 前：期望产出（Y）<br>`:` 后：非期望产出（CO2） |
| | | | |
| `gy` | list[int] | `[0]` | 不调整期望产出 |
| | list[int] | `[1]` | 沿期望产出方向扩张 |
| | | | |
| `gx` | list[int] | `[0,0,0]` | 不调整投入 |
| | list[int] | `[1,0,0]` | 沿第一个投入方向缩减 |
| | | | |
| `gb` | list[int] | `[0]` | 不调整非期望产出 |
| | list[int] | `[1]` | 沿非期望产出方向缩减（碳排放导向） |
| | | | |
| `fun` | str | `FUN_PROD` | 生产前沿 |
| | str | `FUN_COST` | 成本前沿 |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬（相同缩减因子） |
| | str | `RTS_VRS2` | 可变规模报酬（不同缩减因子） |
| | | | |
| **CNLSSDweak 专属** | | | |
| `cet` | str | `CET_MULT` | 乘法复合误差项（**必须**，设为其他值会报错） |
| `solver` | str | `"knitro"` | 非线性求解器（CNLSSDweak 必须用 knitro） |
| | | | |
| **CNLSDDFweak 专属** | | | |
| `solver` | str | `"mosek"` | 线性求解器 |



### 返回值

CNLS 模型需要配合 StoNED 模块进行残差分解，才能获得效率值。

#### StoNED 工作流

```python
model.optimize(solver="...")   # 1. 求解 CNLS 模型
rd = StoNED.StoNED(model)      # 2. 创建 StoNED 对象
res = rd.get_technical_efficiency(method)  # 3. 获取效率
```

#### StoNED 方法

| 方法 | 返回值 | 适用模型 | 说明 |
|------|--------|---------|------|
| `get_technical_efficiency(method)` | numpy array | CNLSSDweak | 技术效率 TE = exp(-u) |
| `get_technical_efficiency_ratio(method)` | numpy array | CNLSDDFweak | 效率比率 TE（按导向公式） |
| `get_technical_inefficiency(method)` | numpy array | 两者 | 技术无效率值 u |
| `get_mean_of_inefficiency(method)` | float | 两者 | 无效率均值 μ |

#### method 参数

| 值 | 全称 | 含义 |
|----|------|------|
| `RED_QLE` | Quasi-Likelihood Estimation | 准似然估计（推荐） |
| `RED_MOM` | Method of Moments | 矩方法 |
| `RED_KDE` | Kernel Deconvolution Estimation | 核密度反卷积估计 |

#### TE 计算方式

| 模型 | 公式 | 说明 |
|------|------|------|
| CNLSSDweak | `TE = exp(-u)` | 指数变换 |
| CNLSDDFweak (gb≥1) | `TE = (-basexy + ε) / (-basexy + ε + ddfhat)` | 非期望产出导向 |
| CNLSDDFweak (gy≥1) | `TE = basexy / (basexy + ddfhat)` | 产出导向 |
| CNLSDDFweak (gx≥1) | `TE = (-basexy + ε) / (-basexy + ε + ddfhat)` | 投入导向 |

---

## CNLSSDweak vs CNLSDDFweak

| 特征 | CNLSSDweak | CNLSDDFweak |
|------|-----------|-------------|
| 类型 | Shephard 距离函数 | 方向距离函数 |
| cet 参数 | **必须** `CET_MULT` | 无（内部加法） |
| 求解器 | `"knitro"`（非线性） | `"mosek"`（线性） |
| StoNED 方法 | `get_technical_efficiency()` | `get_technical_efficiency_ratio()` |
| TE 计算 | `exp(-u)` | ratio 公式（按导向） |

---

## 模型 1：CNLSSDweak — CRS

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

---

## 模型 2：CNLSSDweak — VRS1

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

---

## 模型 3：CNLSSDweak — VRS2

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

---

## 模型 4：CNLSDDFweak — CRS

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

---

## 模型 5：CNLSDDFweak — VRS1

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

---

## 模型 6：CNLSDDFweak — VRS2

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
