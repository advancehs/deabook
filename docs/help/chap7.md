# 第七章：包含非期望产出的生产率

**主题**：使用弱可处置性 Malmquist 指数分解含非期望产出的动态环境效率变动。

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
from deabook import MQDEAweak
from deabook.constant import RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_LOCAL, TOTAL, CONTEMPORARY
import pandas as pd
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表。




---

### 参数详解

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | 面板数据，需含 id 和 year 列 |
| `id` | str | 列名（如 `dmuname`） | DMU 标识列名 |
| `year` | str | `"year"` | 时间列名 |
| `sent` | str | `"K L E=Y:CO2"` | `=` 左边：投入（K, L, E）<br>`=` 右边 `:` 前：期望产出（Y）<br>`:` 后：非期望产出（CO2） |
| | | | |
| `gy` | list[int] | `[0]` | 不调整期望产出 |
| | list[int] | `[1]` | 沿期望产出方向扩大 |
| | | | |
| `gx` | list[int] | `[0,0,0]` | 不调整投入 |
| | list[int] | `[1,0,0]` | 沿第一个投入缩减 |
| | | | |
| `gb` | list[int] | `[0]` | 不调整非期望产出 |
| | list[int] | `[1]` | 沿非期望产出方向缩减 |
| | | | |
| **方向导向** | **判定条件** | | |
| | `input_oriented` | gx≥1, gy\==0, gb==0 | 投入导向 |
| | `output_oriented` | gy≥1, gx\==0, gb==0 | 产出导向 |
| | `undesirable_oriented` | gb≥1, gx\==0, gy==0 | 非期望产出导向 |
| | `hyper_orientedyx` | gx≥1, gy≥1, gb==0 | 投入+产出超效率 |
| | `hyper_orientedyb` | gb≥1, gy≥1, gx==0 | 产出+非期望超效率 |
| | `hyper_orientedxb` | gx≥1, gb≥1, gy==0 | 投入+非期望超效率 |
| | `hyper_orientedyxb` | gx≥1, gb≥1, gy≥1 | 全方向超效率 |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬（相同缩减因子） |
| | str | `RTS_VRS2` | 可变规模报酬（不同缩减因子） |
| | | | |
| `tech` | str | `CONTEMPORARY` | 当期技术（每期仅用当期 DMU 做参考） |
| | str | `TOTAL` | 全局技术（所有时期 DMU 构成统一参考集） |
| | | | |
| `email` | str | `OPT_LOCAL` | 本地求解 |
| | str | 邮箱地址 | 远程求解通知 |
| | | | |
| `solver` | str | `"mosek"` | 商业求解器 |
| | str | `"glpk"` | 开源求解器 |

### 返回值

`optimize()` 返回 DataFrame：

#### 当期技术（CONTEMPORARY）

| 返回列 | 含义 |
|--------|------|
| `D11` | 当期效率距离值 |
| `MQ` | Malmquist 指数（>1 表示环境生产率提升） |
| `MEFFCH` | 效率变化（技术效率追赶效应） |
| `MTECHCH` | 技术变化（前沿移动效应） |

分解关系：`MQ = MEFFCH × MTECHCH`

#### 全局技术（TOTAL）

| 返回列 | 含义 |
|--------|------|
| `D11` | 全局效率距离值 |
| `mqpi` | 全局 Malmquist 指数 |

#### 返回列与方向/RTS 的关系

| 方向 + RTS | 底层返回 | D11 列 |
|-----------|---------|--------|
| 单一导向 + CRS | `te` | `D11 = te` |
| 单一导向 + VRS1/VRS2 | `te` | `D11 = te` |
| hyper_orientedyb + CRS | `te` | `D11 = te` |
| hyper_orientedyb + VRS1/VRS2 | `teuo`, `teo` | `D11_teuo`, `D11_teo` |
| hyper_orientedyx + VRS2 | `tei`, `teo` | `D11_tei`, `D11_teo` |
| hyper_orientedxb + CRS | `te` | `D11 = te` |
| hyper_orientedxb + VRS2 | `tei`, `teuo` | `D11_tei`, `D11_teuo` |
| hyper_orientedyxb + CRS | `te` | `D11 = te` |
| hyper_orientedyxb + VRS2 | `tei`, `teuo`, `teo` | `D11_tei`, `D11_teuo`, `D11_teo` |

---

## 模型 1：MQDEAweak 非期望产出导向 — CRS + 当期技术

```python
model = MQDEAweak.MQDEAweak(
    data, id=dmuname, year='year',
    sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    rts=RTS_CRS,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

---

## 模型 2：MQDEAweak 非期望产出导向 — VRS1 + 当期技术

```python
model = MQDEAweak.MQDEAweak(
    data, id=dmuname, year='year',
    sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    rts=RTS_VRS1,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

---

## 模型 3：MQDEAweak 非期望产出导向 — VRS2 + 当期技术

```python
model = MQDEAweak.MQDEAweak(
    data, id=dmuname, year='year',
    sent="K L E=Y:CO2",
    gy=[0], gx=[0,0,0], gb=[1],
    rts=RTS_VRS2,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

---

## MQDDFweak（方向距离函数版本）

**功能**：基于方向距离函数（DDF）的弱可处置性 Malmquist 生产率指数，支持 Luenberger 指标分解。

### 导入

MQDDFweak 与 MQDEAweak 位于同一模块：

```python
from deabook import MQDEAweak
# 调用方式：MQDEAweak.MQDDFweak(...)
```

### 额外常量导入

```python
from deabook.constant import MAL, LUE, OPT_DEFAULT
```

| 常量 | 含义 |
|------|------|
| `MAL` | Malmquist-Luenberger 分解方式 |
| `LUE` | Luenberger 指标分解方式 |
| `OPT_DEFAULT` | 默认优化选项 |

### 额外参数

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `dynamic` | str | `MAL` | 使用 Malmquist-Luenberger 分解（动态部分） |
| | | `LUE` | 使用 Luenberger 指标分解 |

其余参数（data, id, year, sent, gy, gx, gb, rts, tech, solver）与 MQDEAweak 相同。

### 代码示例

```python
from deabook import MQDEAweak
from deabook.constant import RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_DEFAULT, OPT_LOCAL, TOTAL, CONTEMPORARY, MAL, LUE
import pandas as pd

# MQDDFweak 方向距离函数版本
model = MQDEAweak.MQDDFweak(
    data=data, id=dmuname, year='year',
    sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1],
    rts=RTS_CRS, tech=CONTEMPORARY, dynamic=MAL,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
# 返回值列：year, DMU, MQ, MEFFCH, MTECHCH 等
```

### VRS1 示例

```python
model = MQDEAweak.MQDDFweak(
    data=data, id=dmuname, year='year',
    sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1],
    rts=RTS_VRS1, tech=CONTEMPORARY, dynamic=MAL,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

### VRS2 示例

```python
model = MQDEAweak.MQDDFweak(
    data=data, id=dmuname, year='year',
    sent="K L E=Y:CO2", gy=[0], gx=[0,0,0], gb=[1],
    rts=RTS_VRS2, tech=CONTEMPORARY, dynamic=MAL,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```
