# 第六章：动态环境效率（含非期望产出的 Malmquist 指数）

**主题**：使用弱可处置性 Malmquist 指数分解含非期望产出的动态环境效率变动。

**数据文件**：
- `oecd data.xlsx`（33个OECD国家，2016-2019）
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP), `CO2`：碳排放（非期望产出）

---

## 导入

```python
from deabook import MQDEAweak
from deabook.constant import RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_LOCAL, TOTAL, CONTEMPORARY
import pandas as pd
```

**导入说明**：
- `MQDEAweak`：含非期望产出的 Malmquist 指数模块
- `MQDEAweak.MQDEAweak`：基于 DEAweak2 的 Malmquist 指数
- `CONTEMPORARY`：当期生产技术
- `TOTAL`：全局生产技术

---

## sent 公式格式

与弱可处置性模型相同，使用 `:` 分隔非期望产出：

```
"K L E=Y:CO2"
```

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

## 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `data` | DataFrame | 面板数据 |
| `id` | `"省份"` (dmuname) | DMU 标识列名 |
| `year` | `"year"` | 时间列名 |
| `sent` | `"K L E=Y:CO2"` | 含非期望产出的公式 |
| `gy` | `[0]` | 不调整期望产出 |
| `gx` | `[0,0,0]` | 不调整投入 |
| `gb` | `[1]` | 沿非期望产出方向调整 |
| `rts` | `RTS_CRS`/`RTS_VRS1`/`RTS_VRS2` | 规模报酬 |
| `tech` | `CONTEMPORARY`/`TOTAL` | 生产技术 |
| `email` | `OPT_LOCAL` | 本地求解 |
| `solver` | `"mosek"` | 求解器 |

### tech 参数说明

| 值 | 含义 | 参考集 |
|----|------|--------|
| `CONTEMPORARY` | 当期技术 | 每期仅用当期 DMU 做参考，计算 D11, D12, D21 |
| `TOTAL` | 全局技术 | 所有时期 DMU 构成统一参考集，仅计算 D11 |

---

## 方向导向判定

MQDEAweak 根据方向向量自动判定导向：

| 导向 | 条件 | 说明 |
|------|------|------|
| `input_oriented` | gx≥1, gy==0, gb==0 | 投入导向 |
| `output_oriented` | gy≥1, gx==0, gb==0 | 产出导向 |
| `undesirable_oriented` | gb≥1, gx==0, gy==0 | 非期望产出导向 |
| `hyper_orientedyx` | gx≥1, gy≥1, gb==0 | 投入+产出超效率 |
| `hyper_orientedyb` | gb≥1, gy≥1, gx==0 | 产出+非期望产出超效率 |

---

## 返回值

`res` 为 DataFrame，包含列：

**当期技术（CONTEMPORARY）**：
- `D11`：当期效率距离值
- `MQ`：Malmquist 指数（>1 表示环境生产率提升）
- `MEFFCH`：效率变化（技术效率追赶效应）
- `MTECHCH`：技术变化（前沿移动效应）
- **分解关系**：`MQ = MEFFCH × MTECHCH`

**全局技术（TOTAL）**：
- `D11`：全局效率距离值
- `mqpi`：全局 Malmquist 指数

### 返回列与方向/RTS 的关系

底层使用 DEAweak2 求解，返回列取决于方向和 RTS：

| 方向 + RTS | 底层返回 | 计算 D11 |
|-----------|---------|---------|
| 单一导向 + CRS | `te` | `D11 = te` |
| 单一导向 + VRS1/VRS2 | `te` | `D11 = te` |
| hyper_orientedyb + CRS | `te` | `D11 = te` |
| hyper_orientedyb + VRS1/VRS2 | `teuo`, `teo` | 分别保留 `D11_teuo`, `D11_teo` |
| hyper_orientedyx + VRS2 | `tei`, `teo` | 分别保留 `D11_tei`, `D11_teo` |

---

## MQDEAweak 完整参数列表

```python
MQDEAweak.MQDEAweak(data, id, year, sent, gy=[1], gx=[0], gb=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | 必填 | 面板数据 |
| `id` | str | 必填 | DMU 标识列名 |
| `year` | str | 必填 | 时间列名 |
| `sent` | str | `"K L E=Y:CO2"` | 含非期望产出的公式 |
| `gy` | list | `[1]` | 期望产出方向向量 |
| `gx` | list | `[0]` | 投入方向向量 |
| `gb` | list | `[0]` | 非期望产出方向向量 |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_VRS1`, `RTS_VRS2`, `RTS_CRS` |
| `tech` | str | `TOTAL` | 生产技术：`TOTAL` 或 `CONTEMPORARY` |
| `email` | str | `OPT_LOCAL` | 求解邮箱 |
| `solver` | str | `OPT_DEFAULT` | 求解器 |

### 方法

| 方法 | 返回值 | 说明 |
|------|--------|------|
| `optimize()` | DataFrame | 返回 Malmquist 指数及分解 |

注意：`optimize()` 无参数（求解在 `__init__` 中完成）。

---

## MQDEAweak vs MQDEA 区别

| 特征 | MQDEA（第三章） | MQDEAweak（本章） |
|------|----------------|-------------------|
| 底层模型 | DEA2/DDF2 | DEAweak2 |
| 非期望产出 | 不支持 | 支持（`:` 语法） |
| 方向向量 | gy, gx | gy, gx, **gb** |
| 导向类型 | 2 种 | 5 种（含非期望产出导向） |
