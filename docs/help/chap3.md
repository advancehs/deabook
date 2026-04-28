# 第三章：Malmquist 生产率指数（动态产出效率）

**主题**：使用 Malmquist 指数分解生产率变动（MEFFCH + MTECHCH）。

**数据文件**：
- `oecd data.xlsx`（33个OECD国家，2016-2019）
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP)

---

## 导入

```python
from deabook import MQDEA
from deabook.constant import RTS_VRS1, RTS_CRS, OPT_LOCAL, TOTAL, CONTEMPORARY
import pandas as pd
```

**导入说明**：
- `MQDEA`：Malmquist 指数模块，包含 `MQDEA`（基于 DEA2）和 `MQDDF`（基于 DDF2）
- `CONTEMPORARY`：当期生产技术（每期效率基于当期参考集）
- `TOTAL`：全局生产技术（所有 DMU 构成统一参考集）

---

## 模型 1：MQDEA 产出导向 — VRS + 当期技术

```python
model = MQDEA.MQDEA(
    data, id=dmuname, year='year',
    sent="K L E=Y",
    gx=[0,0,0], gy=[1],
    rts=RTS_VRS1,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

### 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `data` | DataFrame | 面板数据 |
| `id` | `"country"` (dmuname) | DMU 标识列名 |
| `year` | `"year"` | 时间列名 |
| `sent` | `"K L E=Y"` | 投入=产出 |
| `gx` | `[0,0,0]` | 投入方向全0 → 不缩减投入（产出导向） |
| `gy` | `[1]` | 产出方向为1 → 扩大产出 |
| `rts` | `RTS_VRS1` | 可变规模报酬 |
| `tech` | `CONTEMPORARY` | 当期生产技术 |
| `email` | `OPT_LOCAL` | 本地求解 |
| `solver` | `"mosek"` | 求解器 |

### 返回值

`res` 为 DataFrame，包含列：
- `D11`：当期效率距离值
- `MQ`：Malmquist 指数（>1 表示生产率提升）
- `MEFFCH`：效率变化（技术效率追赶效应）
- `MTECHCH`：技术变化（前沿移动效应）
- **分解关系**：`MQ = MEFFCH × MTECHCH`

### tech 参数说明

| 值 | 含义 | 参考集 |
|----|------|--------|
| `CONTEMPORARY` | 当期技术 | 每期仅用当期 DMU 做参考 |
| `TOTAL` | 全局技术 | 所有时期 DMU 构成统一参考集，额外返回 `mqpi` 列 |

### gx/gy 方向与产出/投入导向

| 导向 | gx | gy | 说明 |
|------|----|----|------|
| 产出导向 | `[0,0,0]` | `[1]` | 扩大产出 |
| 投入导向 | `[1,0,0]` | `[0]` | 缩减资本 |

---

## 模型 2：MQDDF 产出导向 — VRS + 当期技术

```python
model = MQDEA.MQDDF(
    data, id=dmuname, year='year',
    sent="K L E=Y",
    gx=[0,0,0], gy=[1],
    rts=RTS_VRS1,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

### MQDEA vs MQDDF 区别

| 特征 | MQDEA | MQDDF |
|------|-------|-------|
| 底层模型 | DEA2（径向） | DDF2（方向距离函数） |
| 距离函数 | 径向比率 | DDF 距离 |
| 效率指标 | 基于 DEA 的 te | 基于 DDF 的 tei/teo |

### 返回值

产出导向时：`D11`, `MQ`, `MEFFCH`, `MTECHCH`（使用 DDF 的 `teo` 效率值）

---

## 模型 3：MQDEA 投入导向 — VRS + 当期技术

```python
model = MQDEA.MQDEA(
    data, id=dmuname, year='year',
    sent="K L E=Y",
    gx=[1,0,0], gy=[0],
    rts=RTS_VRS1,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```

投入导向时使用 `tei` 效率值计算 Malmquist 指数。返回列结构相同。

---

## MQDEA 完整参数列表

```python
MQDEA.MQDEA(data, id, year, sent, gy=[1], gx=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | 必填 | 面板数据 |
| `id` | str | 必填 | DMU 标识列名 |
| `year` | str | 必填 | 时间列名 |
| `sent` | str | `"inputvar=outputvar"` | 投入产出公式 |
| `gy` | list | `[1]` | 产出方向向量 |
| `gx` | list | `[0]` | 投入方向向量 |
| `rts` | str | `RTS_VRS1` | 规模报酬：`RTS_VRS1`, `RTS_CRS` |
| `tech` | str | `TOTAL` | 生产技术：`TOTAL` 或 `CONTEMPORARY` |
| `email` | str | `OPT_LOCAL` | 求解邮箱 |
| `solver` | str | `OPT_DEFAULT` | 求解器 |

## MQDDF 完整参数列表

```python
MQDEA.MQDDF(data, id, year, sent, gy=[1], gx=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)
```

参数与 MQDEA 完全相同。

### 方法

| 方法 | 返回值 | 说明 |
|------|--------|------|
| `optimize()` | DataFrame | 返回 Malmquist 指数及分解 |

注意：MQDEA/MQDDF 的 `optimize()` 无参数（求解在 `__init__` 中完成）。
