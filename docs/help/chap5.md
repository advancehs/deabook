# 第五章：能源效率与能源生产率测度方法

**主题**：使用 DDF 测量静态能源效率，使用 MQDEA 测量动态能源效率变动。

**数据文件**：
- `oecd data.xlsx`（33个OECD国家，2016-2019）
- `china data.xlsx`（30个中国省份，2016-2019）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP)

---

## deabook 包介绍与安装

deabook 是一个用于数据包络分析（DEA）的 Python 包，提供径向 DEA、方向距离函数（DDF）、Malmquist 指数、弱可处置性模型、物质平衡分解、CNLS 随机前沿等功能。

### 安装

```bash
pip install deabook
```

### 本章需要的导入（全部）

```python
from deabook import DEA, MQDEA
from deabook.constant import RTS_VRS1, RTS_CRS, OPT_LOCAL, TOTAL, CONTEMPORARY
import pandas as pd
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表。

---

### 参数详解

#### 第一部分：静态能源效率（DDF2）

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | 包含 sent 中指定列的数据，如K, L, E, Y   |
| `sent` | str | `"投入1 投入2 ...=产出1 产出2 ..."` | `=` 左边为投入变量，`=`右边为产出变量 |
| | | | |
| `gy` | list[int] | `[0]` | 不调整产出（投入导向） |
| | | | |
| `gx` | list[int] | `[0,0,1]` | 只沿能源方向缩减（能源效率） |
| | list[int] | `[1,0,0]` | 只沿资本方向缩减（资本效率） |
| | list[int] | `[0,1,0]` | 只沿劳动方向缩减 |
| | list[int] | `[1,1,1]` | 全要素效率 |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬 |
| | | | |
| `email` | str | `OPT_LOCAL` | 本地求解 |
| | str | 邮箱地址 | 远程求解通知 |
| | | | |
| `solver` | str | `"mosek"` | 商业求解器 |
| | str | `"glpk"` | 开源求解器 |
| | | | |
| `baseindex` | str | `None` | 被评价 DMU 筛选条件，`None` 表示全部 |
| `refindex` | str | `None` | 参考技术筛选条件，`None` 表示全部 |

#### 第二部分：动态能源效率（MQDEA）

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | 面板数据，需含 id 和 year 列 |
| `id` | str | 列名（如 `dmuname`） | DMU 个体标识列名 |
| `year` | str | `"year"` | 时间列名 |
| `sent` | str | `"K L E=Y"` | 投入/产出公式 |
| | | | |
| `gx` | list[int] | `[0,0,1]` | 沿能源方向缩减（能源 Malmquist） |
| | | | |
| `gy` | list[int] | `[0]` | 不调整产出（投入导向） |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬 |
| | | | |
| `tech` | str | `CONTEMPORARY` | 当期生产技术 |
| | str | `TOTAL` | 全局生产技术 |
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

#### 静态（DDF2）

| 返回列 | 含义 |
|--------|------|
| `optimization_status` | 求解状态 |
| `rho` | 距离函数值 |
| `tei` | 能源技术效率（`tei = 1 - rho`，范围 0-1） |

#### 动态（MQDEA）

| 返回列 | 含义 |
|--------|------|
| `D11` | 当期能源距离值 |
| `MQ` | Malmquist 指数 |
| `MEFFCH` | 效率变化分量 |
| `MTECHCH` | 技术变化分量 |

#### 静态 vs 动态对比

| 特征 | 静态（DDF2） | 动态（MQDEA） |
|------|-------------|--------------|
| 时间维度 | 单截面 | 面板跨期 |
| 输入要求 | 仅需 data + sent | 需 data + id + year |
| 输出指标 | tei（单期效率） | MQ, MEFFCH, MTECHCH |
| 效率含义 | 相对当期前沿 | 跨期生产率变动 |

---

## 第一部分：静态能源效率（DDF）

### 模型：DDF 投入导向（能源效率）— VRS

```python
model = DEA.DDF2(data, sent="K L E=Y", gy=[0], gx=[0,0,1], rts=RTS_VRS1)
res = model.optimize(solver="mosek")
```

---

## 第二部分：动态能源效率（MQDEA Malmquist）

### 模型：MQDEA 投入导向（能源）— VRS + 当期技术

```python
model = MQDEA.MQDEA(
    data, id=dmuname, year='year',
    sent="K L E=Y",
    gx=[0,0,1], gy=[0],
    rts=RTS_VRS1,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```
