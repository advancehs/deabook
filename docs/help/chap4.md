# 第四章：Malmquist 生产率指数（动态产出效率）

**主题**：使用 Malmquist 指数分解生产率变动（MEFFCH + MTECHCH）。

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

### 本章需要的导入

```python
from deabook import MQDEA
from deabook.constant import RTS_VRS1, RTS_CRS, OPT_LOCAL, TOTAL, CONTEMPORARY
import pandas as pd
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表中的 `rts`、`tech`、`email` 行。

---

### 参数详解

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | — | 面板数据，需包含 sent 中指定的列、`id` 列和 `year` 列 |
| `id` | str | 数据中任意列名 | DMU 标识列名，用于分组计算跨期效率 |
| `year` | str | 数据中任意时间列名 | 时间列名，需可排序 |
| `sent` | str | `"投入1 投入2 ...=产出1 产出2 ..."` | `=` 左边为投入变量，右边为产出变量 |
| | | | |
| `gy` | list[int] | `[0]` | 不沿该产出方向调整 |
| | list[int] | `[1]` | 沿该产出方向扩大产出 |
| | | | |
| `gx` | list[int] | `[0]` | 不沿该投入方向调整 |
| | list[int] | `[1]` | 沿该投入方向缩减 |
| | | | |
| **gx/gy 组合** | **判断规则** | **`sum(gx)>=1`** | 至少缩减一个投入 |
| | | **`sum(gy)>=1`** | 至少扩大一个产出 |
| | | | |
| | **组合** | ✗ ✓ → 产出导向 | 不缩减投入，只扩大产出 |
| | | ✓ ✗ → 投入导向 | 缩减投入中 gx=1 对应的变量，产出不变 |
| | | ✓ ✓ → 超效率(Hyper) | 同时缩减投入、扩大产出 |
| | | | |
| | **示例**（3投入K,L,E + 1产出Y） | `gx=[0,0,0] gy=[1]` | 产出导向：扩大 Y，投入不变 |
| | | `gx=[1,0,0] gy=[0]` | 投入导向：只缩减 K |
| | | `gx=[0,1,0] gy=[0]` | 投入导向：只缩减 L |
| | | `gx=[1,1,0] gy=[0]` | 投入导向：同时缩减 K 和 L |
| | | `gx=[1,0,0] gy=[1]` | 超效率：缩减 K 同时扩大 Y |
| | | | |
| `rts` | str | `RTS_CRS` | 不变规模报酬 |
| | str | `RTS_VRS1` | 可变规模报酬 |
| | | | |
| `tech` | str | `CONTEMPORARY` | 当期技术——每期仅用当期 DMU 构成参考集，计算 D11、D12、D21 三个距离值 |
| | str | `TOTAL` | 全局技术——所有时期 DMU 构成统一参考集，只计算 D11 一个距离值，额外返回 `mqpi` 列（相邻期 D11 的比值） |
| | | | |
| `email` | str | `OPT_LOCAL` | 本地求解 |
| | str | 邮箱地址 | 远程求解完成后发送通知 |
| | | | |
| `solver` | str | `"mosek"` | 商业求解器，速度快、精度高 |
| | str | `"glpk"` | 开源求解器 |
| | str | `"cbc"` | 开源求解器 |

### 返回值

`optimize()` 返回一个 DataFrame，结构与 **tech 参数和 gx/gy 导向**有关。

#### 情况一：`tech=CONTEMPORARY`，投入/产出导向（非 Hyper）

| 返回列 | 含义 |
|--------|------|
| `D11` | 当期技术下，t 期 DMU 对 t 期参考集的效率距离值 |
| `D12` | 当期技术下，t 期 DMU 对 t-1 期参考集的效率距离值（从第二年起有值） |
| `D21` | 当期技术下，t 期 DMU 对 t+1 期参考集的效率距离值（到倒数第二年起有值） |
| `MQ` | Malmquist 指数，> 1 表示生产率提升，= 1 不变，< 1 下降 |
| `MEFFCH` | 效率变化（技术效率追赶效应），= D11(t) / D11(t-1) |
| `MTECHCH` | 技术变化（前沿移动效应） |

> **分解关系**：`MQ = MEFFCH × MTECHCH`
>
> **第一年无 MQ 值**：MQ 需要相邻两期对比，因此第一年的 `MQ`、`MEFFCH`、`MTECHCH` 均为 NaN。

#### 情况二：`tech=CONTEMPORARY`，超效率导向（Hyper）

Hyper + CRS 时返回结构与情况一相同（单列 `D11`/`MQ` 等，底层 `te = sqrt(rho)`）。

Hyper + VRS 时底层的 DEA2 会返回 `tei`（投入端）和 `teo`（产出端）两个效率值，因此 Malmquist 结果列拆分为投入/产出方向：

| 返回列 | 含义 |
|--------|------|
| `D11_tei` | 当期技术下，t 期对 t 期参考集的投入方向距离值 |
| `D11_teo` | 当期技术下，t 期对 t 期参考集的产出方向距离值 |
| `D12_tei` / `D12_teo` | t 期 DMU 对 t-1 期参考集的投入/产出距离值 |
| `D21_tei` / `D21_teo` | t 期 DMU 对 t+1 期参考集的投入/产出距离值 |
| `MQ_tei` | 投入方向的 Malmquist 指数 |
| `MQ_teo` | 产出方向的 Malmquist 指数 |
| `MEFFCH_tei` / `MEFFCH_teo` | 投入/产出方向的效率变化 |
| `MTECHCH_tei` / `MTECHCH_teo` | 投入/产出方向的技术变化 |

> **分解关系**：`MQ_tei = MEFFCH_tei × MTECHCH_tei`，`MQ_teo` 同理。

#### 情况三：`tech=TOTAL`，投入/产出导向（非 Hyper）

| 返回列 | 含义 |
|--------|------|
| `D11` | 全局技术下，t 期 DMU 对全时期统一参考集的效率距离值 |
| `mqpi` | Malmquist 生产率指数，= D11(t) / D11(t-1)，即相邻期全局效率的比值 |

> **注意**：全局技术下不区分 MEFFCH 和 MTECHCH，只返回比值 `mqpi`。

#### 情况四：`tech=TOTAL`，超效率导向（Hyper）

与情况三同理，Hyper + CRS 返回结构与情况三相同（单列）。Hyper + VRS 时拆分为投入/产出方向：

| 返回列 | 含义 |
|--------|------|
| `D11_tei` | 全局技术下，投入方向距离值 |
| `D11_teo` | 全局技术下，产出方向距离值 |
| `mqpi_tei` | 投入方向的 Malmquist 指数，= D11_tei(t) / D11_tei(t-1) |
| `mqpi_teo` | 产出方向的 Malmquist 指数，= D11_teo(t) / D11_teo(t-1) |

#### MQDEA vs MQDDF 返回值差异

| 特征 | MQDEA | MQDDF |
|------|-------|-------|
| 底层模型 | DEA2（径向） | DDF2（方向距离函数） |
| 效率列名 | 投入导向用 `te` → `D11` | 投入导向用 `tei` → `D11` |
| | 产出导向用 `te` → `D11` | 产出导向用 `teo` → `D11` |
| 分解公式 | 相同 | 相同 |

> 返回列名完全一致，差异仅在距离值的底层计算方式（径向比率 vs DDF 距离）。

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

MQDDF 使用 DDF 的 `teo` 效率值计算 Malmquist 指数，返回列结构相同。

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
---
## 模型 4：MQDDF 投入导向 — VRS + 当期技术

```python
model = MQDEA.MQDDF(
    data, id=dmuname, year='year',
    sent="K L E=Y",
    gx=[1,0,0], gy=[0],
    rts=RTS_VRS1,
    tech=CONTEMPORARY,
    email=OPT_LOCAL, solver="mosek"
)
res = model.optimize()
```
---
