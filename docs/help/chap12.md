# 第十二章：SBM、DDF、NDDF 与 Malmquist 指数从零实现

**主题**：使用 Pyomo 从零实现 SBM（基于松弛的测度）、方向距离函数（DDF）、非径向方向距离函数（NDDF）、Malmquist 指数（MPI）和 Malmquist-Luenberger 指数（MLPI）。

**数据文件**：
- `Ex4.dta`（30个中国省份，3期，含 CO2）

**变量含义**：
- `K`：资本, `L`：劳动, `E`：能源, `Y`：产出(GDP), `CO2`：碳排放（非期望产出）

---

## 安装与导入

本章所有模型使用 Pyomo 从零构建，不依赖 deabook 包。

```bash
pip install pyomo mosek
```

```python
from pyomo.environ import *
import pandas as pd
import numpy as np
import re
import warnings
warnings.filterwarnings("ignore")
```

---

## formula 语法

| 格式 | 示例 | 适用函数 |
|------|------|---------|
| 无非期望产出 | `"Y=K L"` | dea8 |
| 含非期望产出 | `"Y:CO2=K L"` | sbm2, ddf, nddf, sbmdual3, nddfdual2 |

`"Y:CO2=K L"` 解析：`=` 右边 = 投入，`=` 左边 `:` 前 = 期望产出，`:` 后 = 非期望产出。

---

### 参数详解

#### SBM 系列

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `formula` | str | `"Y=K L"` | 无非期望产出 |
| | str | `"Y:CO2=K L"` | 含非期望产出 |
| `dataframe` | DataFrame | — | DMU 数据 |
| `evaquery` | str | `"t==[1,2,3]"` | 评价 DMU 筛选 |
| `refquery` | str | `"t==[1,2,3]"` | 参考集 DMU 筛选 |

> SBM 求解器固定 `"mosek"`。

#### DDF

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `formula` | str | `"Y:CO2=K L"` | 含非期望产出的公式 |
| `gx` | list | `[-1]*n` | 投入缩减方向（负值 = 缩减） |
| `gy` | list | `[1]*m` | 期望产出扩张方向（正值 = 扩张） |
| `gb` | list | `[-1]*p` | 非期望产出缩减方向（负值 = 缩减） |
| `evaquery` | str | `None` | 评价 DMU 筛选 |
| `refquery` | str | `None` | 参考集 DMU 筛选 |

> DDF 非期望产出约束为等式（弱可处置性）。

#### NDDF

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `formula` | str | `"Y:CO2=K L"` | 含非期望产出的公式 |
| `gx` | list | `[-1]*n` | 投入方向 |
| `gy` | list | `[1]*m` | 产出方向 |
| `gb` | list | `[-1]*p` | 非期望产出方向 |
| `weight` | list | `None` | 权重矩阵，默认均等分配 |
| `evaquery` | str | `None` | 评价 DMU 筛选 |
| `refquery` | str | `None` | 参考集 DMU 筛选 |

> 默认权重：每个活跃方向的权重 = `1 / (活跃方向数 × 该方向变量数)`。

#### MPI / MLPI 系列

| 参数 | 类型 | 可选值 | 含义 |
|------|------|--------|------|
| `formula` | str | `"Y=K L"` | MPI 用无非期望产出 |
| | str | `"Y:CO2=K L"` | MLPI 用含非期望产出 |
| `data` | DataFrame | — | 面板数据，需含 id 和 t 列 |
| `id` | str | `"id"` | 个体变量名 |
| `t` | str | `"t"` | 时间变量名 |
| `tech` | str | `None` / `"com"` | 当期技术：每期仅用当期 DMU |
| | str | `"seq"` | 时序技术：累积到当期的所有 DMU |
| | str | `"window N"` | 视窗技术：前后各 N 期窗口 |
| | str | `"global"` | 全局技术：所有时期 DMU（仅 mpi3） |

---

### 返回值

#### SBM

| 返回列 | 含义 |
|--------|------|
| `obj` | SBM 目标函数值（效率值） |
| `slack K`, `slack L` 等 | 投入松弛 |
| `slack Y` 等 | 产出松弛 |
| `slack CO2` 等 | 非期望产出松弛（仅 sbm2） |

#### DDF

| 返回列 | 含义 |
|--------|------|
| `te` | 方向距离函数值（越大表示越无效率） |

#### NDDF

| 返回列 | 含义 |
|--------|------|
| `obj` | 目标函数值 |
| `scale K`, `scale L` 等 | 投入端缩放因子 |
| `scale Y` 等 | 产出端缩放因子 |
| `scale CO2` 等 | 非期望产出端缩放因子 |

#### MPI / MLPI

| 返回列 | 含义 |
|--------|------|
| `mpi` / `mlpi` | Malmquist 指数，> 1 生产率提升，< 1 下降 |

> MPI 额外在原始数据基础上增加列。

---

## DDF vs NDDF 区别

| 特征 | DDF | NDDF |
|------|-----|------|
| 缩放因子 | 单一 theta | 逐变量 thetax, thetay, thetab |
| 目标函数 | maximize theta | maximize 加权和 |
| 权重 | 无 | weight 参数 |

---

## MPI vs MLPI 区别

| 特征 | MPI | MLPI |
|------|-----|------|
| 底层模型 | dea8（径向 DEA） | ddf（方向距离函数） |
| 非期望产出 | 不支持 | 支持 |
| 公式 | `D12/D11 * D11/D21` | `(1+D11)/(1+D12) * (1+D21)/(1+D22)` |

---

## 模型 1：sbm — 无非期望产出

```python
ex4 = pd.read_stata("Ex4.dta")
sbmte = sbm("Y=K L", ex4, "t==[1,2,3]", "t==[1,2,3]")
```

---

## 模型 2：sbm2 — 含非期望产出

```python
sbmte = sbm2("Y:CO2=K L", ex4, "t==[1,2,3]", "t==[1,2,3]")
```

---

## 模型 3：ddf — 方向距离函数

```python
ddfte = ddf("Y:CO2=K L", ex4)
ddfte = ddf("Y:CO2=K L", ex4, gx=[-1,-1], gy=[1], gb=[-1])
```

---

## 模型 4：nddf — 非径向方向距离函数

```python
nddfte = nddf("Y:CO2=K L", ex4)
```

---

## 模型 5：mpi — Malmquist 指数

```python
data2 = mpi("Y=K L", ex4, "id", "t")
```

> 底层用 dea8 计算四个距离值（D11, D12, D21, D22）后组合。

---

## 模型 6：mpi2 — MPI + 生产技术选择

```python
data2 = mpi2("Y=K L", ex4, "id", "t", tech="seq")
```

> 支持 `tech="seq"`（时序技术）和 `tech="window N"`（视窗技术）。

---

## 模型 7：mpi3 — MPI + 全局技术

```python
data2 = mpi3("Y=K L", ex4, "id", "t", tech="global")
```

> 全局技术仅需 D11，`mpi = D11[t-1] / D11[t]`。

---

## 模型 8：mlpi — Malmquist-Luenberger 指数

```python
data = mlpi("Y:CO2=K L", ex4, "id", "t", tech="com")
```

> 底层用 ddf 计算距离值，支持非期望产出。

---

## 函数速查表

| 函数 | formula 格式 | 特点 |
|------|-------------|------|
| `sbm` | `"Y=K L"` | SBM，无非期望产出 |
| `sbm2` | `"Y:CO2=K L"` | SBM，含非期望产出 |
| `ddf` | `"Y:CO2=K L"` | 方向距离函数 |
| `nddf` | `"Y:CO2=K L"` | 非径向 DDF |
| `dea8` | `"Y=K L"` | 径向 DEA（MPI 底层） |
| `mpi` | `"Y=K L"` | Malmquist 指数 |
| `mpi2` | `"Y=K L"` | MPI + 生产技术选择 |
| `mpi3` | `"Y=K L"` | MPI + 全局技术 |
| `mlpi` | `"Y:CO2=K L"` | Malmquist-Luenberger 指数 |
