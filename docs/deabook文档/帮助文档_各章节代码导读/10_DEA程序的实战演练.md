# chap10：从零手写 SBM、DDF、NDDF、MPI、MLPI

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap10`
- 主要代码文件：
  - `chap10.ipynb`

## 2. 本章主要模型 / 函数

- 手写 `sbm`、`sbm2`、`ddf`、`nddf`、`mpi/mpi2/mpi3`、`mlpi`

## 3. 本章使用的数据

`Ex4.dta`，30 个个体、3 个时期 `t=1,2,3`。

## 4. 本章 import 的逐项含义

- `pyomo`：直接写优化模型。
- `re`：解析 `"Y:CO2~K L"` 一类公式。
- `mosek`：作为 notebook 主求解器。

## 5. 本章参数穷尽说明

### `sbm(formula, dataframe, evaquery, refquery)`
- `formula`：`"产出~投入"`
- `dataframe`
- `evaquery`
- `refquery`

### `sbm2(formula, dataframe, evaquery, refquery)`
- `formula`：`"产出:非期望产出~投入"`

### `ddf(formula, dataframe, gx=None, gy=None, gb=None, evaquery=None, refquery=None)`

| 参数 | 穷尽说明 |
|---|---|
| `gx` | `None` 或长度=投入数的列表；`None` 时默认全 `-1` |
| `gy` | `None` 或长度=期望产出数的列表；`None` 时默认全 `1` |
| `gb` | `None` 或长度=非期望产出数的列表；`None` 时默认全 `-1` |

### `nddf(formula, dataframe, gx=None, gy=None, gb=None, weight=None, evaquery=None, refquery=None)`

新增 `weight`：

- `None`：自动平均分配权重
- 列表：长度必须等于 `len(x)+len(y)+len(b)`

### `mpi3(formula, data, id, t, tech=None)` 与 `mlpi(formula, data, id, t, tech=None)`

#### `tech` 的完整枚举

- `None`
- `"com"`
- `"seq"`
- `"window h"`，例如 `"window 2"`
- `"global"`

## 6. 本章输出文件

- 本章主要输出是手写函数返回的 DataFrame。
- 核心中间列包括 `D11`, `D12`, `D21`, `D22`，以及由此得到的动态指数。

## 7. 阅读与运行提醒

- 这章的“穷尽重点”是：`formula`、`gx/gy/gb`、`weight`、`tech` 四组参数。
- `tech` 是最容易忽视但最关键的动态参数。
