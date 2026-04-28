# chap11：对偶问题与影子价格

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap11`
- 主要代码文件：
  - `chap11.ipynb`

## 2. 本章主要模型 / 函数

- 手写 `sbmdual*`、`nddfdual*`

## 3. 本章使用的数据

`Ex4.dta` 与 `Ex5.dta`。

## 4. 本章 import 的逐项含义

- `pyomo`：直接构造对偶规划。
- `pandas`：读 Stata 数据、整理影子价格表。
- `numpy`：提取对偶变量值。
- `re`：解析公式字符串。

## 5. 本章参数穷尽说明

### `sbmdual(dataframe, varname, P, Q)`

| 参数 | 含义 |
|---|---|
| `dataframe` | 数据框 |
| `varname` | 变量列表，如 `['K','L','Y','CO2']` |
| `P` | 投入个数 |
| `Q` | 期望产出个数 |

非期望产出个数由 `总列数 - P - Q` 自动推出。

### `sbmdual2(dataframe, varname, P, Q, evaquery, refquery)`
- 在 `sbmdual` 上新增 `evaquery/refquery`

### `sbmdual3(formula, dataframe, evaquery, refquery)`
- 用 `formula="产出:非期望产出~投入"` 替代 `varname + P + Q`

### `nddfdual(dataframe, varname, P, Q)`
- 自动设定 `gx=-1`, `gy=1`, `gb=-1`
- 自动生成权重 `weight`

### `nddfdual2(formula, dataframe, gx=None, gy=None, gb=None, weight=None, evaquery=None, refquery=None)`

| 参数 | 穷尽说明 |
|---|---|
| `gx/gy/gb` | `None` 或合法方向向量 |
| `weight` | `None` 或合法长度权重向量 |
| `evaquery/refquery` | `None` 或 `query()` 字符串 |

## 6. 本章输出文件

- 主输出是影子价格列：`shadow price K`、`shadow price L`、`shadow price Y`、`shadow price CO2` 等。
- notebook 重点是返回价格表，而不是出图。

## 7. 阅读与运行提醒

- 这章的参数虽然不多，但 `P/Q` 和 `formula` 两套接口必须都看清。
- 一套是“按列数切块”，一套是“按公式解析切块”。
