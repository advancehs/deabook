# chap9：从零手写 DEA 程序

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap9`
- 主要代码文件：
  - `chap9.ipynb`

## 2. 本章主要模型 / 函数

- 纯 pyomo 手写 `dea1` 到 `dea8`

## 3. 本章使用的数据

`Ex3.dta`，20 个 DMU，变量是 `dmu, K, L, Y`。

## 4. 本章 import 的逐项含义

- `from pyomo.environ import *`：导入 `ConcreteModel`、`Set`、`Var`、`Objective`、`Constraint`、`SolverFactory` 等建模语言。
- `pandas`：读 `.dta`。
- `numpy`：提取解和拼接结果。
- `re`：从 `"Y~K L"` 公式字符串里拆变量。
- `warnings`：屏蔽提示。

## 5. 本章参数穷尽说明

这一章的重点就是参数如何一步一步扩展。

### `dea1(data, dataref, numk)`
- `data`：待评价样本
- `dataref`：参照样本
- `numk`：前 `numk` 列是投入

### `dea2(dataframe, varname, numk)`
- `dataframe`：数据框
- `varname`：变量名列表
- `numk`：投入变量个数

### `dea3(dataframe, varname, numk, evaquery)`
- 新增 `evaquery`

### `dea4(dataframe, varname, numk, evaquery, refquery)`
- 新增 `refquery`

### `dea5(dataframe, varname, numk, evaquery, refquery)`
- 参数集合与 `dea4` 相同，写法更规整

### `dea6(formula, dataframe, evaquery, refquery)`
- `formula`：`"Y~K L"`
- `dataframe`
- `evaquery`
- `refquery`

### `dea7(formula, dataframe, rts, evaquery, refquery)`
- `rts` 的完整枚举：`"crs"`、`"vrs"`

### `dea8(formula, dataframe, rts, orient, evaquery=None, refquery=None)`

| 参数 | 合法值 / 规则 |
|---|---|
| `formula` | `"产出~投入"` |
| `rts` | `"crs"`、`"vrs"` |
| `orient` | `"oo"`、`"io"` |
| `evaquery` | `None` 或 `query()` 字符串 |
| `refquery` | `None` 或 `query()` 字符串 |

## 6. 本章输出文件

- 本章不以导出大量图表为主，重点是逐步构造函数。
- 主要输出是 notebook 中的 DataFrame、模型对象、求解结果。

## 7. 阅读与运行提醒

- 如果你想理解 `deabook` 包为什么后来会有 `sent`、`rts`、`baseindex/refindex`，这一章就是底层原型。
- 这章不是“调包调用”，而是“把包要做的事自己先写一遍”。
