# chap12：DEA 调试与模型检查

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap12`
- 主要代码文件：
  - `chap12.ipynb`

## 2. 本章主要模型 / 函数

- 调试版 `dea8` 与 `model.pprint()`

## 3. 本章使用的数据

`Ex3.dta`。

## 4. 本章 import 的逐项含义

- `pyomo`：看模型内部结构。
- `pandas`、`numpy`、`re`：与 chap12 类似。
- `glpk`：这章用它做调试求解器。

## 5. 本章参数穷尽说明

### `dea8(formula, dataframe, rts, orient, evaquery=None, refquery=None)`

参数全集与 chap12 最终版一致：

| 参数 | 合法情况 |
|---|---|
| `formula` | `"Y~K L"` 这类公式 |
| `dataframe` | 包含公式涉及变量的数据框 |
| `rts` | `"crs"`、`"vrs"` |
| `orient` | `"oo"`、`"io"` |
| `evaquery` | `None` 或合法 `query()` 表达式 |
| `refquery` | `None` 或合法 `query()` 表达式 |

本章不是增加新参数，而是增加新的**调试情形**：

- 改用 `glpk`
- 在约束里插入 `print(...)`
- 直接返回 `model`
- 用 `model.pprint()` 检查模型结构

## 6. 本章输出文件

- 主要输出是 `model.pprint()` 的 Pyomo 结构打印结果。
- 不是正式表图章，而是检查模型构造是否正确。

## 7. 阅读与运行提醒

- 这章里有明显试错代码，例如打印中间对象、故意改约束、只返回 `model`。
- 所以它更像“调试工作台”，而不是“最终结果生产脚本”。
