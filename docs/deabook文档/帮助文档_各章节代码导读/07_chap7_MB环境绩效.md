# chap7：MB（物质平衡）环境绩效

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap7`
- 主要代码文件：
  - `chap 7.ipynb`

## 2. 本章主要模型 / 函数

- MB.MB

## 3. 本章使用的数据

`power plant data.xlsx`，notebook 实际筛到 `year>2014`，即 2015–2017。

## 4. 本章 import 的逐项含义

- `from deabook import MB`：导入 `deabook/MB.py`。
- `RTS_CRS`、`RTS_VRS1`：控制规模报酬。
- `CategoricalDtype`：固定 `OE / TE / PIAE / NPIAE / NPOAE` 顺序。

## 5. 本章参数穷尽说明

### `MB.MB(data, sent="inputvar_np + inputvar_p =outputvar_np + outputvar_p:unoutputvar", sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], level=5, rts=RTS_VRS1, baseindex=None, refindex=None)`

| 参数 | 穷尽说明 |
|---|---|
| `data` | 必须能按 MB 语法拆出非污染投入、污染投入、期望产出、非期望产出。 |
| `sent` | MB 专用语法，本章是 `"K L+F=E:CO2"`。 |
| `sx` | 投入侧物质平衡系数矩阵。 |
| `sy` | 产出侧物质平衡系数矩阵。 |
| `level` | **完整枚举：1,2,3,4,5**。 |
| `rts` | `RTS_VRS1` 或 `RTS_CRS`。 |
| `baseindex/refindex` | 同前。 |

#### `level` 的全部分支

| `level` | 含义 |
|---|---|
| `1` | 只变量化 `b` |
| `2` | 再变量化污染投入 `x_p` |
| `3` | 再变量化非污染投入 `x_np` |
| `4` | 再变量化污染产出 `y_p`；若不存在则变量化 `y_np` |
| `5` | 再变量化 `y_np` |

#### `optimize(solver=OPT_DEFAULT, dmu="1")`

| 参数 | 含义 |
|---|---|
| `solver` | 求解器名，如 `mosek` |
| `dmu` | 调用 `info()` 时查看哪个 DMU 的详细信息 |

#### 本章实际取值

- `sent="K L+F=E:CO2"`
- `sx=[[0,0,sx]]`
- `sy=[[sy]]`
- `level=1,2,3,4`
- `rts=RTS_VRS1` 或 `RTS_CRS`

## 6. 本章输出文件

- `table7.1.md`：描述统计
- `MB_CRS_1.md`：`OE`
- `MB_CRS_23.md`：`TE`, `PIAE`
- `MB_CRS_45.md`：`NPIAE`, `NPOAE`

## 7. 阅读与运行提醒

- 这章最特殊的是 `sent` 里的 `+` 号，它不是普通分隔符，而是 MB 语法的一部分。
- `level` 是必须穷尽理解的参数，因为本章就是靠连续抬高 `level` 做分解。
