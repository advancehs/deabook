# 第九章：基于物质平衡的环境绩效测度方法

**主题**：使用物质平衡约束的 DEA 模型分解污染物排放（OE、TE、PIAE、NPIAE、NPOAE）。

**数据文件**：

- `power plant data.xlsx`（30个省份，2015-2017）

**变量含义**：

- `K`：资本, `L`：劳动, `F`：燃料（污染投入）, `E`：电力（期望产出）, `CO2`：碳排放（非期望产出）

---

## deabook 包介绍与安装

deabook 是一个用于数据包络分析（DEA）的 Python 包，提供径向 DEA、方向距离函数（DDF）、Malmquist 指数、弱可处置性模型、物质平衡分解、CNLS 随机前沿等功能。

### 安装

```bash
pip install deabook
```

### 本章需要的导入

```python
from deabook import MB
from deabook.constant import RTS_CRS, RTS_VRS1
import pandas as pd
```

> 各常量的含义与可选值见下方[参数详解](#参数详解)表。

---


### 参数详解

| 参数        | 类型       | 可选值             | 含义                                                         |
| ----------- | ---------- | ------------------ | ------------------------------------------------------------ |
| `data`      | DataFrame  | —                  | 包含 K, L, F, E, CO2 列的数据                                |
| `sent`      | str        | `"K L+F=E:CO2"`    | `=` 左边 `+` 前：非污染投入（K, L）<br>`=` 左边 `+` 后：污染投入（F）<br>`=` 右边 `:` 前：期望产出（E）<br>`:` 后：非期望产出（CO2） |
|             | str        | `"K L+F=E+Y2:CO2"` | 含污染产出 Y2 的情况                                         |
|             |            |                    |                                                              |
| `sx`        | list[list] | `[[0, 0, sx]]`     | 投入端污染系数（per 非期望产出），长度等于总投入个数<br>`0`：非污染投入（K, L）的系数<br>`sx`：污染投入（F）的系数，`sx = CO2.sum() / F.sum()` |
| `sy`        | list[list] | `[[0]]`            | 产出端污染系数（per 非期望产出），长度等于总产出个数<br>`0`：期望产出（E）的系数为 0，表示该期望产出（E）不含污染物质 |
|             |            |                    |                                                              |
| `rts`       | str        | `RTS_CRS`          | 不变规模报酬                                                 |
|             | str        | `RTS_VRS1`         | 可变规模报酬                                                 |
|             |            |                    |                                                              |
| `level`     | int        | `1`                | 仅调整 b（非期望产出）→ B1                                   |
|             | int        | `2`                | + x_p（污染投入）→ B2                                        |
|             | int        | `3`                | + x_np（非污染投入）→ B3                                     |
|             | int        | `4`                | + y（产出）→ B4                                              |
|             | int        | `5`                | + y_np（非污染产出），完整分解                               |
|             |            |                    |                                                              |
| `baseindex` | str        | `None`             | 被评价 DMU 筛选条件                                          |
|             | str        | `"year=[2017]"`    | 仅评价 2017 年 DMU                                           |
|             |            |                    |                                                              |
| `refindex`  | str        | `None`             | 参考 DMU 筛选条件                                            |
|             | str        | `"year=[2017]"`    | 仅用 2017 年 DMU 做参考                                      |

### 返回值

`optimize()` 返回**元组** `(data3, info)`（注意：与其他模块不同，返回元组而非单独 DataFrame）：

#### data3（DataFrame）

| 返回列                | 含义                                   |
| --------------------- | -------------------------------------- |
| `optimization_status` | 求解状态                               |
| `obj`                 | 目标函数值                             |
| `best of Undesirable` | 最优非期望产出值（对应 level 的 B 值） |

#### info（打印/str）

指定 DMU 的 LP 模型详细信息。

### 后处理：污染物分解指标

将 Level 1-4 的 B 值合并后计算分解指标。**注意**：本案例数据不含"不包含污染物质的期望产出"（Y_NP），因此 Level 5 的结果 $B_5 = B_4$，OE 和 NPOAE 直接用 $B_4$ 计算。若数据中存在非污染产出，则应运行 Level 5 并用 $B_5$ 计算 OE 和 NPOAE（见公式 9.20、9.25）。

```python
MB_CRS['OE'] = abs(MB_CRS['B4']) / abs(MB_CRS['CO2'])
MB_CRS['TE'] = abs(MB_CRS['B1']) / abs(MB_CRS['CO2'])
MB_CRS['PIAE'] = abs(MB_CRS['B2']) / abs(MB_CRS['B1'])
MB_CRS['NPIAE'] = MB_CRS['B3'] / MB_CRS['B2']
MB_CRS['NPOAE'] = MB_CRS['B4'] / MB_CRS['B3']
```

| 指标    | 公式（本案例） | 理论公式（含 Y_NP） | 含义               |
| ------- | -------------- | ------------------- | ------------------ |
| `OE`    | `|B4| / |CO2|` | `|B5| / |CO2|`      | 总体排放效率       |
| `TE`    | `|B1| / |CO2|` | `|B1| / |CO2|`      | 技术效率           |
| `PIAE`  | `|B2| / |B1|`  | `|B2| / |B1|`       | 污染投入配置效率   |
| `NPIAE` | `B3 / B2`      | `B3 / B2`           | 非污染投入配置效率 |
| `NPOAE` | `B4 / B3`      | `B5 / B4`           | 非污染产出配置效率 |

B 值层次关系：

```
CO2 → B1(技术效率) → B2(+污染投入) → B3(+非污染投入) → B4(+产出)
      TE              PIAE             NPIAE               NPOAE
```

---

## 模型 1：MB — Level 1-4

### Level 1

```python
model = MB.MB(data, sent="K L+F=E:CO2", sx=[[0,0,sx]], sy=[[sy]],
              rts=RTS_VRS1, baseindex=None, refindex=None, level=1)
res, info = model.optimize("mosek", "1")
```

### Level 2

```python
model = MB.MB(data, sent="K L+F=E:CO2", sx=[[0,0,sx]], sy=[[sy]],
              rts=RTS_CRS, baseindex=None, refindex=None, level=2)
res, info = model.optimize("mosek", "1")
```

### Level 3

```python
model = MB.MB(data, sent="K L+F=E:CO2", sx=[[0,0,sx]], sy=[[sy]],
              rts=RTS_CRS, baseindex=None, refindex=None, level=3)
res, info = model.optimize("mosek", "1")
```

### Level 4

```python
model = MB.MB(data, sent="K L+F=E:CO2", sx=[[0,0,sx]], sy=[[sy]],
              rts=RTS_CRS, baseindex=None, refindex=None, level=4)
res, info = model.optimize("mosek", "1")
```

| `info(dmu)` | dmu: DMU 索引（默认 "all"） | 打印 | 显示 LP 模型详细信息 |

---

## baseindex / refindex 筛选语法

使用 DataFrame.query() 兼容的语法：

```python
baseindex = "year=[2017]"      # 只评价 2017 年的 DMU
refindex = "year=[2015,2016]"  # 使用 2015 和 2016 年的 DMU 做参考
```
