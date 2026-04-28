# 第七章：物质平衡（Material Balance）污染物分解

**主题**：使用物质平衡约束的 DEA 模型分解污染物排放（OE、TE、PIAE、NPIAE、NPOAE）。

**数据文件**：
- `power plant data.xlsx`（30个省份，2015-2017）

**变量含义**：
- `K`：资本, `L`：劳动, `F`：燃料（污染投入）, `E`：电力（期望产出）, `CO2`：碳排放（非期望产出）

---

## 导入

```python
from deabook import MB
from deabook.constant import RTS_CRS, RTS_VRS1
import pandas as pd
```

**导入说明**：
- `MB`：物质平衡模块，包含 `MB`（基础类）
- `RTS_VRS1`：可变规模报酬
- `RTS_CRS`：不变规模报酬

---

## sent 公式格式

MB 模型使用特殊的公式格式，用 `+` 区分污染和非污染变量：

```
"K L+F=E:CO2"
```

- `=` 左边 `+` 前：非污染投入（K, L）
- `=` 左边 `+` 后：污染投入（F）
- `=` 右边 `:` 前：期望产出（E）
- `=` 右边 `:` 后：非期望产出（CO2）

### sent 变体

| sent | 含义 |
|------|------|
| `"K L+F=E:CO2"` | 非污染投入 K,L + 污染投入 F → 期望产出 E + 非期望产出 CO2 |
| `"K L+F=E+Y2:CO2"` | 含污染产出 Y2 的情况 |
| `"K L+F=E:CO2"` | 仅有期望产出 E 的情况 |

---

## sx / sy 污染物质系数

`sx` 和 `sy` 是污染物质系数矩阵，用于构建物质平衡约束：

```python
sx = data['CO2'].sum() / data['F'].sum()  # 平均排放系数
sy = 0
```

### 参数格式

| 参数 | 格式 | 含义 |
|------|------|------|
| `sx` | `[[0, 0, sx_val]]` | 投入端污染系数列表（per 非期望产出）。每个子列表对应一个非期望产出，长度等于总投入个数 |
| `sy` | `[[0]]` | 产出端污染系数列表（per 非期望产出）。每个子列表对应一个非期望产出，长度等于总产出个数 |

- `sx=[[0, 0, sx]]`：K 和 L 的污染系数为 0，F 的污染系数为 sx
- `sy=[[0]]`：期望产出 E 的污染系数为 0

---

## level 参数（求解层级）

`level` 控制模型逐层求解的深度：

| level | 变量化 | 说明 |
|-------|--------|------|
| `1` | `b`（非期望产出） | B1：技术效率对应的排放 |
| `2` | `+ x_p`（污染投入） | B2：加入污染投入调整后的排放 |
| `3` | `+ x_np`（非污染投入） | B3：加入非污染投入调整后的排放 |
| `4` | `+ y_p`/`y_np`（产出） | B4：加入产出调整后的排放（最终值） |
| `5` | `+ y_np`（非污染产出） | 完整分解（当存在非污染产出时使用） |

---

## 模型 1：MB — VRS1（Level 1-4）

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

---

## 参数详解

| 参数 | 值 | 含义 |
|------|----|------|
| `data` | DataFrame | 数据（含 K, L, F, E, CO2 列） |
| `sent` | `"K L+F=E:CO2"` | 非污染投入 + 污染投入 = 期望产出 : 非期望产出 |
| `sx` | `[[0, 0, sx]]` | 每个非期望产出对应的投入端污染系数 |
| `sy` | `[[0]]` | 每个非期望产出对应的产出端污染系数 |
| `rts` | `RTS_VRS1`/`RTS_CRS` | 规模报酬 |
| `baseindex` | `None`/`"year=[2017]"` | 被评价 DMU 筛选条件 |
| `refindex` | `None`/`"year=[2017]"` | 参考 DMU 筛选条件 |
| `level` | `1`/`2`/`3`/`4`/`5` | 求解层级 |

---

## 返回值

`optimize()` 返回元组 `(data3, info)`：

### data3（DataFrame）

| 列名 | 含义 |
|------|------|
| `optimization_status` | 求解状态 |
| `obj` | 目标函数值 |
| `best of Undesirable` | 最优非期望产出值（对应 level 的 B 值） |

### info（打印/str）

指定 DMU 的 LP 模型详细信息。

---

## 后处理：污染物分解指标

将 Level 1-4 的 B 值合并后计算分解指标：

```python
# 重命名 B 列
MB_CRSL1 = MB_CRSL1.rename(columns={'best of Undesirable': 'B1'})
MB_CRSL2 = MB_CRSL2.rename(columns={'best of Undesirable': 'B2'})
MB_CRSL3 = MB_CRSL3.rename(columns={'best of Undesirable': 'B3'})
MB_CRSL4 = MB_CRSL4.rename(columns={'best of Undesirable': 'B4'})

# 计算分解指标
MB_CRS['OE'] = abs(MB_CRS['B4']) / abs(MB_CRS['CO2'])
MB_CRS['TE'] = abs(MB_CRS['B1']) / abs(MB_CRS['CO2'])
MB_CRS['PIAE'] = abs(MB_CRS['B2']) / abs(MB_CRS['B1'])
MB_CRS['NPIAE'] = MB_CRS['B3'] / MB_CRS['B2']
MB_CRS['NPOAE'] = MB_CRS['B4'] / MB_CRS['B3']
```

### 分解指标含义

| 指标 | 公式 | 含义 |
|------|------|------|
| `OE` | `\|B4\| / \|CO2\|` | 总体排放效率（Overall Emission efficiency） |
| `TE` | `\|B1\| / \|CO2\|` | 技术效率（Technical Efficiency） |
| `PIAE` | `\|B2\| / \|B1\|` | 污染投入配置效率（Polluting Input Allocative Efficiency） |
| `NPIAE` | `B3 / B2` | 非污染投入配置效率（Non-Polluting Input Allocative Efficiency） |
| `NPOAE` | `B4 / B3` | 非污染产出配置效率（Non-Polluting Output Allocative Efficiency） |

### B 值层次关系

```
CO2 ──→ B1（技术效率）──→ B2（+污染投入）──→ B3（+非污染投入）──→ B4（+产出）
         TE                PIAE               NPIAE               NPOAE
```

---

## MB 完整参数列表

```python
MB.MB(data, sent="K L+F=E:CO2", sx=[[1,1,1]], sy=[[1]], level=5, rts=RTS_VRS1, baseindex=None, refindex=None)
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `data` | DataFrame | 必填 | 数据 |
| `sent` | str | 见默认 | 含 `+` 和 `:` 的公式 |
| `sx` | list | `[[1,1,1]]` | 投入端污染系数 |
| `sy` | list | `[[1]]` | 产出端污染系数 |
| `level` | int | `5` | 求解层级（1-5） |
| `rts` | str | `RTS_VRS1` | 规模报酬 |
| `baseindex` | str | `None` | 评价 DMU 筛选，如 `"year=[2017]"` |
| `refindex` | str | `None` | 参考 DMU 筛选，如 `"year=[2017]"` |

### 方法

| 方法 | 参数 | 返回值 | 说明 |
|------|------|--------|------|
| `optimize(solver, dmu)` | solver: 求解器; dmu: DMU 索引 | `(DataFrame, str)` | 求解并返回结果和 LP 信息 |
| `info(dmu)` | dmu: DMU 索引（默认 "all"） | 打印 | 显示 LP 模型详细信息 |

---

## baseindex / refindex 筛选语法

使用 DataFrame.query() 兼容的语法：

```python
baseindex = "year=[2017]"      # 只评价 2017 年的 DMU
refindex = "year=[2015,2016]"  # 使用 2015 和 2016 年的 DMU 做参考
```
