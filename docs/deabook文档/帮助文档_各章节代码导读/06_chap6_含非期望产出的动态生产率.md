# chap6：含非期望产出的动态生产率

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap6`
- 主要代码文件：
  - `chap 6.ipynb`
  - `chap 6_province.ipynb`

## 2. 本章主要模型 / 函数

- MQDEAweak.MQDEAweak、MQDEAweak.MQDDFweak

## 3. 本章使用的数据

与 chap5 同源，但进入动态生产率框架。

## 4. 本章 import 的逐项含义

- `from deabook import MQDEAweak`：动态 + 非期望产出模块。
- `RTS_CRS`、`RTS_VRS1`、`RTS_VRS2`：本章都实际使用。
- `TOTAL`、`CONTEMPORARY`：技术集常量；示例主用 `CONTEMPORARY`。
- `MAL`、`LUE`：动态指数类型常量；示例主用 `MAL`。

## 5. 本章参数穷尽说明

### `MQDEAweak.MQDEAweak(data, id, year, sent="inputvar=outputvar:unoutputvar", gy=[1], gx=[0], gb=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)`

参数集合是在 `MQDEA` 上加了：

- `gb`
- `sent` 支持 `:非期望产出`
- `rts` 多了 `RTS_VRS2`

#### 方向分支穷尽

| 条件 | 结果列族 |
|---|---|
| `gx` 非零 | `MQ`, `MEFFCH`, `MTECHCH` |
| `gy` 非零 | `MQ`, `MEFFCH`, `MTECHCH` |
| `gb` 非零 | `MQ`, `MEFFCH`, `MTECHCH` |
| `gx+gy` | 可拆成 `_tei/_teo` 结果 |
| `gy+gb` | 可拆成 `_teuo/_teo` 结果 |

#### 本章实际取值

- `sent="K L E=Y:CO2"`
- `gy=[0], gx=[0,0,0], gb=[1]`
- `rts=RTS_CRS / RTS_VRS1 / RTS_VRS2`
- `tech=CONTEMPORARY`

### `MQDEAweak.MQDDFweak(..., dynamic=MAL, ...)`

比 `MQDEAweak` 多一个关键参数：

| 参数 | 合法值 |
|---|---|
| `dynamic` | `MAL`、`LUE` |

- `dynamic=MAL`：走 Malmquist / ML 型分支
- `dynamic=LUE`：走 Luenberger 分支
- 其他值：源码明确报错

本章实际使用：`dynamic=MAL`。

## 6. 本章输出文件

- `table6.2_`–`table6.4_`：`MQDEAweak`
- `table6.5_`–`table6.7_`：`MQDDFweak`
- 图：`chap6分...png`

## 7. 阅读与运行提醒

- 这章本质上是在问：当污染进入模型后，动态生产率分解是否仍然稳定输出 `MQ / MEFFCH / MTECHCH`。
- 示例是“污染方向”版本，因此 `gb=[1]` 是关键。
