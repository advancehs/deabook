# chap5：含非期望产出的静态环境效率

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap5`
- 主要代码文件：
  - `chap 5.ipynb`
  - `chap 5_province.ipynb`

## 2. 本章主要模型 / 函数

- DEAweak.DEAweak2、DEAweak.DDFweak2、DEAweak.NDDFweak2

## 3. 本章使用的数据

`china data.xlsx` 与 `oecd data.xlsx`，公式统一进入 `Y:CO2` 框架。

## 4. 本章 import 的逐项含义

- `from deabook import DEAweak`：弱可处置、含非期望产出的主模块。
- `RTS_CRS`、`RTS_VRS1`、`RTS_VRS2`：本章都实际跑到了。
- `matplotlib`、`numpy`：画环境效率年度图、地区图。

## 5. 本章参数穷尽说明

### `DEAweak.DEAweak2(data, sent, gy=[1], gx=[0], gb=[1], rts=RTS_VRS1, baseindex=None, refindex=None)`

| 参数 | 穷尽说明 |
|---|---|
| `sent` | 必须是 `输入=期望产出:非期望产出`，本章是 `"K L E=Y:CO2"`。 |
| `gy/gx/gb` | 长度分别匹配期望产出、投入、非期望产出个数。 |
| `rts` | `RTS_CRS`、`RTS_VRS1`、`RTS_VRS2`。 |

#### 全部分支

| 条件 | 返回列 |
|---|---|
| `gx` 非零 | `te` |
| `gy` 非零 | `te` |
| `gb` 非零 | `te` |
| `gx,gy` 非零且 `RTS_CRS` | `te` |
| `gx,gy` 非零且 `RTS_VRS2` | `tei`, `teo` |
| `gy,gb` 非零且 `RTS_VRS1/VRS2` | `teuo`, `teo` |
| `gx,gy` 非零且 `RTS_VRS1` | 源码明确不支持，报错 |

#### 本章实际参数组合

- `gy=[0], gx=[0,0,0], gb=[1]`
- `gy=[1], gx=[0,0,0], gb=[1]`
- `rts=RTS_CRS / RTS_VRS1 / RTS_VRS2`

### `DEAweak.DDFweak2(...)`

全方向分支：

| 条件 | 返回列 |
|---|---|
| 只输入方向 | `tei` |
| 只产出方向 | `teo` |
| 只非期望产出方向 | `teuo` |
| `gx+gy` | `tei`, `teo` |
| `gy+gb` | `teuo`, `teo` |
| `gx+gb` | `tei`, `teuo` |
| `gx+gy+gb` | `tei`, `teuo`, `teo` |

本章实际只跑：`gy=[1], gx=[0,0,0], gb=[1]`。

### `DEAweak.NDDFweak2(...)`

该类不只返回单一效率，而会先算：

- `rhox`
- `rhoy`
- `rhob`

再按权重聚合出：

- `tei`
- `teo`
- `teuo`
- `teuo2o`
- `tei2o`
- `teiuo2o`
- `teuo2i`

本章实际主用 `gy=[1], gx=[0,0,0], gb=[1]` 这一组，也就是污染和期望产出联合方向。

## 6. 本章输出文件

- `table5.2_`–`table5.4_`：`DEAweak2`
- `table5.5_`–`table5.7_`：HYPER 型环境效率
- `table5.8_`–`table5.10_`：`DDFweak2`
- `table5.11_`–`table5.13_`：`NDDFweak2`
- 图：环境效率按年份、地区、国家/省份的比较图

## 7. 阅读与运行提醒

- 本章最重要的不是某一个 `te`，而是 `RTS_VRS1` 和 `RTS_VRS2` 在污染框架下究竟触发了哪一类结果列。
- `sent="K L E=Y:CO2"`、`gb=[1]`、`rts` 三者一起决定模型结构。
