# chap2：经典 DEA、HYPER、DDF

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap2`
- 主要代码文件：
  - `chap 2_country.ipynb`
  - `chap 2_province.ipynb`

## 2. 本章主要模型 / 函数

- DEA.DEA2、DEA.DDF2

## 3. 本章使用的数据

国家样本来自 `oecd data.xlsx`，省级样本来自 `china data.xlsx`。notebook 实际又筛到 `year>2015`，即 2016–2019。

## 4. 本章 import 的逐项含义

- `from deabook import DEA`：导入 `deabook/DEA.py` 模块。
- `from deabook.constant import RTS_CRS, RTS_VRS1`：控制 CRS / VRS。
- `CET_ADDI`、`OPT_LOCAL`：本章导入了，但主计算并未真正使用。
- `pandas`：读 Excel、拼接结果、透视表、输出 markdown。
- `matplotlib`、`numpy`：画年份折线图、国家/省份柱状图。

## 5. 本章参数穷尽说明

### `DEA.DEA2(data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None)`

| 参数 | 穷尽说明 |
|---|---|
| `data` | 必须是 `DataFrame`，且包含 `sent` 里声明的全部变量。 |
| `sent` | 语法是 `输入=期望产出`，本章是 `"K L E=Y"`。 |
| `gy` | 长度 = 期望产出个数。|
| `gx` | 长度 = 投入个数。|
| `rts` | 只支持 `RTS_CRS`、`RTS_VRS1`。 |
| `baseindex` | `None` 或筛选字符串；`None` 表示全样本被评价。 |
| `refindex` | `None` 或筛选字符串；`None` 表示全样本作参照集。 |

#### `gx/gy` 的全部方向分支

| 条件 | 模型含义 | 返回列 |
|---|---|---|
| `sum(gx)>=1, sum(gy)==0` | 投入导向 | `te=rho` |
| `sum(gy)>=1, sum(gx)==0` | 产出导向 | `te=1/rho` |
| `sum(gx)>=1, sum(gy)>=1, rts=RTS_CRS` | 双向 CRS | `te=sqrt(rho)` |
| `sum(gx)>=1, sum(gy)>=1, rts=RTS_VRS1` | 双向 VRS | `tei=1-rho`, `teo=1/(1+rho)` |
| 两边都无方向 | 非法 | 报错 |

#### 本章实际跑到的 `DEA2` 参数组合

- `gy=[0], gx=[1,0,0], rts=RTS_CRS`
- `gy=[0], gx=[1,0,0], rts=RTS_VRS1`
- `gy=[1], gx=[1,0,0], rts=RTS_CRS`
- `gy=[1], gx=[1,0,0], rts=RTS_VRS1`

### `DEA.DDF2(data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None)`

参数集合与 `DEA2` 同构，但底层是方向距离函数。

#### 返回列的全部分支

| 条件 | 返回列 |
|---|---|
| `sum(gx)>=1, sum(gy)==0` | `tei=1-rho` |
| `sum(gy)>=1, sum(gx)==0` | `teo=1/(1+rho)` |
| `sum(gx)>=1, sum(gy)>=1` | 同时返回 `tei`, `teo` |

#### 本章实际跑到的 `DDF2` 参数组合

- `gy=[0], gx=[1,0,0], rts=RTS_CRS`
- `gy=[0], gx=[1,0,0], rts=RTS_VRS1`

## 6. 本章输出文件

- 表：`DEA_CRS.md`、`DEA_VRS.md`、`HYPER_CRS.md`、`HYPER_VRS.md`、`DDF_CRS.md`、`DDF_VRS.md`、`DDF_VRS_province.md`
- 图：`chap2分年份资本效率-国.png`、`chap2分年份资本效率-省.png`、`chap2分国家资本效率.png`、`chap2分省份资本效率.png`、`chap2分年份地区资本效率-省.png`

## 7. 阅读与运行提醒

- `kind` 只要切成 `"province"` 或 `"country"`，就会切换数据源。
- notebook 中有 `ssss`、`SSS` 单元，不能直接从头无脑运行。
- 这章的核心不是普通“资本效率”三个字，而是：同一份数据被分别放进 DEA、HYPER、DDF 三种结构里对比。
