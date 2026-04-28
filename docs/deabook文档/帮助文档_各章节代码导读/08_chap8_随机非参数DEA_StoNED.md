# chap8：随机非参数 DEA / StoNED

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap8`
- 主要代码文件：
  - `chap 8.ipynb`

## 2. 本章主要模型 / 函数

- CNLSSDFDDFweak.CNLSSDweak、CNLSSDFDDFweak.CNLSDDFweak、StoNED.StoNED

## 3. 本章使用的数据

仍使用 `china data.xlsx` 与 `oecd data.xlsx`。

## 4. 本章 import 的逐项含义

- `from deabook import CNLSSDFDDFweak, StoNED`：随机非参数前沿 + 残差分解。
- `FUN_PROD`：生产前沿。
- `CET_MULT`：`CNLSSDweak` 的关键误差结构。
- `RED_QLE`：StoNED 残差分解方法，notebook 主结果实际使用它。

## 5. 本章参数穷尽说明

### `CNLSSDFDDFweak.CNLSSDweak(data, sent, z=None, gy=[1], gx=[0], gb=[0], cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1)`

| 参数 | 穷尽说明 |
|---|---|
| `z` | `None` 或情境变量矩阵。 |
| `cet` | 在这个类里实际强制要求 `CET_MULT`，否则直接报错。 |
| `fun` | `FUN_PROD` 或 `FUN_COST`。 |
| `rts` | `RTS_VRS1` 或 `RTS_VRS2`。 |

#### `optimize(email=OPT_LOCAL, solver=OPT_DEFAULT)`

- `email`：本地 / NEOS
- `solver`：示例里实际用 `knitro`

### `CNLSSDFDDFweak.CNLSDDFweak(data, sent, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS1, deltap=None, gammap=None, kappap=None, epsilonp=None)`

新增参数：

| 参数 | 穷尽说明 |
|---|---|
| `deltap` | `None` 或 `(lower, upper)`，限制 `delta` 边界 |
| `gammap` | `None` 或 `(lower, upper)`，限制 `gamma` 边界 |
| `kappap` | `None` 或 `(lower, upper)`，限制 `kappa` 边界 |
| `epsilonp` | 签名里有，但源码当前未真正启用自定义边界逻辑 |

### `StoNED.StoNED(model)`

#### `model` 的合法输入类（源码识别）

- `CNLSSDweak`
- `CNLSDDFweak`
- `CNLSSD`
- `CNLSDDF`
- `CNLSSDweakmeta`
- `CNLSDDFweakmeta`

本章实际只用：`CNLSSDweak`、`CNLSDDFweak`。

#### `method` 参数的全部枚举

| 方法 | `method` 合法值 |
|---|---|
| `get_mean_of_inefficiency` | `RED_MOM`, `RED_QLE`, `RED_KDE` |
| `get_technical_inefficiency` | 实现重点是 `RED_QLE`, `RED_KDE` |
| `get_technical_efficiency_ratio` | 同上 |
| `get_technical_efficiency` | 同上 |

本章实际用到：`RED_QLE`。

## 6. 本章输出文件

- `table_CNLSSDb_CRS3.md`
- `table_CNLSSDb_VRS_SAME3.md`
- `table_CNLSSDb_VRS_DIFF3.md`
- `table8.2_CNLSDDFb_CRS3.md`
- `table8.3_CNLSDDFb_VRS_SAME3.md`
- `table8.4_CNLSDDFb_VRS_DIFF3.md`
- 图：`chap8分年份碳排放效率-省.png`、`chap8分省份碳排放效率.png`

## 7. 阅读与运行提醒

- 这章一定要分两步理解：先估前沿，再做 StoNED 残差分解。
- 如果只看 `optimize()` 而不看 `StoNED`，就只看到“前沿拟合”，看不到最终技术效率。
