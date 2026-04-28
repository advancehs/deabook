# chap3：传统生产率指数（MQDEA / MQDDF）

## 1. 对应代码目录

- 目录：`D:\BaiduSyncdisk\software\deabook\deabook\example\chap3`
- 主要代码文件：
  - `chap 3.ipynb`
  - `chap 3_province.ipynb`

## 2. 本章主要模型 / 函数

- MQDEA.MQDEA、MQDEA.MQDDF

## 3. 本章使用的数据

继续使用 `china data.xlsx` 与 `oecd data.xlsx`；notebook 主体还是 `year>2015` 后的 2016–2019。

## 4. 本章 import 的逐项含义

- `from deabook import MQDEA`：导入动态生产率模块。
- `from deabook.constant import RTS_VRS1, TOTAL, CONTEMPORARY`：控制规模报酬和技术集。
- `CategoricalDtype`：保证 `MQ`、`MEFFCH`、`MTECHCH` 的展示顺序固定。
- `matplotlib`、`numpy`：画年份和国家/省份比较图。

## 5. 本章参数穷尽说明

### `MQDEA.MQDEA(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0], rts=RTS_VRS1, tech=TOTAL, email=OPT_LOCAL, solver=OPT_DEFAULT)`

| 参数 | 穷尽说明 |
|---|---|
| `data` | 面板型 `DataFrame`。 |
| `id` | 个体列名，例如 `国家`、`省份`。 |
| `year` | 时间列名，源码会检查是否存在。 |
| `sent` | 语法 `输入=期望产出`。 |
| `gy/gx` | 方向向量，长度必须分别匹配产出和投入数。 |
| `rts` | `RTS_CRS` 或 `RTS_VRS1`。 |
| `tech` | 只支持 `TOTAL`、`CONTEMPORARY`。 |
| `email` | 本地或 NEOS 设定。 |
| `solver` | 求解器名，如 `mosek`。 |

#### `tech` 的全部分支

- `TOTAL`：走 `get_total()`。
- `CONTEMPORARY`：走 `get_contemp()`。
- 其他值：直接报错。

#### `gx/gy` 的全部分支

| 条件 | 结果结构 |
|---|---|
| `sum(gx)>=1, sum(gy)==0` | 单列 `D11`，最终给 `MQ`, `MEFFCH`, `MTECHCH` |
| `sum(gy)>=1, sum(gx)==0` | 单列 `D11`，最终给 `MQ`, `MEFFCH`, `MTECHCH` |
| `sum(gx)>=1, sum(gy)>=1, rts=RTS_CRS` | 单列 `D11`，最终给 `MQ`, `MEFFCH`, `MTECHCH` |
| `sum(gx)>=1, sum(gy)>=1, rts=RTS_VRS1` | 双列 `D11_tei`, `D11_teo`，最终给 `_tei/_teo` 成对结果 |

#### 本章实际用到的参数组合

- `gx=[0,0,0], gy=[1], rts=RTS_VRS1, tech=CONTEMPORARY`
- `gx=[1,0,0], gy=[0], rts=RTS_VRS1, tech=CONTEMPORARY`

### `MQDEA.MQDDF(...)`

参数集合与 `MQDEA` 一致，只是底层静态效率从 `DEA2` 换成 `DDF2`。因此：

- 输入方向时，底层读 `tei`
- 产出方向时，底层读 `teo`
- 双向时，可同时保留 `tei`、`teo`

## 6. 本章输出文件

- `table3.1.md`：描述统计
- `table3.2_MQDEA_VRS_.md`：`MQDEA` 主结果
- `table3.3_MQDDF_VRS.md`：`MQDDF` 主结果
- `table3.4.md`：省级描述统计
- 图：四张 `chap3分...png`

## 7. 阅读与运行提醒

- 这章真正要看的是 `MQ`、`MEFFCH`、`MTECHCH` 三列怎么随参数变化。
- 目录里还有 `MQDEA_CRS.md`、`MQDDF_CRS.md` 等旧输出，但 notebook 主过程主要导出 `table3.x` 这一组。
