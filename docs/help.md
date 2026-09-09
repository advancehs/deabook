# deabook 使用帮助文档

本文档按模块详细说明 `deabook` 包中所有**公开接口**的类、方法、参数与返回值。

> **接口风格说明**：从 0.0.11 起，所有静态效率模型都通过**模块级工厂函数**创建（如
> `DEA.ccr`、`DEA.sddf`、`DEAweak.sddf`），工厂返回**未求解**的模型对象，再调用
> `.optimize(email, solver)` 得到 `pandas.DataFrame`。旧版直接实例化的类
> （`DEA.DEA`、`DEA.DDF`、`MQDEA`、`DEAweak.DDFweak2` 等）已移除，本文档只描述新接口。

---

## 1. 常量模块 (`constant`)

```python
from deabook.constant import *
```

### 1.1 复合误差项类型 (CET)

| 常量 | 值 | 含义 |
|------|----|------|
| `CET_ADDI` | `"addi"` | 加性复合误差项 |
| `CET_MULT` | `"mult"` | 乘性复合误差项 |

### 1.2 前沿类型 (FUN)

| 常量 | 值 | 含义 |
|------|----|------|
| `FUN_PROD` | `"prod"` | 生产前沿 |
| `FUN_COST` | `"cost"` | 成本前沿 |

### 1.3 规模报酬类型 (RTS)

| 常量 | 值 | 含义 |
|------|----|------|
| `RTS_VRS1` | `"vrs1"` | 可变规模报酬（统一缩减因子 theta，theta bounds=[0,1]） |
| `RTS_VRS2` | `"vrs2"` | 可变规模报酬（差异化缩减因子 mu，mu bounds=[0, +inf]，按参考单元索引） |
| `RTS_VRS` | `"vrs"` | 通用可变规模报酬（经典 DEA/DDF 用） |
| `RTS_CRS` | `"crs"` | 不变规模报酬（无 VRS 约束） |

### 1.4 残差分解方法 (RED)

| 常量 | 值 | 含义 |
|------|----|------|
| `RED_MOM` | `"MOM"` | 矩估计法 (Method of Moments) |
| `RED_QLE` | `"QLE"` | 拟似然估计法 (Quasi-Likelihood Estimation) |
| `RED_KDE` | `"KDE"` | 核密度估计法 (Kernel Density Estimation)，使用 Richardson-Lucy 盲反卷积 |

### 1.5 求解器选项 (OPT)

| 常量 | 值 | 含义 |
|------|----|------|
| `OPT_LOCAL` | `"local"` | 本地求解 |
| `OPT_DEFAULT` | `None` | 使用默认求解器 |

### 1.6 技术类型 (TECH)

| 常量 | 值 | 含义 |
|------|----|------|
| `CONTEMPORARY` | `"Contemporary production technology"` | 当期生产技术（仅当期参考集合） |
| `SEQUENTIAL` | `"sequential"` | 序列生产技术（参考集合 = 年份 `<= t`） |
| `WINDOW` | `"window"` | 窗口生产技术（以 `t` 为中心的窗口，宽度 `window`） |
| `GLOBAL` | `"Global production technology"` | 全局生产技术（所有期参考集合并） |

### 1.7 指数 / 导向类型

| 常量 | 值 | 含义 |
|------|----|------|
| `ORIENT_IO` / `ORIENT_OO` / `ORIENT_BO` | `"io"` / `"oo"` / `"bo"` | 投入导向 / 产出导向 / 非期望产出导向 |
| `MAL` | — | Malmquist / Malmquist-Luenberger 指数 |
| `LUE` | — | Luenberger 指数 |

---

## 2. DEA 模块（标准模型，无期望产出）

```python
from deabook import DEA
from deabook.constant import RTS_VRS1, OPT_LOCAL
```

所有工厂的通用参数：

| 参数 | 含义 |
|------|------|
| `data` | `DataFrame`，一行一个 DMU |
| `sent` | `"投入=产出"`，如 `"K L E=Y"` |
| `rts` | `RTS_CRS` 或 `RTS_VRS` / `RTS_VRS1`（按接口而定） |
| `evalindex` / `refindex` | 评价 / 参考样本筛选条件，如 `"Year=[2010]"`；默认全部 |
| `email` / `solver` | 传入 `.optimize(email, solver)` |

### 2.1 `DEA.ccr` — 径向 CCR（不变规模报酬）

```python
DEA.ccr(data, sent, orient="oo", evalindex=None, refindex=None)
```
`orient="oo"`（产出导向）或 `"io"`（投入导向）。返回对象 `.optimize()` 含 `te` 列（∈ (0,1]）。

### 2.2 `DEA.bcc` — 径向 BCC（可变规模报酬）

```python
DEA.bcc(data, sent, orient="oo", evalindex=None, refindex=None)
```
与 `ccr` 相同，但 `rts=RTS_VRS`。同样返回 `te`。

### 2.3 `DEA.hyperbolic` — 双曲距离函数

```python
DEA.hyperbolic(data, sent, rts=RTS_CRS, evalindex=None, refindex=None)
```
同时等比例收缩所有投入并扩张所有产出。`.optimize()` 返回 `rho` 与 `objective_value`。

### 2.4 `DEA.ite` — 单投入技术效率

```python
DEA.ite(data, sent, rts=RTS_CRS, evalindex=None, refindex=None)
```
只收缩 `sent` 中**第一个**投入（如能源），其余投入与产出固定。返回 `rho` / `objective_value`。

### 2.5 `DEA.ddf` — 方向距离函数

```python
DEA.ddf(data, sent, gy=[1], gx=[0], rts=RTS_CRS, evalindex=None, refindex=None)
```
沿显式方向向量 `(gx, gy)` 移动。返回 `rho` / `objective_value`。效率换算：投入导向 `E=1-rho`，产出导向 `E=1/(1+rho)`，双向 `E_tei=1-rho`、`E_teo=1/(1+rho)`。

### 2.6 `DEA.sbm` — 基于松弛的测度（SBM，无坏产出）

```python
DEA.sbm(data, sent="inputvar=outputvar", gx=None, gy=None, gb=None,
        rts=RTS_VRS1, evalindex=None, refindex=None)
```
Tone (2001) SBM。返回 `rho` 及 `slack_x_*` / `slack_y_*` 列。`sent` **不能**含 `:` 分隔的坏产出（否则抛错，改用 `DEAweak.sbmu`）。

### 2.7 `DEA.sddf` — 基于松弛的方向距离函数（SD，Färe & Grosskopf 2010）

```python
DEA.sddf(data, sent="inputvar=outputvar", gx=None, gy=None,
         rts=RTS_VRS1, evalindex=None, refindex=None)
```
最大化各投入、各期望产出的独立绝对松弛之和。返回 `objective_value`（= 松弛和，即 SD 无效率）及 `slack_x_*` / `slack_y_*` 列；技术效率 `te = 1/(1+objective_value)`。`sent` **不能**含坏产出（否则抛错，改用 `DEAweak.sddf`）。

---

## 3. DEAweak 模块（弱可处置性，含非期望产出）

```python
from deabook import DEAweak
```
所有工厂的 `sent` 必须形如 `"投入=期望产出:非期望产出"`，如 `"K L E=Y:CO2"`，实现非期望产出的联合弱可处置性。

### 3.1 `DEAweak.ushephard` — 非期望产出导向谢泼德距离（纯环境效率 PEE）

```python
DEAweak.ushephard(data, sent, rts=RTS_CRS, evalindex=None, refindex=None)
```
方向固定为 `gx=0, gy=0, gb=1`。`.optimize()` 返回 `rho`（= PEE）与 `objective_value`。

### 3.2 `DEAweak.hyper` — 含非期望产出的双曲距离函数（混合环境效率）

```python
DEAweak.hyper(data, sent, gx=None, gy=None, gb=None, rts=RTS_CRS,
              evalindex=None, refindex=None)
```
默认方向 `gx=0, gy=1, gb=1`。返回 `rho` / `objective_value`。

### 3.3 `DEAweak.ddf` — 含非期望产出的方向距离函数

```python
DEAweak.ddf(data, sent, gx=None, gy=None, gb=None, rts=RTS_CRS,
            evalindex=None, refindex=None)
```
默认方向 `gx=0, gy=1, gb=1`。返回 `rho` / `objective_value`。

### 3.4 `DEAweak.nddf` — 非径向方向距离函数（各变量独立调整）

```python
DEAweak.nddf(data, sent, gx=None, gy=None, gb=None, rts=RTS_CRS,
             evalindex=None, refindex=None)
```
返回 `objective_value`（各变量 `beta` 的加权和）及各变量的 `beta_<变量名>` 列。

### 3.5 `DEAweak.sbmu` — 基于松弛的测度（SBM，含坏产出）

```python
DEAweak.sbmu(data, sent="inputvar=outputvar:unoutputvar", gx=None, gy=None, gb=None,
              rts=RTS_VRS1, evalindex=None, refindex=None)
```
Tone (2003/2010) 带坏产出的 SBM。返回 `rho` 及 `slack_x_*` / `slack_y_*` / `slack_b_*` 列。`sent` **必须**含 `:` 分隔的坏产出。

### 3.6 `DEAweak.sddf` — SD 模型（含非期望产出，Pham & Zelenyuk 2018）

```python
DEAweak.sddf(data, sent="inputvar=outputvar:unoutputvar", gx=None, gy=None, gb=None,
             rts=RTS_VRS1, evalindex=None, refindex=None)
```
在联合弱可处置性下追加非期望产出的松弛项。返回 `objective_value` 及 `slack_x_*` / `slack_y_*` / `slack_b_*` 列。`sent` **必须**含坏产出。

> 注意：在联合弱可处置性下，仅产出无效率的 DMU 在弱模型中可能变为 SD 有效（扩张期望产出会同时要求扩张非期望产出，使其松弛最小）；只有各方面均被支配的 DMU 才会被判为无效率。这是**正确行为**，不是 bug。

### 3.7 动态环境生产率指数（面板）

以下三个接口接受 `id` / `year` / `tech`（而非 `evalindex` / `refindex`），返回
`MalmquistModel`，调用 `.optimize(email, solver)` 求解：

```python
from deabook import PANEL
DEAweak.dep(data, id, year, sent="inputvar=outputvar:unoutputvar",
            rts=RTS_CRS, tech=CONTEMPORARY, window=3)
DEAweak.depddf(data, id, year, sent="inputvar=outputvar:unoutputvar",
               rts=RTS_CRS, tech=CONTEMPORARY, window=3)
DEAweak.lddf(data, id, year, sent="inputvar=outputvar:unoutputvar",
             orient="oo", rts=RTS_CRS, tech=CONTEMPORARY, window=3)
```

- `dep`：动态环境绩效指数（MEPI），基于 `ushephard`。返回 `D11`/`D12`/`D21` + `MQ`/`MEFFCH`/`MTECHCH`；`GLOBAL` 下为 `D11`/`mqpi`/`MEFFCH`/`MTECHCH`。
- `depddf`：基于方向距离函数的动态环境绩效指数（ML），基于坏产出导向 DDF（`gb=1`）。
- `lddf`：基于方向距离函数的加法 Luenberger 指标，`orient` 取 `"oo"`/`"io"`/`"both"`；返回 `D11`/`D12`/`D21` + `L`/`LEFFCH`/`LTECHCH`（基准 0）。

---

## 4. panel 模块（面板数据包装器）

```python
from deabook import PANEL
```

`PANEL.tech` 把**任意**静态效率工厂变成面板估计量，按所选参考技术在各期上循环。

```python
from deabook import PANEL
p = PANEL.tech(data, id="country", year="year", sent="K L E=Y",
                method=DEA.ccr, orient="oo", tech="contemporary")
long_df = p.solve()        # 长表：每行一个 DMU-期
wide_df = p.pivot()        # 宽表：行=DMU，列=year，单元格=效率
```

| 参数 | 含义 |
|------|------|
| `data` | 面板 `DataFrame`，每行一个 DMU-期观测 |
| `id` | DMU 标识列名 |
| `year` | 时间/期变量列名 |
| `sent` | 变量声明（同各工厂） |
| `method` | 注入的静态效率工厂，如 `DEA.ccr`、`DEAweak.sddf`、`DEAmb.MB` |
| `tech` | `contemporary` / `sequential` / `window` / `global` |
| `window` | `tech=window` 时的窗口宽度（默认 3） |
| `**method_kwargs` | 透传给 `method` 的额外参数（`orient`、`rts`、`gx`/`gy`/`gb`、`sx`/`sy`/`level` 等） |

`method` 不接受的多余参数会被自动忽略；`solve()` 的结果按原始行索引回挂 `id`/`year` 列。

---

## 5. Malmquist / 生产率指数（标准情形）

位于 `deabook.DEA`，用于面板数据，对应旧版的 `MQDEA` / `MQDDF`：

```python
from deabook import PANEL
DEA.mq(data, id, year, sent="inputvar=outputvar", orient="oo",
       rts=RTS_CRS, tech=CONTEMPORARY, window=3)
DEA.mqddf(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0],
          rts=RTS_CRS, tech=CONTEMPORARY, window=3)
DEA.mqfgnz(data, id, year, sent="inputvar=outputvar", orient="oo", ...)
DEA.mqrd(data, id, year, sent="inputvar=outputvar", orient="oo", ...)
DEA.gmqdecom(data, id, year, sent="inputvar=outputvar", orient="oo", ...)
DEA.ep(data, id, year, sent="E K L=Y", rts=RTS_CRS, ...)        # 能源专项（基于 ite）
DEA.epddf(data, id, year, sent="K L E=Y", gy=[0], gx=[0,0,1], ...)  # 能源专项 DDF
```

均返回 `MalmquistModel`，`.optimize(email, solver)` 后给出 `D11`/`D12`/`D21` 及
`MQ`/`MEFFCH`/`MTECHCH`（`GLOBAL` 下为 `D11`/`mqpi`/`MEFFCH`/`MTECHCH`）。

---

## 6. MB 模块（物质平衡分解）

```python
from deabook import DEAmb
```

```python
DEAmb.MB(data, sent="inputvar_np + inputvar_p = outputvar_np + outputvar_p : unoutputvar",
         sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], level=5,
         rts=RTS_VRS1, evalindex=None, refindex=None)
```
`.optimize(solver=OPT_DEFAULT, dmu="1")` 返回 `(DataFrame, info)` 二元组。

面板变体（共享同一构造签名，额外要求 `year`）：`DEAmb.MBx`、`DEAmb.MBxy`、`DEAmb.MB2`。
它们按投入 / 投入+产出端做物质平衡分解，等价类为 `MBx2`、`MBxy2`。

---

## 7. CNLSSDFDDFweak 模块（随机凸非参数前沿，含坏产出）

```python
from deabook.CNLSSDFDDFweak import CNLSSDweak, CNLSDDFweak
```

仍为**类**式接口（直接实例化）：

```python
m = CNLSSDweak(data, sent="K L E=Y:CO2", z=None, gy=[1], gx=[0], gb=[0],
               cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1)
m.optimize(OPT_LOCAL, "highs")
m.get_residual()
```
`CNLSDDFweak` 为方向距离函数版本，接口相同。

---

## 8. StoNED 模块（随机非参数包络）

```python
from deabook.StoNED import StoNED
```

`StoNED` 包装一个**已求解**的 CNLS 模型并做残差分解：

```python
model = CNLSSDweak(data, sent="K L E=Y", cet=CET_MULT, rts=RTS_VRS1)
model.optimize(OPT_LOCAL, "highs")
st = StoNED(model)
u = st.get_technical_inefficiency(method=RED_MOM)   # 或 RED_QLE / RED_KDE
```
`method` 取值 `RED_MOM` / `RED_QLE` / `RED_KDE`，返回企业层面条件无效率 `u_hat`。

---

## 9. 使用示例

### 9.1 标准 DEA 分析

```python
import pandas as pd
from deabook import DEA
from deabook.constant import OPT_LOCAL

data = pd.read_excel("oecd data.xlsx")

# 产出导向 CCR（CRS）
m = DEA.ccr(data, sent="K L E=Y", orient="oo")
res = m.optimize(OPT_LOCAL, "highs")
print(res[["te"]])

# 可变规模报酬 BCC
m = DEA.bcc(data, sent="K L E=Y", orient="oo")
res = m.optimize(OPT_LOCAL, "highs")
```

### 9.2 方向距离函数

```python
from deabook import DEA

# 投入导向 DDF（仅能源 E 缩减）
m = DEA.ddf(data, sent="K L E=Y", gy=[0], gx=[0, 0, 1], rts=RTS_CRS)
res = m.optimize(OPT_LOCAL, "highs")
print(res[["rho", "objective_value"]])
```

### 9.3 弱处置性 DEA（非期望产出）

```python
from deabook import DEAweak

# SD 模型（Pham & Zelenyuk 2018），含 CO2
m = DEAweak.sddf(data, sent="K L E=Y:CO2", rts=RTS_VRS1)
res = m.optimize(OPT_LOCAL, "highs")
print(res[["objective_value", "slack_x_K", "slack_y_Y", "slack_b_CO2"]])
```

### 9.4 Malmquist 生产率指数

```python
from deabook import DEA

m = DEA.mq(data, id="country", year="year", sent="K L E=Y",
           orient="oo", rts=RTS_CRS, tech="contemporary")
res = m.optimize(OPT_LOCAL, "highs")
print(res[["MQ", "MEFFCH", "MTECHCH"]])
```

### 9.5 弱处置性动态环境绩效

```python
from deabook import DEAweak

m = DEAweak.dep(data, id="country", year="year", sent="K L E=Y:CO2",
                rts=RTS_CRS, tech="contemporary")
res = m.optimize(OPT_LOCAL, "highs")
print(res[["MQ", "MEFFCH", "MTECHCH"]])
```

### 9.6 面板数据

```python
from deabook import PANEL, DEA

p = PANEL.tech(data, id="country", year="year", sent="K L E=Y",
                method=DEA.ccr, orient="oo", tech="contemporary")
wide = p.solve().pivot()
```

### 9.7 CNLS 随机前沿 + StoNED 分解

```python
from deabook.CNLSSDFDDFweak import CNLSSDweak
from deabook.StoNED import StoNED
from deabook.constant import CET_MULT, FUN_PROD, RTS_VRS1, OPT_LOCAL

model = CNLSSDweak(data, sent="K L E=Y", cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1)
model.optimize(OPT_LOCAL, "highs")

st = StoNED(model)
u = st.get_technical_inefficiency(method="MOM")
print(u)
```
