# Usage

deabook is built around **factory functions** plus an optional **panel** wrapper.
Every efficiency factory takes your data and a `sent` declaration, returns an
*unsolved* model, and is solved with `.optimize(email, solver)`.

```python
import deabook                       # top-level package
from deabook import DEA, DEAweak, PANEL
from deabook.constant import RTS_VRS1, OPT_LOCAL
```

## The `sent` declaration

`sent` is a small string that names your variables:

| Family | form | example |
| --- | --- | --- |
| Standard | `"inputs=outputs"` | `"K L E=Y"` |
| Bad outputs | `"inputs=outputs:undesirable"` | `"K L E=Y:CO2"` |

## 1. A standard (radial) efficiency model

```python
import pandas as pd
from deabook import DEA
from deabook.constant import OPT_LOCAL

data = pd.read_excel("oecd data.xlsx")   # columns: K, L, E, Y, ...

m = DEA.ccr(data, sent="K L E=Y", orient="oo")   # CRS, output-oriented
res = m.optimize(OPT_LOCAL, "highs")
print(res[["te"]])                                # te in (0, 1]
```

`DEA.bcc` is the VRS analogue; `DEA.hyperbolic`, `DEA.ite`, `DEA.ddf`,
`DEA.sbm`, `DEA.sddf` cover other static measures.

## 2. Undesirable outputs (weak disposability)

```python
from deabook import DEAweak

# Slack-based directional distance function with CO2 as a bad output
# (Pham & Zelenyuk 2018)
m = DEAweak.sddf(data, sent="K L E=Y:CO2", rts=RTS_VRS1)
res = m.optimize(OPT_LOCAL, "highs")
print(res[["objective_value", "slack_x_K", "slack_y_Y", "slack_b_CO2"]])
```

`DEAweak.ushephard` / `hyper` / `ddf` / `nddf` / `sbmu` follow the same pattern.

## 3. Panel data

```python
from deabook import PANEL

p = PANEL.tech(
    data,
    id="country",                       # DMU identifier column
    year="year",                        # time column
    sent="K L E=Y",
    method=DEA.ccr,                    # any static factory
    orient="oo",                       # forwarded to DEA.ccr
    tech="contemporary",               # or sequential / window / global
)
long_df = p.solve()                    # one row per DMU-period
wide_df = p.pivot()                   # rows=DMU, cols=year, cells=efficiency
```

The same `Panel` drives every static interface in the package
(`DEA.*`, `DEAweak.*`, `DEAmb.MB` / `MBx` / `MBxy` / `MB2`).

## Solvers

deabook builds LP/NLP models with [Pyomo](https://pyomo.readthedocs.io) and
needs a solver installed in your environment. The free
[HiGHS](https://highs.dev) solver (`pip install highspy`, use `solver="highs"`)
is recommended and used throughout the test suite. GLPK/CBC and commercial
solvers (Gurobi, CPLEX, Mosek) are supported wherever Pyomo can find them.
