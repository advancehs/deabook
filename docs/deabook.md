# API Reference — Estimators

deabook is organised around **factory functions**, not classes you instantiate
directly. Every efficiency factory takes your data plus a `sent` declaration
and returns an *unsolved* model object; you then call `.optimize(email, solver)`
on it to obtain a `pandas.DataFrame`.

```python
from deabook import DEA, DEAweak, PANEL
from deabook.constant import RTS_VRS1, OPT_LOCAL

m = DEA.ccr(data, sent="K L E=Y", orient="oo")   # build (no solve yet)
res = m.optimize(OPT_LOCAL, "highs")             # solve -> DataFrame
```

The returned `DataFrame` is always indexed by the original DMU row index and
contains at least an `optimization_status` column; the efficiency column(s)
depend on the model (e.g. `te`, `rho`, `objective_value`, `obj`).

The `sent` (specification) string convention used throughout the package is:

| Model family | `sent` form | Example |
| --- | --- | --- |
| Standard (no bad outputs) | `"inputs=outputs"` | `"K L E=Y"` |
| Weak disposability (bad outputs) | `"inputs=outputs:undesirable"` | `"K L E=Y:CO2"` |
| Material balance | `"in_np + in_p = out_np + out_p : undesirable"` | see `deabook.DEAmb.MB` |

---

## `deabook.DEA` — standard models (no undesirable outputs)

::: deabook.DEA.ccr

::: deabook.DEA.bcc

::: deabook.DEA.hyperbolic

::: deabook.DEA.ite

::: deabook.DEA.ddf

::: deabook.DEA.sbm

::: deabook.DEA.sddf

> **Malmquist / productivity indices** are also exposed from this module for
> panel data: `DEA.mq`, `DEA.mqddf`, `DEA.mqfgnz`, `DEA.mqrd`, `DEA.gmqdecom`,
> `DEA.ep`, `DEA.epddf`. They share the `id` / `year` / `tech` interface
> described for the panel wrapper below.

---

## `deabook.DEAweak` — weak-disposability models (with undesirable outputs)

All factories below require a `sent` that declares the undesirable part after
a `:` (e.g. `"K L E=Y:CO2"`). They implement joint weak disposability of the
bad outputs.

::: deabook.DEAweak.ushephard

::: deabook.DEAweak.hyper

::: deabook.DEAweak.ddf

::: deabook.DEAweak.nddf

::: deabook.DEAweak.sbmu

::: deabook.DEAweak.sddf

The following three are *dynamic* (panel) environmental-productivity indices;
they take `id` / `year` / `tech` instead of `evalindex` / `refindex`:

::: deabook.DEAweak.dep

::: deabook.DEAweak.depddf

::: deabook.DEAweak.lddf

---

## `deabook.DEAmb` — material-balance models

::: deabook.DEAmb.MB

---

## `deabook.panel` — panel-data wrapper

`Panel` turns **any** static (single-period) efficiency factory into a panel
estimator by looping it over the observed periods under a chosen reference
technology.

::: PANEL.Panel
