# API Reference — Constants, solvers & technologies

The `deabook.constant` module holds the string constants you pass to the
factories (returns-to-scale, solver mode, reference technology, orientation) and
the `validate_*` helpers used internally.

```python
from deabook.constant import RTS_VRS1, OPT_LOCAL
```

::: deabook.constant

## Quick reference

| Constant | Value | Meaning |
| --- | --- | --- |
| `RTS_CRS` | `"crs"` | Constant returns to scale |
| `RTS_VRS` | `"vrs"` | Variable returns to scale (generic) |
| `RTS_VRS1` | `"vrs1"` | VRS, identical reduction factor across bad outputs |
| `RTS_VRS2` | `"vrs2"` | VRS, different reduction factor per bad output |
| `OPT_LOCAL` | `"local"` | Solve locally (pass to `email` of `.optimize`) |
| `OPT_DEFAULT` | `None` | Use the package default solver (pass to `solver`) |
| `CONTEMPORARY` | — | Reference set = current period only |
| `SEQUENTIAL` | `"sequential"` | Reference set = periods `<= t` |
| `WINDOW` | `"window"` | Centred window of width `window` around `t` |
| `GLOBAL` | `"Global production technology"` | Reference set = all periods |
| `ORIENT_IO` / `ORIENT_OO` / `ORIENT_BO` | `"io"` / `"oo"` / `"bo"` | Input / output / bad-output orientation |
