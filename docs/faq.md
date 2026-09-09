# FAQ

### How do I declare my variables?

Use the `sent` string. Standard models use `"inputs=outputs"` (e.g.
`"K L E=Y"`). Models with undesirable outputs append them after a colon:
`"inputs=outputs:undesirable"` (e.g. `"K L E=Y:CO2"`). The material-balance
models use a `"in_np + in_p = out_np + out_p : undesirable"` form — see
`deabook.DEAmb.MB`.

### Which functions return an efficiency score directly, and which return `rho`?

- `DEA.ccr` / `DEA.bcc` return a ready-made efficiency column `te` (in `(0, 1]`).
- `DEA.hyperbolic`, `DEA.ite`, `DEA.ddf` and the `DEAweak` directional models
  return the raw LP value `rho` (plus `objective_value`); you convert `rho` to an
  efficiency according to the chosen direction (e.g. output direction
  `E = 1/(1+rho)`, input direction `E = 1-rho`).
- `DEA.sbm` / `DEA.sddf` / `DEAweak.sbmu` / `DEAweak.sddf` return `rho` /
  `objective_value` together with per-variable `slack_*` columns.

### How do I choose the solver?

Pass it as the second argument to `.optimize(email, solver)`, e.g.
`m.optimize(OPT_LOCAL, "highs")`. HiGHS (`highs`) is recommended and is what the
test suite uses. See [Installation](installation.md).

### My `.optimize()` raises "solver not found".

deabook does not bundle a solver. Install one (e.g. `pip install highspy`) and
make sure it is on your `PATH` / importable. See
[Installation](installation.md).

### How do I run a model on panel data?

Use `PANEL.Panel`, injecting any static factory as `method`. It loops the
model over periods under a chosen reference technology (`contemporary`,
`sequential`, `window`, `global`) and returns a long `DataFrame` (`.solve()`) or
a wide DMU×year table (`.pivot()`).

### What is the difference between `DEA.sddf` and `DEAweak.sddf`?

`DEA.sddf` is the Färe & Grosskopf (2010) slack-based directional distance
function **without** undesirable outputs. `DEAweak.sddf` is the Pham & Zelenyuk
(2018) extension that adds a bad-output slack under joint weak disposability.
Their `sent` must match: `DEA.sddf` rejects a `:`-separated bad part, and
`DEAweak.sddf` requires one.

### Where do I report bugs or request features?

Open an issue at <https://github.com/advancehs/deabook/issues>.
