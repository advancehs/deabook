# Installation

## Stable release

To install deabook, run this command in your terminal:

```
pip install deabook
```

This is the preferred method to install deabook, as it will always install the
most recent stable release.

If you don't have [pip](https://pip.pip.io) installed, this
[Python installation guide](http://docs.python-guide.org/en/latest/starting/installation/)
can guide you through the process.

## From sources

To install deabook from sources, run this command in your terminal:

```
pip install git+https://github.com/advancehs/deabook
```

## Solver dependency

deabook builds its optimisation models with
[Pyomo](https://pyomo.readthedocs.io) and therefore needs a mathematical
programming solver. Install one of the following:

- **HiGHS** (recommended, free): `pip install highspy` — then pass
  `solver="highs"` to `.optimize()`. This is the solver used by the test suite.
- **GLPK / CBC** (free): install `glpk` / `cbc` and use `solver="glpk"` /
  `solver="cbc"`.
- **Commercial** (Gurobi, CPLEX, Mosek): install the solver and its Python
  bindings; deabook picks them up through Pyomo automatically.

No solver is bundled with deabook — `.optimize()` will raise a clear error until
one is available.
