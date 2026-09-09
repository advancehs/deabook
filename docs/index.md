# Welcome to deabook

[![PyPI version](https://img.shields.io/pypi/v/deabook.svg)](https://pypi.python.org/pypi/deabook)

**deabook** is a Python package for **Data Envelopment Analysis (DEA)** and
related efficiency / productivity measurement, including undesirable (bad)
outputs, panel data, and stochastic frontier methods. It is the companion code
to the book *Data Envelopment Analysis* (in Chinese).

- Free software: MIT License
- Documentation: <https://advancehs.github.io/deabook>
- Source & issues: <https://github.com/advancehs/deabook>

## Features

- **Standard efficiency models** (no undesirable outputs) in `deabook.DEA`:
  radial CCR/BCC, hyperbolic, single-input technical efficiency, directional
  distance function (DDF), Slack-Based Measure (`sbm`), and the
  slack-based directional distance function / SD model (`sddf`,
  Färe & Grosskopf 2010).
- **Weak-disposability models** (with undesirable outputs) in `deabook.DEAweak`:
  undesirable-output Shephard distance, hyperbolic / DDF / non-radial DDF with
  bad outputs, SBM with bads (`sbmu`), the Pham & Zelenyuk (2018) SD model
  (`sddf`), and dynamic environmental-productivity indices (`dep`, `depddf`,
  `lddf`).
- **Material-balance models** in `deabook.DEAmb` (`MB` and its panel variants).
- **Panel-data wrapper** `PANEL.Panel`: turn any static efficiency
  factory into a panel estimator under contemporary / sequential / window /
  global reference technologies.
- **Malmquist & Luenberger productivity indices** for both standard and
  undesirable-output settings.
- **Stochastic / convex frontier methods** (`StoNED`, `CNLSSDFDDF`,
  `CNLSSDFDDFweak`).
- A uniform, factory-based API: every model is built with a `sent` variable
  declaration and solved with `.optimize(email, solver)`, returning a
  `pandas.DataFrame`.

## Design in one sentence

> `data` → `sent` → factory → `.optimize()` → `DataFrame`.

See [Usage](usage.md) for a worked example and the
[API Reference](deabook.md) for every function.
