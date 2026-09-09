# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/).

## [0.0.11]

### Added

- **SD (slack-based directional distance function) models:**
  - `DEA.sddf` — Färe & Grosskopf (2010) SD measure without undesirable outputs.
  - `DEAweak.sddf` — Pham & Zelenyuk (2018) SD measure with undesirable outputs
    under joint weak disposability.
  Both return an `objective_value` (sum of slacks) plus per-variable
  `slack_*` columns and share the `sbm`/`sbmu` engine and return contract.

### Changed

- **Unified factory-function API.** Static efficiency models are now created
  through module-level factory functions (`DEA.ccr`, `DEA.bcc`, `DEA.ddf`,
  `DEA.sbm`, `DEA.sddf`, `DEAweak.ddf`, `DEAweak.sbmu`, `DEAmb.MB`, ...) that
  return an unsolved model; call `.optimize(email, solver)` to solve. The legacy
  directly-instantiated estimator classes (`DEA.DEA`, `DEA.DDF`, `MQDEA`,
  `DEAweak.DDFweak2`, ...) have been removed.
- **Panel support via `PANEL.Panel`.** A single method-agnostic wrapper
  drives every static efficiency factory over panel data under contemporary /
  sequential / window / global reference technologies.

### Fixed

- `DEA.sbm` now rejects `sent` declarations that contain undesirable outputs
  (use `DEAweak.sbmu` instead), and `DEAweak.sbmu` / `DEAweak.sddf` reject
  `sent` without a bad-output part.
