# Changelog

All notable changes are recorded here.

## [0.3.6] - 2026-08-14

### Changed

- Type the public API with `numpy.typing.ArrayLike` inputs and precise return
  types; the stable core now type-checks cleanly under mypy (the experimental
  conformal module is exempt, and a `[tool.mypy]` configuration is included).
  Runtime behaviour is unchanged.

## [0.3.5] - 2026-08-14

### Changed

- The experimental routines (`block_conformal`, `split_conformal`,
  `ConformalBand`, `twoscale_trend`, `wasserstein_decomposition`,
  `TwoScaleResult`) are collected under a new `nsevt.experimental` namespace and
  documented as outside the public-API stability guarantee. They remain
  importable from the top level for backward compatibility.

## [0.3.4] - 2026-08-14

### Documentation

- The README now documents the interval-censored (grouped) fit
  (`gpd_pot_grouped`, `nsevt.grouped`) and the grid-independent minimum-
  detectable effect, so the project page reflects the full 0.3.x feature set.

## [0.3.3] - 2026-08-14

### Added

- A project logo, shown at the top of the README.
- `py.typed` marker (PEP 561): the package now ships type information for
  downstream type checkers and declares the "Typing :: Typed" classifier.

### Changed

- The permutation trend test and the Monte Carlo power / minimum-detectable-
  effect analysis warm-start each refit from the null (or generating) estimate,
  with a multi-start fallback, giving identical results about five times faster.

## [0.3.2] - 2026-08-14

### Changed

- The endpoint bootstrap in `gpd_pot` / `upper_endpoint` warm-starts each
  resample from the point estimate (with a multi-start fallback), giving
  identical intervals about ten times faster.
- The grouped shape-profile interval profiles the scale with a bounded 1-D
  search instead of a simplex over a length-one vector (identical results,
  fewer evaluations), and `interval_cells` accepts a 0-d array grid.

## [0.3.1] - 2026-08-14

### Changed

- Author name recorded as the compound surname `Salomón-Guzmán` in the package
  metadata and citation files, for consistent academic citation.

### Fixed

- The continuous-integration smoke check no longer asserts a hardcoded version,
  so version bumps do not fail the build.

## [0.3.0] - 2026-08-13

### Added

- `nsevt.grouped`: interval-censored (grouped) GPD tail inference for
  discretised exceedances (`gpd_pot_grouped`, `fit_gpd_grouped`,
  `interval_cells`, `profile_ci_xi_grouped`, `profile_endpoint_ci`,
  `GroupedGPDFit`). Fitting the continuous GPD to rounded values biases the
  shape and the finite endpoint; the interval-censored likelihood removes that
  bias, and the endpoint interval is a profile-likelihood interval on the
  reparameterised endpoint rather than a percentile bootstrap.
- `min_detectable_effect` now reports a grid-independent detectable effect from a
  monotone interpolation of the power curve (`emd_positive`, `emd_negative`,
  `emd_per_decade`) with Monte Carlo uncertainty intervals (`emd_*_ci95`),
  alongside the existing grid-based `mde_*` fields.

### Removed

- The maintainer publishing checklist and the PyPI publish workflow are no
  longer tracked in the public repository.

## [0.2.0] - 2026-08-09

### Changed

- `GPDFit.bounded` now means that the full 95% profile interval supports
  `xi < 0`; `bounded_estimate` retains point-estimate semantics.
- Replaced the fixed shape-profile grid with adaptive likelihood-ratio root
  bracketing, explicit numerical-boundary reporting, and the conventional
  likelihood-existence restriction `xi > -1`.
- Clarified endpoint intervals as conditional, model-based summaries rather
  than physical ceilings.
- Calibrated the trend test with complete-block labels, exposed Monte Carlo
  p-value error, and made positive/negative MDEs explicit.
- Replaced binary transportability claims with power- and direction-aware
  multi-source robustness statuses.
- Marked block-conformal and distribution-valued trend routines experimental
  and narrowed their documented guarantees.
- Corrected conformal demo coverage to use raw exceedance values.

### Added

- Statistical assumption, API, and validation documentation.
- Input validation, regression tests, linting, coverage, wheel smoke testing,
  and Trusted Publishing workflow.

## [0.1.0] - 2026-08-09

- Initial PyPI and Zenodo release.
