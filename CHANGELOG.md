# Changelog

All notable changes are recorded here.

## [0.3.0] - Unreleased

### Added

- `nsevt.grouped`: interval-censored (grouped) GPD tail inference for
  discretised exceedances (`gpd_pot_grouped`, `fit_gpd_grouped`,
  `interval_cells`, `profile_ci_xi_grouped`, `profile_endpoint_ci`,
  `GroupedGPDFit`). Fitting the continuous GPD to rounded values biases the
  shape and the finite endpoint; the interval-censored likelihood removes that
  bias, and the endpoint interval is a profile-likelihood interval on the
  reparameterised endpoint rather than a percentile bootstrap.

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
