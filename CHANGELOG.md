# Changelog

All notable changes are recorded here.

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
