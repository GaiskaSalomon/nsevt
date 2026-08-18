# Changelog

All notable changes are recorded here.

## [Unreleased]

## [1.0.2] - 2026-08-16

### Fixed

- Type-annotate `design._validated_inputs` so the validated arrays flow back to
  their callers as `ndarray`; the stable core type-checks cleanly under mypy
  again (`nsevt.design` had regressed).

### Changed

- The continuous-integration lint job now runs `mypy src/nsevt`, so the
  type-clean guarantee is enforced rather than checked only by hand.

### Tested

- Cover the input-validation paths added in 1.0.1 for `nsevt.design`
  (malformed censoring cells, invalid `profile_ci_coef` / `profile_ci_return_level`
  controls, a non-finite shape in `return_level`, and the heavy-tail case where
  `profile_ci_return_level` returns `None`), confirming the hardening triggers as
  intended. No runtime or API change.

## [1.0.1] - 2026-08-15

### Fixed

- Validate the grouped-design sample, design rank, censoring cells, profile
  controls, return-level domain, and optimizer success before returning a fit.
- Validate sequential Monte Carlo kinds, quantiles, tolerances and replicate
  budgets; reject short `draw` results and prevent an initial block from
  exceeding `r_max`.
- Use block quantiles, rather than block means, in the independent-block
  diagnostic for quantile runs.
- Stabilise proportion MCSEs at observed proportions of exactly zero or one so
  a finite Monte Carlo run cannot report zero simulation error.
- Emit a `RuntimeWarning` when permutation refits fail and the reported p-value
  is therefore conditional on the successful refits.
- Reject invalid pseudo-true simulation budgets and non-finite proxy estimates.
- Run a public-API smoke test against clean installations of both the wheel and
  source distribution in CI and before the release workflow can publish.

### Documentation

- Clarify that the interpolated EMD is not restricted to effect-grid nodes but
  still depends on the chosen grid and interpolation; its interval is a
  pointwise-normal approximation that does not model common-random-number
  covariance across the power curve.
- Describe `pseudo_true` as a large-sample simulation proxy whose sensitivity
  to sample size and seed must be checked, rather than as proof of an
  estimator's limiting target.

## [1.0.0] - 2026-08-15

**The public API is now stable.** From this release, the documented public
surface — the `nsevt.gpd`, `nsevt.grouped`, `nsevt.trend`, `nsevt.transportability`,
`nsevt.mc`, `nsevt.calibration` and `nsevt.design` names, their signatures, and
their return schemas (`docs/return-schemas.md`) — is frozen: a backward-incompatible
change to it will require a 2.0.0. Additions (new functions, new optional
arguments, new return-dict keys) remain minor (1.x) changes.

No code changes relative to 0.4.1; this release marks the stability commitment.
The routines under `nsevt.experimental` (conformal aggregation, distribution-valued
trend) remain outside the stability guarantee, as documented.

### Changed

- Development status classifier raised to `5 - Production/Stable`.

## [0.4.1] - 2026-08-14

### Added

- `nsevt.design`: grouped GPD regression and return levels. `fit_grouped_design`
  fits the interval-censored GPD with a log-linear scale design matrix (a trend,
  group-specific scales, or any combination); `profile_ci_coef` is a
  profile-likelihood interval for any coefficient of that design (the interval
  counterpart of the permutation trend test); `return_level` gives the level
  exceeded once per `m` observations at an exceedance rate, and
  `profile_ci_return_level` its profile-likelihood interval, obtained by
  profiling the level itself rather than pushing a profiled shape through the
  return-level formula. The design, covariate coding and exceedance rate are
  caller-supplied; the module is NumPy/SciPy-only.

## [0.4.0] - 2026-08-14

### Added

- `nsevt.mc`: a sequential Monte Carlo precision protocol. Monte Carlo standard
  errors for proportions, means and quantiles (`mcse_proportion`, `mcse_mean`,
  `mcse_quantile`); the replicate budget for a target precision
  (`required_replicates`); a floor-aware permutation/bootstrap p-value
  (`permutation_pvalue`); and `SequentialRun` / `run_sequential`, which grow a
  run in blocks and stop only when the MCSE target, estimate stability and every
  registered qualitative decision have all settled, reporting a full trace,
  an independent-block diagnostic and a `not_stabilised` status when a run
  exhausts its budget. Reproducible, non-interfering random substreams
  (`substream`, `block_streams`) and a cross-seed audit (`multiseed_summary`)
  support extending a run without perturbing the replicates already drawn.
  The module depends only on NumPy and is estimator-agnostic: it decides how
  many replicates any power, coverage, p-value or bootstrap analysis needs.
- `nsevt.calibration`: finite-sample calibration of any estimator or test by
  Monte Carlo, from user-supplied `simulate` / `estimator` / `test` callables.
  `rejection_rate` measures the empirical type-I error (under a null DGP) or
  power (under an alternative); `coverage` measures an interval estimator's
  empirical coverage against a nominal or `pseudo_true` target; `bias_rmse`
  reports the bias, standard deviation and RMSE of a point estimator; and
  `pseudo_true` gives a large-sample simulation proxy for a pseudo-true target
  under a misspecified DGP. The two proportion analyses run on the `nsevt.mc`
  sequential protocol, so each carries a Monte Carlo standard error and a
  stopping decision rather than a fixed replicate budget.

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
  (`gpd_pot_grouped`, `nsevt.grouped`) and the interpolated minimum-detectable
  effect, so the project page reflects the full 0.3.x feature set.

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
- `min_detectable_effect` now reports an interpolated detectable effect from a
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
