# Changelog

All notable changes are recorded here.

## [Unreleased]

### Changed

- Fix the declared build-backend floor: PEP 639 (`license` SPDX string,
  `project.license-files`) needs `setuptools>=77.0.1`, not `>=61.0`. A
  `build-min` CI job now builds the distributions with exactly that pin so the
  floor stays honest.
- Broaden the test matrix and coverage. CI runs Python 3.9 through 3.14 (adds
  3.10 and 3.12), measures branch coverage (`--cov-branch`), adds a
  minimum-dependency job on the declared `numpy>=1.22` / `scipy>=1.7` floors,
  and marks the full-sequential-MC tests `slow` so `pytest -m "not slow"` is a
  fast inner loop. Add `tools/validation_campaign.py`, which writes a
  provenance-stamped JSON archive of per-cell coverage / bias / RMSE with their
  MCSE and replicate counts, marking a cell `"resolved": false` when it does
  not reach the MCSE tolerance.

### Documentation

- Reconcile the stability labels: `split_conformal` is listed as experimental
  (conditional guarantee) in the README, matching `docs/stability.md`. Add an
  implementation-progress table to the review document and link it from
  `docs/index.md`.
- Record the immutable Zenodo DOI assigned to version 1.2.0.

## [1.2.0] - 2026-09-06

### Changed

- Strengthen the type gate. mypy now runs with `check_untyped_defs` (the body
  of an unannotated function is checked, not skipped) and `nsevt.conformal` is
  no longer exempt, so the whole package is type-checked. Add `nsevt._types`
  with `TypedDict` schemas for the main stable-core return dictionaries
  (`permutation_pvalue`, the grouped profile functions, the calibration
  functions, `trend_power` rows, `fit_gpd_grouped`) and annotate those return
  types and the previously-unannotated grouped profile signatures; the runtime
  still returns plain `dict`. Add `tests/test_api_contracts.py`, which asserts
  the frozen 1.x dataclass fields, dict keys, and name aliases are still
  present.
- Pin the lint/type/build tooling (Ruff, mypy, build, twine) in
  `constraints/ci.txt` and apply it in the single-Python `lint` and
  `build`/release jobs and the documented local gate, so a developer machine
  and CI reach the same lint and type verdicts. The multi-version `test` matrix
  and NumPy/SciPy stay unpinned (each interpreter resolves a compatible build;
  pytest 9 has no Python 3.9 distribution); an exact runtime pin for a
  simulation campaign is reserved for a separate `constraints/science.txt`.

### Fixed

- Carry the simulation uncertainty of a power estimate through to the decision
  that uses it. `trend_power` now reports `power` and `power_mcse` at full
  precision (they were rounded to 3 and 4 places before `min_detectable_effect`
  and `multisource_robustness` consumed them) and computes `power_mcse` with the
  Jeffreys-stabilised proportion formula, so a run of 0% or 100% rejections
  still reports a positive error; rows also carry `n_rep`/`n_failed`.
  `min_detectable_effect` drops non-finite power rows from the interpolation and
  adds `emd_*_resolved` (the crossing survives pulling every estimate down two
  MCSE), `emd_*_reps_without_crossing` and `n_power_failed`.
  `multisource_robustness` requires a non-significant source to clear
  `power_threshold` by two MCSE before it is classified
  `not_reproduced_with_power` (otherwise `not_resolved`), and reports
  `power_mcse_for_reference` on each `SourceResult`.
- Make `nsevt.calibration` decisions and failure handling consistent.
  `rejection_rate` now uses one threshold (`alpha + anticonservative_margin`,
  returned as `anticonservative_threshold`) for both the stopping rule and the
  reported `anticonservative` verdict, instead of stopping on the margin and
  reporting on `alpha + 2*MCSE`. `coverage` treats an interval with an infinite
  limit that brackets the target as covering (with `n_infinite`) and a `NaN`
  limit, `lo > hi`, or a raising `estimator` as a failed replicate
  (`n_failed`, `n_effective`), not as non-coverage; `bias_rmse` and
  `pseudo_true` likewise treat a raising estimator as a failure rather than
  aborting the campaign. All four validate `n`, and `coverage` also validates
  `level` and a finite `target`, before simulating.
- Compute the grouped (interval-censored) cell probabilities in log space. The
  cell was formed as `S(a) - S(b)` in ordinary scale and cancelled to zero — so
  the fit reported the penalty value — when the two survivals were close: a wide
  or covariate-dependent scale, a narrow cell, or a shape near zero (which was
  also nudged off zero by `1e-10`). `nsevt.grouped` and `nsevt.design` now share
  `_grouped_loglik`, which uses `log1p`/`expm1` and the exact `xi -> 0` limit
  `-z / sigma`, so the likelihood is continuous through zero shape and matches a
  direct high-precision evaluation.
- Give `nsevt.design.return_level` and `nsevt.GPDFit.return_level` one policy
  outside the tail. For `m * rate <= 1` the target is a non-exceedance quantile
  at or below the threshold; `design.return_level` extrapolated below it
  (returning 36.58 for a period the model cannot represent) while
  `GPDFit.return_level` clipped to the threshold. Both now return the threshold
  there, agree on the common domain, and emit a `RuntimeWarning` when strictly
  below the boundary (`m * rate < 1`).

### Documentation

- Add `docs/stability.md` as the single authoritative list of which symbols the
  semantic-versioning guarantee covers, the stable-core vs `nsevt.experimental`
  split, and the minor/major rule. It supersedes any other stability statement,
  including the manuscript, when they disagree. `docs/index.md` and
  `docs/api.md` now point to it and no longer cite a stale contract version.
- Record the immutable Zenodo DOI assigned to version 1.1.0.

## [1.1.0] - 2026-09-06

### Changed

- Update artifact transfer actions to audited Node.js 24 revisions, removing
  the deprecation annotations observed while publishing version 1.0.3.
- Update the checkout and Python setup actions to their audited Node.js 24
  revisions after full-matrix validation.

### Fixed

- Return an infinite one-sided bound from `split_conformal` / `block_conformal`
  when the calibration sample is too small for the requested `alpha`
  (`ceil((1 - alpha)(m + 1)) > m`). The quantile level was clipped to one and
  the sample maximum returned, whose coverage is only `m / (m + 1)` (0.80 for
  `n=4, alpha=0.1`, not 0.90). `q_standardized` is now `inf` with
  `ConformalBand.underpowered = True`, and the finite-`n` bound is the exact
  `ceil((1 - alpha)(m + 1))`-th order statistic.
- Keep the time labels in `block_bootstrap_trend_ci`. Each resample reassigned
  its sampled blocks to fresh ordered time positions, which mixed scale levels
  from different years and pulled every interval toward zero: a fit of
  `+0.265`/decade returned a bootstrap interval of `[-0.10, 0.08]`. The
  resample now keeps each block's own time (a cluster bootstrap); resamples
  spanning fewer than two distinct times carry no trend information and are
  reported in `n_unidentified` / `n_failed` rather than silently dropped or
  fitted. The result adds `n_boot_requested`, `n_unidentified`, `n_failed`.
- Make the grouped shape and endpoint profile intervals consistent with the
  fit they describe. `profile_ci_xi_grouped` searches outward from `xi_hat`
  instead of inverting a fixed `[-0.95, 0.60]` bracket, so the interval always
  contains the estimate; the grouped kernel now stays in the GPD existence
  region `xi > -1`, matching the continuous fit. `profile_endpoint_ci` grows
  its upper search adaptively and reports `inf` (with `upper_at_bound`) when the
  endpoint is not identified from above at the requested level or the shape is
  not negative, instead of a finite number that is only the search boundary.
  `GroupedGPDFit` gains `level`, `xi_ci`, `endpoint_ci`, `xi_ci_at_bound` and
  `endpoint_ci_at_bound`; `xi_ci95` / `endpoint_ci95` stay as read-only aliases,
  and the summary text and `bounded_supported` now report the level actually
  requested rather than a hard-coded 95%.
- Stop the sequential Monte Carlo protocol from declaring an invalid run
  converged. `SequentialRun` now records failed replicates as `NaN` and stops
  on `n_effective` (finite outcomes) rather than the attempted count, so a run
  padded with failures no longer satisfies `r_min`; a `proportion` run rejects
  any non-0/1 outcome instead of letting the MCSE formula clip it into range;
  `permutation_pvalue` raises on a non-finite `t_obs` or an empty finite null
  instead of returning the resolution floor; and `rejection_rate` counts a
  test that raises or returns an out-of-range p-value as a failed replicate
  rather than a silent non-rejection. `summary()` and `rejection_rate` now
  report `n_attempted` / `n_effective` / `n_failed`.
- Derive the grouped left-truncation point from the threshold and the grid
  instead of assuming the threshold falls on the grid. A mark enters the
  interval-censored sample because it is a grid point above the threshold, so
  the conditioned excess is the lower rounding edge of the smallest such mark,
  not half the grid width. With a shifted threshold the old assumption made the
  conditioning denominator too small and a cell's conditional probability could
  exceed one; aligned thresholds are unaffected.
- Reject non-finite input to the grouped GPD fit before it reaches the
  optimizer. A NaN or infinite grid width, a non-finite threshold, or a
  malformed user-supplied `cells` triple now raises `ValueError` instead of
  letting the Nelder-Mead search "converge" on the constant penalty region and
  return that penalty as a fit. `fit_gpd_grouped` also raises `RuntimeError`
  when no start reaches a feasible optimum, so a penalty objective value is
  never reported as a log-likelihood.

### Documentation

- Record the immutable Zenodo DOI assigned to version 1.0.3.
- Serve the monthly-download badge directly from PePy to avoid intermittent
  upstream rate-limit responses from the previous badge proxy.
- Clarify the single-maintainer authorship policy: external reports and
  proposals remain welcome, while accepted bot or external changes are
  independently reproduced in maintainer-authored commits.

## [1.0.3] - 2026-08-26

### Changed

- Make `src/nsevt/_version.py` the package version source and verify that the
  README, changelog, citation metadata, Zenodo metadata, Git tag, and built
  distributions agree with it before a release can be published.
- Prepare every tagged release as a draft containing the tested wheel, source
  distribution, and SHA-256 manifest; publish those exact artifacts to PyPI
  through Trusted Publishing after the draft is released.
- Add a separate manual TestPyPI workflow and a maintained release checklist.
- Pin release-critical GitHub Actions to audited Node.js 24 revisions so CI no
  longer relies on deprecated Node.js 20 action runtimes.

### Fixed

- Remove `skip-existing` from the PyPI upload. A repeated version now fails
  visibly instead of silently accepting files whose identity was not checked.

### Documentation

- Document the distinction between the all-versions Zenodo concept DOI and the
  immutable DOI assigned to each software version.

There is no change to statistical calculations or the stable public API.

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
