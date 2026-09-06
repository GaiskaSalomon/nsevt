# Return schemas (stable core)

This is the return contract of the stable public API (`nsevt.gpd`,
`nsevt.grouped`, `nsevt.trend`, `nsevt.transportability`, `nsevt.mc`,
`nsevt.calibration`, and `nsevt.design`). Under the 1.0 stability guarantee,
**removing or renaming** a documented key or field is a
breaking change (2.0.0); **adding** keys or fields is a minor change (1.x).
Keys beginning with an underscore are internal and not part of the contract.

## GPD (`nsevt.gpd`)

### `fit_gpd(z) -> dict`
`{"xi": float, "sigma": float, "nll": float}`.

### `profile_ci_xi(z, level=0.95, grid=None) -> tuple`
`(xi_hat: float, sigma_hat: float, (lo: float, hi: float))`. A `RuntimeWarning`
is emitted when `xi_hat <= -0.5` or a limit reaches the search boundary.

### `upper_endpoint(z, threshold, n_boot=2000, seed=...) -> dict`
- `xi`, `sigma`, `endpoint`: float (`endpoint = inf` when `xi >= 0`).
- `endpoint_ci`: `[lo, hi]` percentile interval among bootstrap fits with
  `xi < 0` (`[nan, nan]` if none).
- `bootstrap_fraction_xi_negative`: float, `n_boot_successful`: int.

### `gpd_pot(values, threshold, n_boot=2000, seed=...) -> GPDFit`
Dataclass fields: `threshold, n_exceedances, xi, sigma, xi_ci95 (lo, hi),
endpoint, endpoint_ci95 [lo, hi], bootstrap_fraction_xi_negative,
n_boot_successful, xi_ci95_truncated (bool, bool)`.
Properties: `bounded_estimate` (point estimate `xi < 0`), `bounded_supported`
(whole `xi_ci95` below zero), `bounded` (alias of `bounded_supported`).
Methods: `return_level(return_period, rate=None) -> float` (returns `threshold`
and warns when `return_period * rate < 1`, a quantile below the POT domain),
`summary() -> str`.

## Grouped (`nsevt.grouped`)

### `interval_cells(values, threshold, grid=5.0, tol=1e-6) -> tuple`
`(a, b, trunc)` arrays on the excess scale (cell lower/upper edge and
left-truncation point). `trunc` is the excess of the lower rounding edge of the
smallest grid point strictly above the threshold (``g/2`` when the threshold is
on the grid, smaller otherwise). Raises `ValueError` on non-finite values, a
non-finite threshold, or a non-finite/non-positive grid width.

### `fit_gpd_grouped(values, threshold, grid=5.0, cells=None, starts=...) -> dict`
`{"xi", "sigma", "endpoint", "loglik", "n"}`. A user-supplied `cells` triple is
validated (`ValueError` unless finite, aligned with the exceedances, and
`0 <= a < b`, `0 <= trunc <= b`). `RuntimeError` when no start reaches a
feasible optimum, so a penalty objective is never reported as a fit.

### `profile_ci_xi_grouped(..., level=0.95, xi_floor=-0.999, xi_ceil=5.0) -> dict`
`{"xi_hat", "ci": (lo, hi), "lo_at_bound": bool, "hi_at_bound": bool, "level"}`.
The interval is found by searching outward from `xi_hat`, so it always contains
`xi_hat`; `*_at_bound` marks a limit that only reached `xi_floor` / `xi_ceil`.

### `profile_endpoint_ci(..., level=0.95, gap_init=50.0, gap_cap=1e7) -> dict`
`{"endpoint", "ci": (lo, hi), "upper_at_bound": bool, "lower_at_bound": bool,
"xi", "sigma", "level", "method": str}`. `endpoint` and `ci[1]` are `inf` when
the endpoint is not identified from above at `level` (unbounded shape or the
upper search reached `gap_cap`); `ci[0]` stays finite.

### `gpd_pot_grouped(values, threshold, grid=5.0, level=0.95) -> GroupedGPDFit`
Dataclass fields: `threshold, n_exceedances, xi, sigma, xi_ci (lo, hi),
endpoint, endpoint_ci (lo, hi), loglik, level, xi_ci_at_bound (lo, hi),
endpoint_ci_at_bound (lo, hi)`. `xi_ci95` / `endpoint_ci95` remain as read-only
aliases of `xi_ci` / `endpoint_ci`. The intervals are at `level` (not fixed
0.95); `endpoint_ci` upper limit is `inf` when unidentified.
Property: `bounded_supported`. Method: `summary() -> str`.

## Trend, power, MDE (`nsevt.trend`)

### `trend_permutation(z, block, n_perm=3000, seed=..., ref_block=None) -> dict`
- `trend_per_decade`, `xi`, `log_sigma0`: fitted varying model.
- `xi_null`, `log_sigma_null`: fitted constant-scale null.
- `sigma_change_pct`: percent scale change across the observed span.
- `LR`: two-sided likelihood-ratio statistic.
- `p_asymptotic`: chi-square reference p-value.
- `p_permutation`: plus-one permutation p-value; `p_permutation_mcse`: its
  Monte Carlo standard error.
- `n_permutations`: successful permutations; `permutation_unit`: str.
- A `RuntimeWarning` is emitted if fewer than the requested permutations
  succeed; the returned p-value is then conditional on those successful refits.
- `_null`: internal null-statistic array (not part of the contract).

### `trend_power(z, block, trends, ...) -> list[dict]`
One row per trend: `{"trend_per_decade", "sigma_change_pct" (display, rounded),
"power", "power_mcse", "n_rep", "n_successful", "n_failed"}`. `power` and
`power_mcse` are full precision (not rounded); `power_mcse` is the
Jeffreys-stabilised proportion MCSE, so it is positive even at 0% or 100%.

### `min_detectable_effect(z, block, ...) -> dict`
- Grid-based: `mde_per_decade`, `mde_absolute`, `mde_negative`, `mde_positive`.
- Interpolated: `emd_per_decade`, `emd_negative`, `emd_positive`, each with a
  Monte Carlo interval `emd_negative_ci95` / `emd_positive_ci95` (`[lo, hi]` or
  `None`). The interval is a pointwise-normal approximation; it does not model
  covariance across curve points generated with common random numbers.
- `emd_negative_resolved` / `emd_positive_resolved`: bool — the target is still
  reached with every power estimate pulled down two MCSE (the crossing is not an
  artefact of simulation noise).
- `emd_negative_reps_without_crossing` / `emd_positive_reps_without_crossing`:
  perturbed power curves that never reached the target.
- `n_power_failed`: total failed replicates across the power curve.
- `direction`, `target_power`, and `power_curve` (a `trend_power` list). Rows
  with a non-finite `power` are dropped from the interpolation.

### `block_bootstrap_trend_ci(z, block, n_boot=1000, seed=..., ...) -> dict`
`{"ci95": [lo, hi] (or [None, None]), "n_boot": int, "n_boot_requested": int,
"n_unidentified": int, "n_failed": int}`. Cluster bootstrap: whole blocks are
resampled keeping their own time label; resamples spanning fewer than two
distinct times (or with no variation) are counted in `n_unidentified`, refit
failures in `n_failed`, and `n_boot` is the usable count the interval rests on.

## Sequential Monte Carlo precision (`nsevt.mc`)

### `mcse_proportion(p_hat, R) -> float`, `mcse_mean(values) -> float`, `mcse_quantile(values, q) -> float`
Scalar Monte Carlo standard error; `mcse_proportion` uses a Jeffreys half-count
at observed proportions of exactly zero or one, so finite runs do not report
zero simulation error. `mcse_mean`/`mcse_quantile` return `nan` below their
minimum sample size (2 and 30 finite values).

### `required_replicates(p_hat, epsilon) -> int`
Replicate budget for `MCSE <= epsilon`: `ceil(p (1 - p) / epsilon^2)` away
from zero and one, and the exact inversion of the boundary-stabilised MCSE at
an observed boundary.

### `permutation_pvalue(t_obs, t_null, plus_one=True) -> dict`
- `p`: plus-one Monte Carlo p-value (never zero); `n_exceed`, `B`.
- `n_null_dropped`: non-finite entries removed from `t_null`.
- `floor`: resolution floor `1 / (B + 1)`; `at_floor`: bool.
- `mcse`: Monte Carlo standard error of `p`.
- Raises `ValueError` on a non-finite `t_obs` or an empty finite `t_null`.

### `SequentialRun.summary() -> dict` (also the result of `run_sequential`)
- `analysis`, `kind`, `R_star`, `estimate`, `mcse`, `tolerance`.
- `n_attempted` (= `R_star`), `n_effective` (finite outcomes, drive the
  estimate, MCSE and stopping rule), `n_failed` (non-finite outcomes).
- `last_batch_change`, `stability_tolerance`, `n_stable_blocks`,
  `min_stable_blocks`, `decision_stable`, `seed`.
- `status`: `"converged"` or `"not_stabilised"`.
- `trace`: list of `{"R", "value", "mcse"}` at the pre-specified checkpoints.
- `batch_diagnostic`: `{"n_blocks", "block_size", "observed_sd_between_blocks",
  "theoretical_mcse_per_block", "ratio", "block_means"}` (or
  `{"n_blocks", "ratio": None}` when there are too few blocks). The compatibility
  key `block_means` contains block quantiles when `kind="quantile"`.

### `substream(seed, *tags) -> numpy.random.Generator`, `block_streams(seed, n_blocks, *tags) -> list`
Reproducible, order-independent generators; `block_streams` returns one per
sequential block.

### `multiseed_summary(values_by_seed) -> dict`
`{"n_seeds", "labels", "values", "mean", "sd_across_seeds", "range"}`.

## Finite-sample calibration (`nsevt.calibration`)

### `rejection_rate(test, simulate, n, alpha=0.05, ...) -> dict`
- `rate`, `mcse`, `alpha`, `n`, `R`, `status`.
- `anticonservative_threshold`: `alpha + anticonservative_margin`; the run stops
  on this decision and reports it, so the two agree.
- `n_effective`, `n_failed`: replicates the rate rests on, and replicates whose
  `test` raised or returned a value that was non-finite or outside `[0, 1]`.
- `anticonservative`: bool (`rate > anticonservative_threshold`).
- `stopping`: the `SequentialRun.summary()`.

### `coverage(estimator, simulate, n, target, level=0.95, ...) -> dict`
- `coverage`, `mcse`, `nominal`, `target`, `n`, `R`, `status`.
- `n_effective`, `n_failed`: replicates the coverage rests on, and replicates
  where a limit was `NaN`, `lo > hi`, or `estimator` raised (failures, not
  non-coverage).
- `n_infinite`: replicates whose interval had an infinite limit but still
  bracketed `target` (counted as covering, but uninformative).
- `miscalibration`: `coverage - nominal`.
- `stopping`: the `SequentialRun.summary()`.
- Rejects a non-positive `n`, `level` outside `(0, 1)`, or a non-finite `target`
  before simulating.

### `bias_rmse(estimator, simulate, n, truth, n_rep=5000, ...) -> dict`
`{"bias", "bias_mcse", "sd", "rmse", "rmse_mcse", "mean_estimate", "truth",
"n", "n_rep", "n_failed"}`. An `estimator` that raises or returns a non-finite
value is a failed replicate (`n_failed`), not a campaign abort.

### `pseudo_true(estimator, simulate, R=20000, ...) -> dict`
`{"pseudo_true": float, "R": int}`.

## Grouped regression and return levels (`nsevt.design`)

### `fit_grouped_design(values, threshold, design, grid=5.0, cells=None) -> dict`
`{"xi", "coef" (list), "sigma0", "loglik", "n", "p", "truncation"}`.

### `profile_ci_coef(values, threshold, design, coef=1, ...) -> dict`
`{"coef" (int), "estimate", "ci": [lo, hi], "at_bound", "level", "loglik"}`.

### `return_level(xi, sigma, threshold, rate, m) -> ndarray`
The return level(s) for scalar or array `m`. Entries with `m * rate <= 1` are a
non-exceedance quantile outside the peaks-over-threshold model and are returned
as `threshold` (matching `GPDFit.return_level`), with a `RuntimeWarning` when
any entry is strictly below the domain (`m * rate < 1`).

### `profile_ci_return_level(values, threshold, rate, m, ...) -> dict | None`
`{"m", "return_level", "ci": [lo, hi], "upper_at_bound", "xi_at_max", "loglik",
"level"}`, or `None` when the point return level is not finite and positive.

## Multi-source robustness (`nsevt.transportability`)

### `multisource_robustness(sources, threshold, reference=None, ...) -> ArenaResult`
`ArenaResult` fields: `sources (list[SourceResult]), alpha, shape_bounded_all,
trend_reproduces, reference_source, trend_status, verdict`.
Method: `table() -> str`. `trend_status` is one of `reproduced`,
`inconsistent_direction`, `not_reproduced_with_power`, `not_resolved`,
`no_reference_signal`, `single_source_only`. `not_reproduced_with_power`
requires every non-significant source to clear `power_threshold` by two MCSE;
a power that only sits near the threshold falls to `not_resolved`.

`SourceResult` fields: `name, n, xi, xi_ci95, bounded_estimate,
bounded_supported, endpoint, trend_per_decade, p_permutation,
trend_significant, trend_direction, power_for_reference,
power_mcse_for_reference`. Property: `bounded`.

The name `transportability` is a backward-compatible alias of
`multisource_robustness`.
