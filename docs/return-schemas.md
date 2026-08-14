# Return schemas (stable core)

This is the return contract of the stable public API (`nsevt.gpd`,
`nsevt.grouped`, `nsevt.trend`, `nsevt.transportability`). Under the planned 1.0
stability guarantee, **removing or renaming** a documented key or field is a
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
Methods: `return_level(return_period, rate=None) -> float`, `summary() -> str`.

## Grouped (`nsevt.grouped`)

### `interval_cells(values, threshold, grid=5.0, tol=1e-6) -> tuple`
`(a, b, trunc)` arrays on the excess scale (cell lower/upper edge and
per-observation left-truncation point).

### `fit_gpd_grouped(values, threshold, grid=5.0, cells=None, starts=...) -> dict`
`{"xi", "sigma", "endpoint", "loglik", "n"}`.

### `profile_ci_xi_grouped(...) -> dict`
`{"xi_hat", "ci": (lo, hi), "lo_at_bound": bool, "hi_at_bound": bool, "level"}`.

### `profile_endpoint_ci(...) -> dict`
`{"endpoint", "ci": (lo, hi), "upper_at_bound": bool, "xi", "sigma", "level",
"method": str}`.

### `gpd_pot_grouped(values, threshold, grid=5.0, level=0.95) -> GroupedGPDFit`
Dataclass fields: `threshold, n_exceedances, xi, sigma, xi_ci95 (lo, hi),
endpoint, endpoint_ci95 (lo, hi), loglik`.
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
- `_null`: internal null-statistic array (not part of the contract).

### `trend_power(z, block, trends, ...) -> list[dict]`
One row per trend: `{"trend_per_decade", "sigma_change_pct", "power",
"power_mcse", "n_successful"}`.

### `min_detectable_effect(z, block, ...) -> dict`
- Grid-based: `mde_per_decade`, `mde_absolute`, `mde_negative`, `mde_positive`.
- Interpolated: `emd_per_decade`, `emd_negative`, `emd_positive`, each with a
  Monte Carlo interval `emd_negative_ci95` / `emd_positive_ci95` (`[lo, hi]` or
  `None`).
- `direction`, `target_power`, and `power_curve` (a `trend_power` list).

### `block_bootstrap_trend_ci(z, block, n_boot=1000, seed=..., ...) -> dict`
`{"ci95": [lo, hi] (or [None, None]), "n_boot": int}`.

## Sequential Monte Carlo precision (`nsevt.mc`)

### `mcse_proportion(p_hat, R) -> float`, `mcse_mean(values) -> float`, `mcse_quantile(values, q) -> float`
Scalar Monte Carlo standard error; `mcse_mean`/`mcse_quantile` return `nan`
below their minimum sample size (2 and 30 finite values).

### `required_replicates(p_hat, epsilon) -> int`
Replicate budget `ceil(p (1 - p) / epsilon^2)` for `MCSE <= epsilon`.

### `permutation_pvalue(t_obs, t_null, plus_one=True) -> dict`
- `p`: plus-one Monte Carlo p-value (never zero); `n_exceed`, `B`.
- `floor`: resolution floor `1 / (B + 1)`; `at_floor`: bool.
- `mcse`: Monte Carlo standard error of `p`.

### `SequentialRun.summary() -> dict` (also the result of `run_sequential`)
- `analysis`, `kind`, `R_star`, `estimate`, `mcse`, `tolerance`.
- `last_batch_change`, `stability_tolerance`, `n_stable_blocks`,
  `min_stable_blocks`, `decision_stable`, `seed`.
- `status`: `"converged"` or `"not_stabilised"`.
- `trace`: list of `{"R", "value", "mcse"}` at the pre-specified checkpoints.
- `batch_diagnostic`: `{"n_blocks", "block_size", "observed_sd_between_blocks",
  "theoretical_mcse_per_block", "ratio", "block_means"}` (or
  `{"n_blocks", "ratio": None}` when there are too few blocks).

### `substream(seed, *tags) -> numpy.random.Generator`, `block_streams(seed, n_blocks, *tags) -> list`
Reproducible, order-independent generators; `block_streams` returns one per
sequential block.

### `multiseed_summary(values_by_seed) -> dict`
`{"n_seeds", "labels", "values", "mean", "sd_across_seeds", "range"}`.

## Multi-source robustness (`nsevt.transportability`)

### `multisource_robustness(sources, threshold, reference=None, ...) -> ArenaResult`
`ArenaResult` fields: `sources (list[SourceResult]), alpha, shape_bounded_all,
trend_reproduces, reference_source, trend_status, verdict`.
Method: `table() -> str`. `trend_status` is one of `reproduced`,
`inconsistent_direction`, `not_reproduced_with_power`, `not_resolved`,
`no_reference_signal`, `single_source_only`.

`SourceResult` fields: `name, n, xi, xi_ci95, bounded_estimate,
bounded_supported, endpoint, trend_per_decade, p_permutation,
trend_significant, trend_direction, power_for_reference`. Property: `bounded`.

The name `transportability` is a backward-compatible alias of
`multisource_robustness`.
