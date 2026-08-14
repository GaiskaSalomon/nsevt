# API guide

## Stable core

### `gpd_pot(values, threshold, n_boot=2000, seed=...)`

Accepts raw values and fits positive excesses. The returned `GPDFit` separates
`bounded_estimate` from `bounded_supported` and provides profile and conditional
bootstrap summaries.

### `trend_permutation(z, block, n_perm=3000, seed=...)`

Accepts excesses and matching numeric block times. It returns the signed
log-scale trend per decade, LR statistic, asymptotic comparison p-value,
plus-one block-label permutation p-value, Monte Carlo standard error, and null
draws.

### `min_detectable_effect(z, block, direction="both", ...)`

Computes conditional simulation power on a user-supplied or default magnitude
grid. `mde_positive` and `mde_negative` are the smallest grid points reaching the
target power; `emd_positive` and `emd_negative` are the crossings of a monotone
interpolant of the power curve, which do not depend on the grid spacing and come
with Monte Carlo uncertainty intervals (`emd_positive_ci95`, `emd_negative_ci95`).
`mde_per_decade` is retained as the smallest signed grid effect by absolute
magnitude for compatibility.

### `multisource_robustness(sources, threshold, reference=..., ...)`

Each source is `(name, raw_values, block_times)`. The same threshold and
procedures are applied to every source. The older function name
`transportability` remains a compatibility alias.

### `gpd_pot_grouped(values, threshold, grid=5.0, level=0.95)`

Fits an interval-censored GPD to exceedances recorded on a discrete grid (a
single precision, or several widths for a mixed-precision record). The returned
`GroupedGPDFit` carries the shape and the finite endpoint, each with a
profile-likelihood interval; `interval_cells`, `fit_gpd_grouped`,
`profile_ci_xi_grouped`, and `profile_endpoint_ci` expose the individual pieces.
Prefer this over `gpd_pot` when the values are rounded, because fitting the
continuous GPD to rounded data biases the shape and the endpoint. The endpoint
interval is a profile interval on the reparameterised endpoint, not a percentile
bootstrap.

### `nsevt.mc` — sequential Monte Carlo precision

`mc.mcse_proportion(p_hat, R)`, `mc.mcse_mean(values)` and
`mc.mcse_quantile(values, q)` give the Monte Carlo standard error of a
proportion, a mean and a sample quantile; `mc.required_replicates(p_hat,
epsilon)` inverts the proportion formula to a replicate budget.
`mc.permutation_pvalue(t_obs, t_null)` returns the floor-aware
`(1 + #exceed) / (B + 1)` p-value with `at_floor` and its MCSE.

`mc.run_sequential(name, draw, *, kind, epsilon, ...)` drives a
`mc.SequentialRun`: it seeds the run with `r0` replicate outcomes from
`draw(k, block_index)`, then grows it in blocks of `block` (never below `r_min`,
never past `r_max`) and stops only when the MCSE is at or below `epsilon`, the
estimate has been stable for `min_stable_blocks` checkpoints, and every callable
in `decision_rules` has returned the same qualitative decision across that
window. `kind` selects the MCSE formula (`"proportion"`, `"mean"`,
`"quantile"`). The run exposes `summary()` (with `R_star`, `status`
—`converged` or `not_stabilised`— `trace()`, and an independent-block
`batch_diagnostic()`). `mc.substream(seed, *tags)` and `mc.block_streams(seed,
n_blocks, *tags)` build reproducible, non-interfering random streams so that
extending a run never perturbs the replicates already drawn, and
`mc.multiseed_summary(values_by_seed)` audits an estimate across independent
streams.

## Experimental APIs (`nsevt.experimental`)

`split_conformal`, `block_conformal`, `ConformalBand`, `twoscale_trend`,
`wasserstein_decomposition`, and `TwoScaleResult` are collected under
`nsevt.experimental`. `split_conformal` is stable only under its documented
exchangeability and score-fitting assumptions; the others make narrower claims
whose finite-sample guarantees are not established. They remain importable from
the top level for backward compatibility, but are excluded from the public-API
stability guarantee and may change without a major version bump.

All randomized public functions accept an explicit seed. See function
docstrings for argument-level validation and returned fields, and
[`return-schemas.md`](return-schemas.md) for the full return contract of the
stable core (the surface that the 1.0 stability guarantee will cover).
