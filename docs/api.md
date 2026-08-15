# API guide

## Stable core

### `gpd_pot(values, threshold, n_boot=2000, seed=...)`

Accepts raw values and fits positive excesses. The returned `GPDFit` separates
`bounded_estimate` from `bounded_supported` and provides profile and conditional
bootstrap summaries.

### `trend_permutation(z, block, n_perm=3000, seed=...)`

Accepts excesses and matching numeric block times. It returns the signed
log-scale trend per decade, LR statistic, asymptotic comparison p-value,
plus-one block-label permutation p-value, Monte Carlo standard error, successful
refit count, and null draws. If any requested refit fails, a `RuntimeWarning`
makes explicit that the p-value is conditional on the successful refits.

### `min_detectable_effect(z, block, direction="both", ...)`

Computes conditional simulation power on a user-supplied or default magnitude
grid. `mde_positive` and `mde_negative` are the smallest grid points reaching the
target power; `emd_positive` and `emd_negative` are crossings of a monotone
interpolant of the power curve. They are not restricted to grid nodes but still
depend on the supplied grid, the interpolation, and Monte Carlo error;
uncertainty intervals are returned as `emd_positive_ci95` and
`emd_negative_ci95`. Those intervals perturb pointwise power estimates by their
simulation errors; they are approximate and do not represent the covariance
induced by common random numbers across the curve.
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
At an observed proportion of exactly zero or one, `mcse_proportion` uses a
Jeffreys half-count so a finite run does not claim zero simulation error.
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

Invalid kinds, quantiles, tolerances and replicate budgets are rejected before
the first draw. Each `draw(k, block_index)` result must contain exactly `k`
values, and quantile runs compare block quantiles in their batch diagnostic.

### `nsevt.calibration` — finite-sample calibration

Given a `simulate(rng, n)` callable that draws one data set from your DGP and
your own estimator or test, this module measures finite-sample properties.
`calibration.rejection_rate(test, simulate, n, alpha=...)` returns the empirical
rejection probability of `test(sample) -> p` (the type-I error under a null DGP,
power under an alternative). `calibration.coverage(estimator, simulate, n,
target, level=...)` returns the empirical coverage of an interval
`estimator(sample) -> (lo, hi)` against `target` (a float or a zero-argument
callable). `calibration.bias_rmse(estimator, simulate, n, truth, n_rep=...)`
returns the bias, standard deviation and RMSE of a point estimator.
`calibration.pseudo_true(estimator, simulate, R=...)` computes a reproducible
large-sample proxy for the estimator's limiting target under a misspecified
DGP. Convergence of that proxy should be checked over increasing `R` and more
than one seed before it is used as a coverage or bias target. The two proportion analyses accept and
forward the `nsevt.mc` sequential parameters (`epsilon`, `r0`, `r_min`, `r_max`,
`block`, `min_stable_blocks`, `seed`); set them small when the estimator is
expensive, since each replicate refits the model.

### `nsevt.design` — grouped regression and return levels

`design.fit_grouped_design(values, threshold, design, grid=5.0, cells=None)` fits
the interval-censored GPD with a log-linear scale `design[i] @ coef` (intercept in
column 0); `values` are the exceedance marks and `design` the aligned `(n, p)`
matrix. `design.profile_ci_coef(values, threshold, design, coef, ...)` is a
profile-likelihood interval for column `coef` (the shape and other coefficients
maximised out). `design.return_level(xi, sigma, threshold, rate, m)` is the level
exceeded once per `m` observations at exceedance rate `rate`, and
`design.profile_ci_return_level(values, threshold, rate, m, ...)` its profile
interval (bounded tails only; `None` when the point level is not finite and
positive). The covariate coding and the exceedance rate are caller-supplied.

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
