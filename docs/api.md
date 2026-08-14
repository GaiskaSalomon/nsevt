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
grid. Use `mde_positive` and `mde_negative`; `mde_per_decade` is retained as the
smallest signed effect by absolute magnitude for compatibility.

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

## Additional APIs

`split_conformal` is stable only under its documented exchangeability and
score-fitting assumptions. `block_conformal`, `twoscale_trend`, and
`wasserstein_decomposition` are marked experimental in their results or
documentation.

All randomized public functions accept an explicit seed. See function
docstrings for argument-level validation and returned fields.
