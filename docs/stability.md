# API stability

This is the single source of truth for which symbols are covered by the
stability guarantee and which are not. It supersedes any stability statement
elsewhere in the documentation or the manuscript when they disagree.

## The guarantee

The **stable core** has followed semantic versioning since `1.0.0`
(published 2026-08-15):

- **patch** (`1.1.0` → `1.1.1`): a defect is corrected relative to the
  documented behaviour. Numeric outputs may change for the affected case; the
  change and its claim boundary are recorded in `CHANGELOG.md`.
- **minor** (`1.1.0` → `1.2.0`): backward-compatible additions — a new function,
  a new optional argument, a new key in a returned `dict`, a new dataclass
  field. Existing names, signatures, keys and fields keep working.
- **major** (`1.x` → `2.0.0`): a documented name, signature, return key or
  dataclass field is removed or renamed, or the meaning of a documented result
  is deliberately changed in an incompatible way.

The frozen return contract of the stable core is `docs/return-schemas.md`; the
assumptions and claim boundaries are `docs/assumptions.md`.

## Stable core

Covered by the guarantee above. Importable from the top level (`nsevt.<name>`)
and from the module shown.

| Module | Public symbols |
| --- | --- |
| `nsevt.gpd` | `gpd_pot`, `GPDFit`, `fit_gpd`, `profile_ci_xi`, `upper_endpoint` |
| `nsevt.grouped` | `gpd_pot_grouped`, `GroupedGPDFit`, `fit_gpd_grouped`, `interval_cells`, `profile_ci_xi_grouped`, `profile_endpoint_ci` |
| `nsevt.trend` | `trend_permutation`, `trend_power`, `min_detectable_effect`, `block_bootstrap_trend_ci` |
| `nsevt.transportability` | `multisource_robustness`, `transportability` (permanent alias), `ArenaResult`, `SourceResult` |
| `nsevt.mc` | `mcse_proportion`, `mcse_mean`, `mcse_quantile`, `required_replicates`, `permutation_pvalue`, `SequentialRun`, `Checkpoint`, `run_sequential`, `substream`, `block_streams`, `multiseed_summary` |
| `nsevt.calibration` | `rejection_rate`, `coverage`, `bias_rmse`, `pseudo_true` |
| `nsevt.design` | `fit_grouped_design`, `profile_ci_coef`, `return_level`, `profile_ci_return_level` |

`GPDFit` fields/methods: `threshold, n_exceedances, xi, sigma, xi_ci95,
endpoint, endpoint_ci95, bootstrap_fraction_xi_negative, n_boot_successful,
xi_ci95_truncated`; `bounded_estimate`, `bounded_supported`, `bounded`,
`return_level(...)`, `summary()`.

`GroupedGPDFit` fields/methods: `threshold, n_exceedances, xi, sigma, xi_ci,
endpoint, endpoint_ci, loglik, level, xi_ci_at_bound, endpoint_ci_at_bound`;
`xi_ci95` / `endpoint_ci95` are permanent read-only aliases of `xi_ci` /
`endpoint_ci`; `bounded_supported`, `summary()`.

Additions since `1.0.0` (all minor, backward compatible): `n_attempted`,
`n_effective`, `n_failed` in `SequentialRun.summary()` and `rejection_rate`;
`n_null_dropped` in `permutation_pvalue`; `n_boot_requested`, `n_unidentified`,
`n_failed` in `block_bootstrap_trend_ci`; `lower_at_bound` in
`profile_endpoint_ci`; `level`, `xi_ci`, `endpoint_ci`, `xi_ci_at_bound`,
`endpoint_ci_at_bound` on `GroupedGPDFit`; `underpowered` on `ConformalBand`.

## Experimental

**Not** covered by the guarantee. Collected under `nsevt.experimental`; the
same names stay importable from the top level as silent backward-compatible
aliases (no `DeprecationWarning`). These may change without a major version
bump.

| Module | Public symbols | Status |
| --- | --- | --- |
| `nsevt.conformal` | `split_conformal`, `ConformalBand` | Marginal finite-sample guarantee **only** under the documented exchangeability and independent-score-construction assumptions. |
| `nsevt.conformal` | `block_conformal` | Block-aggregation diagnostic for ordered dependent exceedances. No general finite-sample or beta-mixing guarantee. |
| `nsevt.twoscale` | `twoscale_trend`, `wasserstein_decomposition`, `TwoScaleResult` | Residual circular moving-block bootstrap diagnostic for empirical quantile functions. Guarantees not established. |

`nsevt.conformal` is also exempt from the mypy gate (`[tool.mypy]` override).

## Notes

- The manuscript in `paper/` is a point-in-time description and is not
  authoritative for the current API; this file is.
- `ROADMAP.md` and `API_FREEZE_1.0.md` are maintainer working notes kept out of
  version control on purpose; the durable decisions they record that affect
  users are reflected here and in `CHANGELOG.md`.
