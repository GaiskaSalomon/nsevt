"""nsevt: non-stationary extreme-value tail inference and robustness.

A small, dependency-light toolkit for calibrated analysis of trends in
environmental extremes:

* :func:`gpd_pot` -- peaks-over-threshold GPD fit with a profile-likelihood
  interval for the shape and a bootstrap of the finite upper endpoint;
* :func:`trend_permutation` / :func:`min_detectable_effect` -- a
  permutation-calibrated test of a tail-scale trend and a Monte-Carlo
  power/MDE analysis (what can this record actually resolve?);
* :func:`multisource_robustness` -- apply the same specification across sources
  and distinguish disagreement from limited power;
* :func:`gpd_pot_grouped` -- interval-censored (grouped) GPD fit for discretised
  exceedances, with profile-likelihood intervals for the shape and the endpoint;
* :mod:`nsevt.mc` -- a sequential Monte Carlo precision protocol (MCSE, an
  MCSE/stability/decision stopping rule and reproducible traces) that decides
  how many replicates a power, coverage, p-value or bootstrap actually needs.

The conformal block aggregation and distribution-valued trend routines are
experimental and are grouped under :mod:`nsevt.experimental`; they are excluded
from the public-API stability guarantee. They remain importable from the top
level for backward compatibility.
"""
from .conformal import ConformalBand, block_conformal, split_conformal
from .gpd import GPDFit, fit_gpd, gpd_pot, profile_ci_xi, upper_endpoint
from .grouped import (
    GroupedGPDFit,
    fit_gpd_grouped,
    gpd_pot_grouped,
    interval_cells,
    profile_ci_xi_grouped,
    profile_endpoint_ci,
)
from .mc import (
    Checkpoint,
    SequentialRun,
    block_streams,
    mcse_mean,
    mcse_proportion,
    mcse_quantile,
    multiseed_summary,
    permutation_pvalue,
    required_replicates,
    run_sequential,
    substream,
)
from .transportability import ArenaResult, SourceResult, multisource_robustness, transportability
from .trend import block_bootstrap_trend_ci, min_detectable_effect, trend_permutation, trend_power
from .twoscale import TwoScaleResult, twoscale_trend, wasserstein_decomposition

__version__ = "0.4.0"

__all__ = [
    "fit_gpd", "profile_ci_xi", "upper_endpoint", "gpd_pot", "GPDFit",
    "gpd_pot_grouped", "fit_gpd_grouped", "interval_cells",
    "profile_ci_xi_grouped", "profile_endpoint_ci", "GroupedGPDFit",
    "trend_permutation", "trend_power", "min_detectable_effect",
    "block_bootstrap_trend_ci",
    "mcse_proportion", "mcse_mean", "mcse_quantile", "required_replicates",
    "permutation_pvalue", "SequentialRun", "Checkpoint", "run_sequential",
    "substream", "block_streams", "multiseed_summary",
    "multisource_robustness", "transportability", "ArenaResult", "SourceResult",
    "block_conformal", "split_conformal", "ConformalBand",
    "twoscale_trend", "wasserstein_decomposition", "TwoScaleResult",
    "__version__",
]
