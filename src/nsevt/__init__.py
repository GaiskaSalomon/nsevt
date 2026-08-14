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
  exceedances, with profile-likelihood intervals for the shape and the endpoint.

The conformal block aggregation and distribution-valued trend modules remain
experimental and make narrower claims than the stable inferential core.
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
from .transportability import ArenaResult, SourceResult, multisource_robustness, transportability
from .trend import block_bootstrap_trend_ci, min_detectable_effect, trend_permutation, trend_power
from .twoscale import TwoScaleResult, twoscale_trend, wasserstein_decomposition

__version__ = "0.2.0"

__all__ = [
    "fit_gpd", "profile_ci_xi", "upper_endpoint", "gpd_pot", "GPDFit",
    "gpd_pot_grouped", "fit_gpd_grouped", "interval_cells",
    "profile_ci_xi_grouped", "profile_endpoint_ci", "GroupedGPDFit",
    "trend_permutation", "trend_power", "min_detectable_effect",
    "block_bootstrap_trend_ci",
    "multisource_robustness", "transportability", "ArenaResult", "SourceResult",
    "block_conformal", "split_conformal", "ConformalBand",
    "twoscale_trend", "wasserstein_decomposition", "TwoScaleResult",
    "__version__",
]
