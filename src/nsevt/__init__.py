"""nsevt: non-stationary extreme-value tail risk with honest uncertainty.

A small, dependency-light toolkit for the *honest* analysis of trends in
environmental extremes:

* :func:`gpd_pot` -- peaks-over-threshold GPD fit with a profile-likelihood
  interval for the shape and a bootstrap of the finite upper endpoint;
* :func:`trend_permutation` / :func:`min_detectable_effect` -- a
  permutation-calibrated test of a tail-scale trend and a Monte-Carlo
  power/MDE analysis (what can this record actually resolve?);
* :func:`transportability` -- the multi-source "evidence arena": does a trend
  survive changing the data source, or is it an artifact?
* :func:`block_conformal` -- distribution-free prediction bands for extreme
  quantiles with coverage guarantees under temporal dependence;
* :func:`twoscale_trend` -- a Wasserstein/Frechet-mean trend test for a series
  of distributions estimated from small per-period samples.

The numerics of the GPD and permutation machinery are ported verbatim from the
frozen, unit-tested code of the accompanying research, so results are identical.
"""
from .gpd import fit_gpd, profile_ci_xi, upper_endpoint, gpd_pot, GPDFit
from .trend import (trend_permutation, trend_power, min_detectable_effect,
                    block_bootstrap_trend_ci)
from .transportability import transportability, ArenaResult, SourceResult
from .conformal import block_conformal, split_conformal, ConformalBand
from .twoscale import twoscale_trend, wasserstein_decomposition, TwoScaleResult

__version__ = "0.1.0"

__all__ = [
    "fit_gpd", "profile_ci_xi", "upper_endpoint", "gpd_pot", "GPDFit",
    "trend_permutation", "trend_power", "min_detectable_effect",
    "block_bootstrap_trend_ci",
    "transportability", "ArenaResult", "SourceResult",
    "block_conformal", "split_conformal", "ConformalBand",
    "twoscale_trend", "wasserstein_decomposition", "TwoScaleResult",
    "__version__",
]
