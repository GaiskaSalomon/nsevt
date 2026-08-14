"""Experimental and assumptions-dependent APIs, outside the stable core.

These routines are provided for evaluation. Their finite-sample guarantees are
either not established (`block_conformal`, `twoscale_trend`,
`wasserstein_decomposition`) or hold only under stated assumptions
(`split_conformal`). They may change without a major version bump and are **not**
covered by the package's public-API stability guarantee.

Prefer the stable core for load-bearing inference:
`nsevt.gpd`, `nsevt.grouped`, `nsevt.trend`, and `nsevt.transportability`.

The same names remain importable from the top level for backward compatibility,
but `nsevt.experimental` is their documented home.
"""
from ..conformal import ConformalBand, block_conformal, split_conformal
from ..twoscale import TwoScaleResult, twoscale_trend, wasserstein_decomposition

__all__ = [
    "block_conformal",
    "split_conformal",
    "ConformalBand",
    "twoscale_trend",
    "wasserstein_decomposition",
    "TwoScaleResult",
]
