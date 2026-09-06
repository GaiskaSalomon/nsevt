"""Typed schemas for the return dictionaries of the stable core.

These mirror the frozen contract in ``docs/return-schemas.md``. They are hints
only -- the functions still return plain ``dict`` at runtime -- so a typed client
that reads a wrong key or assumes a wrong type is caught by a type checker
without changing any behaviour. ``total=False`` marks keys that are only present
in some branches.
"""
from __future__ import annotations

from typing import TypedDict


class PermutationPValue(TypedDict):
    """Return of :func:`nsevt.mc.permutation_pvalue`."""

    p: float
    n_exceed: int
    B: int
    n_null_dropped: int
    floor: float
    at_floor: bool
    mcse: float


class GroupedShapeCI(TypedDict):
    """Return of :func:`nsevt.grouped.profile_ci_xi_grouped`."""

    xi_hat: float
    ci: tuple[float, float]
    lo_at_bound: bool
    hi_at_bound: bool
    level: float


class GroupedEndpointCI(TypedDict):
    """Return of :func:`nsevt.grouped.profile_endpoint_ci`."""

    endpoint: float
    ci: tuple[float, float]
    upper_at_bound: bool
    lower_at_bound: bool
    xi: float
    sigma: float
    level: float
    method: str


class RejectionRate(TypedDict):
    """Return of :func:`nsevt.calibration.rejection_rate`."""

    rate: float
    mcse: float
    alpha: float
    anticonservative_threshold: float
    n: int
    R: int
    n_effective: int
    n_failed: int
    status: str
    anticonservative: bool
    stopping: dict


class Coverage(TypedDict):
    """Return of :func:`nsevt.calibration.coverage`."""

    coverage: float
    mcse: float
    nominal: float
    target: float
    n: int
    R: int
    n_effective: int
    n_failed: int
    n_infinite: int
    status: str
    miscalibration: float
    stopping: dict


class BiasRMSE(TypedDict):
    """Return of :func:`nsevt.calibration.bias_rmse`."""

    bias: float
    bias_mcse: float
    sd: float
    rmse: float
    rmse_mcse: float
    mean_estimate: float
    truth: float
    n: int
    n_rep: int
    n_failed: int


class PseudoTrue(TypedDict):
    """Return of :func:`nsevt.calibration.pseudo_true`."""

    pseudo_true: float
    R: int


class PowerRow(TypedDict):
    """One row of :func:`nsevt.trend.trend_power`."""

    trend_per_decade: float
    sigma_change_pct: float
    power: float
    power_mcse: float
    n_rep: int
    n_successful: int
    n_failed: int


class GPDGroupedFitDict(TypedDict):
    """Return of :func:`nsevt.grouped.fit_gpd_grouped`."""

    xi: float
    sigma: float
    endpoint: float
    loglik: float
    n: int


__all__ = [
    "PermutationPValue",
    "GroupedShapeCI",
    "GroupedEndpointCI",
    "RejectionRate",
    "Coverage",
    "BiasRMSE",
    "PseudoTrue",
    "PowerRow",
    "GPDGroupedFitDict",
]
