"""Experimental distribution-valued trend diagnostics.

This module estimates trends in empirical quantile functions and calibrates a
descriptive statistic with a residual circular moving-block bootstrap.  The
procedure is useful for exploratory sensitivity analysis, but nsevt does not
claim an exact finite-sample test or the previously stated two-scale validity
frontier.  It is intentionally outside the package's stable inferential core.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import numpy.typing as npt


def _validate_grid(grid):
    grid = np.asarray(grid, dtype=float)
    if (
        grid.ndim != 1
        or grid.size < 3
        or np.any(~np.isfinite(grid))
        or np.any((grid <= 0) | (grid >= 1))
        or np.any(np.diff(grid) <= 0)
    ):
        raise ValueError("grid must be strictly increasing finite levels in (0, 1)")
    return grid


def _integration_weights(grid, weight):
    # Trapezoidal quadrature weights, normalized so constants integrate to one.
    delta = np.diff(grid)
    quadrature = np.empty_like(grid)
    quadrature[0], quadrature[-1] = delta[0] / 2, delta[-1] / 2
    quadrature[1:-1] = (delta[:-1] + delta[1:]) / 2
    if weight is not None:
        weight = np.asarray(weight, dtype=float)
        if weight.shape != grid.shape or np.any(~np.isfinite(weight)) or np.any(weight < 0):
            raise ValueError("weight must be finite, non-negative, and match grid")
        quadrature *= weight
    total = float(np.sum(quadrature))
    if total <= 0:
        raise ValueError("integration weights must have positive total")
    return quadrature / total


def _quantile_grid(samples, grid):
    rows = []
    for index, sample in enumerate(samples):
        sample = np.asarray(sample, dtype=float)
        if sample.ndim != 1 or sample.size < 2 or np.any(~np.isfinite(sample)):
            raise ValueError(f"samples[{index}] must contain >=2 finite observations")
        rows.append(np.quantile(sample, grid))
    return np.asarray(rows), np.asarray([len(sample) for sample in samples], dtype=float)


@dataclass
class TwoScaleResult:
    statistic: float
    p_value: float
    n_periods: int
    contamination_ratio: float
    slope_l2_norm: float
    n_boot: int
    experimental: bool = True

    def summary(self) -> str:
        return (
            f"experimental quantile-function trend diagnostic: S_W={self.statistic:.4g}, "
            f"p_boot={self.p_value:.4f} (residual circular MBB, B={self.n_boot}); "
            f"T={self.n_periods}; sampling-to-between variation ratio="
            f"{self.contamination_ratio:.2f}"
        )


def twoscale_trend(
    samples: Sequence[npt.ArrayLike],
    grid: Optional[np.ndarray] = None,
    weight=None,
    n_boot: int = 2000,
    block_length: Optional[int] = None,
    seed: int = 20260722,
) -> TwoScaleResult:
    """Compute an experimental trend diagnostic for per-period distributions.

    The squared slope norm is integrated over quantile level with trapezoidal
    weights.  A circular moving-block bootstrap resamples detrended residual
    quantile functions under a constant-mean null.  The p-value is exploratory
    and depends on the residual stationarity and chosen block length.
    """
    periods = len(samples)
    if periods < 4:
        raise ValueError("need >= 4 periods")
    if not isinstance(n_boot, (int, np.integer)) or n_boot < 1:
        raise ValueError("n_boot must be a positive integer")
    grid = _validate_grid(np.linspace(0.05, 0.95, 41) if grid is None else grid)
    integration = _integration_weights(grid, weight)
    quantiles, sizes = _quantile_grid(samples, grid)
    x = (np.arange(1, periods + 1) - (periods + 1) / 2.0) / periods
    denominator = float(np.sum(x**2))

    def statistic(q_matrix):
        slope = (x @ q_matrix) / denominator
        value = float(np.sum(slope**2 * integration))
        return value, slope

    observed, slope = statistic(quantiles)

    derivative = np.gradient(quantiles, grid, axis=1)
    density = 1.0 / np.maximum(np.abs(derivative), 1e-6)
    quantile_variance = grid * (1 - grid) / (sizes[:, None] * density**2)
    sampling_variance = (
        np.sum(x[:, None] ** 2 * quantile_variance, axis=0) / denominator**2
    )
    between_variance = np.var(quantiles, axis=0, ddof=1)
    between = float(np.sum(between_variance * integration))
    contamination = (
        float(np.sqrt(np.sum(sampling_variance * integration) / between))
        if between > 0
        else np.inf
    )

    if block_length is None:
        block_length = max(2, int(round(np.sqrt(periods))))
    if (
        not isinstance(block_length, (int, np.integer))
        or block_length < 1
        or block_length > periods
    ):
        raise ValueError("block_length must be an integer between 1 and n_periods")

    intercept = np.mean(quantiles, axis=0)
    residuals = quantiles - intercept - x[:, None] * slope
    residuals -= np.mean(residuals, axis=0)
    rng = np.random.default_rng(seed)
    blocks_needed = int(np.ceil(periods / block_length))
    null = np.empty(n_boot)
    for index in range(n_boot):
        starts = rng.integers(0, periods, size=blocks_needed)
        selected = np.concatenate(
            [(np.arange(start, start + block_length) % periods) for start in starts]
        )[:periods]
        null[index] = statistic(intercept + residuals[selected])[0]
    p_value = float((1 + np.sum(null >= observed)) / (n_boot + 1))
    return TwoScaleResult(
        statistic=observed,
        p_value=p_value,
        n_periods=periods,
        contamination_ratio=contamination,
        slope_l2_norm=float(np.sqrt(observed)),
        n_boot=n_boot,
    )


def wasserstein_decomposition(
    sample1: npt.ArrayLike,
    sample2: npt.ArrayLike,
    grid: Optional[np.ndarray] = None,
) -> dict:
    """Numerical location/scale/shape decomposition on a quantile grid.

    The decomposition is exact for the discrete quadrature definition returned
    here.  It approximates continuous ``W2^2`` according to the supplied grid.
    """
    grid = _validate_grid(np.linspace(0.005, 0.995, 400) if grid is None else grid)
    integration = _integration_weights(grid, None)
    sample1 = np.asarray(sample1, dtype=float)
    sample2 = np.asarray(sample2, dtype=float)
    if (
        sample1.ndim != 1
        or sample2.ndim != 1
        or sample1.size < 2
        or sample2.size < 2
        or np.any(~np.isfinite(sample1))
        or np.any(~np.isfinite(sample2))
    ):
        raise ValueError("each sample must be a finite 1-D array with >=2 values")
    q1, q2 = np.quantile(sample1, grid), np.quantile(sample2, grid)
    mean1, mean2 = float(np.sum(q1 * integration)), float(np.sum(q2 * integration))
    centered1, centered2 = q1 - mean1, q2 - mean2
    sd1 = float(np.sqrt(np.sum(centered1**2 * integration)))
    sd2 = float(np.sqrt(np.sum(centered2**2 * integration)))
    if sd1 > 0 and sd2 > 0:
        correlation = float(np.sum(centered1 * centered2 * integration) / (sd1 * sd2))
        correlation = float(np.clip(correlation, -1.0, 1.0))
    else:
        correlation = 1.0
    location = float((mean2 - mean1) ** 2)
    scale = float((sd1 - sd2) ** 2)
    shape = float(max(0.0, 2 * sd1 * sd2 * (1 - correlation)))
    return {
        "w2_squared": location + scale + shape,
        "location": location,
        "scale": scale,
        "shape": shape,
        "correlation": correlation,
        "mean_shift": mean2 - mean1,
        "sd1": sd1,
        "sd2": sd2,
    }
