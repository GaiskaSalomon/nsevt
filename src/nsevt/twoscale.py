"""Two-scale trend test for a time series of distributions (Wasserstein/Frechet).

When each period's distribution is not observed but *estimated* from a finite
within-period sample, a trend in the Frechet (Wasserstein-2 barycenter) mean is
tested by the slope of the empirical quantile functions in L2, with a
**moving-block bootstrap** null.  The two-scale decomposition separates a
between-period signal from within-period sampling noise: **validity** is
governed by the smallest within-period sample size, **power** by the number of
periods, and the two decouple.  A contamination ratio ``r_T`` locates a record
relative to the validity frontier.

Implements the two-scale-trend manuscript.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np


def _quantile_grid(samples, grid):
    """Empirical quantile functions of each period's sample on ``grid``."""
    return np.array([np.quantile(np.asarray(s, float), grid) for s in samples])


@dataclass
class TwoScaleResult:
    statistic: float
    p_value: float
    n_periods: int
    contamination_ratio: float
    slope_l2_norm: float
    n_boot: int

    def summary(self) -> str:
        return (f"two-scale Frechet trend test: S_W={self.statistic:.4g}, "
                f"p={self.p_value:.4f} (block bootstrap, B={self.n_boot}); "
                f"T={self.n_periods} periods; contamination r_T={self.contamination_ratio:.2f} "
                f"({'noise-dominated' if self.contamination_ratio > 1 else 'signal-resolvable'})")


def twoscale_trend(samples: Sequence[Sequence[float]],
                   grid: Optional[np.ndarray] = None,
                   weight=None, n_boot: int = 2000, block_length: Optional[int] = None,
                   seed: int = 20260722) -> TwoScaleResult:
    """Test for a trend in the Frechet mean of ``T`` per-period samples.

    Parameters
    ----------
    samples : sequence of arrays
        ``samples[t]`` is the within-period sample of period ``t`` (chronological
        order); sizes may differ and be small.
    grid : array, optional
        Quantile levels ``u`` on a compact ``[u_L, u_U] subset (0,1)`` where the
        slope is evaluated (default 41 levels on ``[0.05, 0.95]``).
    weight : array, optional
        Weight ``w(u)`` over the grid (default uniform).
    block_length : int, optional
        Moving-block length ``ell`` (default ``~ sqrt(T)``).
    """
    T = len(samples)
    if T < 4:
        raise ValueError("need >= 4 periods")
    if grid is None:
        grid = np.linspace(0.05, 0.95, 41)
    grid = np.asarray(grid, float)
    w = np.ones_like(grid) if weight is None else np.asarray(weight, float)
    Q = _quantile_grid(samples, grid)                     # (T, m)
    ns = np.array([len(s) for s in samples], float)

    x = (np.arange(1, T + 1) - (T + 1) / 2.0) / T          # centered, rescaled
    D_T = float(np.sum(x ** 2))

    def s_w(Qmat):
        beta = (x @ Qmat) / D_T                             # slope per level
        return float(np.mean(beta ** 2 * w)), beta

    obs, beta = s_w(Q)

    # --- contamination ratio r_T (validity-frontier diagnostic) ---
    # within-period quantile variance ~ u(1-u)/(n f^2), f from finite diffs
    du = np.gradient(grid)
    dens = np.where(np.abs(np.gradient(Q, grid, axis=1)) > 1e-9,
                    1.0 / np.abs(np.gradient(Q, grid, axis=1)), 0.0)  # f_t(u)
    vt = grid * (1 - grid) / (ns[:, None] * np.maximum(dens, 1e-6) ** 2)
    vbar = (D_T ** -1) * np.sum(x[:, None] ** 2 * vt, axis=0)         # ~ within noise
    sinf2 = np.var(Q, axis=0, ddof=1)                                  # between-period
    denom = np.sum(sinf2 * w)
    r_T = float(np.sqrt(np.sum(vbar * w) / denom)) if denom > 0 else np.inf

    # --- moving-block bootstrap null (reassign the time ordering) ---
    if block_length is None:
        block_length = max(2, int(round(np.sqrt(T))))
    rng = np.random.default_rng(seed)
    n_blocks = int(np.ceil(T / block_length))
    null = np.empty(n_boot)
    for b in range(n_boot):
        starts = rng.integers(0, T, size=n_blocks)
        idx = np.concatenate([(np.arange(s, s + block_length) % T) for s in starts])[:T]
        null[b] = s_w(Q[idx])[0]
    p = (1 + np.sum(null >= obs)) / (1 + n_boot)
    return TwoScaleResult(statistic=obs, p_value=float(p), n_periods=T,
                          contamination_ratio=r_T,
                          slope_l2_norm=float(np.sqrt(np.sum(beta ** 2 * w))),
                          n_boot=n_boot)


def wasserstein_decomposition(sample1: Sequence[float], sample2: Sequence[float],
                              grid: Optional[np.ndarray] = None) -> dict:
    """Exact energy split of ``W2^2`` between two samples into location/scale/shape.

    ``W2^2 = (mu2-mu1)^2 + (sigma1-sigma2)^2 + 2 sigma1 sigma2 (1-rho)``, with the
    location term L2-orthogonal to the rest.  The scale/shape split is an exact,
    non-negative additive energy budget (equal to a projection only when
    ``rho=1`` or ``sigma1=sigma2``).
    """
    if grid is None:
        grid = np.linspace(0.005, 0.995, 400)
    q1 = np.quantile(np.asarray(sample1, float), grid)
    q2 = np.quantile(np.asarray(sample2, float), grid)
    mu1, mu2 = q1.mean(), q2.mean()
    c1, c2 = q1 - mu1, q2 - mu2
    s1, s2 = np.sqrt(np.mean(c1 ** 2)), np.sqrt(np.mean(c2 ** 2))
    rho = float(np.mean(c1 * c2) / (s1 * s2)) if s1 > 0 and s2 > 0 else 1.0
    loc = float((mu2 - mu1) ** 2)
    scale = float((s1 - s2) ** 2)
    shape = float(2 * s1 * s2 * (1 - rho))
    return dict(w2_squared=loc + scale + shape, location=loc, scale=scale,
                shape=shape, correlation=rho,
                mean_shift=float(mu2 - mu1), sd1=float(s1), sd2=float(s2))
