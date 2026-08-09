"""Block-conformal prediction bands for non-stationary extreme quantiles.

Classical split-conformal prediction gives finite-sample coverage under
*exchangeability*, which temporally dependent extremes violate.  Block conformal
calibrates over contiguous temporal **blocks** treated as the approximately
exchangeable units, yielding

    P(X in C_alpha) >= 1 - alpha - o(1)

under beta-mixing, with the o(1) term controlled by the block length.  The
non-conformity score is the **standardized excess** ``s = (x - u)/sigma``: by the
pivotal property of the GPD, dividing by the (possibly covariate-dependent)
scale removes the covariate/non-stationarity, leaving a common score law across
covariate values.

Implements the method of the block-conformal-EVT manuscript.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .gpd import fit_gpd


@dataclass
class ConformalBand:
    """A one-sided upper conformal prediction band for a tail."""
    threshold: float
    alpha: float
    q_standardized: float
    method: str
    n_blocks: int
    block_length: int

    def predict_upper(self, scale: float) -> float:
        """Upper prediction bound ``u + q * sigma`` at scale ``sigma``."""
        return float(self.threshold + self.q_standardized * scale)

    def coverage(self, x_new: Sequence[float], scale) -> float:
        """Empirical coverage of the band on new exceedances ``x_new``.

        ``scale`` is a scalar or per-observation array of ``sigma`` values used
        to build the band for each new point.
        """
        x_new = np.asarray(x_new, float)
        scale = np.broadcast_to(np.asarray(scale, float), x_new.shape)
        upper = self.threshold + self.q_standardized * scale
        return float(np.mean(x_new <= upper))


def _blocks(n: int, block_length: Optional[int], n_blocks: Optional[int]):
    if block_length is None and n_blocks is None:
        block_length = max(2, int(round(n ** 0.5)))  # ell ~ sqrt(T)
    if block_length is None:
        block_length = max(2, n // int(n_blocks))
    edges = list(range(0, n, block_length))
    return [(edges[i], min(edges[i] + block_length, n)) for i in range(len(edges))]


def block_conformal(x: Sequence[float], threshold: float, alpha: float = 0.1,
                    scale=None, order: Optional[Sequence] = None,
                    block_length: Optional[int] = None,
                    n_blocks: Optional[int] = None) -> ConformalBand:
    """One-sided block-conformal upper band for peaks over ``threshold``.

    Parameters
    ----------
    x : array
        Raw observations; only ``x > threshold`` are used as calibration
        exceedances.
    threshold : float
        Tail threshold ``u``.
    alpha : float
        Miscoverage level; the band targets coverage ``>= 1 - alpha``.
    scale : float or array, optional
        Per-observation GPD scale ``sigma`` used to standardize excesses.  If
        ``None``, a constant scale is estimated by GPD MLE on the exceedances.
        For non-stationary/covariate-dependent tails, pass ``sigma(z_i)``.
    order : array, optional
        Temporal ordering key for the exceedances (e.g. time); blocks are formed
        on the sorted order.  If ``None``, input order is assumed temporal.
    block_length, n_blocks : int, optional
        Block size / count.  Default block length ``~ sqrt(#exceedances)``.
    """
    x = np.asarray(x, float)
    mask = x > threshold
    z = x[mask] - threshold
    if z.size < 4:
        raise ValueError(f"only {z.size} exceedances; need >= 4 for calibration")
    if scale is None:
        sigma = fit_gpd(z)["sigma"]
        s = z / sigma
    else:
        scale = np.asarray(scale, float)
        sc = scale[mask] if scale.shape == x.shape else scale
        s = z / sc
    if order is not None:
        order = np.asarray(order)[mask] if np.asarray(order).shape == x.shape else np.asarray(order)
        s = s[np.argsort(order, kind="stable")]

    blocks = _blocks(len(s), block_length, n_blocks)
    # per-block aggregate score: within-block (1-alpha) quantile
    agg = np.array([np.quantile(s[a:b], 1 - alpha) for a, b in blocks if b > a])
    K = len(agg)
    if K < 2:
        raise ValueError("need >= 2 calibration blocks; reduce block_length")
    # finite-sample conformal quantile of the block aggregates
    lvl = min(1.0, np.ceil((1 - alpha) * (K + 1)) / K)
    q = float(np.quantile(agg, lvl))
    bl = blocks[0][1] - blocks[0][0]
    return ConformalBand(threshold=float(threshold), alpha=float(alpha),
                         q_standardized=q, method="block-conformal",
                         n_blocks=K, block_length=int(bl))


def split_conformal(x: Sequence[float], threshold: float, alpha: float = 0.1,
                    scale=None) -> ConformalBand:
    """Marginal split-conformal upper band (ignores temporal dependence).

    Provided as the baseline whose coverage the block version corrects under
    dependence; use :func:`block_conformal` for dependent extremes.
    """
    x = np.asarray(x, float)
    mask = x > threshold
    z = x[mask] - threshold
    if z.size < 4:
        raise ValueError(f"only {z.size} exceedances; need >= 4")
    sigma = fit_gpd(z)["sigma"] if scale is None else None
    s = z / (sigma if scale is None else np.asarray(scale, float)[mask])
    n = len(s)
    lvl = min(1.0, np.ceil((1 - alpha) * (n + 1)) / n)
    q = float(np.quantile(s, lvl))
    return ConformalBand(threshold=float(threshold), alpha=float(alpha),
                         q_standardized=q, method="split-conformal",
                         n_blocks=1, block_length=n)
