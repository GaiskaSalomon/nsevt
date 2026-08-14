"""Conformal upper bounds for threshold exceedances.

``split_conformal`` provides the standard finite-sample marginal guarantee when
the calibration scores and future tail score are exchangeable.  The target is
a future observation conditional on exceeding the same fixed threshold.

``block_conformal`` is an experimental block-aggregation diagnostic for ordered
dependent exceedances.  It may be useful in sensitivity analyses, but this
implementation does not claim a general finite-sample or beta-mixing coverage
guarantee.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import numpy.typing as npt

from .gpd import fit_gpd


@dataclass
class ConformalBand:
    """One-sided upper conformal bound for a future tail observation."""

    threshold: float
    alpha: float
    q_standardized: float
    method: str
    n_blocks: int
    block_length: int
    experimental: bool = False

    def predict_upper(self, scale: float) -> float:
        """Return ``u + q * sigma`` for a positive finite scale."""
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("scale must be positive and finite")
        return float(self.threshold + self.q_standardized * scale)

    def coverage(self, x_new: npt.ArrayLike, scale) -> float:
        """Empirical conditional-tail coverage on raw values ``x_new > u``."""
        x_new = np.asarray(x_new, dtype=float)
        if x_new.ndim != 1 or x_new.size == 0 or np.any(~np.isfinite(x_new)):
            raise ValueError("x_new must be a non-empty finite 1-D array")
        if np.any(x_new <= self.threshold):
            raise ValueError("x_new must contain raw exceedances above threshold, not excesses")
        scale = np.broadcast_to(np.asarray(scale, dtype=float), x_new.shape)
        if np.any(~np.isfinite(scale)) or np.any(scale <= 0):
            raise ValueError("scale must be positive, finite, and broadcastable to x_new")
        upper = self.threshold + self.q_standardized * scale
        return float(np.mean(x_new <= upper))


def _validate_inputs(x, threshold, alpha):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or np.any(~np.isfinite(x)):
        raise ValueError("x must be a finite 1-D array")
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")
    mask = x > threshold
    z = x[mask] - threshold
    if z.size < 4:
        raise ValueError(f"only {z.size} exceedances; need >= 4 for calibration")
    return x, mask, z


def _standardized_scores(z, mask, x_shape, scale):
    if scale is None:
        sigma = fit_gpd(z)["sigma"]
        return z / sigma
    scale = np.asarray(scale, dtype=float)
    if scale.ndim == 0:
        selected = np.full(z.shape, float(scale))
    elif scale.shape == x_shape:
        selected = scale[mask]
    elif scale.shape == z.shape:
        selected = scale
    else:
        raise ValueError("scale must be scalar or match x or its exceedances")
    if np.any(~np.isfinite(selected)) or np.any(selected <= 0):
        raise ValueError("scale values must be positive and finite")
    return z / selected


def _blocks(n: int, block_length: Optional[int], n_blocks: Optional[int]):
    if block_length is not None and n_blocks is not None:
        raise ValueError("specify block_length or n_blocks, not both")
    if block_length is not None and (
        not isinstance(block_length, (int, np.integer)) or block_length < 1
    ):
        raise ValueError("block_length must be a positive integer")
    if n_blocks is not None and (
        not isinstance(n_blocks, (int, np.integer)) or n_blocks < 2
    ):
        raise ValueError("n_blocks must be an integer >= 2")
    if block_length is None and n_blocks is None:
        block_length = max(2, int(round(n**0.5)))
    if block_length is None:
        block_length = max(1, int(np.ceil(n / n_blocks)))
    return [(a, min(a + block_length, n)) for a in range(0, n, block_length)]


def block_conformal(
    x: npt.ArrayLike,
    threshold: float,
    alpha: float = 0.1,
    scale=None,
    order: Optional[Sequence] = None,
    block_length: Optional[int] = None,
    n_blocks: Optional[int] = None,
) -> ConformalBand:
    """Experimental block-aggregated upper bound for ordered exceedances.

    Contiguous calibration blocks are summarized by within-block upper
    quantiles, followed by a conformal order statistic across blocks.  This is
    exposed for dependence sensitivity analysis and is not part of nsevt's
    stable inferential core.
    """
    x, mask, z = _validate_inputs(x, threshold, alpha)
    scores = _standardized_scores(z, mask, x.shape, scale)
    if order is not None:
        order = np.asarray(order)
        if order.shape == x.shape:
            order = order[mask]
        if order.shape != z.shape:
            raise ValueError("order must match x or its exceedances")
        scores = scores[np.argsort(order, kind="stable")]

    blocks = _blocks(len(scores), block_length, n_blocks)
    aggregates = np.array(
        [np.quantile(scores[a:b], 1 - alpha, method="higher") for a, b in blocks]
    )
    if aggregates.size < 2:
        raise ValueError("need >= 2 calibration blocks; reduce block_length")
    level = min(1.0, np.ceil((1 - alpha) * (aggregates.size + 1)) / aggregates.size)
    q = float(np.quantile(aggregates, level, method="higher"))
    return ConformalBand(
        threshold=float(threshold),
        alpha=float(alpha),
        q_standardized=q,
        method="experimental-block-aggregate",
        n_blocks=int(aggregates.size),
        block_length=int(blocks[0][1] - blocks[0][0]),
        experimental=True,
    )


def split_conformal(
    x: npt.ArrayLike, threshold: float, alpha: float = 0.1, scale=None
) -> ConformalBand:
    """Marginal split-conformal upper bound for exchangeable tail scores.

    Scale values must be fitted independently of these calibration scores (or
    otherwise fixed) for the usual split-conformal guarantee to apply.  If
    ``scale=None``, nsevt estimates a common GPD scale on the same observations;
    that convenience mode is model-assisted and does not retain the exact
    distribution-free split guarantee.
    """
    x, mask, z = _validate_inputs(x, threshold, alpha)
    scores = _standardized_scores(z, mask, x.shape, scale)
    n = len(scores)
    level = min(1.0, np.ceil((1 - alpha) * (n + 1)) / n)
    q = float(np.quantile(scores, level, method="higher"))
    return ConformalBand(
        threshold=float(threshold),
        alpha=float(alpha),
        q_standardized=q,
        method="split-conformal",
        n_blocks=1,
        block_length=n,
        experimental=False,
    )
