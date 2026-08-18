"""Grouped GPD regression: a covariate-dependent scale, and return levels.

This module generalises the interval-censored (grouped) fit of :mod:`nsevt.grouped`
from a single scale to a **log-linear scale design**, and adds the two
profile-likelihood intervals that a tail regression needs:

* :func:`fit_grouped_design` -- the grouped GPD maximum-likelihood fit with an
  arbitrary log-scale design matrix ``X`` (intercept in column 0), so the same
  routine covers a stationary tail, a trend, group-specific scales, or any
  combination;
* :func:`profile_ci_coef` -- a profile-likelihood interval for any coefficient of
  that design (the shape and the remaining coefficients are maximised out at each
  trial value), the natural interval counterpart to the permutation trend test in
  :mod:`nsevt.trend`;
* :func:`return_level` and :func:`profile_ci_return_level` -- the level exceeded
  once per ``m`` observations at exceedance rate ``rate``, and a profile-likelihood
  interval for that level. Profiling the level itself (rather than pushing a
  profiled shape through the formula, or resampling) inherits the likelihood's
  reparameterisation invariance and respects the non-linearity of the level as
  the shape approaches zero.

The design matrix, the covariate coding and the exceedance rate are supplied by
the caller: this module knows nothing about what the covariate means. The
numerics are ported from unit-tested research code and use NumPy and SciPy only.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2

from .grouped import _gpd_surv, _grouped_nll_nu, interval_cells


# --------------------------------------------------------------------------
# grouped likelihood with a log-scale design
# --------------------------------------------------------------------------
def _grouped_nll_design(par, a, b, X, trunc):
    """Grouped negative log-likelihood with a general log-scale design.

    ``X`` is ``(n, p)`` with an intercept in column 0, so the same routine covers
    the stationary model, a trend, group-specific scales, or any combination.
    """
    xi = par[0]
    if abs(xi) < 1e-10:
        xi = 1e-10 if xi >= 0 else -1e-10
    s = np.exp(X @ par[1:])
    if not np.all(np.isfinite(s)) or np.any(s <= 0):
        return 1e10
    if np.any(1.0 + xi * a / s <= 0):
        return 1e10
    cell = _gpd_surv(a, xi, s) - _gpd_surv(b, xi, s)
    if np.any(~np.isfinite(cell)) or np.any(cell <= 1e-300):
        return 1e10
    ll = float(np.sum(np.log(cell)))
    trunc = np.asarray(trunc, dtype=float)
    if np.any(trunc > 0):
        surv = np.broadcast_to(_gpd_surv(trunc, xi, s), np.shape(a))
        if np.any(surv <= 1e-300):
            return 1e10
        ll -= float(np.sum(np.log(surv)))
    return -ll


def _cells_of(values, threshold, grid, cells):
    if cells is not None:
        if not isinstance(cells, (tuple, list)) or len(cells) != 3:
            raise ValueError("cells must be an (a, b, trunc) triple")
        out = tuple(np.asarray(c, dtype=float) for c in cells)
        if any(c.ndim != 1 or c.size != values.size for c in out):
            raise ValueError("each cells array must be 1-D and aligned with values")
        if any(np.any(~np.isfinite(c)) for c in out):
            raise ValueError("cells must contain only finite values")
        a, b, trunc = out
        if np.any(a < 0) or np.any(b <= a) or np.any(trunc < 0) or np.any(trunc > b):
            raise ValueError("cells must satisfy 0 <= a < b and 0 <= trunc <= b")
        return out
    return interval_cells(values, threshold, grid)


def _validated_inputs(
    values: npt.ArrayLike,
    threshold: float,
    design: npt.ArrayLike,
    grid: float | npt.ArrayLike,
    cells: Optional[tuple],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    X = np.asarray(design, dtype=float)
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if values.ndim != 1 or np.any(~np.isfinite(values)):
        raise ValueError("values must be a finite 1-D array")
    if values.size < 3 or np.any(values <= threshold):
        raise ValueError("values must contain at least 3 marks strictly above threshold")
    if X.ndim != 2 or X.shape[0] != values.size or X.shape[1] < 1:
        raise ValueError("design must be a 2-D matrix aligned with values")
    if np.any(~np.isfinite(X)):
        raise ValueError("design must contain only finite values")
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError("design must have full column rank")
    a, b, trunc = _cells_of(values, threshold, grid, cells)
    return values, X, a, b, trunc


def fit_grouped_design(
    values: npt.ArrayLike,
    threshold: float,
    design: npt.ArrayLike,
    grid: float | npt.ArrayLike = 5.0,
    cells: Optional[tuple] = None,
    starts: tuple = (-0.4, -0.25, -0.1, 0.05),
) -> dict:
    """Grouped GPD MLE with a log-linear scale design.

    ``values`` are the exceedance marks (all above ``threshold``) and ``design``
    is the aligned ``(n, p)`` matrix with an intercept in column 0; the log-scale
    of observation ``i`` is ``design[i] @ coef``. Pass ``cells`` (the
    ``(a, b, trunc)`` triple from :func:`~nsevt.grouped.interval_cells`) to reuse
    mixed-precision cells; otherwise they are built from ``grid``.

    Returns ``{"xi", "coef", "sigma0", "loglik", "n", "p", "truncation"}``, where
    ``coef`` are the log-scale coefficients and ``sigma0 = exp(coef[0])`` is the
    baseline scale.
    """
    values, X, a, b, trunc = _validated_inputs(
        values, threshold, design, grid, cells
    )
    starts_arr = np.asarray(starts, dtype=float)
    if starts_arr.ndim != 1 or starts_arr.size == 0 or np.any(~np.isfinite(starts_arr)):
        raise ValueError("starts must be a non-empty finite 1-D sequence")
    z = values - float(threshold)
    best = (np.inf, None)
    for x0 in starts_arr:
        init = np.concatenate([[x0, np.log(max(float(np.mean(z)), 1e-3))],
                               np.zeros(X.shape[1] - 1)])
        r = minimize(_grouped_nll_design, init, args=(a, b, X, trunc),
                     method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-10,
                              "maxiter": 200000, "maxfev": 200000})
        if r.success and np.isfinite(r.fun) and r.fun < best[0]:
            best = (float(r.fun), r.x)
    par = best[1]
    if par is None:
        raise RuntimeError("grouped design fit did not converge from any start")
    tr = np.unique(np.asarray(trunc, dtype=float)).tolist()
    return {"xi": float(par[0]), "coef": par[1:].tolist(),
            "sigma0": float(np.exp(par[1])), "loglik": -best[0],
            "n": int(z.size), "p": int(X.shape[1]),
            "truncation": tr[0] if len(tr) == 1 else tr}


def profile_ci_coef(
    values: npt.ArrayLike,
    threshold: float,
    design: npt.ArrayLike,
    coef: int = 1,
    grid: float | npt.ArrayLike = 5.0,
    level: float = 0.95,
    cells: Optional[tuple] = None,
    fit: Optional[dict] = None,
    n_bisect: int = 32,
    span: float = 0.6,
) -> dict:
    """Profile-likelihood interval for one coefficient of the scale design.

    ``coef`` selects the column of ``design`` to profile; the shape and the
    remaining coefficients are maximised out at each trial value. The interval
    refers to the same grouped estimator whose point value is reported, and (as
    a profile interval) it is invariant to reparameterisation. Use it, for
    example, as the interval counterpart of a permutation trend test on a time
    column.

    Returns ``{"coef", "estimate", "ci", "at_bound", "level", "loglik"}``.
    """
    values, X, a, b, trunc = _validated_inputs(
        values, threshold, design, grid, cells
    )
    if not isinstance(coef, (int, np.integer)) or not 0 <= coef < X.shape[1]:
        raise ValueError("coef must index a design column")
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between 0 and 1")
    if not isinstance(n_bisect, (int, np.integer)) or n_bisect < 1:
        raise ValueError("n_bisect must be a positive integer")
    if not np.isfinite(span) or span <= 0:
        raise ValueError("span must be positive and finite")
    if fit is None:
        fit = fit_grouped_design(values, threshold, X, cells=(a, b, trunc))
    fit_coef = np.asarray(fit.get("coef", []), dtype=float)
    if (fit_coef.shape != (X.shape[1],) or np.any(~np.isfinite(fit_coef))
            or not np.isfinite(fit.get("xi", np.nan))
            or not np.isfinite(fit.get("loglik", np.nan))):
        raise ValueError("fit is not a finite grouped-design fit for this design")
    b_hat = float(fit_coef[coef])
    ll_max = float(fit["loglik"])
    target = ll_max - float(chi2.ppf(level, 1)) / 2.0

    keep = [j for j in range(X.shape[1]) if j != coef]
    x_col, X_rest = X[:, coef], X[:, keep]
    start = np.concatenate([[fit["xi"]], np.asarray(fit["coef"], float)[keep]])

    def prof(bv):
        # the fixed coefficient enters as a known offset of the log-scale
        off = bv * x_col
        Xa = np.column_stack([X_rest, off]) if X_rest.size else off[:, None]

        def nll(par):
            return _grouped_nll_design(np.concatenate([par, [1.0]]), a, b, Xa, trunc)

        r = minimize(nll, start, method="Nelder-Mead",
                     options={"xatol": 1e-9, "fatol": 1e-9,
                              "maxiter": 40000, "maxfev": 40000})
        return -float(r.fun) if r.success and np.isfinite(r.fun) else -np.inf

    scale = max(abs(b_hat), 0.05)

    def hunt(direction):
        step = span * scale
        inside, outside = b_hat, b_hat + direction * step
        for _ in range(8):
            if prof(outside) < target:
                break
            inside, outside = outside, outside + direction * step
        else:
            return float(outside), True          # never dropped: open bound
        for _ in range(n_bisect):
            mid = 0.5 * (inside + outside)
            if prof(mid) >= target:
                inside = mid
            else:
                outside = mid
        return float(0.5 * (inside + outside)), False

    lo, lo_open = hunt(-1.0)
    hi, hi_open = hunt(+1.0)
    return {"coef": int(coef), "estimate": b_hat, "ci": [lo, hi],
            "at_bound": bool(lo_open or hi_open), "level": level,
            "loglik": ll_max}


# --------------------------------------------------------------------------
# return levels
# --------------------------------------------------------------------------
def return_level(xi: float, sigma: float, threshold: float, rate: float,
                 m: npt.ArrayLike) -> np.ndarray:
    """Level exceeded once per ``m`` observations at exceedance ``rate``.

    ``rate`` is the probability that an observation exceeds ``threshold`` (the
    peaks-over-threshold rate ``zeta``); ``m`` is the return period in the same
    observation unit. Returns ``threshold + (sigma / xi)[(m rate)^xi - 1]``, or
    the ``xi -> 0`` limit ``threshold + sigma log(m rate)``.
    """
    m = np.asarray(m, dtype=float)
    xi = float(xi)
    sigma = float(sigma)
    threshold = float(threshold)
    rate = float(rate)
    if not all(np.isfinite(v) for v in (xi, sigma, threshold, rate)):
        raise ValueError("xi, sigma, threshold and rate must be finite")
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if not 0 < rate <= 1:
        raise ValueError("rate must lie in (0, 1]")
    if m.size == 0 or np.any(~np.isfinite(m)) or np.any(m <= 0):
        raise ValueError("m must contain positive finite return periods")
    with np.errstate(over="ignore", invalid="ignore"):
        if abs(xi) < 1e-8:
            return threshold + sigma * np.log(m * rate)
        return threshold + (sigma / xi) * ((m * rate) ** xi - 1.0)


def profile_ci_return_level(
    values: npt.ArrayLike,
    threshold: float,
    rate: float,
    m: float,
    grid: float | npt.ArrayLike = 5.0,
    level: float = 0.95,
    cells: Optional[tuple] = None,
    n_bisect: int = 40,
    span: tuple = (0.2, 6.0),
) -> Optional[dict]:
    """Profile-likelihood interval for the ``m``-observation return level.

    The return level is a non-linear function of ``(xi, sigma)``, so an interval
    built on the shape and pushed through the formula is not the profile interval
    for the level. Here the level is made a parameter and the shape maximised out
    at each trial value, through the ``(xi, nu)`` reparameterisation the endpoint
    interval uses. Restricted to bounded tails (``xi < 0``); returns ``None`` when
    the point return level is not finite and positive.

    Returns ``{"m", "return_level", "ci", "upper_at_bound", "xi_at_max",
    "loglik", "level"}``.
    """
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if not np.isfinite(rate) or not 0 < rate <= 1:
        raise ValueError("rate must lie in (0, 1]")
    if not np.isfinite(m) or m <= 0:
        raise ValueError("m must be positive and finite")
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between 0 and 1")
    if not isinstance(n_bisect, (int, np.integer)) or n_bisect < 1:
        raise ValueError("n_bisect must be a positive integer")
    span_arr = np.asarray(span, dtype=float)
    if (span_arr.shape != (2,) or np.any(~np.isfinite(span_arr))
            or not 0 < span_arr[0] < 1 < span_arr[1]):
        raise ValueError("span must be a finite (lower, upper) pair around 1")
    values = np.asarray(values, dtype=float)
    X = np.ones((values.size, 1), dtype=float)
    values, _, a, b, trunc = _validated_inputs(
        values, threshold, X, grid, cells
    )
    a_max = float(a.max())
    mz = float(m) * float(rate)

    def nll(xi, y):
        c = mz ** xi - 1.0
        if not np.isfinite(c) or abs(c) < 1e-12:
            return 1e10
        nu = -y / c
        return _grouped_nll_nu(xi, nu, a, b, trunc)

    def prof(y):
        r = minimize_scalar(lambda x: nll(x, y), bounds=(-0.95, -1e-4),
                            method="bounded", options={"xatol": 1e-10})
        if not r.success or not np.isfinite(r.fun):
            return -np.inf, float("nan")
        return -float(r.fun), float(r.x)

    f = fit_grouped_design(values, threshold, X,
                           cells=(a, b, trunc))
    if f["xi"] >= 0:
        return None
    y_hat = float(return_level(f["xi"], f["sigma0"], threshold, rate, m) - threshold)
    if not np.isfinite(y_hat) or y_hat <= 0:
        return None

    mesh = y_hat * np.linspace(span[0], span[1], 120)
    vals = np.array([prof(v)[0] for v in mesh])
    if np.all(~np.isfinite(vals)):
        raise RuntimeError("return-level profile optimization failed")
    k = int(np.argmax(vals))
    r = minimize_scalar(lambda v: -prof(v)[0],
                        bounds=(mesh[max(k - 1, 0)], mesh[min(k + 1, mesh.size - 1)]),
                        method="bounded", options={"xatol": 1e-8})
    if not r.success or not np.isfinite(r.fun):
        raise RuntimeError("return-level profile optimization failed")
    y_max = float(r.x)
    ll_max = -float(r.fun)
    target = ll_max - float(chi2.ppf(level, 1)) / 2.0

    def bisect(inside, outside):
        for _ in range(n_bisect):
            mid = 0.5 * (inside + outside)
            if prof(mid)[0] >= target:
                inside = mid
            else:
                outside = mid
        return 0.5 * (inside + outside)

    lo_out = max(a_max * 1e-6, y_max * 1e-3)
    y_lo = bisect(y_max, lo_out) if prof(lo_out)[0] < target else lo_out
    hi_out = y_max * span[1]
    at_bound = prof(hi_out)[0] >= target
    y_hi = hi_out if at_bound else bisect(y_max, hi_out)

    return {"m": float(m), "return_level": float(threshold + y_max),
            "ci": [float(threshold + y_lo), float(threshold + y_hi)],
            "upper_at_bound": bool(at_bound), "xi_at_max": prof(y_max)[1],
            "loglik": ll_max, "level": level}
