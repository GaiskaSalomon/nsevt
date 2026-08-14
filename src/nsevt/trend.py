"""Permutation-calibrated trend inference for a GPD tail scale.

The model is ``log sigma(t) = beta0 + beta1 t``, with ``t`` measured in
decades.  The primary statistic is the likelihood-ratio (LR) statistic for a
constant versus linearly varying scale.  Calibration permutes complete block
labels: observations within a block remain together while the block's time
label is reassigned.  The resulting test is exact only under exchangeability
of those block labels under the null; otherwise it is a design-based
sensitivity analysis, not a universal dependence correction.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize
from scipy.stats import chi2

_XI0 = (-0.4, -0.25, -0.1, 0.05, 0.2)


def _validate_design(z, block) -> tuple[np.ndarray, np.ndarray]:
    z = np.asarray(z, dtype=float)
    block = np.asarray(block)
    if z.ndim != 1 or block.ndim != 1 or z.size != block.size:
        raise ValueError("z and block must be matching 1-D arrays")
    if z.size < 3:
        raise ValueError("need at least 3 excesses")
    if np.any(~np.isfinite(z)) or np.any(z < 0) or np.ptp(z) == 0:
        raise ValueError("z must contain finite, non-negative, varying excesses")
    try:
        block_float = block.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("block labels must be numeric time values") from exc
    if np.any(~np.isfinite(block_float)) or np.unique(block_float).size < 2:
        raise ValueError("block must contain at least 2 distinct finite time labels")
    return z, block_float


def _nll_ns(par, z, t, varying):
    xi = float(par[0])
    log_s = par[1] + par[2] * t if varying else np.full_like(t, par[1])
    if not np.isfinite(xi) or xi <= -0.999 or np.any(~np.isfinite(log_s)):
        return 1e100
    sigma = np.exp(log_s)
    if np.any(~np.isfinite(sigma)):
        return 1e100
    y = 1.0 + xi * z / sigma
    if np.any(y <= 0):
        return 1e100
    if abs(xi) < 1e-8:
        return float(np.sum(log_s + z / sigma))
    return float(np.sum(log_s) + (1.0 + 1.0 / xi) * np.sum(np.log(y)))


def _fit_ns(z, t, varying):
    best, bp = np.inf, None
    mean_z = max(float(np.mean(z)), 1e-8)
    for xi0 in _XI0:
        p0 = [xi0, np.log(mean_z), 0.0] if varying else [xi0, np.log(mean_z)]
        r = minimize(
            _nll_ns,
            p0,
            args=(z, t, varying),
            method="Nelder-Mead",
            options={"maxiter": 30000, "maxfev": 30000, "xatol": 1e-8, "fatol": 1e-8},
        )
        if r.success and np.isfinite(r.fun) and r.fun < best:
            best, bp = float(r.fun), np.asarray(r.x, dtype=float)
    if bp is None:
        raise RuntimeError("non-stationary GPD optimization failed")
    return bp, best


def _fit_ns_from(z, t, varying, start):
    """Single-start non-stationary GPD fit from a given ``start`` guess.

    Warm-starts a permutation from the null maximum-likelihood estimate (the
    constant-scale fit with zero trend), which sits near every permuted-label
    optimum; ``None`` signals a fall-back to the multi-start :func:`_fit_ns`.
    """
    r = minimize(
        _nll_ns, start, args=(z, t, varying), method="Nelder-Mead",
        options={"maxiter": 30000, "maxfev": 30000, "xatol": 1e-8, "fatol": 1e-8},
    )
    if r.success and np.isfinite(r.fun):
        return np.asarray(r.x, dtype=float), float(r.fun)
    return None


def _decades(block, ref: Optional[float]) -> np.ndarray:
    block = np.asarray(block, dtype=float)
    ref = float(np.min(block)) if ref is None else float(ref)
    if not np.isfinite(ref):
        raise ValueError("ref_block must be finite")
    return (block - ref) / 10.0


def trend_permutation(
    z: npt.ArrayLike,
    block: npt.ArrayLike,
    n_perm: int = 3000,
    seed: int = 20260722,
    ref_block: Optional[float] = None,
) -> dict:
    """Test a linear trend in GPD log-scale by complete-block permutation.

    ``block`` identifies the exchangeable unit (for example, season or year).
    Every permutation reassigns time labels to complete blocks and retains all
    observations within their original block.  The reported Monte Carlo
    p-value uses the plus-one correction; ``p_permutation_mcse`` quantifies its
    simulation error.  The two-sided LR statistic is primary.
    """
    z, block = _validate_design(z, block)
    if not isinstance(n_perm, (int, np.integer)) or n_perm < 1:
        raise ValueError("n_perm must be a positive integer")
    t = _decades(block, ref_block)
    p1, l1 = _fit_ns(z, t, True)
    p0, l0 = _fit_ns(z, t, False)
    lr = max(0.0, 2.0 * (l0 - l1))
    rng = np.random.default_rng(seed)
    uniq = np.unique(block)
    warm = [p0[0], p0[1], 0.0]          # null MLE: constant scale, zero trend
    null = []
    for _ in range(n_perm):
        permuted = rng.permutation(uniq)
        mapping = dict(zip(uniq, permuted))
        tb = _decades(np.array([mapping[b] for b in block]), ref_block)
        fit = _fit_ns_from(z, tb, True, warm)
        if fit is None:
            try:
                fit = _fit_ns(z, tb, True)
            except RuntimeError:
                continue
        null.append(max(0.0, 2.0 * (l0 - fit[1])))
    if not null:
        raise RuntimeError("all permutation fits failed")
    null_arr = np.asarray(null, dtype=float)
    p_perm = float((1 + np.sum(null_arr >= lr)) / (1 + len(null_arr)))
    mcse = float(np.sqrt(p_perm * (1.0 - p_perm) / (len(null_arr) + 1)))
    span = float(t.max() - t.min())
    return {
        "trend_per_decade": float(p1[2]),
        "xi": float(p1[0]),
        "log_sigma0": float(p1[1]),
        "xi_null": float(p0[0]),
        "log_sigma_null": float(p0[1]),
        "sigma_change_pct": float(100.0 * (np.exp(p1[2] * span) - 1.0)),
        "LR": lr,
        "p_asymptotic": float(1.0 - chi2.cdf(lr, 1)),
        "p_permutation": p_perm,
        "p_permutation_mcse": mcse,
        "n_permutations": int(len(null_arr)),
        "permutation_unit": "complete block label",
        "_null": null_arr,
    }


def trend_power(
    z: npt.ArrayLike,
    block: npt.ArrayLike,
    trends: npt.ArrayLike,
    xi_true: Optional[float] = None,
    log_sigma_true: Optional[float] = None,
    crit: Optional[float] = None,
    n_rep: int = 300,
    seed: int = 20260722,
    ref_block: Optional[float] = None,
    n_perm_calibration: int = 499,
    alpha: float = 0.05,
) -> list:
    """Estimate power for signed log-scale trends at the observed design.

    Simulations hold the number and temporal layout of excesses fixed.  Thus
    the curve describes detectability under the fitted GPD simulation model;
    it does not establish power against arbitrary misspecification.
    """
    z, block = _validate_design(z, block)
    trends = np.asarray(trends, dtype=float)
    if trends.ndim != 1 or trends.size == 0 or np.any(~np.isfinite(trends)):
        raise ValueError("trends must be a non-empty finite 1-D sequence")
    if not isinstance(n_rep, (int, np.integer)) or n_rep < 1:
        raise ValueError("n_rep must be a positive integer")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")
    t = _decades(block, ref_block)
    span = float(t.max() - t.min())
    if xi_true is None or log_sigma_true is None or crit is None:
        res = trend_permutation(
            z,
            block,
            n_perm=n_perm_calibration,
            seed=seed,
            ref_block=ref_block,
        )
        xi_true = res["xi_null"] if xi_true is None else float(xi_true)
        log_sigma_true = (
            res["log_sigma_null"] if log_sigma_true is None else float(log_sigma_true)
        )
        crit = (
            float(np.quantile(res["_null"], 1 - alpha, method="higher"))
            if crit is None
            else float(crit)
        )
    if not all(np.isfinite(v) for v in (xi_true, log_sigma_true, crit)):
        raise ValueError("xi_true, log_sigma_true and crit must be finite")
    if xi_true <= -0.999 or crit < 0:
        raise ValueError("xi_true must exceed -0.999 and crit must be non-negative")

    t_mean = float(np.mean(t))
    out = []
    for trend in trends:
        # Common random numbers make signed effect comparisons less noisy.
        rng = np.random.default_rng(seed)
        # warm-start each replicate's fits from the generating parameters
        warm_v = [xi_true, log_sigma_true, trend]
        warm_c = [xi_true, log_sigma_true + trend * t_mean]
        rejected = successful = 0
        for _ in range(n_rep):
            sigma = np.exp(log_sigma_true + trend * t)
            u = rng.uniform(size=len(z))
            if abs(xi_true) < 1e-8:
                simulated = -sigma * np.log1p(-u)
            else:
                simulated = sigma / xi_true * ((1.0 - u) ** (-xi_true) - 1.0)
            fv = _fit_ns_from(simulated, t, True, warm_v)
            if fv is None:
                try:
                    fv = _fit_ns(simulated, t, True)
                except RuntimeError:
                    continue
            fc = _fit_ns_from(simulated, t, False, warm_c)
            if fc is None:
                try:
                    fc = _fit_ns(simulated, t, False)
                except RuntimeError:
                    continue
            successful += 1
            rejected += max(0.0, 2.0 * (fc[1] - fv[1])) >= crit
        power = rejected / successful if successful else np.nan
        mcse = np.sqrt(power * (1.0 - power) / successful) if successful else np.nan
        out.append(
            {
                "trend_per_decade": float(trend),
                "sigma_change_pct": round(100.0 * (np.exp(trend * span) - 1.0), 1),
                "power": round(float(power), 3),
                "power_mcse": round(float(mcse), 4),
                "n_successful": successful,
            }
        )
    return out


def _pchip_root(x, y, target):
    """Monotone-interpolated crossing of a power curve at ``target``.

    The detectable effect is the root of the monotone interpolant, not the first
    grid point above ``target``, so it does not depend on the arbitrary mesh.
    Returns ``None`` when the curve never reaches ``target``.
    """
    from scipy.interpolate import PchipInterpolator

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    y = np.maximum.accumulate(y)          # enforce monotonicity in the mesh
    if x.size < 2 or y[-1] < target:
        return None
    if y[0] >= target:
        return float(x[0])
    f = PchipInterpolator(x, y)
    lo, hi = float(x[0]), float(x[-1])
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi))


def _emd_interp(curve, sign, target, rng, reps):
    """Signed interpolated EMD and its Monte-Carlo uncertainty for one direction.

    The crossing is resampled by perturbing each power estimate within its
    simulation standard error, so the EMD is reported as an interval rather than
    an exact constant read off an arbitrary grid.
    """
    rows = [r for r in curve if np.sign(r["trend_per_decade"]) == sign]
    if len(rows) < 2:
        return None, None
    x = np.array([abs(r["trend_per_decade"]) for r in rows])
    y = np.array([r["power"] for r in rows])
    se = np.array([r.get("power_mcse", 0.0) for r in rows])
    root = _pchip_root(x, y, target)
    if root is None:
        return None, None
    draws = []
    for _ in range(reps):
        yj = np.clip(y + rng.standard_normal(y.size) * se, 0.0, 1.0)
        rj = _pchip_root(x, yj, target)
        if rj is not None:
            draws.append(sign * rj)
    ci = ([float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))]
          if draws else None)
    return float(sign * root), ci


def min_detectable_effect(
    z: npt.ArrayLike,
    block: npt.ArrayLike,
    target_power: float = 0.80,
    grid: Optional[npt.ArrayLike] = None,
    n_rep: int = 300,
    seed: int = 20260722,
    ref_block: Optional[float] = None,
    direction: str = "both",
    n_perm_calibration: int = 499,
    alpha: float = 0.05,
    emd_uncertainty_reps: int = 2000,
) -> dict:
    """Return detectable positive and/or negative per-decade effects.

    ``direction`` is ``"positive"``, ``"negative"``, or ``"both"``.  Two readings
    of the detectable effect are returned. ``mde_*`` is the smallest grid point
    reaching ``target_power`` (kept for compatibility). ``emd_*`` is the crossing
    of a monotone interpolant of the power curve, which does not depend on the
    grid spacing and carries a Monte-Carlo uncertainty interval (``emd_*_ci95``)
    obtained by resampling each power estimate within its simulation error. The
    compatibility field ``mde_per_decade`` is the signed grid effect with the
    smallest absolute magnitude reaching ``target_power``.
    """
    if not 0 < target_power < 1:
        raise ValueError("target_power must lie strictly between 0 and 1")
    if direction not in {"positive", "negative", "both"}:
        raise ValueError("direction must be 'positive', 'negative', or 'both'")
    if grid is None:
        magnitudes = np.array([0.02, 0.05, 0.08, 0.10, 0.13, 0.15, 0.20, 0.25, 0.30])
    else:
        magnitudes = np.abs(np.asarray(grid, dtype=float))
        if magnitudes.ndim != 1 or magnitudes.size == 0 or np.any(~np.isfinite(magnitudes)):
            raise ValueError("grid must be a non-empty finite 1-D sequence")
        magnitudes = np.unique(magnitudes[magnitudes > 0])
        if magnitudes.size == 0:
            raise ValueError("grid must contain at least one non-zero effect")
    effects = []
    if direction in {"negative", "both"}:
        effects.extend((-magnitudes[::-1]).tolist())
    if direction in {"positive", "both"}:
        effects.extend(magnitudes.tolist())
    curve = trend_power(
        z,
        block,
        effects,
        alpha=alpha,
        n_rep=n_rep,
        n_perm_calibration=n_perm_calibration,
        seed=seed,
        ref_block=ref_block,
    )

    def reached(sign):
        candidates = [
            row["trend_per_decade"]
            for row in curve
            if row["power"] >= target_power and np.sign(row["trend_per_decade"]) == sign
        ]
        return min(candidates, key=abs) if candidates else None

    neg = reached(-1)
    pos = reached(1)
    available = [x for x in (neg, pos) if x is not None]
    mde = min(available, key=abs) if available else None

    rng_u = np.random.default_rng(seed + 7919)
    emd_neg, emd_neg_ci = _emd_interp(curve, -1, target_power, rng_u, emd_uncertainty_reps)
    emd_pos, emd_pos_ci = _emd_interp(curve, 1, target_power, rng_u, emd_uncertainty_reps)
    emd_available = [e for e in (emd_neg, emd_pos) if e is not None]
    emd = min(emd_available, key=abs) if emd_available else None
    return {
        "mde_per_decade": mde,
        "mde_absolute": abs(mde) if mde is not None else None,
        "mde_negative": neg,
        "mde_positive": pos,
        "emd_per_decade": emd,
        "emd_negative": emd_neg,
        "emd_positive": emd_pos,
        "emd_negative_ci95": emd_neg_ci,
        "emd_positive_ci95": emd_pos_ci,
        "direction": direction,
        "target_power": target_power,
        "power_curve": curve,
    }


def block_bootstrap_trend_ci(
    z: npt.ArrayLike,
    block: npt.ArrayLike,
    n_boot: int = 1000,
    seed: int = 20260722,
    ref_block: Optional[float] = None,
) -> dict:
    """Descriptive percentile interval from a pairs block bootstrap.

    Complete observed blocks are sampled with replacement and assigned to the
    corresponding ordered bootstrap positions.  This complements rather than
    replaces the block-label permutation test.
    """
    z, block = _validate_design(z, block)
    if not isinstance(n_boot, (int, np.integer)) or n_boot < 1:
        raise ValueError("n_boot must be a positive integer")
    blocks = np.unique(block)
    by_block = {b: np.flatnonzero(block == b) for b in blocks}
    ordered_times = np.sort(blocks)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        sampled = rng.choice(blocks, size=len(blocks), replace=True)
        z_parts, t_parts = [], []
        for target_time, source_block in zip(ordered_times, sampled):
            idx = by_block[source_block]
            z_parts.append(z[idx])
            t_parts.append(np.full(idx.size, target_time))
        zb = np.concatenate(z_parts)
        tb = np.concatenate(t_parts)
        try:
            par, _ = _fit_ns(zb, _decades(tb, ref_block), True)
        except RuntimeError:
            continue
        draws.append(float(par[2]))
    if not draws:
        return {"ci95": [None, None], "n_boot": 0}
    return {
        "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
        "n_boot": len(draws),
    }
