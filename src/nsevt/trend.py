"""Non-stationary GPD tail scale: calibrated trend test and detectability.

A time-varying log-scale ``log sigma(t) = beta0 + beta1 t`` (with ``t`` in
decades from the first block) lets the tail widen or narrow over time.  The
inferential question is treated honestly:

* the trend is tested by a **season/block permutation** test that is exact in
  finite samples under exchangeability of the block labels, avoiding the
  asymptotic chi-square that is unreliable for a boundary-adjacent parameter on
  a few hundred exceedances;
* a **Monte-Carlo power / minimum-detectable-effect (MDE)** analysis turns any
  non-rejection into a quantitative statement of what the record can resolve.

Ported from the frozen thesis implementation.
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from scipy.stats import chi2

_XI0 = (-0.4, -0.25, -0.1, 0.05, 0.2)


def _nll_ns(par, z, t, varying):
    from scipy.optimize import minimize  # noqa: F401  (kept local; see _fit_ns)
    xi = par[0]
    s = np.exp(par[1] + par[2] * t) if varying else np.exp(np.full_like(t, par[1]))
    if np.any(~np.isfinite(s)):
        return 1e10
    if abs(xi) < 1e-8:
        return np.sum(np.log(s) + z / s)
    y = 1.0 + xi * z / s
    if np.any(y <= 1e-12):
        return 1e10
    return np.sum(np.log(s)) + (1.0 + 1.0 / xi) * np.sum(np.log(y))


def _fit_ns(z, t, varying):
    from scipy.optimize import minimize
    best, bp = np.inf, None
    for xi0 in _XI0:
        p0 = [xi0, np.log(max(np.mean(z), 1e-3)), 0.0]
        r = minimize(_nll_ns, p0, args=(z, t, varying), method="Nelder-Mead",
                     options={"maxiter": 60000, "maxfev": 60000,
                              "xatol": 1e-10, "fatol": 1e-10})
        if r.fun < best:
            best, bp = r.fun, r.x
    return bp, float(best)


def _decades(block, ref: Optional[float]) -> np.ndarray:
    block = np.asarray(block, float)
    ref = float(np.min(block)) if ref is None else float(ref)
    return (block - ref) / 10.0


def trend_permutation(z: Sequence[float], block: Sequence[float],
                      n_perm: int = 3000, seed: int = 20260722,
                      ref_block: Optional[float] = None) -> dict:
    """Permutation-calibrated test of a time trend in the GPD log-scale.

    Parameters
    ----------
    z : array
        Excesses over the threshold (``x - u``).
    block : array
        The block label (e.g. year/season) of each excess; blocks are the
        exchangeable unit reshuffled under the null.
    n_perm : int
        Number of block-label permutations.

    Returns a dict with ``trend_per_decade``, ``sigma_change_pct`` over the
    record, the likelihood-ratio statistic ``LR``, the asymptotic and the
    (primary) permutation p-values, and ``_null`` (the permutation LR draws,
    reused by :func:`min_detectable_effect`).
    """
    z = np.asarray(z, float)
    block = np.asarray(block)
    t = _decades(block, ref_block)
    p1, l1 = _fit_ns(z, t, True)
    _, l0 = _fit_ns(z, t, False)
    lr = 2 * (l0 - l1)
    rng = np.random.default_rng(seed)
    uniq = np.unique(block)
    null = []
    for _ in range(n_perm):
        mp = dict(zip(uniq, rng.permutation(uniq)))
        tb = _decades(np.array([mp[b] for b in block]), ref_block)
        try:
            _, la = _fit_ns(z, tb, True)
            null.append(2 * (l0 - la))  # stationary fit is permutation-invariant
        except Exception:
            pass
    null = np.asarray(null)
    p_perm = (1 + np.sum(null >= lr)) / (1 + len(null))
    span = t.max() - t.min()
    return dict(trend_per_decade=float(p1[2]), xi=float(p1[0]),
                log_sigma0=float(p1[1]),
                sigma_change_pct=float(100 * (np.exp(p1[2] * span) - 1)),
                LR=float(lr), p_asymptotic=float(1 - chi2.cdf(lr, 1)),
                p_permutation=float(p_perm), n_permutations=int(len(null)),
                _null=null)


def trend_power(z: Sequence[float], block: Sequence[float], trends: Sequence[float],
                xi_true: Optional[float] = None, log_sigma_true: Optional[float] = None,
                crit: Optional[float] = None, n_rep: int = 300,
                seed: int = 20260722, ref_block: Optional[float] = None) -> list:
    """Monte-Carlo power to detect each trend in ``trends`` at the observed design.

    Data are simulated under the fitted null shape/scale with an injected trend;
    each replicate is tested against ``crit`` (the 95th percentile of the
    permutation null by default).  Holds ``n`` and the block layout fixed at the
    observed values, so power reflects the *record*, not an idealized design.
    """
    z = np.asarray(z, float)
    block = np.asarray(block)
    t = _decades(block, ref_block)
    span = t.max() - t.min()
    n = len(z)
    if xi_true is None or log_sigma_true is None or crit is None:
        res = trend_permutation(z, block, seed=seed, ref_block=ref_block)
        xi_true = res["xi"] if xi_true is None else xi_true
        log_sigma_true = res["log_sigma0"] if log_sigma_true is None else log_sigma_true
        crit = float(np.percentile(res["_null"], 95)) if crit is None else crit
    rng = np.random.default_rng(seed)
    out = []
    for tr in trends:
        rej, R = 0, n_rep
        for _ in range(R):
            s = np.exp(log_sigma_true + tr * t)
            u = rng.uniform(size=n)
            if abs(xi_true) < 1e-8:
                zs = -s * np.log(1 - u)
            else:
                zs = s / xi_true * ((1 - u) ** (-xi_true) - 1)
            try:
                _, la = _fit_ns(zs, t, True)
                _, lb = _fit_ns(zs, t, False)
                if 2 * (lb - la) >= crit:
                    rej += 1
            except Exception:
                R -= 1
        out.append(dict(trend_per_decade=float(tr),
                        sigma_change_pct=round(100 * (np.exp(tr * span) - 1), 1),
                        power=round(rej / max(R, 1), 3)))
    return out


def min_detectable_effect(z: Sequence[float], block: Sequence[float],
                          target_power: float = 0.80,
                          grid: Optional[Sequence[float]] = None,
                          n_rep: int = 300, seed: int = 20260722,
                          ref_block: Optional[float] = None) -> dict:
    """Smallest per-decade trend detectable at ``target_power`` (default 0.80)."""
    if grid is None:
        grid = [0.02, 0.05, 0.08, 0.10, 0.13, 0.15, 0.20, 0.25, 0.30]
    curve = trend_power(z, block, grid, n_rep=n_rep, seed=seed, ref_block=ref_block)
    ok = [c["trend_per_decade"] for c in curve if c["power"] >= target_power]
    return dict(mde_per_decade=(min(ok) if ok else None),
                target_power=target_power, power_curve=curve)


def block_bootstrap_trend_ci(z: Sequence[float], block: Sequence[float],
                             n_boot: int = 1000, seed: int = 20260722,
                             ref_block: Optional[float] = None) -> dict:
    """Block bootstrap interval for the log-scale trend coefficient.

    Entire blocks are resampled so shared within-block conditions are not
    treated as independent replication.  Descriptive; complements the
    permutation test rather than replacing it.
    """
    z = np.asarray(z, float)
    block = np.asarray(block)
    t = _decades(block, ref_block)
    blocks = np.unique(block)
    by_block = {b: np.flatnonzero(block == b) for b in blocks}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        sampled = rng.choice(blocks, size=len(blocks), replace=True)
        idx = np.concatenate([by_block[b] for b in sampled])
        try:
            par, _ = _fit_ns(z[idx], t[idx], True)
            draws.append(float(par[2]))
        except Exception:
            pass
    if not draws:
        return {"ci95": [None, None], "n_boot": 0}
    return {"ci95": [float(np.percentile(draws, 2.5)),
                     float(np.percentile(draws, 97.5))],
            "n_boot": int(len(draws))}
