"""Peaks-over-threshold generalized Pareto (GPD) tail estimation.

The excess ``Z = X - u`` of a variable over a high threshold ``u`` is modeled
by a two-parameter GPD with survival function

    P(Z > z | Z > 0) = (1 + xi z / sigma)_+^(-1/xi),   sigma > 0,

estimated by maximum likelihood.  A negative shape (``xi < 0``) implies a
*bounded* tail with finite upper endpoint ``x* = u - sigma/xi``; this is the
finite-ceiling regime central to the package.

Numerics (multi-start Nelder--Mead, profile-likelihood interval for ``xi``,
storm/observation bootstrap of the endpoint) are ported verbatim from the
frozen, unit-tested implementation used in the thesis, so results are identical.

References
----------
Coles (2001), *An Introduction to Statistical Modeling of Extreme Values*;
Smith (1985), Biometrika 72:67-90 (regularity for ``xi > -1/2``).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

_XI0 = (-0.4, -0.25, -0.1, 0.05, 0.2)


# --------------------------------------------------------------------------- #
# constant-scale GPD log-likelihood and MLE
# --------------------------------------------------------------------------- #
def _nll_const(par, z):
    xi, log_s = par
    s = np.exp(log_s)
    y = 1.0 + xi * z / s
    if np.any(y <= 1e-12):
        return 1e10
    if abs(xi) < 1e-8:
        return len(z) * log_s + np.sum(z) / s
    return len(z) * log_s + (1.0 + 1.0 / xi) * np.sum(np.log(y))


def fit_gpd(z: Sequence[float]) -> dict:
    """Two-parameter GPD MLE of the excesses ``z >= 0``.

    Returns ``dict(xi, sigma, nll)``.
    """
    z = np.asarray(z, float)
    if z.ndim != 1 or z.size < 3:
        raise ValueError("z must be a 1-D array of at least 3 excesses")
    if np.any(z < 0):
        raise ValueError("excesses z must be non-negative (z = x - u for x > u)")
    best, bp = np.inf, None
    for xi0 in _XI0:
        for s0 in (np.mean(z), np.std(z) + 1e-6):
            if s0 <= 0:
                continue
            r = minimize(_nll_const, [xi0, np.log(s0)], args=(z,),
                         method="Nelder-Mead",
                         options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 20000})
            if r.fun < best:
                best, bp = r.fun, r.x
    return dict(xi=float(bp[0]), sigma=float(np.exp(bp[1])), nll=float(best))


def profile_ci_xi(z: Sequence[float], level: float = 0.95,
                  grid: Optional[np.ndarray] = None) -> tuple:
    """Profile-likelihood confidence interval for the shape ``xi``.

    Inverts the profile deviance against ``chi2_1`` (Wilks); the interval
    respects the asymmetry of the likelihood and, unlike a Wald interval, does
    not impose constant curvature near the ``xi = 0`` boundary.  Returns
    ``(xi_hat, sigma_hat, (lo, hi))``.
    """
    z = np.asarray(z, float)
    fit = fit_gpd(z)
    xih, sh, l0 = fit["xi"], fit["sigma"], fit["nll"]
    crit = chi2.ppf(level, 1) / 2.0
    if grid is None:
        grid = np.linspace(-0.9, 0.9, 361)

    def pnll(xi):
        best = np.inf
        for s0 in (np.mean(z), np.std(z) + 1e-6, sh):
            r = minimize(lambda ls: _nll_const([xi, ls[0]], z),
                         [np.log(max(s0, 1e-6))], method="Nelder-Mead",
                         options={"xatol": 1e-10, "fatol": 1e-10})
            best = min(best, r.fun)
        return best

    dev = np.array([pnll(x) - l0 for x in grid])
    ok = grid[dev <= crit]
    ci = (float(ok.min()), float(ok.max())) if len(ok) else (np.nan, np.nan)
    return xih, sh, ci


def upper_endpoint(z: Sequence[float], threshold: float,
                   n_boot: int = 2000, seed: int = 20260722) -> dict:
    """Finite upper endpoint ``x* = u - sigma/xi`` with a percentile bootstrap.

    The endpoint is finite iff ``xi < 0``.  ``bootstrap_fraction_xi_negative``
    is the share of resamples with a bounded tail: a direct read of how firmly
    the finite-ceiling conclusion holds.
    """
    z = np.asarray(z, float)
    fit = fit_gpd(z)
    xi, sig = fit["xi"], fit["sigma"]
    endpoint = threshold - sig / xi if xi < 0 else np.inf
    rng = np.random.default_rng(seed)
    ends, xis = [], []
    for _ in range(n_boot):
        zb = rng.choice(z, size=len(z), replace=True)
        try:
            f = fit_gpd(zb)
            xis.append(f["xi"])
            if f["xi"] < -1e-3:
                ends.append(threshold - f["sigma"] / f["xi"])
        except Exception:
            pass
    ends = np.asarray(ends)
    ci = ([float(np.percentile(ends, 2.5)), float(np.percentile(ends, 97.5))]
          if ends.size else [np.nan, np.nan])
    return dict(xi=xi, sigma=sig, endpoint=float(endpoint), endpoint_ci=ci,
                bootstrap_fraction_xi_negative=float(np.mean(np.asarray(xis) < 0))
                if xis else np.nan)


# --------------------------------------------------------------------------- #
# high-level container
# --------------------------------------------------------------------------- #
@dataclass
class GPDFit:
    """Result of a peaks-over-threshold GPD fit."""
    threshold: float
    n_exceedances: int
    xi: float
    sigma: float
    xi_ci95: tuple
    endpoint: float
    endpoint_ci95: list
    bootstrap_fraction_xi_negative: float
    bounded: bool = field(init=False)

    def __post_init__(self):
        self.bounded = self.xi < 0

    def return_level(self, return_period: float, rate: Optional[float] = None) -> float:
        """Return level for a 1-in-``return_period`` event.

        ``rate`` is P(X > u) (the exceedance probability of the threshold); if
        given, the return level is on the scale of *all* observations, otherwise
        it is conditional on exceeding ``u``.
        """
        p = 1.0 / return_period
        q = p / rate if rate is not None else p
        if q >= 1:
            return float(self.threshold)
        if abs(self.xi) < 1e-8:
            return float(self.threshold - self.sigma * np.log(q))
        return float(self.threshold + self.sigma / self.xi * (q ** (-self.xi) - 1))

    def summary(self) -> str:
        tail = "bounded (finite ceiling)" if self.bounded else "unbounded"
        end = f"{self.endpoint:.2f}" if np.isfinite(self.endpoint) else "inf"
        return (f"GPD POT fit (u={self.threshold:g}, n={self.n_exceedances}): "
                f"xi={self.xi:.3f} {self.xi_ci95}, sigma={self.sigma:.3f}; "
                f"tail {tail}, endpoint={end} {self.endpoint_ci95}; "
                f"P(xi<0) bootstrap={self.bootstrap_fraction_xi_negative:.3f}")


def gpd_pot(values: Sequence[float], threshold: float,
            n_boot: int = 2000, seed: int = 20260722) -> GPDFit:
    """Fit a GPD to the peaks over ``threshold`` of ``values``.

    Parameters
    ----------
    values : array
        Raw observations (not excesses).
    threshold : float
        The tail threshold ``u``; excesses are ``z = values[values > u] - u``.
    """
    values = np.asarray(values, float)
    z = values[values > threshold] - threshold
    if z.size < 3:
        raise ValueError(f"only {z.size} exceedances over u={threshold}; need >= 3")
    xih, sh, ci = profile_ci_xi(z)
    ep = upper_endpoint(z, threshold, n_boot=n_boot, seed=seed)
    return GPDFit(threshold=float(threshold), n_exceedances=int(z.size),
                  xi=xih, sigma=sh, xi_ci95=(round(ci[0], 4), round(ci[1], 4)),
                  endpoint=ep["endpoint"], endpoint_ci95=[round(x, 3) for x in ep["endpoint_ci"]],
                  bootstrap_fraction_xi_negative=ep["bootstrap_fraction_xi_negative"])
