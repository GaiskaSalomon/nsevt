"""Peaks-over-threshold generalized Pareto (GPD) tail estimation.

For a threshold ``u``, positive excesses ``Z = X - u | X > u`` are modeled
with a two-parameter GPD.  A negative fitted shape implies a finite *statistical*
endpoint under the selected threshold and GPD model.  It is not, by itself,
evidence of a physical upper bound.

The module deliberately separates a negative point estimate from inferential
support for a negative shape.  The latter requires the entire profile-
likelihood confidence interval to lie below zero.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.stats import chi2

_XI0 = (-0.4, -0.25, -0.1, 0.05, 0.2)


def _validate_excesses(z: Sequence[float], minimum: int = 3) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size < minimum:
        raise ValueError(f"z must be a 1-D array of at least {minimum} excesses")
    if np.any(~np.isfinite(z)):
        raise ValueError("excesses z must all be finite")
    if np.any(z < 0):
        raise ValueError("excesses z must be non-negative (z = x - u for x > u)")
    if not np.any(z > 0) or np.ptp(z) == 0:
        raise ValueError("excesses z must contain positive variation")
    return z


def _nll_const(par, z):
    xi, log_s = par
    # For xi <= -1 the likelihood may be unbounded at the sample maximum.
    # Keep optimization in the conventional existence region xi > -1.
    if not np.isfinite(xi) or xi <= -0.999 or not np.isfinite(log_s):
        return 1e100
    s = np.exp(log_s)
    y = 1.0 + xi * z / s
    if not np.isfinite(s) or np.any(y <= 0):
        return 1e100
    if abs(xi) < 1e-8:
        return len(z) * log_s + np.sum(z) / s
    return len(z) * log_s + (1.0 + 1.0 / xi) * np.sum(np.log(y))


def fit_gpd(z: Sequence[float]) -> dict:
    """Return the two-parameter GPD maximum-likelihood fit of ``z >= 0``.

    Multiple starts reduce sensitivity to local numerical failures.  A
    ``RuntimeError`` is raised when no optimizer run yields a finite feasible
    solution; invalid input is never converted silently into a fit.
    """
    z = _validate_excesses(z)
    starts = {float(np.mean(z)), float(np.std(z, ddof=0))}
    best, bp = np.inf, None
    for xi0 in _XI0:
        for s0 in starts:
            if not np.isfinite(s0) or s0 <= 0:
                continue
            r = minimize(
                _nll_const,
                [xi0, np.log(s0)],
                args=(z,),
                method="Nelder-Mead",
                options={"xatol": 1e-9, "fatol": 1e-9, "maxiter": 20000},
            )
            if r.success and np.isfinite(r.fun) and r.fun < best:
                best, bp = float(r.fun), np.asarray(r.x, dtype=float)
    if bp is None:
        raise RuntimeError("GPD maximum-likelihood optimization failed")
    return {"xi": float(bp[0]), "sigma": float(np.exp(bp[1])), "nll": best}


def _fit_gpd_from(z, xi0, log_s0):
    """Single-start GPD fit from a given ``(xi, log sigma)`` guess.

    Used to warm-start bootstrap resamples from the point estimate, which lies
    close to every resample's optimum, so one start suffices; ``None`` signals
    that the caller should fall back to the multi-start :func:`fit_gpd`.
    """
    r = minimize(
        _nll_const, [xi0, log_s0], args=(z,), method="Nelder-Mead",
        options={"xatol": 1e-9, "fatol": 1e-9, "maxiter": 20000},
    )
    if r.success and np.isfinite(r.fun):
        return {"xi": float(r.x[0]), "sigma": float(np.exp(r.x[1])), "nll": float(r.fun)}
    return None


def _profile_nll(xi: float, z: np.ndarray, sigma_hint: float) -> float:
    """Profile negative log-likelihood at fixed ``xi``."""
    lower = max(0.0, -float(xi) * float(np.max(z)))
    scale_ref = max(float(sigma_hint), float(np.mean(z)), 1e-8)

    # sigma = lower + exp(eta) enforces the support constraint for xi < 0.
    def objective(eta):
        sigma = lower + np.exp(float(eta))
        return _nll_const((xi, np.log(sigma)), z)

    centre = np.log(max(scale_ref - lower, scale_ref * 1e-3, 1e-10))
    r = minimize_scalar(
        objective,
        bounds=(centre - 18.0, centre + 18.0),
        method="bounded",
        options={"xatol": 1e-9, "maxiter": 2000},
    )
    return float(r.fun) if r.success and np.isfinite(r.fun) else np.inf


def _profile_ci_details(
    z: Sequence[float],
    level: float = 0.95,
    grid: Optional[np.ndarray] = None,
) -> tuple[float, float, tuple[float, float], tuple[bool, bool]]:
    z = _validate_excesses(z)
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between 0 and 1")
    fit = fit_gpd(z)
    xih, sh, l0 = fit["xi"], fit["sigma"], fit["nll"]
    crit = float(chi2.ppf(level, 1))

    def equation(xi):
        return 2.0 * (_profile_nll(float(xi), z, sh) - l0) - crit

    if grid is not None:
        grid = np.asarray(grid, dtype=float)
        if grid.ndim != 1 or grid.size < 3 or np.any(~np.isfinite(grid)):
            raise ValueError("grid must be a finite 1-D array with at least 3 values")
        grid = np.unique(np.append(grid, xih))
        vals = np.array([equation(x) for x in grid])
        accepted = grid[vals <= 0]
        if accepted.size == 0:
            return xih, sh, (np.nan, np.nan), (True, True)
        lo, hi = float(accepted.min()), float(accepted.max())
        return xih, sh, (lo, hi), (lo == grid.min(), hi == grid.max())

    # Bracket each likelihood-ratio root adaptively.  These wide numerical
    # bounds avoid the old silent clipping at [-0.9, 0.9].
    bounds = (-0.999, 5.0)

    def find_side(direction: int) -> tuple[float, bool]:
        inner = float(xih)
        step = 0.025
        boundary = bounds[0] if direction < 0 else bounds[1]
        while True:
            outer = inner + direction * step
            if (direction < 0 and outer <= boundary) or (direction > 0 and outer >= boundary):
                outer = boundary
            value = equation(outer)
            if np.isfinite(value) and value >= 0:
                a, b = sorted((inner, outer))
                return float(brentq(equation, a, b, xtol=1e-8, rtol=1e-8)), False
            if outer == boundary:
                return float(boundary), True
            inner = outer
            step *= 1.6

    lo, lo_truncated = find_side(-1)
    hi, hi_truncated = find_side(1)
    return xih, sh, (lo, hi), (lo_truncated, hi_truncated)


def profile_ci_xi(
    z: Sequence[float], level: float = 0.95, grid: Optional[np.ndarray] = None
) -> tuple:
    """Profile-likelihood confidence interval for GPD shape ``xi``.

    The interval inverts the one-degree-of-freedom likelihood-ratio statistic.
    This is an asymptotic interval; the usual GPD likelihood regularity becomes
    problematic at and below ``xi = -1/2``.  A warning is emitted when the fit
    enters that regime.  The return value is ``(xi_hat, sigma_hat, (lo, hi))``.

    ``grid`` is retained for reproducibility experiments.  The default uses
    adaptive root bracketing and does not silently truncate at the old fixed
    grid limits.
    """
    xih, sh, ci, truncated = _profile_ci_details(z, level=level, grid=grid)
    if xih <= -0.5:
        warnings.warn(
            "xi_hat <= -0.5: standard profile-likelihood asymptotics may be unreliable",
            RuntimeWarning,
            stacklevel=2,
        )
    if any(truncated):
        warnings.warn(
            "profile-likelihood confidence limit reached the numerical search boundary",
            RuntimeWarning,
            stacklevel=2,
        )
    return xih, sh, ci


def upper_endpoint(
    z: Sequence[float],
    threshold: float,
    n_boot: int = 2000,
    seed: int = 20260722,
) -> dict:
    """Estimate the model-conditional endpoint and its conditional bootstrap CI.

    The endpoint ``u - sigma/xi`` exists for a negative fitted shape.  The
    percentile interval is computed only among bootstrap fits with ``xi < 0``;
    ``bootstrap_fraction_xi_negative`` must therefore be reported alongside it.
    It is not a confidence interval for a physical limit.
    """
    z = _validate_excesses(z)
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if not isinstance(n_boot, (int, np.integer)) or n_boot < 0:
        raise ValueError("n_boot must be a non-negative integer")
    fit = fit_gpd(z)
    xi, sig = fit["xi"], fit["sigma"]
    log_sig = float(np.log(sig))
    endpoint = threshold - sig / xi if xi < 0 else np.inf
    rng = np.random.default_rng(seed)
    ends, xis = [], []
    for _ in range(n_boot):
        zb = rng.choice(z, size=len(z), replace=True)
        if np.ptp(zb) == 0:
            continue
        # warm-start each resample from the point estimate; fall back to the
        # multi-start fit only if the single start fails to converge
        f = _fit_gpd_from(zb, xi, log_sig)
        if f is None:
            try:
                f = fit_gpd(zb)
            except (RuntimeError, ValueError):
                continue
        xis.append(f["xi"])
        if f["xi"] < 0:
            ends.append(threshold - f["sigma"] / f["xi"])
    ends = np.asarray(ends, dtype=float)
    ci = (
        [float(np.percentile(ends, 2.5)), float(np.percentile(ends, 97.5))]
        if ends.size
        else [np.nan, np.nan]
    )
    return {
        "xi": xi,
        "sigma": sig,
        "endpoint": float(endpoint),
        "endpoint_ci": ci,
        "bootstrap_fraction_xi_negative": (
            float(np.mean(np.asarray(xis) < 0)) if xis else np.nan
        ),
        "n_boot_successful": len(xis),
    }


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
    n_boot_successful: int = 0
    xi_ci95_truncated: tuple = (False, False)

    @property
    def bounded_estimate(self) -> bool:
        """Whether the point estimate implies a finite statistical endpoint."""
        return bool(self.xi < 0)

    @property
    def bounded_supported(self) -> bool:
        """Whether the 95% profile interval lies wholly below zero."""
        return bool(np.isfinite(self.xi_ci95[1]) and self.xi_ci95[1] < 0)

    @property
    def bounded(self) -> bool:
        """Compatibility alias for :attr:`bounded_supported` (since v0.2)."""
        return self.bounded_supported

    def return_level(self, return_period: float, rate: Optional[float] = None) -> float:
        """Return the model-based level for a positive return period.

        ``rate`` is the threshold exceedance probability on the scale of all
        observations.  Without it, the return period is conditional on an
        exceedance.  Extrapolation far beyond the observed record remains
        model-dependent even when a finite endpoint is estimated.
        """
        if not np.isfinite(return_period) or return_period <= 1:
            raise ValueError("return_period must be finite and greater than 1")
        if rate is not None and (not np.isfinite(rate) or not 0 < rate <= 1):
            raise ValueError("rate must lie in (0, 1]")
        p = 1.0 / return_period
        q = p / rate if rate is not None else p
        if q >= 1:
            return float(self.threshold)
        if abs(self.xi) < 1e-8:
            return float(self.threshold - self.sigma * np.log(q))
        value = self.threshold + self.sigma / self.xi * (q ** (-self.xi) - 1)
        if self.bounded_estimate:
            value = min(value, self.endpoint)
        return float(value)

    def summary(self) -> str:
        if self.bounded_supported:
            tail = "95% profile interval supports xi < 0"
        elif self.bounded_estimate:
            tail = "xi_hat < 0, but its 95% profile interval includes zero"
        else:
            tail = "xi_hat >= 0"
        endpoint = f"{self.endpoint:.2f}" if np.isfinite(self.endpoint) else "inf"
        shape_ci = f"({self.xi_ci95[0]:.4f}, {self.xi_ci95[1]:.4f})"
        endpoint_ci = (
            f"[{self.endpoint_ci95[0]:.3f}, {self.endpoint_ci95[1]:.3f}]"
            if np.all(np.isfinite(self.endpoint_ci95))
            else "[not available]"
        )
        return (
            f"GPD POT fit (u={self.threshold:g}, n={self.n_exceedances}): "
            f"xi={self.xi:.3f} {shape_ci}, sigma={self.sigma:.3f}; "
            f"{tail}; model-conditional endpoint={endpoint} {endpoint_ci}; "
            f"bootstrap fraction xi<0={self.bootstrap_fraction_xi_negative:.3f} "
            f"(successful B={self.n_boot_successful})"
        )


def gpd_pot(
    values: Sequence[float],
    threshold: float,
    n_boot: int = 2000,
    seed: int = 20260722,
) -> GPDFit:
    """Fit a GPD to raw observations strictly above ``threshold``."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)):
        raise ValueError("values must be a finite 1-D array")
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    z = values[values > threshold] - threshold
    if z.size < 3:
        raise ValueError(f"only {z.size} exceedances over u={threshold}; need >= 3")
    xih, sh, ci, truncated = _profile_ci_details(z)
    if xih <= -0.5:
        warnings.warn(
            "xi_hat <= -0.5: standard profile-likelihood asymptotics may be unreliable",
            RuntimeWarning,
            stacklevel=2,
        )
    if any(truncated):
        warnings.warn(
            "profile-likelihood confidence limit reached the numerical search boundary",
            RuntimeWarning,
            stacklevel=2,
        )
    ep = upper_endpoint(z, threshold, n_boot=n_boot, seed=seed)
    return GPDFit(
        threshold=float(threshold),
        n_exceedances=int(z.size),
        xi=xih,
        sigma=sh,
        xi_ci95=(float(ci[0]), float(ci[1])),
        endpoint=ep["endpoint"],
        endpoint_ci95=[float(x) for x in ep["endpoint_ci"]],
        bootstrap_fraction_xi_negative=ep["bootstrap_fraction_xi_negative"],
        n_boot_successful=ep["n_boot_successful"],
        xi_ci95_truncated=truncated,
    )
