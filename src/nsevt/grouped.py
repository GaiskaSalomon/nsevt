"""Interval-censored (grouped) generalized Pareto tail inference.

When exceedances are recorded on a discrete grid (for example wind speeds in
5 kt steps, or any rounded measurement), fitting the *continuous* GPD to the
rounded values biases the shape toward zero and can move the estimated finite
endpoint by a large margin. This module fits the GPD to the interval each
recorded value implies (its rounding cell), and reports profile-likelihood
intervals for the shape and for the endpoint that respect the discretisation.

Two honesty properties motivate the module:

* the interval-censored likelihood removes the discretisation bias, so the
  reported shape and endpoint refer to the estimator whose interval is quoted;
* the endpoint interval is obtained by *profiling the endpoint itself* (a
  reparameterisation of the GPD), not by substituting a profiled shape into
  ``M* = u - sigma/xi`` and not by a percentile bootstrap; the profile interval
  is invariant under reparameterisation and respects the non-linearity of the
  endpoint as ``xi`` approaches zero.

The numerics are ported from unit-tested research code and use NumPy and SciPy
only.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2

from ._types import GPDGroupedFitDict, GroupedEndpointCI, GroupedShapeCI

#: Finite value returned by the grouped negative log-likelihood on an infeasible
#: parameter (non-positive scale, support violated, an underflowed cell). A fit
#: whose best objective is not strictly below this never left the penalty region
#: and must be reported as a failure rather than accepted.
_PENALTY = 1e10


def _log_gpd_surv(z: np.ndarray, xi: float, sigma) -> np.ndarray:
    """Log GPD survival ``(-1/xi) log(1 + xi z / sigma)``.

    Uses the exact ``xi -> 0`` limit ``-z / sigma`` (so the likelihood is
    continuous through zero shape) and ``log1p`` for the general term, and
    returns ``-inf`` beyond the support. ``sigma`` may be a scalar or aligned
    with ``z``.
    """
    z = np.asarray(z, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    # the optimiser transiently visits extreme scales; a non-finite intermediate
    # there is expected and handled by the caller's feasibility checks
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        if abs(float(xi)) < 1e-8:
            return -z / sigma
        t = xi * z / sigma
        return np.where(t > -1.0, -np.log1p(t) / xi, -np.inf)


def _grouped_loglik(a, b, trunc, xi: float, sigma) -> float | None:
    """Summed ``log P(cell_i | X > trunc_i)`` for the grouped GPD, in log space.

    The cell probability is a difference of survivals; forming it as
    ``exp(log S(a)) - exp(log S(b))`` cancels catastrophically when the two
    survivals are close (a wide scale, a narrow cell, shape near zero). This
    keeps the whole computation in logs:
    ``log[S(a) - S(b)] = log S(a) + log1p(-exp(log S(b) - log S(a)))``.
    Returns ``None`` when a cell has no representable mass or the support fails.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    lsa = _log_gpd_surv(a, xi, sigma)
    lsb = _log_gpd_surv(b, xi, sigma)
    if np.any(~np.isfinite(lsa)):
        return None
    d = np.minimum(lsb - lsa, -1e-15)          # b > a, so log S(b) < log S(a)
    with np.errstate(over="ignore"):
        log_cell = lsa + np.log1p(-np.exp(d))
    if np.any(~np.isfinite(log_cell)):
        return None
    ll = float(np.sum(log_cell))
    trunc = np.asarray(trunc, dtype=float)
    if np.any(trunc > 0):
        lst = _log_gpd_surv(trunc, xi, sigma)
        if np.any(~np.isfinite(lst)):
            return None
        ll -= float(np.sum(np.broadcast_to(lst, np.shape(log_cell))))
    return ll


def interval_cells(
    values: npt.ArrayLike,
    threshold: float,
    grid: float | npt.ArrayLike = 5.0,
    tol: float = 1e-6,
) -> tuple:
    """Rounding cell and truncation point of each exceedance, in excess scale.

    Each recorded value ``x`` above ``threshold`` is treated as censored to the
    interval ``[x - g/2, x + g/2)`` of its recording precision ``g``, expressed
    as an excess over the threshold and clipped at zero. ``grid`` may be a single
    width (one precision for the whole sample) or several widths, in which case
    each value is assigned the *coarsest* grid it lies on (a value that is a
    multiple of 5 is also a multiple of 1, and the coarser cell is kept).

    A value enters the sample because its *recorded* mark is a grid point
    strictly above the threshold. The event the likelihood conditions on is
    therefore "the true excess is at least the lower rounding edge of the
    smallest admissible mark", ``max(ceil-grid(threshold) - g/2 - threshold, 0)``,
    returned as ``trunc``. When the threshold sits on the grid this equals
    ``g/2``; when it does not, assuming ``g/2`` conditions on too large an excess
    and can drive a cell's conditional probability above one.

    Returns ``(a, b, trunc)`` arrays on the excess scale.
    """
    m = np.asarray(values, dtype=float)
    if m.ndim != 1:
        raise ValueError("values must be a 1-D array")
    if m.size and np.any(~np.isfinite(m)):
        raise ValueError("values must be finite")
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    grid_arr = np.atleast_1d(np.asarray(grid, dtype=float))
    if grid_arr.size == 0 or np.any(~np.isfinite(grid_arr)) or np.any(grid_arr <= 0):
        raise ValueError("grid widths must be positive and finite")
    grids = tuple(float(g) for g in grid_arr)
    half = np.full(m.size, min(grids) / 2.0)
    for g in sorted(grids):
        on = np.abs(m / g - np.round(m / g)) < tol
        half = np.where(on, g / 2.0, half)
    z = m - float(threshold)
    a = np.maximum(z - half, 0.0)
    b = z + half
    # Left-truncation is the excess of the lower rounding edge of the smallest
    # grid point strictly above the threshold, at each observation's precision.
    # With an aligned threshold this is g/2 (the previous assumption); with a
    # shifted threshold it is smaller, and using g/2 there made the conditioning
    # denominator too small so a cell probability could exceed one.
    step = 2.0 * half
    first_mark = (np.floor(float(threshold) / step + tol) + 1.0) * step
    trunc = np.maximum(first_mark - half - float(threshold), 0.0)
    return a, b, trunc


def _validate_cells(cells: tuple, n: int) -> tuple:
    """Check a user-supplied ``(a, b, trunc)`` triple before it reaches the fit."""
    if not isinstance(cells, (tuple, list)) or len(cells) != 3:
        raise ValueError("cells must be an (a, b, trunc) triple")
    out = tuple(np.asarray(c, dtype=float) for c in cells)
    if any(c.ndim != 1 or c.size != n for c in out):
        raise ValueError("each cells array must be 1-D and aligned with the exceedances")
    if any(np.any(~np.isfinite(c)) for c in out):
        raise ValueError("cells must contain only finite values")
    a, b, trunc = out
    if np.any(a < 0) or np.any(b <= a) or np.any(trunc < 0) or np.any(trunc > b):
        raise ValueError("cells must satisfy 0 <= a < b and 0 <= trunc <= b")
    return out


def _grouped_nll(par, a, b, trunc):
    """Stationary grouped negative log-likelihood in ``(xi, log sigma)``."""
    xi = float(par[0])
    # keep the search in the GPD existence region xi > -1, matching the
    # continuous kernel; below it the interval-censored MLE is not regular and
    # the profile interval would not refer to the same estimator
    if not np.isfinite(xi) or xi <= -0.999:
        return _PENALTY
    sigma = np.exp(par[1])
    if not np.isfinite(sigma) or sigma <= 0:
        return _PENALTY
    ll = _grouped_loglik(a, b, trunc, xi, float(sigma))
    return _PENALTY if ll is None or not np.isfinite(ll) else -ll


def fit_gpd_grouped(
    values: npt.ArrayLike,
    threshold: float,
    grid: float | npt.ArrayLike = 5.0,
    cells: tuple | None = None,
    starts: npt.ArrayLike = (-0.4, -0.25, -0.1, 0.05),
) -> GPDGroupedFitDict:
    """Interval-censored GPD maximum-likelihood fit above ``threshold``.

    ``cells`` accepts the ``(a, b, trunc)`` triple from :func:`interval_cells`;
    otherwise the cells are built from ``grid``. Returns a dictionary with the
    shape, scale, endpoint and log-likelihood.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)):
        raise ValueError("values must be a finite 1-D array")
    if not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    z = values[values > threshold] - float(threshold)
    if z.size < 3:
        raise ValueError(f"only {z.size} exceedances over u={threshold}; need >= 3")
    if cells is None:
        cells = interval_cells(values[values > threshold], threshold, grid)
    else:
        cells = _validate_cells(cells, z.size)
    a, b, trunc = cells
    best = (np.inf, None)
    for x0 in np.asarray(starts, dtype=float):
        init = np.array([x0, np.log(max(float(np.mean(z)), 1e-3))])
        r = minimize(
            _grouped_nll, init, args=(a, b, trunc), method="Nelder-Mead",
            options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 200000, "maxfev": 200000},
        )
        if r.success and np.isfinite(r.fun) and r.fun < best[0]:
            best = (float(r.fun), r.x)
    if best[1] is None:
        raise RuntimeError("grouped GPD maximum-likelihood optimization failed")
    if best[0] >= _PENALTY:
        raise RuntimeError(
            "grouped GPD fit did not reach a feasible optimum: the likelihood "
            "stayed in the penalty region (invalid scale, support violation or "
            "an underflowed cell probability) from every start"
        )
    xi, sigma = float(best[1][0]), float(np.exp(best[1][1]))
    endpoint = float(threshold) - sigma / xi if xi < 0 else np.inf
    return {"xi": xi, "sigma": sigma, "endpoint": endpoint, "loglik": -best[0],
            "n": int(z.size)}


def profile_ci_xi_grouped(
    values: npt.ArrayLike,
    threshold: float,
    grid: float | npt.ArrayLike = 5.0,
    cells: tuple | None = None,
    level: float = 0.95,
    n_bisect: int = 40,
    xi_floor: float = -0.999,
    xi_ceil: float = 5.0,
    fit: GPDGroupedFitDict | None = None,
) -> GroupedShapeCI:
    """Profile-likelihood interval for the shape under the grouped likelihood.

    Inverts the one-degree-of-freedom likelihood-ratio statistic,
    ``2[ell_p(xi_hat) - ell_p(xi)] = chi2_{1,level}``, by searching outward from
    ``xi_hat`` on each side, so the interval always contains the estimator whose
    point value is reported. ``xi_floor`` (the GPD existence region) and
    ``xi_ceil`` bound the numerical search; a limit that reaches one of them is
    returned with ``lo_at_bound`` / ``hi_at_bound`` set, to distinguish a
    numerical search boundary from a statistical limit.
    """
    values = np.asarray(values, dtype=float)
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between 0 and 1")
    if cells is None:
        cells = interval_cells(values[values > threshold], threshold, grid)
    a, b, trunc = cells
    if fit is None:
        fit = fit_gpd_grouped(values, threshold, grid=grid, cells=cells)
    xi_hat = float(fit["xi"])
    log_sigma_hat = float(np.log(fit["sigma"]))

    def profile(xi):
        # profiling out the single scale is a smooth 1-D problem: a bounded
        # scalar search is faster and steadier than simplex over a 1-vector
        r = minimize_scalar(
            lambda ls: _grouped_nll(np.array([xi, ls]), a, b, trunc),
            bounds=(log_sigma_hat - 8.0, log_sigma_hat + 8.0),
            method="bounded", options={"xatol": 1e-9})
        return -float(r.fun)

    # reference the target at the re-profiled value at xi_hat so the interval is
    # guaranteed to contain xi_hat even when the fit optimiser stopped a hair off
    target = profile(xi_hat) - float(chi2.ppf(level, 1)) / 2.0

    def inside(xi):
        return profile(xi) >= target

    def find_side(direction: int) -> tuple:
        boundary = xi_floor if direction < 0 else xi_ceil
        if inside(boundary):
            return float(boundary), True
        inner, step = xi_hat, 0.02
        while True:
            outer = inner + direction * step
            if (direction < 0 and outer <= boundary) or (direction > 0 and outer >= boundary):
                outer = boundary
            if not inside(outer):
                break
            inner = outer
            step *= 1.6
        t_end, f_end = inner, outer
        for _ in range(n_bisect):
            mid = 0.5 * (t_end + f_end)
            if inside(mid):
                t_end = mid
            else:
                f_end = mid
        return 0.5 * (t_end + f_end), False

    lo, lo_at_bound = find_side(-1)
    hi, hi_at_bound = find_side(+1)
    return {"xi_hat": xi_hat, "ci": (float(lo), float(hi)),
            "lo_at_bound": bool(lo_at_bound), "hi_at_bound": bool(hi_at_bound),
            "level": level}


def _grouped_nll_nu(xi, nu, a, b, trunc):
    """Grouped negative log-likelihood in the ``(xi, nu)`` parameterisation.

    ``nu`` is the endpoint on the excess scale (``sigma = -xi * nu``), so the
    support condition is simply ``x < nu`` with no cancellation.
    """
    if nu <= float(np.max(a)) or xi >= 0 or xi <= -0.999:
        return 1e10

    def ls(x):
        return -np.log1p(-np.minimum(x, nu * (1 - 1e-12)) / nu) / xi

    la, lb, lc = ls(a), ls(b), ls(trunc)
    d = lb - la
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        log_cell = la + np.log1p(-np.exp(np.minimum(d, -1e-15)))
        out = float(np.sum(log_cell - lc))
    return 1e10 if not np.isfinite(out) else -out


def profile_endpoint_ci(
    values: npt.ArrayLike,
    threshold: float,
    grid: float | npt.ArrayLike = 5.0,
    cells: tuple | None = None,
    level: float = 0.95,
    n_bisect: int = 40,
    gap_init: float = 50.0,
    gap_cap: float = 1.0e7,
    fit: GPDGroupedFitDict | None = None,
) -> GroupedEndpointCI:
    """Profile-likelihood interval for the finite endpoint ``M*`` itself.

    Makes the endpoint a parameter (``nu`` on the excess scale) and profiles the
    shape out, so the interval respects the non-linearity of ``M* = u -
    sigma/xi`` where a percentile bootstrap and a Wald interval do not. Returns
    the interval on the ``M*`` scale.

    The upper search grows adaptively: if the profile is still above the LR
    threshold once the bracket reaches ``gap_cap`` the endpoint is not bounded
    at ``level`` and the upper limit is returned as ``inf`` with
    ``upper_at_bound=True``, instead of a finite number that is only the search
    boundary. When the unconstrained shape fit is not negative the point
    endpoint is ``inf`` as well and only the lower limit is finite.
    """
    values = np.asarray(values, dtype=float)
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between 0 and 1")
    if cells is None:
        cells = interval_cells(values[values > threshold], threshold, grid)
    a, b, trunc = (np.asarray(c, dtype=float) for c in cells)
    a_max = float(a.max())
    if fit is None:
        fit = fit_gpd_grouped(values, threshold, grid=grid, cells=cells)
    crit = float(chi2.ppf(level, 1)) / 2.0

    def prof(nu):
        r = minimize_scalar(lambda x: _grouped_nll_nu(x, nu, a, b, trunc),
                            bounds=(-0.999, -1e-8), method="bounded",
                            options={"xatol": 1e-10})
        return -float(r.fun), float(r.x)

    # locate the endpoint MLE on a log mesh, growing the range while the maximum
    # sits at the far edge (the profile still rising as nu -> inf, xi -> 0-)
    gap = float(gap_init)
    while True:
        mesh = a_max + np.exp(np.linspace(np.log(0.05), np.log(gap), 120))
        vals = np.array([prof(v)[0] for v in mesh])
        k = int(np.argmax(vals))
        if k < mesh.size - 1 or gap >= gap_cap:
            break
        gap *= 8.0
    rising_to_infinity = k == mesh.size - 1 and gap >= gap_cap
    unbounded_point = bool(fit["xi"] >= -1e-6 or rising_to_infinity)
    if k == mesh.size - 1:
        nu_hat, ll_max = float(mesh[-1]), float(vals[-1])
    else:
        r = minimize_scalar(lambda v: -prof(v)[0],
                            bounds=(mesh[max(k - 1, 0)], mesh[min(k + 1, mesh.size - 1)]),
                            method="bounded", options={"xatol": 1e-8})
        nu_hat, ll_max = float(r.x), -float(r.fun)
    target = ll_max - crit

    # lower limit: bisect between the largest observed cell edge and nu_hat
    lo_lim, hi_lim = a_max * (1 + 1e-9), nu_hat
    lower_at_bound = prof(lo_lim)[0] >= target
    if not lower_at_bound:
        for _ in range(n_bisect):
            mid = 0.5 * (lo_lim + hi_lim)
            if prof(mid)[0] >= target:
                hi_lim = mid
            else:
                lo_lim = mid
    nu_lo = lo_lim if lower_at_bound else 0.5 * (lo_lim + hi_lim)

    # upper limit: grow the bracket until the profile drops below target; a cap
    # hit means the endpoint is not identified from above at this level
    upper_at_bound = unbounded_point
    if not upper_at_bound:
        far = max(nu_hat * 2.0, nu_hat + gap_init)
        while prof(far)[0] >= target:
            far *= 8.0
            if far >= a_max + gap_cap:
                upper_at_bound = True
                break
    if upper_at_bound:
        nu_hi = np.inf
    else:
        lo_lim, hi_lim = nu_hat, far
        for _ in range(n_bisect):
            mid = 0.5 * (lo_lim + hi_lim)
            if prof(mid)[0] >= target:
                lo_lim = mid
            else:
                hi_lim = mid
        nu_hi = 0.5 * (lo_lim + hi_lim)

    xi_hat = prof(nu_hat)[1]
    endpoint = np.inf if unbounded_point else float(threshold + nu_hat)
    hi_val = np.inf if not np.isfinite(nu_hi) else float(threshold + nu_hi)
    return {"endpoint": endpoint,
            "ci": (float(threshold + nu_lo), hi_val),
            "upper_at_bound": bool(upper_at_bound),
            "lower_at_bound": bool(lower_at_bound),
            "xi": xi_hat, "sigma": float(-xi_hat * nu_hat),
            "level": level,
            "method": "profile likelihood on the reparameterised endpoint"}


@dataclass
class GroupedGPDFit:
    """Interval-censored GPD fit with profile intervals for shape and endpoint.

    ``xi_ci`` and ``endpoint_ci`` are at ``level`` (default 0.95). The historical
    names ``xi_ci95`` / ``endpoint_ci95`` remain as read-only aliases. A limit
    that only reached the numerical search boundary is marked in
    ``xi_ci_at_bound`` / ``endpoint_ci_at_bound`` (``(lower, upper)``); an
    unidentified upper endpoint limit is reported as ``inf``.
    """

    threshold: float
    n_exceedances: int
    xi: float
    sigma: float
    xi_ci: tuple
    endpoint: float
    endpoint_ci: tuple
    loglik: float
    level: float = 0.95
    xi_ci_at_bound: tuple = (False, False)
    endpoint_ci_at_bound: tuple = (False, False)

    @property
    def xi_ci95(self) -> tuple:
        """Historical alias of :attr:`xi_ci` (the interval is at :attr:`level`)."""
        return self.xi_ci

    @property
    def endpoint_ci95(self) -> tuple:
        """Historical alias of :attr:`endpoint_ci`."""
        return self.endpoint_ci

    @property
    def bounded_supported(self) -> bool:
        """Whether the profile interval for the shape lies wholly below zero."""
        return bool(np.isfinite(self.xi_ci[1]) and self.xi_ci[1] < 0)

    def summary(self) -> str:
        pct = f"{self.level:.0%}"
        shape_ci = f"({self.xi_ci[0]:.4f}, {self.xi_ci[1]:.4f})"
        if any(self.xi_ci_at_bound):
            shape_ci += " [search bound reached]"
        ep = f"{self.endpoint:.2f}" if np.isfinite(self.endpoint) else "inf"
        hi = self.endpoint_ci[1]
        ep_hi = "inf" if not np.isfinite(hi) else f"{hi:.2f}"
        ep_ci = f"[{self.endpoint_ci[0]:.2f}, {ep_hi}]"
        if self.endpoint_ci_at_bound[1]:
            ep_ci += " (upper limit not identified)"
        tail = "supported" if self.bounded_supported else f"not supported by the {pct} interval"
        return (f"grouped GPD fit (u={self.threshold:g}, n={self.n_exceedances}): "
                f"xi={self.xi:.4f} {shape_ci}, sigma={self.sigma:.3f}; "
                f"bounded tail {tail}; endpoint={ep} {ep_ci} (profile, {pct})")


def gpd_pot_grouped(
    values: npt.ArrayLike,
    threshold: float,
    grid: float | npt.ArrayLike = 5.0,
    level: float = 0.95,
) -> GroupedGPDFit:
    """Fit an interval-censored GPD to discretised exceedances over ``threshold``.

    ``grid`` is the recording precision (a single width, or several widths for a
    mixed-precision record). ``level`` sets the profile-interval level (default
    0.95). Returns a :class:`GroupedGPDFit` with the shape and endpoint each
    carrying a profile-likelihood interval at ``level`` and a boundary flag.
    """
    if not 0 < level < 1:
        raise ValueError("level must lie strictly between 0 and 1")
    values = np.asarray(values, dtype=float)
    cells = interval_cells(values[values > threshold], threshold, grid)
    fit = fit_gpd_grouped(values, threshold, grid=grid, cells=cells)
    xi_prof = profile_ci_xi_grouped(values, threshold, grid=grid, cells=cells,
                                    level=level, fit=fit)
    ep = profile_endpoint_ci(values, threshold, grid=grid, cells=cells,
                             level=level, fit=fit)
    return GroupedGPDFit(
        threshold=float(threshold), n_exceedances=fit["n"], xi=fit["xi"],
        sigma=fit["sigma"],
        xi_ci=(float(xi_prof["ci"][0]), float(xi_prof["ci"][1])),
        endpoint=fit["endpoint"], endpoint_ci=tuple(ep["ci"]), loglik=fit["loglik"],
        level=float(level),
        xi_ci_at_bound=(bool(xi_prof["lo_at_bound"]), bool(xi_prof["hi_at_bound"])),
        endpoint_ci_at_bound=(bool(ep["lower_at_bound"]), bool(ep["upper_at_bound"])),
    )
