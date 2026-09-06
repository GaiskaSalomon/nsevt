"""Tests for the interval-censored (grouped) GPD module."""
import numpy as np
import pytest

import nsevt
from nsevt import grouped
from nsevt.grouped import interval_cells


def _discretised_bounded(rng, xi=-0.20, sigma=15.0, u=40.0, n=6000, grid=5.0):
    """True GPD excesses over ``u`` with marks rounded to ``grid``."""
    uu = rng.uniform(size=n)
    z = sigma / xi * ((1 - uu) ** (-xi) - 1.0)
    marks = np.round((u + z) / grid) * grid
    return marks[marks > u]


def test_interval_cells_single_grid():
    marks = np.array([45.0, 50.0, 55.0])
    a, b, trunc = interval_cells(marks, threshold=40.0, grid=5.0)
    # excess of 45 kt over u=40 is 5; cell [2.5, 7.5), truncation 2.5
    assert np.allclose(a, [2.5, 7.5, 12.5])
    assert np.allclose(b, [7.5, 12.5, 17.5])
    assert np.allclose(trunc, 2.5)


def test_interval_cells_shifted_threshold_keeps_conditional_probability_valid():
    # threshold 42 is not on the 5-unit grid; the first admissible mark is 45
    # with lower rounding edge 42.5, so the conditioned excess is 0.5, not g/2.
    a, b, trunc = interval_cells([45.0, 50.0], threshold=42.0, grid=5.0)
    assert np.allclose(trunc, 0.5)
    # unit-scale exponential tail: P(cell | excess >= trunc) must lie in [0, 1]
    ratio = (np.exp(-a) - np.exp(-b)) / np.exp(-trunc)
    assert np.all((ratio >= 0.0) & (ratio <= 1.0))
    # aligned thresholds are unchanged: the truncation is still g/2
    _, _, trunc_aligned = interval_cells([45.0, 50.0], threshold=40.0, grid=5.0)
    assert np.allclose(trunc_aligned, 2.5)


def test_interval_cells_mixed_precision():
    # 50 is a multiple of 5 (coarse cell), 48 is only a multiple of 1 (fine cell)
    a, b, trunc = interval_cells(np.array([50.0, 48.0]), threshold=40.0, grid=(5.0, 1.0))
    assert np.isclose(trunc[0], 2.5) and np.isclose(trunc[1], 0.5)
    assert np.isclose(b[0] - a[0], 5.0) and np.isclose(b[1] - a[1], 1.0)


def test_grouped_recovers_shape_where_continuous_is_biased():
    rng = np.random.default_rng(0)
    xi_true, sigma_true, u = -0.20, 15.0, 40.0
    marks = _discretised_bounded(rng, xi_true, sigma_true, u, n=8000)

    grouped = nsevt.gpd_pot_grouped(marks, threshold=u, grid=5.0)
    continuous = nsevt.gpd_pot(marks, threshold=u, n_boot=0)

    # the grouped fit recovers the true shape; the continuous fit is more biased
    assert abs(grouped.xi - xi_true) < 0.03
    assert abs(grouped.xi - xi_true) < abs(continuous.xi - xi_true)


def test_grouped_endpoint_profile_ci_covers_and_orders():
    rng = np.random.default_rng(1)
    xi_true, sigma_true, u = -0.20, 15.0, 40.0
    m_star_true = u - sigma_true / xi_true  # 115.0
    marks = _discretised_bounded(rng, xi_true, sigma_true, u, n=8000)

    fit = nsevt.gpd_pot_grouped(marks, threshold=u, grid=5.0)
    lo, hi = fit.endpoint_ci95
    assert lo < fit.endpoint < hi
    assert lo <= m_star_true <= hi
    assert fit.bounded_supported


def test_gpd_pot_grouped_summary_and_fields():
    rng = np.random.default_rng(2)
    marks = _discretised_bounded(rng, n=4000)
    fit = nsevt.gpd_pot_grouped(marks, threshold=40.0, grid=5.0)
    assert fit.xi < 0 and np.isfinite(fit.endpoint)
    assert fit.xi_ci95[0] < fit.xi_ci95[1]
    assert isinstance(fit.summary(), str) and "grouped GPD fit" in fit.summary()


def test_grouped_requires_enough_exceedances():
    with pytest.raises(ValueError):
        nsevt.gpd_pot_grouped(np.array([41.0, 42.0]), threshold=40.0, grid=5.0)


# -- H01: invalid input and infeasible optima are never returned as a fit ------
def test_fit_gpd_grouped_rejects_nonfinite_grid():
    # a NaN grid width used to reach the optimizer, which "converged" on the
    # constant penalty region and returned xi=-0.4, loglik=-1e10 as a fit
    with pytest.raises(ValueError):
        nsevt.fit_gpd_grouped([45.0, 50.0, 55.0], 40.0, grid=np.nan)


def test_fit_gpd_grouped_rejects_nonfinite_threshold():
    with pytest.raises(ValueError):
        nsevt.fit_gpd_grouped([45.0, 50.0, 55.0], np.inf)


def test_interval_cells_rejects_nonfinite_inputs():
    with pytest.raises(ValueError):
        interval_cells([45.0, np.nan, 55.0], threshold=40.0, grid=5.0)
    with pytest.raises(ValueError):
        interval_cells([45.0, 50.0, 55.0], threshold=40.0, grid=np.inf)


def test_fit_gpd_grouped_validates_user_supplied_cells():
    with pytest.raises(ValueError):
        nsevt.fit_gpd_grouped(
            [45.0, 50.0, 55.0], 40.0,
            cells=(np.array([5.0, 6.0, 7.0]),      # upper edge below lower edge
                   np.array([1.0, 2.0, 3.0]),
                   np.array([0.5, 0.5, 0.5])),
        )
    with pytest.raises(ValueError):
        nsevt.fit_gpd_grouped(
            [45.0, 50.0, 55.0], 40.0,
            cells=(np.array([2.5, 7.5]), np.array([7.5, 12.5]),  # wrong length
                   np.array([2.5, 2.5])),
        )


def test_fit_gpd_grouped_reports_penalty_optimum_as_failure(monkeypatch):
    # if every likelihood evaluation stays in the penalty region the run is a
    # failure, not a fit whose objective is the penalty constant
    monkeypatch.setattr(grouped, "_grouped_nll",
                        lambda *args, **kwargs: grouped._PENALTY)
    with pytest.raises(RuntimeError):
        nsevt.fit_gpd_grouped([45.0, 50.0, 55.0], 40.0)


# -- H03: profile intervals consistent with the fit, open bounds exposed -------
def test_profile_ci_xi_grouped_contains_the_estimate():
    # a fixed [-0.95, 0.60] bracket used to exclude xi_hat for small samples
    out = nsevt.profile_ci_xi_grouped([45.0, 50.0, 55.0], 40.0)
    lo, hi = out["ci"]
    assert lo <= out["xi_hat"] <= hi
    assert "lo_at_bound" in out and "hi_at_bound" in out


def test_profile_ci_xi_grouped_rejects_bad_level():
    with pytest.raises(ValueError):
        nsevt.profile_ci_xi_grouped([45.0, 50.0, 55.0], 40.0, level=1.0)


def test_profile_endpoint_ci_reports_unbounded_tail_as_infinite():
    # heavy tail (xi > 0): the endpoint is not finite, and the upper limit must
    # be inf with a flag, not a number that is only the search boundary
    rng = np.random.default_rng(13)
    u = rng.uniform(size=300)
    z = 8.0 / 0.25 * ((1.0 - u) ** (-0.25) - 1.0)
    marks = np.round((40.0 + z) / 5.0) * 5.0
    marks = marks[marks > 40.0]
    ep = nsevt.profile_endpoint_ci(marks, 40.0)
    assert not np.isfinite(ep["endpoint"])
    assert not np.isfinite(ep["ci"][1])
    assert ep["upper_at_bound"] is True
    assert np.isfinite(ep["ci"][0])          # data still bound it from below


def test_gpd_pot_grouped_honours_requested_level():
    rng = np.random.default_rng(0)
    marks = _discretised_bounded(rng, n=6000)
    fit80 = nsevt.gpd_pot_grouped(marks, 40.0, grid=5.0, level=0.80)
    fit95 = nsevt.gpd_pot_grouped(marks, 40.0, grid=5.0, level=0.95)
    assert fit80.level == 0.80
    assert "80%" in fit80.summary()
    # historical aliases still resolve
    assert fit80.xi_ci95 == fit80.xi_ci
    assert fit80.endpoint_ci95 == fit80.endpoint_ci
    # a lower level gives a narrower shape interval
    assert (fit80.xi_ci[1] - fit80.xi_ci[0]) < (fit95.xi_ci[1] - fit95.xi_ci[0])
    assert fit80.xi_ci_at_bound == (False, False)
