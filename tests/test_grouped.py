"""Tests for the interval-censored (grouped) GPD module."""
import numpy as np
import pytest

import nsevt
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
