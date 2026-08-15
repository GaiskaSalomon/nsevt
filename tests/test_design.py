"""Tests for grouped GPD regression and return levels (nsevt.design)."""
import numpy as np
import pytest

import nsevt
from nsevt import design


def _grouped_trend(rng, n=4000, xi=-0.2, log_sigma0=2.7, slope=0.15,
                   u=40.0, grid=5.0):
    """Bounded GPD excesses with a log-linear trend, rounded to ``grid``."""
    t = np.linspace(0.0, 4.0, n)
    sigma = np.exp(log_sigma0 + slope * t)
    uu = rng.uniform(size=n)
    z = sigma / xi * ((1 - uu) ** (-xi) - 1.0)
    marks = np.round((u + z) / grid) * grid
    keep = marks > u
    return marks[keep], np.column_stack([np.ones(keep.sum()), t[keep]])


# -- fit_grouped_design ----------------------------------------------------
def test_intercept_only_matches_gpd_pot_grouped():
    rng = np.random.default_rng(0)
    marks, _ = _grouped_trend(rng, slope=0.0)
    X = np.ones(marks.size)[:, None]
    f = design.fit_grouped_design(marks, 40.0, X, grid=5.0)
    g = nsevt.gpd_pot_grouped(marks, 40.0, grid=5.0)
    assert f["xi"] == pytest.approx(g.xi, abs=1e-6)
    assert f["sigma0"] == pytest.approx(g.sigma, abs=1e-4)
    assert f["loglik"] == pytest.approx(g.loglik, abs=1e-4)
    assert f["p"] == 1 and f["n"] == g.n_exceedances


def test_fit_recovers_a_positive_trend():
    rng = np.random.default_rng(1)
    marks, X = _grouped_trend(rng, slope=0.15)
    f = design.fit_grouped_design(marks, 40.0, X, grid=5.0)
    assert f["p"] == 2
    assert f["coef"][1] > 0.0                    # positive log-scale trend
    assert f["xi"] < 0.0                         # bounded tail preserved


@pytest.mark.parametrize(
    "values, design, message",
    [
        ([45.0, 50.0, 55.0], np.ones((2, 1)), "aligned"),
        ([45.0, 50.0, 55.0], np.ones((3, 2)), "full column rank"),
        ([40.0, 50.0, 55.0], np.ones((3, 1)), "strictly above"),
        ([45.0, np.nan, 55.0], np.ones((3, 1)), "finite 1-D"),
    ],
)
def test_fit_rejects_invalid_or_unidentified_design(values, design, message):
    with pytest.raises(ValueError, match=message):
        nsevt.fit_grouped_design(values, 40.0, design)


def test_fit_rejects_misaligned_cells():
    cells = (np.zeros(2), np.ones(2), np.zeros(2))
    with pytest.raises(ValueError, match="aligned"):
        design.fit_grouped_design([45.0, 50.0, 55.0], 40.0,
                                  np.ones((3, 1)), cells=cells)


# -- profile_ci_coef -------------------------------------------------------
def test_coef_ci_brackets_the_estimate():
    rng = np.random.default_rng(2)
    marks, X = _grouped_trend(rng, slope=0.15)
    out = design.profile_ci_coef(marks, 40.0, X, coef=1, grid=5.0)
    lo, hi = out["ci"]
    assert lo < out["estimate"] < hi
    assert out["coef"] == 1
    assert set(out) >= {"estimate", "ci", "at_bound", "level", "loglik"}


def test_flat_design_ci_contains_zero():
    rng = np.random.default_rng(3)
    marks, X = _grouped_trend(rng, slope=0.0)
    out = design.profile_ci_coef(marks, 40.0, X, coef=1, grid=5.0)
    lo, hi = out["ci"]
    assert lo < 0.0 < hi                         # no trend -> interval covers 0


def test_coef_ci_reuses_a_supplied_fit():
    rng = np.random.default_rng(4)
    marks, X = _grouped_trend(rng, slope=0.15)
    f = design.fit_grouped_design(marks, 40.0, X, grid=5.0)
    out = design.profile_ci_coef(marks, 40.0, X, coef=1, grid=5.0, fit=f)
    assert out["estimate"] == pytest.approx(f["coef"][1])


# -- return_level ----------------------------------------------------------
def test_return_level_matches_formula_and_is_monotone():
    xi, sigma, u, rate = -0.2, 15.0, 40.0, 0.4
    m = np.array([10.0, 100.0, 1000.0])
    rl = design.return_level(xi, sigma, u, rate, m)
    expected = u + (sigma / xi) * ((m * rate) ** xi - 1.0)
    assert np.allclose(rl, expected)
    assert rl[0] < rl[1] < rl[2]                 # rises with the return period
    # xi -> 0 limit is the exponential/Gumbel level
    rl0 = design.return_level(0.0, sigma, u, rate, m)
    assert np.allclose(rl0, u + sigma * np.log(m * rate))


def test_bounded_return_level_stays_below_the_endpoint():
    xi, sigma, u, rate = -0.25, 15.0, 40.0, 0.4
    endpoint = u - sigma / xi
    rl = design.return_level(xi, sigma, u, rate, 1e6)
    assert rl < endpoint                         # can never exceed the ceiling


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"sigma": 0.0}, "sigma"),
        ({"rate": 0.0}, "rate"),
        ({"rate": 1.1}, "rate"),
        ({"m": 0.0}, "return periods"),
    ],
)
def test_return_level_rejects_invalid_domain(kwargs, message):
    args = {"xi": -0.2, "sigma": 15.0, "threshold": 40.0,
            "rate": 0.4, "m": 100.0}
    args.update(kwargs)
    with pytest.raises(ValueError, match=message):
        design.return_level(**args)


# -- profile_ci_return_level ----------------------------------------------
def test_return_level_profile_interval_brackets_the_level():
    rng = np.random.default_rng(5)
    marks, _ = _grouped_trend(rng, slope=0.0)
    out = design.profile_ci_return_level(marks, 40.0, rate=0.4, m=50.0, grid=5.0)
    assert out is not None
    lo, hi = out["ci"]
    assert lo <= out["return_level"] <= hi
    assert out["xi_at_max"] < 0.0
    assert set(out) >= {"m", "return_level", "ci", "upper_at_bound", "loglik"}


# -- public surface --------------------------------------------------------
def test_public_reexports():
    for name in ("fit_grouped_design", "profile_ci_coef", "return_level",
                 "profile_ci_return_level"):
        assert hasattr(nsevt, name)
