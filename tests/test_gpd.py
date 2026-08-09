import numpy as np
import pytest

import nsevt


def _bounded(n, seed=0, hi=25.0):
    rng = np.random.default_rng(seed)
    return 40.0 + hi * rng.beta(1.2, 2.5, size=n)


def test_bounded_tail_detected():
    x = _bounded(600)
    fit = nsevt.gpd_pot(x, threshold=40, n_boot=50)
    assert fit.xi < 0
    assert fit.bounded_estimate
    assert fit.bounded_supported
    assert fit.bounded
    assert np.isfinite(fit.endpoint)
    # endpoint below the true bound (40+25=65) but above the sample max minus tiny
    assert 40 < fit.endpoint < 75
    assert fit.bootstrap_fraction_xi_negative > 0.9


def test_bounded_requires_interval_support_not_only_point_estimate():
    fit = nsevt.GPDFit(
        threshold=40,
        n_exceedances=100,
        xi=-0.08,
        sigma=8.0,
        xi_ci95=(-0.25, 0.12),
        endpoint=140.0,
        endpoint_ci95=[80.0, 300.0],
        bootstrap_fraction_xi_negative=0.7,
    )
    assert fit.bounded_estimate
    assert not fit.bounded_supported
    assert not fit.bounded


def test_exponential_tail_near_zero_shape():
    rng = np.random.default_rng(1)
    x = 40.0 + rng.exponential(8.0, size=800)   # xi ~ 0
    fit = nsevt.gpd_pot(x, threshold=40, n_boot=30)
    assert abs(fit.xi) < 0.2
    lo, hi = fit.xi_ci95
    assert lo <= 0.05  # interval consistent with xi <= 0-ish


def test_return_level_monotone_and_capped():
    x = _bounded(600)
    fit = nsevt.gpd_pot(x, threshold=40, n_boot=100)
    rate = (x > 40).mean()
    rls = [fit.return_level(rp, rate=rate) for rp in (10, 100, 1000, 10000)]
    assert all(b >= a for a, b in zip(rls, rls[1:]))       # non-decreasing
    assert rls[-1] <= fit.endpoint + 1e-6                  # capped by endpoint


def test_too_few_exceedances_raises():
    with pytest.raises(ValueError):
        nsevt.gpd_pot(np.array([1.0, 2.0, 41.0]), threshold=40)


@pytest.mark.parametrize("bad", [[1.0, np.nan, 2.0], [1.0, -0.1, 2.0], [1.0, 1.0, 1.0]])
def test_invalid_excesses_raise(bad):
    with pytest.raises(ValueError):
        nsevt.fit_gpd(bad)


def test_small_sample_fit_stays_in_likelihood_existence_region():
    fit = nsevt.fit_gpd([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert fit["xi"] > -0.999


def test_return_level_validates_inputs():
    x = _bounded(600)
    fit = nsevt.gpd_pot(x, threshold=40, n_boot=0)
    with pytest.raises(ValueError):
        fit.return_level(1)
    with pytest.raises(ValueError):
        fit.return_level(10, rate=0)


def test_profile_ci_brackets_point_estimate():
    x = _bounded(500)
    xih, sh, (lo, hi) = nsevt.profile_ci_xi(x[x > 40] - 40)
    assert lo <= xih <= hi
    assert lo < 0  # bounded
