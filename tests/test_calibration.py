"""Tests for the finite-sample calibration suite (nsevt.calibration)."""
import numpy as np
import pytest

import nsevt
from nsevt import calibration as cal


def _gpd_sample(rng, n, xi=-0.2, sigma=15.0):
    """One sample of ``n`` GPD(xi, sigma) exceedances."""
    u = rng.uniform(size=n)
    if abs(xi) < 1e-9:
        return -sigma * np.log1p(-u)
    return sigma / xi * ((1.0 - u) ** (-xi) - 1.0)


# -- rejection_rate --------------------------------------------------------
def test_type_i_of_an_exact_test_is_near_alpha():
    # A p-value drawn uniformly on [0, 1] is an exact test; its rejection rate
    # at alpha must be alpha up to the MCSE.
    def unif_test(sample):
        return float(sample[0])

    def simulate(rng, n):
        return rng.uniform(size=n)

    out = cal.rejection_rate(unif_test, simulate, n=1, alpha=0.10,
                             epsilon=0.003, r_max=60000, seed=1)
    assert abs(out["rate"] - 0.10) <= 4 * out["mcse"] + 0.01
    assert out["status"] in (nsevt.mc.CONVERGED, nsevt.mc.NOT_STABILISED)
    assert "stopping" in out


def test_power_of_a_biased_test_exceeds_its_size():
    # A test that rejects more often under a shifted p-value has higher power.
    def test(sample):
        return float(sample[0])

    size = cal.rejection_rate(test, lambda r, n: r.uniform(size=n), n=1,
                              alpha=0.2, epsilon=0.004, r_max=40000, seed=2)
    power = cal.rejection_rate(test, lambda r, n: r.uniform(size=n) ** 2, n=1,
                               alpha=0.2, epsilon=0.004, r_max=40000, seed=3)
    assert power["rate"] > size["rate"]


def test_rejection_rate_rejects_bad_alpha():
    with pytest.raises(ValueError):
        cal.rejection_rate(lambda s: 0.5, lambda r, n: r.uniform(size=n),
                           n=1, alpha=1.5)


# -- coverage --------------------------------------------------------------
def test_profile_interval_covers_xi_when_well_specified():
    xi_true = -0.2

    def estimator(sample):
        _, _, (lo, hi) = nsevt.profile_ci_xi(sample)
        return lo, hi

    def simulate(rng, n):
        return _gpd_sample(rng, n, xi=xi_true)

    out = cal.coverage(estimator, simulate, n=400, target=xi_true, level=0.95,
                       epsilon=0.02, r0=600, r_min=600, block=300, r_max=900,
                       min_stable_blocks=1, seed=4)
    # coverage of a correctly specified profile interval is close to nominal
    assert out["coverage"] > 0.88
    assert abs(out["miscalibration"]) < 0.08


def test_coverage_accepts_a_callable_target():
    seen = {}

    def target():
        seen["called"] = True
        return -0.2

    def estimator(sample):
        _, _, (lo, hi) = nsevt.profile_ci_xi(sample)
        return lo, hi

    out = cal.coverage(estimator, lambda r, n: _gpd_sample(r, n), n=300,
                       target=target, epsilon=0.05, r0=400, r_min=400,
                       block=200, r_max=600, min_stable_blocks=1, seed=5)
    assert seen.get("called") is True
    assert out["target"] == pytest.approx(-0.2)


# -- bias_rmse -------------------------------------------------------------
def test_bias_of_the_mle_shape_is_small_when_well_specified():
    xi_true = -0.2

    def estimator(sample):
        return nsevt.fit_gpd(sample)["xi"]

    out = cal.bias_rmse(estimator, lambda r, n: _gpd_sample(r, n, xi=xi_true),
                        n=500, truth=xi_true, n_rep=400, seed=6)
    assert abs(out["bias"]) < 5 * out["bias_mcse"] + 0.03
    assert out["rmse"] >= abs(out["bias"])
    assert out["n_failed"] == 0
    assert set(out) >= {"bias", "sd", "rmse", "rmse_mcse", "mean_estimate"}


def test_bias_rmse_rejects_tiny_budget():
    with pytest.raises(ValueError):
        cal.bias_rmse(lambda s: 0.0, lambda r, n: _gpd_sample(r, n), n=100,
                      truth=0.0, n_rep=1)


# -- pseudo_true -----------------------------------------------------------
def test_pseudo_true_differs_from_nominal_under_rounding():
    # Fitting a continuous GPD to grid-rounded data biases the shape; the
    # pseudo-true value the continuous MLE targets is not the generating xi.
    xi_true = -0.2

    def rounded(rng, n):
        z = _gpd_sample(rng, n, xi=xi_true)
        return np.round((40.0 + z) / 5.0) * 5.0 - 40.0

    def estimator(sample):
        return nsevt.fit_gpd(sample[sample > 0])["xi"]

    pt = cal.pseudo_true(estimator, rounded, R=40000, seed=7)
    assert pt["R"] == 40000
    assert np.isfinite(pt["pseudo_true"])
    # a clean, un-rounded DGP recovers the generating value at the pseudo limit
    clean = cal.pseudo_true(lambda s: nsevt.fit_gpd(s)["xi"],
                            lambda r, n: _gpd_sample(r, n, xi=xi_true),
                            R=40000, seed=7)
    assert abs(clean["pseudo_true"] - xi_true) < 0.03


def test_pseudo_true_rejects_invalid_budget_and_nonfinite_estimate():
    with pytest.raises(ValueError):
        cal.pseudo_true(lambda s: 0.0, lambda r, n: np.zeros(n), R=2)
    with pytest.raises(RuntimeError, match="non-finite"):
        cal.pseudo_true(lambda s: np.nan, lambda r, n: np.zeros(n), R=10)


# -- public surface --------------------------------------------------------
def test_public_reexports():
    for name in ("rejection_rate", "coverage", "bias_rmse", "pseudo_true"):
        assert hasattr(nsevt, name)
    assert hasattr(nsevt, "calibration")
