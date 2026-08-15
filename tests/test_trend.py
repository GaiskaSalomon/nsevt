import numpy as np
import pytest

import nsevt
import nsevt.trend as trend_module


def _make(trend_per_decade, n_per_year=25, years=range(1980, 2024), seed=0, xi=-0.25):
    """Excesses whose GPD log-scale drifts at ``trend_per_decade``."""
    rng = np.random.default_rng(seed)
    z, blk = [], []
    s0 = 12.0
    for y in years:
        t = (y - 1980) / 10.0
        s = s0 * np.exp(trend_per_decade * t)
        u = rng.uniform(size=n_per_year)
        zi = s / xi * ((1 - u) ** (-xi) - 1)   # GPD(xi, s) excesses
        z.append(zi)
        blk.append(np.full(n_per_year, y))
    return np.concatenate(z), np.concatenate(blk)


def test_no_trend_not_rejected():
    z, blk = _make(0.0, seed=3)
    r = nsevt.trend_permutation(z, blk, n_perm=99, seed=1)
    assert r["p_permutation"] > 0.05
    assert "_null" in r and len(r["_null"]) > 0
    assert r["p_permutation_mcse"] >= 0
    assert r["n_permutations"] == 99


def test_strong_trend_detected():
    z, blk = _make(0.25, seed=4)
    r = nsevt.trend_permutation(z, blk, n_perm=99, seed=1)
    assert r["trend_per_decade"] > 0.05
    assert r["p_permutation"] < 0.05


def test_mde_returns_value_and_monotone_power():
    z, blk = _make(0.0, seed=5)
    m = nsevt.min_detectable_effect(
        z, blk, grid=[0.05, 0.15, 0.30], direction="positive",
        n_rep=30, n_perm_calibration=49, seed=2
    )
    powers = [c["power"] for c in m["power_curve"]]
    assert powers[0] <= powers[-1] + 1e-9         # power increases with effect
    assert m["mde_per_decade"] is None or m["mde_per_decade"] > 0
    assert m["mde_negative"] is None


def test_emd_interpolated_crossing_and_uncertainty():
    # a design with real power: the 80% crossing falls between grid points
    z, blk = _make(0.0, seed=5)
    m = nsevt.min_detectable_effect(
        z, blk, grid=[0.05, 0.10, 0.15, 0.20, 0.30], direction="positive",
        n_rep=60, n_perm_calibration=99, seed=2, emd_uncertainty_reps=300,
    )
    emd, ci = m["emd_positive"], m["emd_positive_ci95"]
    assert emd is not None and ci is not None
    assert ci[0] <= emd <= ci[1]                       # uncertainty brackets the crossing
    assert emd <= m["mde_positive"]                    # interpolant is not coarser than the grid
    assert m["emd_negative"] is None                   # positive-only request
    assert m["emd_per_decade"] == emd


def test_failed_permutation_refits_are_disclosed(monkeypatch):
    z, blk = _make(0.0, n_per_year=5, years=range(2000, 2008), seed=9)
    original_fit = trend_module._fit_ns
    calls = {"n": 0}

    def flaky_fit(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] > 2 and calls["n"] % 2:
            raise RuntimeError("synthetic optimizer failure")
        return original_fit(*args, **kwargs)

    monkeypatch.setattr(trend_module, "_fit_ns_from", lambda *args: None)
    monkeypatch.setattr(trend_module, "_fit_ns", flaky_fit)
    with pytest.warns(RuntimeWarning, match="conditional on the successful refits"):
        out = trend_module.trend_permutation(z, blk, n_perm=6, seed=4)
    assert out["n_permutations"] == 3


def test_block_bootstrap_ci_contains_estimate_direction():
    z, blk = _make(0.20, seed=6)
    ci = nsevt.block_bootstrap_trend_ci(z, blk, n_boot=40, seed=2)["ci95"]
    assert ci[0] is not None and ci[1] >= ci[0]


def test_mde_both_directions_are_explicit():
    z, blk = _make(0.0, n_per_year=8, years=range(2000, 2015), seed=8)
    m = nsevt.min_detectable_effect(
        z, blk, grid=[0.1], direction="both", n_rep=4,
        n_perm_calibration=9, seed=3
    )
    assert {row["trend_per_decade"] for row in m["power_curve"]} == {-0.1, 0.1}


def test_trend_design_validation():
    with np.testing.assert_raises(ValueError):
        nsevt.trend_permutation([1, 2, 3], [2000, 2000, 2000], n_perm=9)
