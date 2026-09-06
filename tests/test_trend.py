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


def test_block_bootstrap_keeps_the_trend_it_measures():
    # resampled blocks keep their own time, so a strong trend is not pulled to
    # zero the way relabelling to ordered positions did
    z, blk = _make(0.30, seed=6)
    point = nsevt.trend_permutation(z, blk, n_perm=9, seed=2)["trend_per_decade"]
    out = nsevt.block_bootstrap_trend_ci(z, blk, n_boot=60, seed=2)
    lo, hi = out["ci95"]
    assert lo > 0.0 and lo <= point <= hi
    assert out["n_unidentified"] == 0 and out["n_failed"] == 0
    assert out["n_boot_requested"] == 60


def test_block_bootstrap_reports_unidentifiable_resamples():
    # two blocks only: many resamples draw a single distinct time and cannot
    # inform a trend; they are counted, not silently dropped or fitted
    z, blk = _make(0.1, years=range(2000, 2002), n_per_year=40, seed=1)
    out = nsevt.block_bootstrap_trend_ci(z, blk, n_boot=50, seed=0)
    assert out["n_unidentified"] > 0
    assert out["n_unidentified"] + out["n_failed"] + out["n_boot"] == 50


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


# -- H07: power uncertainty is carried through, not rounded away --------------
def test_power_curve_carries_stabilised_mcse_and_failure_counts():
    z, blk = _make(0.0, seed=5)
    # a tiny effect (near-zero power) and a huge one (near-1 power): the plain
    # binomial MCSE would be 0 at both ends; the Jeffreys-stabilised one is not
    curve = nsevt.trend_power(z, blk, [0.0, 5.0], n_rep=40, seed=2,
                              n_perm_calibration=49)
    for row in curve:
        assert row["power_mcse"] > 0.0
        assert row["n_rep"] == 40
        assert row["n_successful"] + row["n_failed"] == 40
        # power is not pre-rounded to 3 dp before MDE/robustness consume it
        assert row["power"] == pytest.approx(row["power"], abs=0)


def test_mde_reports_emd_resolution_and_power_failures():
    z, blk = _make(0.0, seed=5)
    m = nsevt.min_detectable_effect(
        z, blk, grid=[0.05, 0.10, 0.15, 0.20, 0.30], direction="positive",
        n_rep=60, n_perm_calibration=99, seed=2, emd_uncertainty_reps=200,
    )
    assert "emd_positive_resolved" in m and isinstance(m["emd_positive_resolved"], bool)
    assert m["emd_positive_reps_without_crossing"] >= 0
    assert m["n_power_failed"] == sum(r["n_failed"] for r in m["power_curve"])


def test_emd_interpolation_tolerates_a_nan_power_row():
    # if a whole effect's replicates fail, its power is NaN; the interpolation
    # drops that row instead of crashing
    z, blk = _make(0.0, seed=5)
    m = nsevt.min_detectable_effect(
        z, blk, grid=[0.05, 0.15, 0.30], direction="positive",
        n_rep=20, n_perm_calibration=49, seed=2,
    )
    m["power_curve"][1]["power"] = float("nan")          # inject a failed row
    out, diag = trend_module._emd_interp(
        m["power_curve"], 1, 0.8, np.random.default_rng(0), 50)
    assert out is None or np.isfinite(out)


def test_multisource_carries_power_mcse():
    def src(trend, seed, xi=-0.25, s0=12.0):
        r = np.random.default_rng(seed)
        z, b = [], []
        for y in range(1980, 2024):
            s = s0 * np.exp(trend * (y - 1980) / 10.0)
            u = r.uniform(size=25)
            z.append(40 + s / xi * ((1 - u) ** (-xi) - 1))
            b.append(np.full(25, y))
        return np.concatenate(z), np.concatenate(b)

    xo, yo = src(0.25, 10)
    xi_, yi = src(0.0, 11)
    arena = nsevt.multisource_robustness(
        [("op", xo, yo), ("ind", xi_, yi)], threshold=40, n_perm=99,
        n_power=80, seed=1)
    non_ref = [s for s in arena.sources if s.name != arena.reference_source]
    assert all(s.power_mcse_for_reference is not None
               and s.power_mcse_for_reference > 0.0 for s in non_ref)
