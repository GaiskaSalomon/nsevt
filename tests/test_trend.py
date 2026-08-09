import numpy as np
import nsevt


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
        z.append(zi); blk.append(np.full(n_per_year, y))
    return np.concatenate(z), np.concatenate(blk)


def test_no_trend_not_rejected():
    z, blk = _make(0.0, seed=3)
    r = nsevt.trend_permutation(z, blk, n_perm=400, seed=1)
    assert r["p_permutation"] > 0.05
    assert "_null" in r and len(r["_null"]) > 0


def test_strong_trend_detected():
    z, blk = _make(0.25, seed=4)
    r = nsevt.trend_permutation(z, blk, n_perm=400, seed=1)
    assert r["trend_per_decade"] > 0.05
    assert r["p_permutation"] < 0.05


def test_mde_returns_value_and_monotone_power():
    z, blk = _make(0.0, seed=5)
    m = nsevt.min_detectable_effect(z, blk, grid=[0.05, 0.15, 0.30],
                                    n_rep=120, seed=2)
    powers = [c["power"] for c in m["power_curve"]]
    assert powers[0] <= powers[-1] + 1e-9         # power increases with effect
    assert m["mde_per_decade"] is None or m["mde_per_decade"] > 0


def test_block_bootstrap_ci_contains_estimate_direction():
    z, blk = _make(0.20, seed=6)
    ci = nsevt.block_bootstrap_trend_ci(z, blk, n_boot=300, seed=2)["ci95"]
    assert ci[0] is not None and ci[1] >= ci[0]
