import numpy as np
import nsevt


# --------------------------- conformal --------------------------------------
def test_block_conformal_covers_holdout():
    rng = np.random.default_rng(0)
    x = 40.0 + rng.exponential(9.0, size=1500)
    cal, test = x[:1000], x[1000:]
    band = nsevt.block_conformal(cal, threshold=40, alpha=0.10)
    sigma = nsevt.fit_gpd(cal[cal > 40] - 40)["sigma"]
    cov = band.coverage(test[test > 40], sigma)
    assert cov >= 0.85            # target 0.90, allow Monte-Carlo slack


def test_split_and_block_both_run():
    rng = np.random.default_rng(1)
    x = 40.0 + rng.exponential(9.0, size=800)
    b = nsevt.block_conformal(x, threshold=40, alpha=0.2)
    s = nsevt.split_conformal(x, threshold=40, alpha=0.2)
    assert b.q_standardized > 0 and s.q_standardized > 0
    assert b.n_blocks >= 2 and s.n_blocks == 1


# --------------------------- two-scale --------------------------------------
def test_twoscale_no_trend_not_rejected():
    rng = np.random.default_rng(2)
    samples = [rng.normal(0, 1, size=int(rng.integers(6, 14))) for _ in range(40)]
    r = nsevt.twoscale_trend(samples, n_boot=300, seed=1)
    assert r.p_value > 0.05


def test_twoscale_trend_detected():
    rng = np.random.default_rng(3)
    samples = [rng.normal(0.06 * t, 1, size=40) for t in range(40)]  # clear drift
    r = nsevt.twoscale_trend(samples, n_boot=300, seed=1)
    assert r.p_value < 0.05


def test_wasserstein_decomposition_pure_shift():
    rng = np.random.default_rng(4)
    a = rng.normal(0, 1, size=6000)
    b = a + 2.0                      # pure location shift, same shape/scale
    d = nsevt.wasserstein_decomposition(a, b)
    assert abs(d["location"] - 4.0) < 0.05      # (2.0)^2
    assert d["scale"] < 0.02 and d["shape"] < 0.02
    assert abs(d["w2_squared"] - (d["location"] + d["scale"] + d["shape"])) < 1e-9


# --------------------------- arena ------------------------------------------
def test_transportability_verdict_and_table():
    rng = np.random.default_rng(5)

    def src(trend, seed, xi=-0.25, s0=12.0):
        r = np.random.default_rng(seed)
        z, blk = [], []
        for y in range(1980, 2024):
            s = s0 * np.exp(trend * (y - 1980) / 10.0)
            u = r.uniform(size=25)
            z.append(40 + s / xi * ((1 - u) ** (-xi) - 1)); blk.append(np.full(25, y))
        return np.concatenate(z), np.concatenate(blk)

    xo, yo = src(0.25, 10)     # operational: strong trend
    xi_, yi = src(0.0, 11)     # independent: no trend
    arena = nsevt.transportability(
        [("operational", xo, yo), ("independent", xi_, yi)],
        threshold=40, n_perm=400, n_boot=200, check_power=True)
    assert arena.shape_bounded_all          # both bounded (xi<0)
    assert not arena.trend_reproduces       # trend does not survive
    assert "does NOT reproduce" in arena.verdict
    assert "VERDICT" in arena.table()
