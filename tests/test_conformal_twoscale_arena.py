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
    assert b.experimental and not s.experimental


def test_conformal_coverage_rejects_excess_scale_input():
    rng = np.random.default_rng(12)
    x = 40.0 + rng.exponential(5.0, size=100)
    band = nsevt.split_conformal(x, threshold=40, scale=5.0)
    with np.testing.assert_raises(ValueError):
        band.coverage(x[:10] - 40.0, 5.0)


# --------------------------- two-scale --------------------------------------
def test_twoscale_no_trend_not_rejected():
    rng = np.random.default_rng(2)
    samples = [rng.normal(0, 1, size=int(rng.integers(6, 14))) for _ in range(40)]
    r = nsevt.twoscale_trend(samples, n_boot=199, seed=1)
    assert r.p_value > 0.05
    assert r.experimental


def test_twoscale_trend_detected():
    rng = np.random.default_rng(3)
    samples = [rng.normal(0.06 * t, 1, size=40) for t in range(40)]  # clear drift
    r = nsevt.twoscale_trend(samples, n_boot=199, seed=1)
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
    def src(trend, seed, xi=-0.25, s0=12.0):
        r = np.random.default_rng(seed)
        z, blk = [], []
        for y in range(1980, 2024):
            s = s0 * np.exp(trend * (y - 1980) / 10.0)
            u = r.uniform(size=25)
            z.append(40 + s / xi * ((1 - u) ** (-xi) - 1))
            blk.append(np.full(25, y))
        return np.concatenate(z), np.concatenate(blk)

    xo, yo = src(0.25, 10)     # operational: strong trend
    xi_, yi = src(0.0, 11)     # independent: no trend
    arena = nsevt.multisource_robustness(
        [("operational", xo, yo), ("independent", xi_, yi)],
        threshold=40, n_perm=99, n_boot=30, check_power=False)
    assert arena.shape_bounded_all == all(s.bounded_supported for s in arena.sources)
    assert not arena.trend_reproduces       # trend does not survive
    assert arena.trend_status == "not_resolved"
    assert "not significant in every source" in arena.verdict
    assert "VERDICT" in arena.table()


def test_multisource_validates_names_and_reference():
    values = np.array([41.0, 42.0, 44.0, 45.0])
    blocks = np.array([2000, 2001, 2002, 2003])
    with np.testing.assert_raises(ValueError):
        nsevt.multisource_robustness(
            [("a", values, blocks), ("a", values, blocks)], threshold=40
        )
    with np.testing.assert_raises(ValueError):
        nsevt.multisource_robustness(
            [("a", values, blocks)], threshold=40, reference="missing"
        )


def test_opposite_significant_trends_are_not_reproduction():
    def source(effect, seed):
        rng = np.random.default_rng(seed)
        values, blocks = [], []
        xi = -0.25
        for year in range(1990, 2020):
            sigma = 8.0 * np.exp(effect * (year - 1990) / 10.0)
            u = rng.uniform(size=18)
            values.append(40 + sigma / xi * ((1 - u) ** (-xi) - 1))
            blocks.append(np.full(18, year))
        return np.concatenate(values), np.concatenate(blocks)

    positive, years_positive = source(0.45, 31)
    negative, years_negative = source(-0.45, 32)
    result = nsevt.multisource_robustness(
        [("positive", positive, years_positive),
         ("negative", negative, years_negative)],
        threshold=40,
        reference="positive",
        n_perm=49,
        n_boot=0,
        check_power=False,
    )
    assert not result.trend_reproduces
    assert result.trend_status == "inconsistent_direction"
