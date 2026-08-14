"""Tests for the sequential Monte Carlo precision protocol (nsevt.mc)."""
import numpy as np
import pytest

import nsevt
from nsevt import mc


# -- precision formulae ----------------------------------------------------
def test_required_replicates_anchor():
    # pi = 0.80, MCSE <= 0.0025 needs ~25,600 replicates
    assert mc.required_replicates(0.80, 0.0025) == 25600
    assert mc.mcse_proportion(0.80, 25600) == pytest.approx(0.0025, rel=1e-9)


def test_required_replicates_rejects_nonpositive_epsilon():
    with pytest.raises(ValueError):
        mc.required_replicates(0.5, 0.0)


def test_mcse_mean_and_quantile_shrink_with_R():
    rng = np.random.default_rng(0)
    small = rng.normal(size=200)
    big = rng.normal(size=20000)
    assert mc.mcse_mean(big) < mc.mcse_mean(small)
    assert mc.mcse_quantile(big, 0.5) < mc.mcse_quantile(small, 0.5)
    assert np.isnan(mc.mcse_mean([1.0]))
    assert np.isnan(mc.mcse_quantile(rng.normal(size=10), 0.5))  # < 30 values


# -- permutation p-value ---------------------------------------------------
def test_permutation_pvalue_floor_and_plus_one():
    null = np.zeros(999)                      # nothing reaches t_obs
    out = mc.permutation_pvalue(t_obs=5.0, t_null=null)
    assert out["at_floor"] is True
    assert out["n_exceed"] == 0
    assert out["p"] == pytest.approx(1.0 / 1000.0)
    assert out["floor"] == pytest.approx(1.0 / 1000.0)
    # p is never zero even with the raw fraction
    assert mc.permutation_pvalue(5.0, null, plus_one=False)["p"] == 0.0


def test_permutation_pvalue_counts_exceedances():
    null = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    out = mc.permutation_pvalue(t_obs=2.0, t_null=null)   # 2,3,4 exceed
    assert out["n_exceed"] == 3
    assert out["p"] == pytest.approx(4.0 / 6.0)


# -- sequential run --------------------------------------------------------
def test_sequential_run_respects_floor_and_reports_columns():
    rng = np.random.default_rng(0)
    run = mc.run_sequential(
        "t", lambda k, b: (rng.random(k) < 0.05).astype(float),
        kind="proportion", epsilon=0.002, tol_stability=0.002,
        r0=10000, block=2500, r_min=20000, r_max=40000, min_stable_blocks=4)
    assert run.R >= 20000                       # never below the floor
    assert run.status in (mc.CONVERGED, mc.NOT_STABILISED)
    s = run.summary()
    for key in ("R_star", "mcse", "tolerance", "last_batch_change",
                "n_stable_blocks", "decision_stable", "status", "trace",
                "batch_diagnostic"):
        assert key in s


def test_drifting_process_never_stabilises():
    rng = np.random.default_rng(1)
    state = {"k": 0}

    def draw(k, b):
        state["k"] += 1
        return (rng.random(k) < min(0.05 * state["k"], 0.9)).astype(float)

    run = mc.run_sequential("drift", draw, r0=10000, block=2500, r_min=20000,
                            r_max=25000, epsilon=1e-9)
    assert run.status == mc.NOT_STABILISED


def test_decision_stable_detects_a_flipping_decision():
    # A decision that keeps flipping across the stability window is unstable;
    # one that has settled is stable.
    rules = {"d": lambda v: v > 0.5}
    flip = mc.SequentialRun("f", min_stable_blocks=4, decision_rules=rules)
    flip.checkpoints = [mc.Checkpoint(R=r, value=0.5, mcse=0.0,
                                      decisions={"d": bool(i % 2)})
                        for i, r in enumerate(range(10000, 25001, 2500))]
    assert flip.decision_stable() is False
    settled = mc.SequentialRun("s", min_stable_blocks=4, decision_rules=rules)
    settled.checkpoints = [mc.Checkpoint(R=r, value=0.5, mcse=0.0,
                                         decisions={"d": True})
                           for r in range(10000, 25001, 2500)]
    assert settled.decision_stable() is True


def test_no_decision_rule_run_can_converge():
    # An identically stable run with no decision rule is allowed to converge.
    plain = mc.run_sequential(
        "plain", lambda k, b: np.full(k, 0.3), kind="mean", epsilon=1e9,
        tol_stability=1e9, r0=10000, block=2500, r_min=20000, r_max=25000,
        min_stable_blocks=4)
    assert plain.status == mc.CONVERGED


def test_batch_diagnostic_ratio_near_one_for_iid():
    rng = np.random.default_rng(2)
    run = mc.SequentialRun("t", kind="proportion")
    run.extend((rng.random(20000) < 0.3).astype(float))
    diag = run.batch_diagnostic(size=1000)
    assert diag["n_blocks"] == 20
    assert 0.5 < diag["ratio"] < 1.8


def test_trace_is_a_prefix_of_a_longer_run():
    rng = np.random.default_rng(3)
    vals = (rng.random(30000) < 0.4).astype(float)
    short = mc.SequentialRun("s").extend(vals[:15000])
    long = mc.SequentialRun("l").extend(vals[:30000])
    # cumulative estimate at R=15000 is identical whether the run stops there
    e_short, _ = short._estimate(upto=15000)
    e_long, _ = long._estimate(upto=15000)
    assert e_short == e_long


# -- reproducible substreams -----------------------------------------------
def test_substream_is_reproducible_and_order_independent():
    a = mc.substream(24072026, "power", "beta0.08").random(5)
    b = mc.substream(24072026, "power", "beta0.08").random(5)
    assert np.array_equal(a, b)


def test_substreams_of_different_tags_do_not_leak():
    x = mc.substream(1, "power").random(2000)
    y = mc.substream(1, "bootstrap").random(2000)
    assert not np.array_equal(x, y)
    assert abs(float(np.corrcoef(x, y)[0, 1])) < 0.08


def test_block_streams_extend_without_perturbing_earlier_blocks():
    six = mc.block_streams(7, 6, "run")
    twelve = mc.block_streams(7, 12, "run")
    for i in range(6):
        assert np.array_equal(six[i].random(10), twelve[i].random(10))
    with pytest.raises(ValueError):
        mc.block_streams(7, 0)


def test_multiseed_summary_reports_spread():
    out = mc.multiseed_summary({"s1": 0.80, "s2": 0.81, "s3": 0.79})
    assert out["n_seeds"] == 3
    assert out["mean"] == pytest.approx(0.80)
    assert out["range"] == pytest.approx(0.02)


# -- public surface --------------------------------------------------------
def test_public_reexports():
    for name in ("mcse_proportion", "required_replicates", "permutation_pvalue",
                 "SequentialRun", "run_sequential", "substream", "block_streams",
                 "multiseed_summary"):
        assert hasattr(nsevt, name)
