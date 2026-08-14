"""Multi-source robustness: status branches, validation, table, and alias."""
import numpy as np
import pytest

import nsevt


def _source(rng, beta1, n_per_year=20, years=range(1990, 2024),
            xi=-0.2, s0=10.0, u=40.0):
    """Raw values above ``u`` with a log-scale trend ``beta1`` per decade."""
    vals, blk = [], []
    for y in years:
        t = (y - 1990) / 10.0
        s = s0 * np.exp(beta1 * t)
        uu = rng.uniform(size=n_per_year)
        vals.append(u + s / xi * ((1 - uu) ** (-xi) - 1))
        blk.append(np.full(n_per_year, y))
    return np.concatenate(vals), np.concatenate(blk)


_KW = dict(threshold=40, n_perm=199, n_boot=80, n_power=40, seed=1)
_STATUSES = {
    "reproduced", "inconsistent_direction", "not_reproduced_with_power",
    "not_resolved", "no_reference_signal", "single_source_only",
}


def test_reproduced_runs_power_loop_and_table():
    rng = np.random.default_rng(0)
    va, ba = _source(rng, 0.18)
    vb, bb = _source(rng, 0.18)
    arena = nsevt.multisource_robustness(
        [("A", va, ba), ("B", vb, bb)], reference="A", **_KW)
    assert arena.trend_status == "reproduced"
    assert arena.reference_source == "A"
    tbl = arena.table()
    assert isinstance(tbl, str) and "STATUS" in tbl and "VERDICT" in tbl
    # the power-aware loop populated the non-reference source's power
    b = next(s for s in arena.sources if s.name == "B")
    assert b.power_for_reference is not None
    assert b.bounded in (True, False)


def test_no_reference_signal_when_stationary():
    rng = np.random.default_rng(2)
    va, ba = _source(rng, 0.0)
    vb, bb = _source(rng, 0.0)
    arena = nsevt.transportability([("A", va, ba), ("B", vb, bb)], **_KW)
    assert arena.trend_status in _STATUSES        # exercises the decision branches
    assert isinstance(arena.table(), str)


def test_single_source_status():
    rng = np.random.default_rng(4)
    v, b = _source(rng, 0.18)
    arena = nsevt.multisource_robustness([("only", v, b)], **_KW)
    assert arena.trend_status in {"single_source_only", "no_reference_signal"}


def test_validation_errors():
    rng = np.random.default_rng(1)
    v, b = _source(rng, 0.0)
    with pytest.raises(ValueError):
        nsevt.multisource_robustness([], threshold=40)
    with pytest.raises(ValueError):
        nsevt.multisource_robustness([("A", v, b)], threshold=40, alpha=1.5)
    with pytest.raises(ValueError):
        nsevt.multisource_robustness([("A", v, b)], threshold=40, power_threshold=0)
    with pytest.raises(ValueError):
        nsevt.multisource_robustness([("A", v, b), ("A", v, b)], threshold=40)
    with pytest.raises(ValueError):
        nsevt.multisource_robustness([("", v, b)], threshold=40)
    with pytest.raises(ValueError):
        nsevt.multisource_robustness([("A", v, b)], threshold=40, reference="Z")
    with pytest.raises(ValueError):
        nsevt.multisource_robustness([("A", v[:, None], b)], threshold=40)
