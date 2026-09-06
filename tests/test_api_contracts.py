"""1.x contract guard: frozen dataclass fields, dict keys, and aliases.

A field, key, or alias listed here shipped in a 1.x release and must not
disappear or be renamed without a major version bump (see `docs/stability.md`).
New keys/fields are allowed; these are subset checks.
"""
import dataclasses

import numpy as np
import pytest

import nsevt
import nsevt.experimental
from nsevt import _types, design, mc
from nsevt import calibration as cal


# --------------------------------------------------------------------------
# fixtures: the smallest inputs that exercise each stable entry point
# --------------------------------------------------------------------------
def _bounded(n, seed=0, xi=-0.25, sigma=12.0, u=40.0):
    rng = np.random.default_rng(seed)
    uu = rng.uniform(size=n)
    z = sigma / xi * ((1 - uu) ** (-xi) - 1.0)
    return u + z


def _blocks(seed=0):
    rng = np.random.default_rng(seed)
    z, blk = [], []
    for y in range(1990, 2020):
        s = 10.0 * np.exp(0.15 * (y - 1990) / 10.0)
        u = rng.uniform(size=12)
        z.append(s / -0.25 * ((1 - u) ** 0.25 - 1))
        blk.append(np.full(12, y))
    return np.concatenate(z), np.concatenate(blk)


# --------------------------------------------------------------------------
# dataclasses
# --------------------------------------------------------------------------
DATACLASS_FIELDS = {
    "GPDFit": {"threshold", "n_exceedances", "xi", "sigma", "xi_ci95", "endpoint",
              "endpoint_ci95", "bootstrap_fraction_xi_negative", "n_boot_successful",
              "xi_ci95_truncated"},
    "GroupedGPDFit": {"threshold", "n_exceedances", "xi", "sigma", "xi_ci",
                     "endpoint", "endpoint_ci", "loglik"},
    "SourceResult": {"name", "n", "xi", "xi_ci95", "bounded_estimate",
                    "bounded_supported", "endpoint", "trend_per_decade",
                    "p_permutation", "trend_significant", "trend_direction",
                    "power_for_reference"},
    "ArenaResult": {"sources", "alpha", "shape_bounded_all", "trend_reproduces",
                   "reference_source", "trend_status", "verdict"},
    "TwoScaleResult": None,     # existence only
    "ConformalBand": {"threshold", "alpha", "q_standardized", "method",
                     "n_blocks", "block_length", "experimental"},
    "Checkpoint": None,
    "SequentialRun": None,
}


@pytest.mark.parametrize("name, required", DATACLASS_FIELDS.items())
def test_stable_dataclass_fields_are_preserved(name, required):
    cls = getattr(nsevt, name)
    assert dataclasses.is_dataclass(cls)
    if required is not None:
        have = {f.name for f in dataclasses.fields(cls)}
        assert required <= have, f"{name} lost {required - have}"


def test_dataclass_aliases():
    fit = nsevt.gpd_pot(_bounded(400), threshold=40.0, n_boot=0)
    assert fit.bounded == fit.bounded_supported
    g = nsevt.gpd_pot_grouped(_bounded(4000), threshold=40.0, grid=5.0)
    assert g.xi_ci95 == g.xi_ci and g.endpoint_ci95 == g.endpoint_ci


# --------------------------------------------------------------------------
# dict return contracts
# --------------------------------------------------------------------------
def test_dict_return_keys_are_preserved():
    raw = _bounded(4000)
    z = raw[raw > 40] - 40.0
    assert {"xi", "sigma", "nll"} <= set(nsevt.fit_gpd(z))
    marks = np.round(raw / 5.0) * 5.0
    marks = marks[marks > 40.0]              # design fits require strict exceedances
    assert {"xi", "sigma", "endpoint", "loglik", "n"} <= set(
        nsevt.fit_gpd_grouped(marks, 40.0))
    assert {"xi_hat", "ci", "lo_at_bound", "hi_at_bound", "level"} <= set(
        nsevt.profile_ci_xi_grouped(marks, 40.0))
    assert {"endpoint", "ci", "upper_at_bound", "xi", "sigma", "level",
            "method"} <= set(nsevt.profile_endpoint_ci(marks, 40.0))

    zb, blk = _blocks()
    tp = nsevt.trend_permutation(zb, blk, n_perm=99, seed=1)
    assert {"trend_per_decade", "xi_null", "log_sigma_null", "LR",
            "p_asymptotic", "p_permutation", "p_permutation_mcse"} <= set(tp)
    row = nsevt.trend_power(zb, blk, [0.1], n_rep=20, seed=1, n_perm_calibration=49)[0]
    assert {"trend_per_decade", "sigma_change_pct", "power", "power_mcse",
            "n_successful"} <= set(row)
    m = nsevt.min_detectable_effect(zb, blk, grid=[0.1, 0.3], direction="both",
                                    n_rep=8, n_perm_calibration=19, seed=1)
    assert {"mde_per_decade", "mde_negative", "mde_positive", "emd_per_decade",
            "emd_negative", "emd_positive", "emd_negative_ci95",
            "emd_positive_ci95", "direction", "target_power",
            "power_curve"} <= set(m)

    out = nsevt.block_bootstrap_trend_ci(zb, blk, n_boot=20, seed=1)
    assert {"ci95", "n_boot"} <= set(out)

    pv = mc.permutation_pvalue(1.0, [0.0, 2.0, 3.0])
    assert {"p", "n_exceed", "B", "floor", "at_floor", "mcse"} <= set(pv)
    run = mc.run_sequential("c", lambda k, b: np.zeros(k), kind="mean",
                            r0=10, r_min=10, r_max=10, block=5, epsilon=1e9,
                            tol_stability=1e9, min_stable_blocks=1)
    assert {"analysis", "kind", "R_star", "estimate", "mcse", "tolerance",
            "status", "trace", "batch_diagnostic"} <= set(run.summary())

    ctrl = dict(r0=20, r_min=20, r_max=40, block=20, min_stable_blocks=1, epsilon=0.2)
    rr = cal.rejection_rate(lambda s: 0.5, lambda r, n: r.uniform(size=n), n=1, **ctrl)
    assert {"rate", "mcse", "alpha", "n", "R", "status", "anticonservative",
            "stopping"} <= set(rr)
    cov = cal.coverage(lambda s: (-1.0, 1.0), lambda r, n: r.normal(size=n),
                       n=5, target=0.0, **ctrl)
    assert {"coverage", "mcse", "nominal", "target", "n", "R", "status",
            "miscalibration", "stopping"} <= set(cov)
    br = cal.bias_rmse(lambda s: 0.1, lambda r, n: r.normal(size=n), n=5,
                       truth=0.1, n_rep=10)
    assert {"bias", "bias_mcse", "sd", "rmse", "rmse_mcse", "mean_estimate",
            "truth", "n", "n_rep", "n_failed"} <= set(br)
    pt = cal.pseudo_true(lambda s: 0.1, lambda r, n: r.normal(size=n), R=10)
    assert {"pseudo_true", "R"} <= set(pt)

    X = np.ones((marks.size, 1))
    fd = design.fit_grouped_design(marks, 40.0, X, grid=5.0)
    assert {"xi", "coef", "sigma0", "loglik", "n", "p", "truncation"} <= set(fd)


# --------------------------------------------------------------------------
# name-level aliases and re-exports
# --------------------------------------------------------------------------
def test_name_aliases_and_experimental_reexports():
    # `transportability` is a permanent backward-compatible wrapper (not the
    # same object) that returns the same ArenaResult
    assert callable(nsevt.transportability)
    zb, blk = _blocks()
    src = [("x", 40.0 + zb, blk)]            # multisource takes raw values
    a = nsevt.multisource_robustness(src, threshold=40, n_perm=49,
                                     check_power=False, seed=1)
    b = nsevt.transportability(src, threshold=40, n_perm=49,
                               check_power=False, seed=1)
    assert type(a) is type(b) and a.trend_status == b.trend_status
    for name in ("split_conformal", "block_conformal", "ConformalBand",
                 "twoscale_trend", "wasserstein_decomposition", "TwoScaleResult"):
        assert getattr(nsevt.experimental, name) is getattr(nsevt, name)


def test_typed_schemas_match_the_runtime_dicts():
    # the TypedDicts in nsevt._types are a subset of what the functions return
    zb, blk = _blocks()
    checks = [
        (_types.PermutationPValue, mc.permutation_pvalue(1.0, [0.0, 2.0])),
        (_types.PowerRow, nsevt.trend_power(zb, blk, [0.1], n_rep=15, seed=1,
                                            n_perm_calibration=19)[0]),
    ]
    for schema, value in checks:
        assert set(schema.__annotations__) <= set(value)
