"""Fast end-to-end smoke test for an installed nsevt distribution.

This script intentionally uses public entry points only.  Release workflows run
it outside the source tree against both the wheel and the source distribution;
the statistical test suite remains responsible for numerical calibration and
branch-level coverage.
"""
from __future__ import annotations

import importlib.metadata as metadata

import numpy as np

import nsevt
import nsevt.experimental as experimental


def tail_sample(seed: int = 7, n: int = 240) -> tuple[np.ndarray, np.ndarray]:
    """A small bounded-tail sample with repeated block labels."""
    rng = np.random.default_rng(seed)
    xi, sigma = -0.2, 12.0
    u = rng.uniform(size=n)
    excess = sigma / xi * ((1.0 - u) ** (-xi) - 1.0)
    block = np.repeat(np.arange(2000, 2012), n // 12)
    return excess, block


def main() -> None:
    installed = metadata.version("nsevt")
    assert nsevt.__version__ == installed

    excess, block = tail_sample()
    values = 40.0 + excess

    # Continuous and grouped POT-GPD paths.
    fit = nsevt.gpd_pot(values, threshold=40.0, n_boot=8, seed=1)
    assert fit.n_exceedances == values.size and np.isfinite(fit.xi)
    marks = np.round(values / 5.0) * 5.0
    marks = marks[marks > 40.0]
    grouped = nsevt.gpd_pot_grouped(marks, threshold=40.0, grid=5.0)
    assert grouped.n_exceedances == marks.size and np.isfinite(grouped.loglik)

    # Trend, grouped scale design, and return-level paths.
    trend = nsevt.trend_permutation(excess, block, n_perm=9, seed=2)
    assert trend["n_permutations"] == 9 and 0.0 < trend["p_permutation"] <= 1.0
    design = np.ones((marks.size, 1))
    regression = nsevt.fit_grouped_design(marks, 40.0, design, grid=5.0)
    assert regression["n"] == marks.size and regression["p"] == 1
    level = nsevt.return_level(-0.2, 12.0, 40.0, 0.4, np.array([10.0, 100.0]))
    assert np.all(np.isfinite(level)) and level[1] > level[0]

    # Sequential Monte Carlo and finite-sample calibration paths.
    rng = np.random.default_rng(3)
    run = nsevt.run_sequential(
        "smoke",
        lambda k, b: (rng.random(k) < 0.3).astype(float),
        r0=20,
        block=10,
        r_min=20,
        r_max=40,
        min_stable_blocks=1,
        epsilon=0.2,
        tol_stability=0.2,
    )
    assert run.R <= 40 and run.summary()["status"] in {
        nsevt.mc.CONVERGED,
        nsevt.mc.NOT_STABILISED,
    }
    proxy = nsevt.pseudo_true(
        np.mean, lambda r, n: r.normal(loc=2.0, size=n), R=100, seed=4
    )
    assert np.isfinite(proxy["pseudo_true"]) and proxy["R"] == 100

    # Multi-source status and experimental namespace imports remain usable.
    arena = nsevt.multisource_robustness(
        [("smoke", values, block)],
        threshold=40.0,
        n_perm=9,
        n_boot=5,
        check_power=False,
        seed=5,
    )
    assert arena.reference_source == "smoke" and "VERDICT" in arena.table()
    conformal = experimental.split_conformal(values, threshold=40.0, alpha=0.2)
    assert conformal.predict_upper(fit.sigma) > 40.0 and not conformal.experimental
    samples = [np.random.default_rng(i).normal(size=10) for i in range(12)]
    two_scale = experimental.twoscale_trend(samples, n_boot=9, seed=6)
    assert 0.0 < two_scale.p_value <= 1.0 and two_scale.experimental

    print(f"nsevt {installed} installed-distribution smoke test: PASS")


if __name__ == "__main__":
    main()
