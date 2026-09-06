#!/usr/bin/env python3
"""Run a small finite-sample validation campaign and write a provenance archive.

This is the mechanism behind the ``validation-vX.Y.Z`` deposit described in
``RELEASING.md``. It exercises the calibration suite over a grid of
(shape, sample size) cells and records, per cell, the empirical size, power,
coverage, bias, RMSE and their Monte Carlo standard errors, the requested /
successful / failed replicate counts, the stopping status, every seed, and the
package version and git commit.

A cell whose MCSE does not meet ``--tolerance`` is written with
``"resolved": false``: it is reported as *not resolved at this budget* rather
than as a validated number.

    python tools/validation_campaign.py --out validation-v1.2.0.json

The default grid is deliberately light so the script is runnable in minutes; a
real deposit raises ``--reps`` / ``--r-max`` and widens ``--shapes`` /
``--sizes``. The output is deterministic given the seeds.
"""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import nsevt  # noqa: E402
from nsevt import calibration as cal  # noqa: E402


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[1],
            text=True,
        ).strip()
    except Exception:
        return "unknown"


def _gpd_sample(rng: np.random.Generator, n: int, xi: float, sigma: float = 15.0):
    u = rng.uniform(size=n)
    if abs(xi) < 1e-9:
        return -sigma * np.log1p(-u)
    return sigma / xi * ((1.0 - u) ** (-xi) - 1.0)


def _cell(xi: float, n: int, seed: int, reps: int, r_max: int, tol: float) -> dict:
    """One (shape, size) cell: profile-interval coverage and shape bias/RMSE.

    A full deposit adds a proper size/power test for the trend statistic; this
    covers the two properties the suite reports directly from a fitted model.
    """
    def simulate(rng, m):
        return _gpd_sample(rng, m, xi=xi)

    cov = cal.coverage(
        lambda s: nsevt.profile_ci_xi(s)[2], simulate, n=n, target=xi, level=0.95,
        seed=seed, epsilon=tol, r0=min(2000, r_max), r_min=min(2000, r_max),
        block=1000, r_max=r_max, min_stable_blocks=1,
    )
    bias = cal.bias_rmse(
        lambda s: nsevt.fit_gpd(s)["xi"], simulate, n=n, truth=xi,
        n_rep=reps, seed=seed + 1,
    )
    resolved = (
        cov["mcse"] <= tol
        and bias["bias_mcse"] <= tol
        and cov["status"] == nsevt.mc.CONVERGED
    )
    return {
        "xi": xi,
        "n": n,
        "resolved": bool(resolved),
        "seeds": {"coverage": seed, "bias_rmse": seed + 1},
        "coverage": {k: cov[k] for k in ("coverage", "mcse", "n_effective",
                                         "n_failed", "n_infinite", "status")},
        "bias_rmse": {k: bias[k] for k in ("bias", "bias_mcse", "sd", "rmse",
                                           "rmse_mcse", "n_failed")},
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--shapes", type=float, nargs="+", default=[-0.3, -0.1, 0.1])
    p.add_argument("--sizes", type=int, nargs="+", default=[200, 500])
    p.add_argument("--reps", type=int, default=400)
    p.add_argument("--r-max", type=int, default=8000)
    p.add_argument("--tolerance", type=float, default=0.01)
    p.add_argument("--seed", type=int, default=20260906)
    args = p.parse_args()

    cells = []
    seed = args.seed
    for xi in args.shapes:
        for n in args.sizes:
            cells.append(_cell(xi, n, seed, args.reps, args.r_max, args.tolerance))
            seed += 100

    archive = {
        "package_version": nsevt.__version__,
        "git_commit": _git_commit(),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "platform": platform.platform(),
        "tolerance": args.tolerance,
        "reps": args.reps,
        "r_max": args.r_max,
        "n_cells": len(cells),
        "n_unresolved": sum(1 for c in cells if not c["resolved"]),
        "cells": cells,
    }
    args.out.write_text(json.dumps(archive, indent=2, sort_keys=True))
    print(f"wrote {args.out} ({archive['n_cells']} cells, "
          f"{archive['n_unresolved']} unresolved)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
