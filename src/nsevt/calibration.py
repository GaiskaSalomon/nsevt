"""Finite-sample calibration of an estimator or test by Monte Carlo.

A method's asymptotic guarantees (a nominal 0.05 level, a 95% interval, a
consistent point estimate) need not hold at the sample size and under the data
generating process you actually face.  This module measures three finite-sample
properties directly, from user-supplied callables, so the result is a property
of *your* estimator on *your* DGP rather than of an idealised model:

* :func:`rejection_rate` -- the empirical rejection probability of a test.  Run
  it under a null DGP to get the type-I error, or under an alternative to get
  power.
* :func:`coverage` -- the empirical probability that an interval estimator
  covers a target value.  Under a misspecified DGP the estimator may converge
  to a *pseudo-true* value rather than the generating parameter;
  :func:`pseudo_true` supplies a large-sample simulation proxy for that target.
* :func:`bias_rmse` -- the bias, variance and RMSE of a point estimator.

Each proportion (rejection rate, coverage) is driven by the sequential protocol
of :mod:`nsevt.mc`, so it comes with a Monte Carlo standard error and a stopping
decision rather than a fixed, arbitrary replicate budget.  The module is
estimator-agnostic and depends only on NumPy and :mod:`nsevt.mc`.

The user supplies two callables:

``simulate(rng, n) -> sample``
    Draw one data set of size ``n`` from the DGP using the given NumPy
    ``Generator``.  The returned ``sample`` is whatever the estimator/test
    consumes (typically a 1-D array of exceedances).

``estimator(sample) -> value`` or ``(lo, hi)``; ``test(sample) -> p``
    Map a sample to a point estimate, an interval, or a p-value.
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from . import mc

Sample = object
Simulate = Callable[[np.random.Generator, int], Sample]


def _block_draw(seed: int, tag: str, per_replicate: Callable[[np.random.Generator], float]):
    """Build a ``run_sequential`` draw that gives each block its own substream."""

    def draw(k: int, block_index: int) -> np.ndarray:
        rng = mc.substream(seed, tag, "block", block_index)
        return np.array([per_replicate(rng) for _ in range(int(k))], dtype=float)

    return draw


def rejection_rate(
    test: Callable[[Sample], float],
    simulate: Simulate,
    n: int,
    *,
    alpha: float = 0.05,
    seed: int = 20260814,
    epsilon: float = 0.0025,
    tag: str = "rejection_rate",
    anticonservative_margin: float = 0.005,
    **seq_kwargs,
) -> dict:
    """Empirical rejection probability of ``test`` at level ``alpha``.

    ``test(sample)`` must return a p-value; the test rejects when ``p < alpha``.
    Run under a null DGP this is the empirical type-I error (compare it with
    ``alpha``); run under an alternative it is power.  A decision rule flags the
    run as ``anticonservative`` while the estimated rate sits more than
    ``anticonservative_margin`` above ``alpha``, so the run cannot stop until
    that verdict has settled.
    """
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")

    def per_replicate(rng: np.random.Generator) -> float:
        return float(test(simulate(rng, n)) < alpha)

    run = mc.run_sequential(
        tag, _block_draw(seed, tag, per_replicate), kind="proportion",
        epsilon=epsilon,
        decision_rules={"anticonservative":
                        lambda v: v > alpha + anticonservative_margin},
        seed_label=f"{tag}/{seed}", **seq_kwargs)
    est, se = run._estimate()
    return {
        "rate": est,
        "mcse": se,
        "alpha": alpha,
        "n": int(n),
        "R": run.R,
        "status": run.status,
        "anticonservative": bool(est > alpha + 2 * se),
        "stopping": run.summary(),
    }


def coverage(
    estimator: Callable[[Sample], tuple],
    simulate: Simulate,
    n: int,
    target: float | Callable[[], float],
    *,
    level: float = 0.95,
    seed: int = 20260814,
    epsilon: float = 0.0025,
    tag: str = "coverage",
    **seq_kwargs,
) -> dict:
    """Empirical probability that ``estimator`` covers ``target``.

    ``estimator(sample)`` must return an interval ``(lo, hi)``.  ``target`` is
    the value coverage is judged against: the generating parameter when the DGP
    is well specified, or the :func:`pseudo_true` value when it is not (pass a
    float, or a zero-argument callable evaluated once).  A non-finite limit
    counts as non-coverage.  The gap between ``level`` and the reported coverage
    is the interval's finite-sample miscalibration.
    """
    tgt = float(target() if callable(target) else target)

    def per_replicate(rng: np.random.Generator) -> float:
        lo, hi = estimator(simulate(rng, n))
        lo, hi = float(lo), float(hi)
        return float(np.isfinite(lo) and np.isfinite(hi) and lo <= tgt <= hi)

    run = mc.run_sequential(
        tag, _block_draw(seed, tag, per_replicate), kind="proportion",
        epsilon=epsilon,
        decision_rules={"covers_nominal": lambda v: v >= level},
        seed_label=f"{tag}/{seed}", **seq_kwargs)
    est, se = run._estimate()
    return {
        "coverage": est,
        "mcse": se,
        "nominal": level,
        "target": tgt,
        "n": int(n),
        "R": run.R,
        "status": run.status,
        "miscalibration": est - level,
        "stopping": run.summary(),
    }


def bias_rmse(
    estimator: Callable[[Sample], float],
    simulate: Simulate,
    n: int,
    truth: float,
    *,
    n_rep: int = 5000,
    seed: int = 20260814,
    tag: str = "bias_rmse",
) -> dict:
    """Bias, standard deviation and RMSE of a point ``estimator``.

    ``estimator(sample)`` must return a scalar.  ``truth`` is the value the
    estimate is compared against (usually the generating parameter, or the
    :func:`pseudo_true` value under misspecification).  Non-finite estimates are
    dropped and counted in ``n_failed``.  This is a fixed-budget summary; its
    Monte Carlo errors (``bias_mcse``, ``rmse_mcse``) say whether ``n_rep`` was
    enough.
    """
    if not isinstance(n_rep, (int, np.integer)) or n_rep < 2:
        raise ValueError("n_rep must be an integer >= 2")
    rng = mc.substream(seed, tag)
    est = np.array([estimator(simulate(rng, n)) for _ in range(int(n_rep))],
                   dtype=float)
    finite = est[np.isfinite(est)]
    R = finite.size
    if R < 2:
        raise RuntimeError("fewer than two finite estimates")
    err = finite - float(truth)
    bias = float(err.mean())
    var = float(finite.var(ddof=1))
    rmse = float(np.sqrt(np.mean(err ** 2)))
    # MCSE of RMSE by the delta method: sd(err^2) / (2 RMSE sqrt(R)).
    sq = err ** 2
    rmse_mcse = (float(sq.std(ddof=1) / (2.0 * rmse * np.sqrt(R)))
                 if rmse > 0 else float("nan"))
    return {
        "bias": bias,
        "bias_mcse": mc.mcse_mean(finite),
        "sd": float(np.sqrt(var)),
        "rmse": rmse,
        "rmse_mcse": rmse_mcse,
        "mean_estimate": float(finite.mean()),
        "truth": float(truth),
        "n": int(n),
        "n_rep": int(n_rep),
        "n_failed": int(n_rep - R),
    }


def pseudo_true(
    estimator: Callable[[Sample], float],
    simulate: Simulate,
    *,
    R: int = 20000,
    seed: int = 20260814,
    tag: str = "pseudo_true",
) -> dict:
    """Large-sample simulation proxy for an estimator's pseudo-true target.

    When the DGP is misspecified for the estimator's model (discretisation, a
    mixture, dependence), the estimator may converge to a target different
    from the nominal parameter.  This function approximates that target by
    fitting one simulated sample of size ``R``.  It does not prove convergence,
    existence, or uniqueness of a pseudo-true value; callers should examine
    sensitivity to increasing ``R`` and independent seeds before using the
    result as a coverage or bias target.
    """
    if not isinstance(R, (int, np.integer)) or R < 3:
        raise ValueError("R must be an integer >= 3")
    rng = mc.substream(seed, tag)
    sample = simulate(rng, int(R))
    value = float(estimator(sample))
    if not np.isfinite(value):
        raise RuntimeError("estimator returned a non-finite pseudo-true proxy")
    return {"pseudo_true": value, "R": int(R)}
