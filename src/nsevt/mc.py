"""Sequential Monte Carlo precision: MCSE, stopping rules and traces.

A Monte Carlo estimate (a power, a coverage, a permutation p-value, a bootstrap
interval limit) is only as trustworthy as its own simulation error.  This module
replaces a fixed replicate budget by an explicit precision rule that decides how
many replicates ``R*`` a run needs:

    R* = min{ R : MCSE_R <= epsilon
                  AND the estimate is stable over the last k blocks
                  AND no qualitative decision changed over those blocks }

Three conditions, not one.  Reaching the MCSE target is not sufficient, and
neither is numerical stability while a p-value is still crossing its threshold
or a power curve is still crossing its target.  A run that exhausts ``r_max``
without satisfying all three is reported as ``not_stabilised`` rather than
quietly truncated.

Everything a run accumulates -- the value at every checkpoint, the per-block
estimates, the stopping decision -- is written to the trace, so precision tables
and convergence figures are generated from data rather than typed by hand.

The module is deliberately estimator-agnostic and depends only on NumPy: the
caller supplies replicate outcomes (0/1 rejections, error magnitudes, interval
limits) and chooses ``kind`` to select the matching MCSE formula.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import numpy.typing as npt

#: Convergence statuses reported by :meth:`SequentialRun.finalise`.
CONVERGED = "converged"
NOT_STABILISED = "not_stabilised"
RUNNING = "running"


# --------------------------------------------------------------------------
# precision
# --------------------------------------------------------------------------
def mcse_proportion(p_hat: float, R: int) -> float:
    """Monte Carlo standard error of a proportion.

    ``MCSE = sqrt(p (1 - p) / R)`` away from the boundaries.  At an observed
    proportion of exactly zero or one, a Jeffreys half-count stabilisation is
    used so a finite run never reports zero simulation error.  Use for power,
    coverage and type-I error.
    """
    R = max(int(R), 1)
    p_hat = float(np.clip(p_hat, 0.0, 1.0))
    if p_hat in {0.0, 1.0}:
        p_hat = (p_hat * R + 0.5) / (R + 1.0)
    return math.sqrt(p_hat * (1.0 - p_hat) / R)


def mcse_mean(values: npt.ArrayLike) -> float:
    """Monte Carlo standard error of a mean (bias, RMSE, interval limits)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return float("nan")
    return float(v.std(ddof=1) / math.sqrt(v.size))


def mcse_quantile(values: npt.ArrayLike, q: float) -> float:
    """MCSE of a sample quantile from the spacing of the order statistics.

    ``se = sqrt(q (1 - q) / R) / f(x_q)`` with the density ``f`` estimated from
    the local spacing of the order statistics.  Use for the interval limits
    (for example 2.5/50/97.5) whose stabilisation you want to require before
    stopping a bootstrap.  Returns ``nan`` for fewer than 30 finite values.
    """
    if not np.isfinite(q) or not 0 < q < 1:
        raise ValueError("q must lie strictly between 0 and 1")
    v = np.sort(np.asarray(values, dtype=float))
    v = v[np.isfinite(v)]
    R = v.size
    if R < 30:
        return float("nan")
    k = int(np.clip(round(q * R), 1, R - 2))
    h = max(int(round(0.5 * R ** 0.5)), 2)
    lo, hi = max(k - h, 0), min(k + h, R - 1)
    dens = (hi - lo) / (R * max(v[hi] - v[lo], 1e-12))
    return float(math.sqrt(q * (1 - q) / R) / max(dens, 1e-12))


def required_replicates(p_hat: float, epsilon: float) -> int:
    """Replicates needed for ``MCSE <= epsilon`` at a proportion ``p_hat``.

    ``R = ceil(p (1 - p) / epsilon^2)`` away from the boundaries. At an
    observed zero or one, the returned budget inverts the same Jeffreys-
    stabilised MCSE as :func:`mcse_proportion`. At ``p = 0.80`` and
    ``epsilon = 0.0025`` this returns ``25600``, a useful sanity anchor when
    budgeting a power or coverage study near its most demanding point.
    """
    if not (epsilon > 0):
        raise ValueError("epsilon must be positive")
    p_hat = float(np.clip(p_hat, 0.0, 1.0))
    if p_hat in {0.0, 1.0}:
        upper = 1
        while mcse_proportion(p_hat, upper) > epsilon:
            upper *= 2
        lower = upper // 2
        while lower + 1 < upper:
            mid = (lower + upper) // 2
            if mcse_proportion(p_hat, mid) <= epsilon:
                upper = mid
            else:
                lower = mid
        return upper
    return int(math.ceil(p_hat * (1.0 - p_hat) / (epsilon ** 2)))


# --------------------------------------------------------------------------
# permutation / bootstrap p-values
# --------------------------------------------------------------------------
def permutation_pvalue(
    t_obs: float, t_null: npt.ArrayLike, plus_one: bool = True
) -> dict:
    """Monte Carlo p-value with the plus-one correction.

    ``p = (1 + #{T_b >= T_obs}) / (B + 1)``.  ``p = 0`` is never returned: when
    no permuted statistic reaches the observed one the result equals the
    resolution floor ``1 / (B + 1)`` and ``at_floor`` says so, so the p-value
    can be reported as "below the resolution of the experiment" rather than as
    an artificially precise zero.  Set ``plus_one=False`` for the raw
    exceedance fraction (diagnostics only).
    """
    t_null = np.asarray(t_null, dtype=float)
    t_null = t_null[np.isfinite(t_null)]
    B = int(t_null.size)
    exceed = int(np.sum(t_null >= t_obs))
    p = (1.0 + exceed) / (B + 1.0) if plus_one else exceed / max(B, 1)
    return {
        "p": float(p),
        "n_exceed": exceed,
        "B": B,
        "floor": float(1.0 / (B + 1.0)),
        "at_floor": bool(exceed == 0),
        "mcse": mcse_proportion(p, B),
    }


# --------------------------------------------------------------------------
# sequential run
# --------------------------------------------------------------------------
@dataclass
class Checkpoint:
    """One recorded point of a :class:`SequentialRun`."""

    R: int
    value: float
    mcse: float
    delta: Optional[float] = None
    decisions: dict = field(default_factory=dict)


@dataclass
class SequentialRun:
    """Accumulate a Monte Carlo run in blocks and apply the stopping rule.

    Parameters
    ----------
    name :
        Identifier used in the precision/stopping summary.
    kind :
        ``"proportion"``, ``"mean"`` or ``"quantile"``; selects the MCSE
        formula and hence the meaning of ``epsilon``.
    epsilon :
        MCSE tolerance, declared before the run.
    tol_stability :
        Maximum absolute change tolerated between consecutive checkpoints.
    r0, block, r_min, r_max, min_stable_blocks :
        The sequential rule: seed the run with ``r0`` replicates, then grow in
        blocks of ``block`` up to ``r_max``, never stopping below ``r_min`` and
        never before ``min_stable_blocks`` consecutive stable checkpoints.
    decision_rules :
        Mapping ``name -> callable(value) -> hashable`` evaluated at every
        checkpoint; if any output changes during the stability window the run
        must continue (guards against stopping while a decision is still
        flipping, even if the MCSE target is already met).

    Notes
    -----
    The object stores the raw per-replicate values, so extending a run never
    invalidates the trace already recorded: the checkpoint at ``R = 15000`` is
    exactly the first 15000 values of the run that later reaches ``30000``.
    """

    name: str
    kind: str = "proportion"
    epsilon: float = 0.0025
    tol_stability: float = 0.002
    r0: int = 10000
    block: int = 2500
    r_min: int = 20000
    r_max: int = 100000
    min_stable_blocks: int = 4
    trace_at: tuple = (1000, 2500, 5000, 10000, 15000, 20000, 30000, 50000)
    decision_rules: dict = field(default_factory=dict)
    quantile: float = 0.5
    seed_label: str = ""

    values: list = field(default_factory=list)
    checkpoints: list = field(default_factory=list)
    status: str = RUNNING

    def __post_init__(self) -> None:
        """Reject an incoherent precision protocol before drawing replicates."""
        if self.kind not in {"proportion", "mean", "quantile"}:
            raise ValueError("kind must be 'proportion', 'mean', or 'quantile'")
        if not np.isfinite(self.epsilon) or self.epsilon <= 0:
            raise ValueError("epsilon must be positive and finite")
        if not np.isfinite(self.tol_stability) or self.tol_stability < 0:
            raise ValueError("tol_stability must be non-negative and finite")
        for name in ("r0", "block", "r_min", "r_max", "min_stable_blocks"):
            value = getattr(self, name)
            if not isinstance(value, (int, np.integer)) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.r0 > self.r_max:
            raise ValueError("r0 must not exceed r_max")
        if self.r_min > self.r_max:
            raise ValueError("r_min must not exceed r_max")
        if not 0 < self.quantile < 1:
            raise ValueError("quantile must lie strictly between 0 and 1")
        if not isinstance(self.decision_rules, dict) or not all(
            callable(fn) for fn in self.decision_rules.values()
        ):
            raise ValueError("decision_rules must map names to callables")
        if any(not isinstance(r, (int, np.integer)) or r < 1 for r in self.trace_at):
            raise ValueError("trace_at must contain positive integers")

    # -- accumulation ------------------------------------------------------
    def extend(self, new_values: npt.ArrayLike) -> SequentialRun:
        """Add a block of replicate outcomes and record a checkpoint."""
        self.values.extend(np.asarray(new_values, dtype=float).ravel().tolist())
        self._checkpoint()
        return self

    @property
    def R(self) -> int:
        return len(self.values)

    def _estimate(self, upto: Optional[int] = None) -> tuple[float, float]:
        v = np.asarray(self.values[:upto] if upto else self.values, dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            return float("nan"), float("nan")
        if self.kind == "proportion":
            est = float(v.mean())
            return est, mcse_proportion(est, v.size)
        if self.kind == "quantile":
            return float(np.quantile(v, self.quantile)), mcse_quantile(v, self.quantile)
        return float(v.mean()), mcse_mean(v)

    def _checkpoint(self) -> None:
        est, se = self._estimate()
        prev = self.checkpoints[-1].value if self.checkpoints else None
        decisions = {k: fn(est) for k, fn in self.decision_rules.items()}
        self.checkpoints.append(Checkpoint(
            R=self.R, value=est, mcse=se,
            delta=None if prev is None else abs(est - prev),
            decisions=decisions))

    # -- stopping rule -----------------------------------------------------
    def stable_blocks(self) -> int:
        """Consecutive trailing checkpoints within ``tol_stability``."""
        k = 0
        for cp in reversed(self.checkpoints):
            if cp.delta is None or cp.delta >= self.tol_stability:
                break
            k += 1
        return k

    def decision_stable(self) -> bool:
        """No qualitative decision changed during the stability window."""
        if not self.decision_rules:
            return True
        window = self.checkpoints[-(self.min_stable_blocks + 1):]
        if len(window) < 2:
            return False
        first = window[0].decisions
        return all(cp.decisions == first for cp in window[1:])

    def should_stop(self) -> bool:
        if self.R < max(self.r0, self.r_min):
            return False
        _, se = self._estimate()
        return (se <= self.epsilon
                and self.stable_blocks() >= self.min_stable_blocks
                and self.decision_stable())

    def finalise(self) -> str:
        """Set and return the convergence status for the summary."""
        self.status = CONVERGED if self.should_stop() else NOT_STABILISED
        return self.status

    # -- diagnostics -------------------------------------------------------
    def batch_diagnostic(self, size: int = 1000) -> dict:
        """Independent-block diagnostic of the reported MCSE.

        Splits the run into blocks of ``size`` and compares the observed
        between-block scatter with the MCSE the theory predicts for one block.
        A ``ratio`` far from 1 signals a dependence or bookkeeping problem
        rather than a merely noisy run.
        """
        if not isinstance(size, (int, np.integer)) or size < 1:
            raise ValueError("size must be a positive integer")
        v = np.asarray(self.values, dtype=float)
        v = v[np.isfinite(v)]
        nb = v.size // size
        if nb < 2:
            return {"n_blocks": int(nb), "ratio": None}
        shaped = v[:nb * size].reshape(nb, size)
        if self.kind == "quantile":
            blocks = np.quantile(shaped, self.quantile, axis=1)
            full_mcse = mcse_quantile(v, self.quantile)
            theoretical = float(full_mcse * math.sqrt(v.size / size))
        else:
            blocks = shaped.mean(axis=1)
        observed = float(blocks.std(ddof=1))
        if self.kind == "proportion":
            theoretical = mcse_proportion(float(v.mean()), size)
        elif self.kind == "mean":
            theoretical = float(v.std(ddof=1) / math.sqrt(size))
        return {
            "n_blocks": int(nb),
            "block_size": int(size),
            "observed_sd_between_blocks": observed,
            "theoretical_mcse_per_block": float(theoretical),
            "ratio": float(observed / theoretical) if theoretical > 0 else None,
            "block_means": blocks.tolist(),
        }

    def trace(self) -> list:
        """Cumulative estimate at the pre-specified checkpoints."""
        out = []
        for R in self.trace_at:
            if R > self.R:
                continue
            est, se = self._estimate(upto=R)
            out.append({"R": int(R), "value": est, "mcse": se})
        est, se = self._estimate()
        if not out or out[-1]["R"] != self.R:
            out.append({"R": int(self.R), "value": est, "mcse": se})
        return out

    # -- reporting ---------------------------------------------------------
    def summary(self) -> dict:
        """One row of a Monte Carlo precision-and-stopping table."""
        est, se = self._estimate()
        last = self.checkpoints[-1] if self.checkpoints else None
        return {
            "analysis": self.name,
            "kind": self.kind,
            "R_star": int(self.R),
            "estimate": est,
            "mcse": se,
            "tolerance": self.epsilon,
            "last_batch_change": None if last is None else last.delta,
            "stability_tolerance": self.tol_stability,
            "n_stable_blocks": int(self.stable_blocks()),
            "min_stable_blocks": int(self.min_stable_blocks),
            "decision_stable": bool(self.decision_stable()),
            "seed": self.seed_label,
            "status": self.status,
            "trace": self.trace(),
            "batch_diagnostic": self.batch_diagnostic(),
        }


def run_sequential(
    name: str,
    draw: Callable[[int, int], npt.ArrayLike],
    *,
    kind: str = "proportion",
    epsilon: float = 0.0025,
    tol_stability: float = 0.002,
    r0: int = 10000,
    block: int = 2500,
    r_min: int = 20000,
    r_max: int = 100000,
    min_stable_blocks: int = 4,
    decision_rules: Optional[dict] = None,
    quantile: float = 0.5,
    seed_label: str = "",
    trace_at: Optional[tuple] = None,
    progress: Optional[Callable[[SequentialRun], None]] = None,
) -> SequentialRun:
    """Drive a sequential run to ``R*``.

    ``draw(k, block_index)`` must return ``k`` fresh replicate outcomes.  Give
    each block index its own random substream (see :func:`block_streams`) so
    the trace is reproducible and extending a run leaves earlier blocks
    untouched.
    """
    run = SequentialRun(
        name=name, kind=kind, epsilon=epsilon, tol_stability=tol_stability,
        r0=r0, block=block, r_min=r_min, r_max=r_max,
        min_stable_blocks=min_stable_blocks,
        decision_rules=decision_rules or {}, quantile=quantile,
        seed_label=seed_label,
        trace_at=tuple(trace_at) if trace_at else SequentialRun.trace_at)

    bi = 0
    def draw_exact(k: int, block_index: int) -> np.ndarray:
        values = np.asarray(draw(k, block_index), dtype=float).ravel()
        if values.size != k:
            raise ValueError(f"draw({k}, {block_index}) returned {values.size} values")
        return values

    run.extend(draw_exact(r0, bi))
    if progress:
        progress(run)
    while not run.should_stop() and run.R < r_max:
        bi += 1
        run.extend(draw_exact(min(block, r_max - run.R), bi))
        if progress:
            progress(run)
    run.finalise()
    return run


# --------------------------------------------------------------------------
# reproducible substreams
# --------------------------------------------------------------------------
def _entropy(seed: int, *tags: object) -> list:
    parts = [int(seed)]
    for t in tags:
        digest = hashlib.sha256(str(t).encode("utf-8")).digest()[:8]
        parts.append(int.from_bytes(digest, "big") >> 1)
    return parts


def substream(seed: int, *tags: object) -> np.random.Generator:
    """A reproducible ``Generator`` for a named ``(seed, *tags)`` combination.

    The generator is a pure function of the master ``seed`` and the string form
    of the tags, so requesting streams in a different order never changes any of
    them and independent analyses (different tags) do not share numbers.  Use it
    to give each Monte Carlo analysis its own stream from one master seed.
    """
    return np.random.default_rng(np.random.SeedSequence(_entropy(seed, *tags)))


def block_streams(seed: int, n_blocks: int, *tags: object) -> list:
    """One reproducible ``Generator`` per sequential block.

    Giving block ``b`` its own substream means the trace at an earlier ``R`` is
    exactly the first blocks of a longer run: extending a run never perturbs the
    replicates already drawn.  Pass ``streams[block_index]`` inside the ``draw``
    callback of :func:`run_sequential`.
    """
    if not isinstance(n_blocks, (int, np.integer)) or n_blocks < 1:
        raise ValueError("n_blocks must be a positive integer")
    children = np.random.SeedSequence(_entropy(seed, *tags)).spawn(int(n_blocks))
    return [np.random.default_rng(c) for c in children]


# --------------------------------------------------------------------------
# multi-seed audit
# --------------------------------------------------------------------------
def multiseed_summary(values_by_seed: dict) -> dict:
    """Compare the same estimate across independent substreams.

    ``values_by_seed`` maps a substream label to its estimate.  The spread
    across seeds should be compatible with the reported MCSE; a larger spread
    means the run stopped too early whatever its nominal MCSE said.
    """
    labels = list(values_by_seed)
    vals = np.array([values_by_seed[k] for k in labels], dtype=float)
    return {
        "n_seeds": int(vals.size),
        "labels": labels,
        "values": vals.tolist(),
        "mean": float(vals.mean()),
        "sd_across_seeds": float(vals.std(ddof=1)) if vals.size > 1 else 0.0,
        "range": float(vals.max() - vals.min()) if vals.size else 0.0,
    }
