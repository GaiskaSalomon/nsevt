"""Pre-specified multi-source robustness analysis.

The same threshold and inferential procedure are applied to each source.  The
analysis can show whether conclusions reproduce under a change of source; it
cannot, without additional identification assumptions, prove that a discrepancy
is an instrumental artifact or that agreement is a physical property.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .gpd import gpd_pot
from .trend import trend_permutation, trend_power


@dataclass
class SourceResult:
    """Results for one data source under the common specification."""

    name: str
    n: int
    xi: float
    xi_ci95: tuple
    bounded_estimate: bool
    bounded_supported: bool
    endpoint: float
    trend_per_decade: float
    p_permutation: float
    trend_significant: bool
    trend_direction: str
    power_for_reference: Optional[float] = None

    @property
    def bounded(self) -> bool:
        """Compatibility alias for inferentially supported negative shape."""
        return self.bounded_supported


@dataclass
class ArenaResult:
    """Cross-source result and calibrated interpretation status."""

    sources: list
    alpha: float
    shape_bounded_all: bool
    trend_reproduces: bool
    reference_source: str
    trend_status: str
    verdict: str

    def table(self) -> str:
        head = (
            f"{'source':<20}{'n':>6}{'xi':>8}{'xi CI95':>16}"
            f"{'xi<0*':>8}{'trend/dec':>11}{'p_perm':>9}{'power':>9}"
        )
        rows = [head, "-" * len(head)]
        for source in self.sources:
            ci = f"[{source.xi_ci95[0]:.2f},{source.xi_ci95[1]:.2f}]"
            power = (
                f"{source.power_for_reference:.2f}"
                if source.power_for_reference is not None
                else "--"
            )
            rows.append(
                f"{source.name:<20}{source.n:>6}{source.xi:>8.3f}{ci:>16}"
                f"{str(source.bounded_supported):>8}{source.trend_per_decade:>11.3f}"
                f"{source.p_permutation:>9.4f}{power:>9}"
            )
        rows.extend(["", "* Entire 95% profile interval below zero."])
        rows.append(f"STATUS: {self.trend_status}")
        rows.append(f"VERDICT: {self.verdict}")
        return "\n".join(rows)


def _direction(effect: float, tolerance: float = 1e-12) -> str:
    if effect > tolerance:
        return "positive"
    if effect < -tolerance:
        return "negative"
    return "zero"


def multisource_robustness(
    sources: Sequence[tuple],
    threshold: float,
    reference: Optional[str] = None,
    alpha: float = 0.05,
    n_perm: int = 3000,
    n_boot: int = 2000,
    check_power: bool = True,
    seed: int = 20260722,
    n_power: int = 200,
    power_threshold: float = 0.80,
) -> ArenaResult:
    """Apply one frozen GPD/trend specification across multiple sources.

    A trend is called reproduced only when the reference source is significant
    and every source has a significant estimate in the same direction.  When a
    source does not reject, simulation power for the signed reference effect
    distinguishes ``not_reproduced_with_power`` from ``not_resolved``.
    """
    if not sources:
        raise ValueError("provide at least one (name, values, block) source")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between 0 and 1")
    if not 0 < power_threshold < 1:
        raise ValueError("power_threshold must lie strictly between 0 and 1")
    names = [item[0] for item in sources]
    if any(not isinstance(name, str) or not name.strip() for name in names):
        raise ValueError("every source name must be a non-empty string")
    if len(set(names)) != len(names):
        raise ValueError("source names must be unique")
    reference = reference or names[0]
    if reference not in names:
        raise ValueError(f"reference source {reference!r} is not present")

    results, internals = [], {}
    for name, values, block in sources:
        values = np.asarray(values, dtype=float)
        block = np.asarray(block)
        if values.ndim != 1 or block.ndim != 1 or values.size != block.size:
            raise ValueError(f"source {name!r}: values and block must be matching 1-D arrays")
        if np.any(~np.isfinite(values)):
            raise ValueError(f"source {name!r}: values must be finite")
        fit = gpd_pot(values, threshold, n_boot=n_boot, seed=seed)
        mask = values > threshold
        z, selected_blocks = values[mask] - threshold, block[mask]
        trend = trend_permutation(z, selected_blocks, n_perm=n_perm, seed=seed)
        internals[name] = (z, selected_blocks, trend)
        results.append(
            SourceResult(
                name=name,
                n=fit.n_exceedances,
                xi=fit.xi,
                xi_ci95=fit.xi_ci95,
                bounded_estimate=fit.bounded_estimate,
                bounded_supported=fit.bounded_supported,
                endpoint=fit.endpoint,
                trend_per_decade=trend["trend_per_decade"],
                p_permutation=trend["p_permutation"],
                trend_significant=trend["p_permutation"] < alpha,
                trend_direction=_direction(trend["trend_per_decade"]),
            )
        )

    reference_result = next(row for row in results if row.name == reference)
    reference_trend = reference_result.trend_per_decade
    if check_power and reference_result.trend_significant:
        for row in results:
            if row.name == reference:
                row.power_for_reference = None
                continue
            z, selected_blocks, trend = internals[row.name]
            critical = float(np.quantile(trend["_null"], 1 - alpha, method="higher"))
            curve = trend_power(
                z,
                selected_blocks,
                [reference_trend],
                xi_true=trend["xi_null"],
                log_sigma_true=trend["log_sigma_null"],
                crit=critical,
                alpha=alpha,
                n_rep=n_power,
                seed=seed,
            )
            row.power_for_reference = curve[0]["power"]

    shape_all = all(row.bounded_supported for row in results)
    if not reference_result.trend_significant:
        trend_status = "no_reference_signal"
        reproduces = False
        trend_text = "the reference source does not reject the no-trend null"
    else:
        nonreference = [row for row in results if row.name != reference]
        opposite = any(
            row.trend_significant and row.trend_direction != reference_result.trend_direction
            for row in nonreference
        )
        reproduces = bool(nonreference) and all(
            row.trend_significant and row.trend_direction == reference_result.trend_direction
            for row in nonreference
        )
        if not nonreference:
            trend_status = "single_source_only"
            trend_text = "cross-source reproduction cannot be assessed with one source"
        elif opposite:
            trend_status = "inconsistent_direction"
            trend_text = "at least one source has a significant trend in the opposite direction"
        elif reproduces:
            trend_status = "reproduced"
            trend_text = "the signed trend is significant in the same direction in every source"
        else:
            unresolved = [row for row in nonreference if not row.trend_significant]
            adequate = bool(unresolved) and all(
                row.power_for_reference is not None
                and row.power_for_reference >= power_threshold
                for row in unresolved
            )
            if check_power and adequate:
                trend_status = "not_reproduced_with_power"
                trend_text = (
                    "the trend is absent in at least one source with simulated power "
                    f">= {power_threshold:.0%} for the signed reference effect"
                )
            else:
                trend_status = "not_resolved"
                trend_text = (
                    "the trend is not significant in every source, but available power "
                    "does not resolve whether this is disagreement or limited information"
                )

    if shape_all:
        shape_text = "all 95% profile intervals support a negative GPD shape"
    else:
        shape_text = "a negative GPD shape is not supported by every 95% profile interval"
    verdict = (
        f"Under the common threshold and model, {shape_text}; {trend_text}. "
        "This is a robustness assessment, not attribution of discrepancies to a cause."
    )
    return ArenaResult(
        sources=results,
        alpha=alpha,
        shape_bounded_all=shape_all,
        trend_reproduces=reproduces,
        reference_source=reference,
        trend_status=trend_status,
        verdict=verdict,
    )


def transportability(*args, **kwargs) -> ArenaResult:
    """Backward-compatible name for :func:`multisource_robustness`."""
    return multisource_robustness(*args, **kwargs)
