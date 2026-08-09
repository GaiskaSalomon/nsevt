"""Multi-source transportability ("evidence arena").

A trend estimated on one product can be an artifact of that product's
heterogeneity rather than a property of the phenomenon.  This module re-runs the
*same* frozen specification (same threshold, same GPD and permutation
procedures) on several intensity/data sources and reports two verdicts:

* is the tail **bounded** (``xi < 0``) in every source?  A shape that reproduces
  across sources is a property of the phenomenon.
* does the **trend survive** the change of source?  A trend present in one
  product but absent in an independent or homogenized one is the signature of an
  instrumental artifact, not a physical change.

This is the design that, in the motivating application, separated a robust
finite ceiling from an apparent trend that did not survive homogenization.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .gpd import gpd_pot
from .trend import trend_permutation, min_detectable_effect


@dataclass
class SourceResult:
    name: str
    n: int
    xi: float
    xi_ci95: tuple
    bounded: bool
    endpoint: float
    trend_per_decade: float
    p_permutation: float
    trend_significant: bool
    power_for_reference: Optional[float] = None


@dataclass
class ArenaResult:
    sources: list
    alpha: float
    shape_bounded_all: bool
    trend_reproduces: bool
    reference_source: str
    verdict: str

    def table(self) -> str:
        head = (f"{'source':<20}{'n':>6}{'xi':>8}{'xi CI95':>16}"
                f"{'endpoint':>10}{'trend/dec':>11}{'p_perm':>9}")
        rows = [head, "-" * len(head)]
        for s in self.sources:
            ci = f"[{s.xi_ci95[0]:.2f},{s.xi_ci95[1]:.2f}]"
            end = f"{s.endpoint:.1f}" if np.isfinite(s.endpoint) else "inf"
            rows.append(f"{s.name:<20}{s.n:>6}{s.xi:>8.3f}{ci:>16}"
                        f"{end:>10}{s.trend_per_decade:>11.3f}{s.p_permutation:>9.4f}")
        rows.append("")
        rows.append(f"VERDICT: {self.verdict}")
        return "\n".join(rows)


def transportability(sources: Sequence[tuple], threshold: float,
                     reference: Optional[str] = None, alpha: float = 0.05,
                     n_perm: int = 3000, n_boot: int = 2000,
                     check_power: bool = True, seed: int = 20260722) -> ArenaResult:
    """Run the frozen GPD + trend specification across intensity sources.

    Parameters
    ----------
    sources : sequence of ``(name, values, block)``
        Each source is a raw-value array with matching block labels (e.g. year).
    threshold : float
        Common tail threshold ``u`` applied identically to every source.
    reference : str, optional
        Name of the source whose estimated trend defines the effect size for the
        power check on the other sources (default: the first source).
    check_power : bool
        If true, evaluate each non-reference source's Monte-Carlo power to detect
        the reference trend, so a null can be judged genuine vs. underpowered.
    """
    if not sources:
        raise ValueError("provide at least one (name, values, block) source")
    reference = reference or sources[0][0]
    results, ref_trend = [], None
    fits = {}
    for name, values, block in sources:
        fit = gpd_pot(values, threshold, n_boot=n_boot, seed=seed)
        z = np.asarray(values, float)
        z = z[z > threshold] - threshold
        tr = trend_permutation(z, np.asarray(block)[np.asarray(values, float) > threshold],
                               n_perm=n_perm, seed=seed)
        fits[name] = (z, np.asarray(block)[np.asarray(values, float) > threshold])
        r = SourceResult(
            name=name, n=fit.n_exceedances, xi=fit.xi, xi_ci95=fit.xi_ci95,
            bounded=fit.bounded, endpoint=fit.endpoint,
            trend_per_decade=tr["trend_per_decade"], p_permutation=tr["p_permutation"],
            trend_significant=tr["p_permutation"] < alpha)
        results.append(r)
        if name == reference:
            ref_trend = tr["trend_per_decade"]

    if check_power and ref_trend is not None:
        for r in results:
            if r.name == reference:
                r.power_for_reference = 1.0
                continue
            z, blk = fits[r.name]
            pc = min_detectable_effect(z, blk, grid=[abs(ref_trend)], n_rep=200,
                                       seed=seed)["power_curve"]
            r.power_for_reference = pc[0]["power"] if pc else None

    shape_all = all(r.bounded for r in results)
    sig = [r for r in results if r.trend_significant]
    # trend "reproduces" only if significant in every source
    trend_repro = all(r.trend_significant for r in results)
    ref = next(r for r in results if r.name == reference)
    if shape_all and not trend_repro and ref.trend_significant:
        verdict = ("bounded tail is ROBUST across all sources; the apparent trend "
                   "does NOT reproduce across sources (present in the reference, "
                   "absent elsewhere) -> not resolvable / likely source artifact")
    elif shape_all and trend_repro:
        verdict = "bounded tail robust AND trend reproduces in every source"
    elif not shape_all:
        verdict = "tail shape is NOT bounded in every source; interpret with care"
    else:
        verdict = "no significant trend in any source"
    return ArenaResult(sources=results, alpha=alpha, shape_bounded_all=shape_all,
                       trend_reproduces=trend_repro, reference_source=reference,
                       verdict=verdict)
