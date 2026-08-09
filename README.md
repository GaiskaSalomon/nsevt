# nsevt — non-stationary extreme-value tail risk with honest uncertainty

[![CI](https://github.com/GaiskaSalomon/nsevt/actions/workflows/ci.yml/badge.svg)](https://github.com/GaiskaSalomon/nsevt/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/nsevt.svg)](https://pypi.org/project/nsevt/)
[![Python versions](https://img.shields.io/pypi/pyversions/nsevt.svg)](https://pypi.org/project/nsevt/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/1328485209.svg)](https://zenodo.org/badge/latestdoi/1328485209)

`nsevt` is a small, dependency-light Python package (NumPy + SciPy) for the
**honest** analysis of trends in environmental extremes. It packages a workflow
that existing EVT tools (`extRemes`, `ismev`, `POT`, `pyextremes`, `texmex`) do
not offer as a single, tested pipeline:

1. **Bounded-tail detection** — a peaks-over-threshold GPD fit with a
   profile-likelihood interval for the shape and a bootstrap of the finite upper
   endpoint (`ξ<0` ⇒ a finite physical ceiling).
2. **A permutation-calibrated trend test** on the tail scale — exact in finite
   samples, avoiding the unreliable asymptotic χ² for a boundary-adjacent
   parameter on a few hundred exceedances.
3. **Power / minimum-detectable-effect** — turns any non-rejection into a
   quantitative statement of what the record can resolve.
4. **Multi-source transportability** ("evidence arena") — does an apparent trend
   *survive* changing the data source, or is it an instrumental artifact?
5. **Block-conformal prediction bands** — distribution-free coverage for extreme
   quantiles **under temporal dependence**.
6. **A two-scale Wasserstein trend test** for a series of distributions each
   estimated from a small per-period sample, with an exact
   location/scale/shape energy decomposition.

The numerics of the GPD and permutation machinery are ported verbatim from the
frozen, unit-tested research code, so results are reproducible and identical.

## Install

```bash
pip install -e .            # from this directory
pip install -e ".[demo]"    # also install the Streamlit demo dependencies
```

## Quick start

```python
import numpy as np, nsevt

rng = np.random.default_rng(0)
# 500 observations; tail excesses above u=40 are bounded (xi<0)
x = 40 + rng.gamma(2.0, 8.0, size=500)

fit = nsevt.gpd_pot(x, threshold=40)
print(fit.summary())
print("1-in-100 return level:", fit.return_level(100, rate=(x > 40).mean()))

# is there a trend in the tail scale? (block = year label per observation)
year = rng.integers(1980, 2024, size=500)
tr = nsevt.trend_permutation(x[x > 40] - 40, year[x > 40])
print("trend/decade:", round(tr["trend_per_decade"], 3),
      "p_perm:", tr["p_permutation"])

# what can the record resolve?
mde = nsevt.min_detectable_effect(x[x > 40] - 40, year[x > 40])
print("MDE per decade:", mde["mde_per_decade"])

# distribution-free 90% prediction band under dependence
band = nsevt.block_conformal(x, threshold=40, alpha=0.10)
print("upper bound at sigma:", band.predict_upper(fit.sigma))
```

## The multi-source arena

```python
arena = nsevt.transportability(
    [("operational", x_op, year_op),
     ("independent", x_ind, year_ind),
     ("homogenized", x_homog, year_homog)],
    threshold=40)
print(arena.table())
```

The verdict distinguishes a **robust bounded tail** (shape reproduces in every
source) from an **apparent trend that does not survive** the change of source —
the core methodological contribution of the accompanying research.

## Modules

| module | purpose |
|---|---|
| `nsevt.gpd` | POT-GPD fit, profile CI for ξ, endpoint bootstrap, return levels |
| `nsevt.trend` | permutation trend test, power/MDE, block-bootstrap trend CI |
| `nsevt.transportability` | multi-source evidence arena |
| `nsevt.conformal` | block / split conformal prediction bands for tails |
| `nsevt.twoscale` | two-scale Wasserstein distributional trend test + energy split |

## Tests

```bash
pip install ".[test]" && pytest
```

## License

MIT. See `LICENSE`.
