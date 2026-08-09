---
title: 'nsevt: honest non-stationary extreme-value tail risk with permutation-calibrated trends, power analysis, multi-source transportability, and block-conformal prediction'
tags:
  - Python
  - extreme value theory
  - generalized Pareto distribution
  - non-stationarity
  - conformal prediction
  - climate extremes
authors:
  - name: Elí Gaiska Salomón Guzmán
    orcid: 0000-0002-3533-3900
    affiliation: 1
affiliations:
  - name: Colegio de Postgraduados, Campus Montecillo, Texcoco, Mexico
    index: 1
date: 8 August 2026
bibliography: paper.bib
---

# Summary

`nsevt` is a small, dependency-light Python package (NumPy + SciPy) for the
*honest* analysis of trends in environmental extremes. It bundles, as one tested
pipeline, a workflow that is scattered across the literature and largely absent
from existing extreme-value software: (i) peaks-over-threshold generalized
Pareto (GPD) estimation with a profile-likelihood interval for the shape and a
bootstrap of the finite upper endpoint, which detects a *bounded* tail
(`xi < 0`, a finite physical ceiling); (ii) a **permutation-calibrated** test of
a trend in the tail scale, exact in finite samples under exchangeability of the
period labels; (iii) a Monte-Carlo **power / minimum-detectable-effect** (MDE)
analysis that converts any non-rejection into a quantitative statement of what a
record can resolve; (iv) a **multi-source transportability** check that re-runs
the identical specification across data sources and reports whether an apparent
trend survives the change of source; (v) **block-conformal** prediction bands
that give distribution-free coverage for extreme quantiles under temporal
dependence; and (vi) a **two-scale Wasserstein / Fréchet-mean** trend test for a
series of distributions each estimated from a small per-period sample, with an
exact location/scale/shape energy decomposition.

The GPD and permutation numerics are ported verbatim from the frozen,
unit-tested code of the accompanying research, so results are reproducible and
identical. The package ships a test suite and a Streamlit browser demo.

# Statement of need

Practitioners routinely fit trends to environmental extremes—heatwave
intensities, flood peaks, wind gusts, tropical-cyclone rapid-intensification
magnitude—and over-interpret them. Three failure modes are pervasive: reporting
a trend without asking whether the record could even *detect* it (no power
analysis); relying on an asymptotic chi-square that is unreliable for a
boundary-adjacent shape parameter on a few hundred exceedances; and treating a
trend found in one heterogeneous data product as a property of the phenomenon
rather than of the instrument.

Established EVT toolkits—`extRemes` [@gilleland2016extremes], `ismev`, `POT`,
`texmex` in R, and `pyextremes` in Python—provide excellent stationary and
non-stationary GPD/GEV fitting, but none packages the *honest* companions:
permutation calibration, a power/MDE analysis, a multi-source transportability
("evidence arena") verdict, and conformal prediction bands that remain valid
under temporal dependence. Conformal prediction [@vovk2005; @barber2023] and
distribution-of-distributions trend testing in Wasserstein space
[@panaretos2019] have likewise lacked an implementation targeted at the
sparse, dependent, tail regime. `nsevt` fills this gap with a single, small API,
built on the Pickands–Balkema–de Haan foundation [@coles2001] and the
GPD-regularity results of @smith1985.

The design was distilled from a program of work on the extreme magnitude of
tropical-cyclone rapid intensification, where the same tools separated a robust
finite ceiling (bounded tail reproduced across every intensity product,
including a homogenized satellite record) from an apparent scale trend that did
not survive homogenization—an artifact of best-track heterogeneity. The
transportability and power machinery that made that distinction defensible is
general, and `nsevt` exposes it for any environmental extreme.

# Functionality

The public API mirrors the workflow:

```python
import nsevt
fit  = nsevt.gpd_pot(x, threshold=u)        # xi, profile CI, endpoint, return levels
tr   = nsevt.trend_permutation(z, block)    # permutation-calibrated trend p-value
mde  = nsevt.min_detectable_effect(z, block)# what can the record resolve?
arena= nsevt.transportability(sources, u)   # does the trend survive the source?
band = nsevt.block_conformal(x, u, alpha)   # coverage-guaranteed tail bound
ts   = nsevt.twoscale_trend(samples)        # distributional trend from small samples
```

Each routine is documented with the estimand and the finite-sample guarantee it
provides, and returns plain dictionaries or small dataclasses for downstream
use. A Streamlit application reproduces the full analysis interactively.

# Acknowledgements

This work was conducted as part of doctoral studies at the Colegio de
Postgraduados. The author acknowledges a doctoral scholarship from Mexico's
SECIHTI.

# References
