---
title: 'nsevt: calibrated GPD tail-trend inference and multi-source robustness in Python'
tags:
  - Python
  - extreme value theory
  - generalized Pareto distribution
  - non-stationarity
  - permutation test
  - power analysis
authors:
  - name: Elí Gaiska Salomón Guzmán
    orcid: 0000-0002-3533-3900
    affiliation: 1
affiliations:
  - name: Colegio de Postgraduados, Campus Montecillo, Texcoco, Mexico
    index: 1
date: 9 August 2026
bibliography: paper.bib
---

# Summary

`nsevt` is a dependency-light Python package for inferential questions that
arise after a peaks-over-threshold generalized Pareto (GPD) model has been
selected. It connects four steps that are often performed separately: GPD shape
estimation with a profile-likelihood interval; a likelihood-ratio (LR) test for
a linear trend in log scale calibrated by complete-block label permutation;
Monte Carlo power and signed minimum-detectable-effect (MDE) analysis at the
observed temporal design; and application of the same frozen specification to
multiple data sources with power-aware interpretation of non-reproduction.

The software separates estimation from inferential support. A negative shape
point estimate implies a finite statistical endpoint under the selected
threshold and GPD model. `nsevt` labels a negative shape as supported only when
the entire 95% profile interval lies below zero, and describes the resulting
endpoint as model-conditional rather than physical. Likewise, a non-significant
trend in a second source is not treated as disagreement unless that source has
adequate simulated power for the signed reference effect.

# Statement of need

Environmental records can be short, heterogeneous, and clustered by season or
year. In that setting, three shortcuts can produce conclusions stronger than
the data warrant. First, the sign of the fitted GPD shape can be reported as if
it were known without uncertainty. Second, a non-stationary fit can be judged
only by a default asymptotic p-value even when the effective number of temporal
blocks is modest. Third, failure to reproduce a trend across products can be
called an artifact without checking whether the comparison product could have
detected the reference effect.

`nsevt` provides one auditable workflow for these questions. Its intended users
are researchers comparing tail behavior or tail-scale trends in environmental
products while retaining explicit control over threshold, block definition,
significance level, simulation count, and random seed. The package does not
choose the threshold, establish independence of exceedances, identify causal
effects, or turn a statistical endpoint into a physical limit. Those decisions
remain part of the scientific design and are documented as assumptions rather
than hidden defaults.

# State of the field

Established extreme-value packages provide broad model-fitting functionality.
In R, `extRemes` supports stationary and non-stationary extreme-value analysis
[@gilleland2016extremes], while `POT` focuses on peaks-over-threshold methods
[@ribatet2007]. In Python, `pyextremes` supplies extraction, fitting,
diagnostics, and return-value workflows [@bocharov2026]. These packages are more
general EVT environments than `nsevt` and should be preferred for many standard
analyses.

The narrower contribution of `nsevt` is integration of profile-shape
interpretation, complete-block permutation calibration, design-specific
power/MDE, and a multi-source status that uses effect direction and power. The
implementation builds on standard GPD theory [@coles2001] and acknowledges the
nonregular likelihood regime described by @smith1985. Plus-one Monte Carlo
p-values avoid zero estimates and preserve randomization-test validity under
the stated exchangeability design [@phipson2010].

The repository also contains explicitly experimental conformal and
distribution-valued diagnostics, motivated by work on conformal prediction
beyond exchangeability [@barber2023] and Wasserstein statistics
[@panaretos2019]. They are labeled experimental because the present
implementations do not establish general dependent-data guarantees. They are
not part of the stable inferential contribution claimed here.

# Software design

The stable core has three modules and plain NumPy inputs. `nsevt.gpd` fits the
two-parameter GPD by multi-start maximum likelihood. The shape interval is
obtained by adaptive inversion of the profile LR statistic rather than a fixed
grid that can silently clip a confidence limit. The API exposes separate
`bounded_estimate` and `bounded_supported` properties. If the point estimate is
negative, endpoint bootstrap intervals are computed conditionally among
negative-shape bootstrap fits and returned together with the fraction of all
successful fits having negative shape. A warning marks estimates at or below
the usual regularity boundary \(\xi=-1/2\); optimization is restricted to
\(\xi>-1\), below which the likelihood may be unbounded at the sample maximum.

`nsevt.trend` models
\(\log\{\sigma(t)\}=\beta_0+\beta_1t\), where \(t\) is measured in decades. The
primary statistic is twice the maximized log-likelihood difference between the
varying- and constant-scale models. Under the null, complete block labels are
permuted: all exceedances in a block stay together while that block receives a
permuted time label. The interpretation is exact under exchangeability of these
block labels; otherwise it is a design-based sensitivity analysis. The output
includes the plus-one permutation p-value, its Monte Carlo standard error, the
asymptotic reference p-value for comparison, and successful permutation count.

Power simulations retain the observed number of excesses and their block
layout, simulate from a fitted GPD with a signed injected trend, and evaluate
the same LR critical value. MDE output is directional because increasing and
decreasing scale need not have identical finite-sample detectability. Simulation
standard errors and successful-replicate counts accompany the curve.

`nsevt.transportability` applies that frozen specification across named sources.
A trend is reproduced only if the reference result is significant and all
other sources are significant in the same direction. Otherwise the result is
classified as opposite-direction evidence, non-reproduction with adequate
simulated power, unresolved, or no reference signal. Wording in the returned
verdict explicitly avoids causal attribution.

Input validation rejects non-finite values, invalid excesses, unmatched source
layouts, duplicate source names, and invalid simulation controls. Randomized
routines use NumPy generators with user-visible seeds. Unit and statistical
regression tests cover shape semantics, profile bracketing, null and alternative
trend behavior, signed MDEs, multi-source decisions, and experimental API claim
boundaries. Continuous integration tests the supported Python versions and builds
both source and wheel distributions.

# Research impact statement

`nsevt` was extracted to make a recurring research workflow inspectable and
reusable outside one analysis pipeline. It can support sensitivity analyses in
which investigators must report not only a trend estimate but also its
randomization calibration, detectable effect scale, and behavior across
measurement products. Because the project is newly public, no independent
community adoption or downstream scientific discoveries are claimed at the
time of writing. The immediate impact is reproducibility: methodological claim
boundaries, deterministic seeds, validation tests, and power-aware verdicts are
available in one installable implementation.

# Acknowledgements

This work was conducted as part of doctoral studies at the Colegio de
Postgraduados. The author acknowledges a doctoral scholarship from Mexico's
SECIHTI.

# References
