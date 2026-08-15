# Validation scope

The test suite is designed to catch implementation and interpretation
regressions.

It checks:

- recovery of a negative GPD shape in controlled bounded-tail simulations;
- the distinction between a negative estimate and interval-supported negativity;
- adaptive profile intervals that contain the MLE;
- monotonic return levels and input validation;
- non-rejection under one seeded no-trend design and rejection under one strong
  seeded alternative;
- explicit positive and negative power grids;
- multi-source direction/status logic and source validation;
- conformal raw-value semantics and experimental labels;
- a null and a clear alternative for the exploratory quantile-function trend
  diagnostic; and
- the numerical location/scale/shape decomposition under a pure shift.

These tests do **not** prove nominal size, power, or coverage uniformly across
data-generating processes. Before scientific use, analysts should run
application-specific simulations over plausible tail shapes, sample sizes,
block dependence, missingness, thresholds, and measurement error.

Release validation consists of:

    ruff check src tests demo
    pytest --cov=nsevt --cov-report=term-missing --cov-fail-under=80
    python -m build
    python -m twine check dist/*

Continuous integration repeats lint, tests, package build, and a clean wheel
and source-distribution installation smoke test. The installed-distribution
smoke test exercises the continuous and grouped fits, block permutation,
grouped design and return levels, sequential Monte Carlo, finite-sample
calibration, multi-source status, and experimental namespace using public entry
points only. The publication workflow runs the same test on both artifacts
before the PyPI job can start.
