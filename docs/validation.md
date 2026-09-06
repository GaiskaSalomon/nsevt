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

Every finding of the 2026-09-06 review (`docs/revision-arquitectura-2026-09-06.md`)
has an automated regression; `tests/test_api_contracts.py` locks the frozen 1.x
return keys, dataclass fields and aliases. Tests that run a full sequential
Monte Carlo or power loop are marked `slow`; `pytest -m "not slow"` is the fast
inner loop, and CI runs the whole suite.

These tests do **not** prove nominal size, power, or coverage uniformly across
data-generating processes. Before scientific use, analysts should run
application-specific simulations over plausible tail shapes, sample sizes,
block dependence, missingness, thresholds, and measurement error.
`tools/validation_campaign.py` runs such a grid and writes a JSON archive with,
per (shape, size) cell, the coverage / bias / RMSE and their MCSE, the
requested/successful/failed replicate counts, the stopping status, every seed,
and the package version and git commit; a cell that does not reach the MCSE
tolerance is written with `"resolved": false`. A full run of it is deposited
separately as `validation-vX.Y.Z` (see `RELEASING.md`).

Release validation consists of:

    python tools/release_check.py
    ruff check src tests demo tools
    mypy src/nsevt
    pytest --cov=nsevt --cov-branch --cov-report=term-missing --cov-fail-under=80
    python -m build
    python -m twine check dist/*
    python tools/write_checksums.py dist
    python tools/check_distributions.py dist

Continuous integration repeats lint, tests (branch coverage, on Python 3.9
through 3.14, plus a minimum-dependency job on the declared NumPy/SciPy floors
and a build against the minimum declared setuptools), package build, and a
clean wheel and source-distribution installation smoke test. The installed-distribution
smoke test exercises the continuous and grouped fits, block permutation,
grouped design and return levels, sequential Monte Carlo, finite-sample
calibration, multi-source status, and experimental namespace using public entry
points only. A tagged candidate is attached to a draft GitHub release only after
all quality gates pass. The publication workflow then downloads those exact
artifacts, verifies their SHA-256 manifest and embedded versions, and repeats
both smoke installations before the PyPI job can start.
