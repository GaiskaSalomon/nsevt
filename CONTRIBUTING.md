# Contributing to nsevt

`nsevt` is maintained as a single-author software project. Bug reports,
reproducible statistical cases, documentation corrections, and feature
requests are welcome, but external code commits and pull requests are not
merged. `GaiskaSalomon` is the sole commit author on `main`.

## Reporting issues

Open an issue describing the problem with a minimal reproducible example
(package version, Python version, and a short script). For statistical
questions, please state the estimand and the expected behavior.

## Proposing a change

Open an issue with the proposed behaviour, estimand, assumptions, and a minimal
example. If the proposal is accepted, the maintainer will reproduce and
implement it in a maintainer-authored commit. Opening an issue does not imply
that a feature will be added or that an inferential claim has been validated.

Automated dependency pull requests are advisory only. Their release notes,
source revision, compatibility, and CI result are reviewed; accepted updates
are then reproduced in a maintainer-authored commit and the automated pull
request is closed without merging its commits.

## Maintainer development

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev,demo]"
ruff check src tests demo
pytest --cov=nsevt
```

## Change requirements

1. Add or update tests under `tests/` for any behavior change; new
   estimators must ship with a test that checks a known statistical property
   (e.g. size under the null, power under an alternative, coverage of a band).
2. Keep the dependency footprint minimal (NumPy + SciPy for the core).
3. Document new public functions with the estimand, assumptions, and exact
   boundary of any finite-sample or asymptotic claim.
4. Update `CHANGELOG.md` for user-visible behavior.
5. Ensure the release check, lint, type checks, tests, and builds pass before
   updating `main`.

## Conduct and support

Participation is governed by [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md). See
[SUPPORT.md](SUPPORT.md) for help and [SECURITY.md](SECURITY.md) for private
vulnerability reports.
