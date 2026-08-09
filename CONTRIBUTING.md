# Contributing to nsevt

Contributions, bug reports, and feature requests are welcome.

## Reporting issues

Open an issue describing the problem with a minimal reproducible example
(package version, Python version, and a short script). For statistical
questions, please state the estimand and the expected behavior.

## Development

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev,demo]"
ruff check src tests demo
pytest --cov=nsevt
```

## Pull requests

1. Fork and create a feature branch.
2. Add or update tests under `tests/` for any behavior you change; new
   estimators must ship with a test that checks a known statistical property
   (e.g. size under the null, power under an alternative, coverage of a band).
3. Keep the dependency footprint minimal (NumPy + SciPy for the core).
4. Document new public functions with the estimand, assumptions, and exact
   boundary of any finite-sample or asymptotic claim.
5. Update `CHANGELOG.md` for user-visible behavior.
6. Ensure lint and tests pass and open the PR against `main`.

## Conduct and support

Participation is governed by [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md). See
[SUPPORT.md](SUPPORT.md) for help and [SECURITY.md](SECURITY.md) for private
vulnerability reports.
