# Contributing to nsevt

Contributions, bug reports, and feature requests are welcome.

## Reporting issues
Open an issue describing the problem with a minimal reproducible example
(package version, Python version, and a short script). For statistical
questions, please state the estimand and the expected behavior.

## Development
```bash
python -m venv .venv && . .venv/bin/activate
pip install -e ".[test,demo]"
pytest            # run the test suite
```

## Pull requests
1. Fork and create a feature branch.
2. Add or update tests under `tests/` for any behavior you change; new
   estimators must ship with a test that checks a known statistical property
   (e.g. size under the null, power under an alternative, coverage of a band).
3. Keep the dependency footprint minimal (NumPy + SciPy for the core).
4. Document new public functions with the estimand and the finite-sample
   guarantee they provide.
5. Ensure `pytest` passes and open the PR against `main`.

## Code of conduct
Be respectful and constructive. By participating you agree to uphold a
harassment-free experience for everyone.
