# constraints/

Version pins for reproducibility, kept separate from the library's dependency
ranges in `pyproject.toml` (which stay deliberately wide).

## `ci.txt`

The development and CI tooling that the `main`-branch quality gate depends on:
Ruff, mypy, pytest, pytest-cov, coverage, build, twine. Install with

```bash
pip install -e ".[dev,demo]" -c constraints/ci.txt
```

so a local run of

```bash
python tools/release_check.py
ruff check src tests demo tools
mypy src/nsevt
pytest -q --cov=nsevt --cov-report=term-missing --cov-fail-under=80
```

matches CI. It does not pin NumPy/SciPy: the library supports a range and each
supported Python resolves a compatible build.

## `science.txt` (not present yet)

An exact runtime pin (Python, NumPy, SciPy and their transitive builds) for a
versioned simulation campaign. Add it alongside the campaign's data deposit
(`validation-vX.Y.Z`), naming the interpreter and platform it was resolved on.
