# constraints/

Version pins for reproducibility, kept separate from the library's dependency
ranges in `pyproject.toml` (which stay deliberately wide).

## `ci.txt`

Pins the tooling whose pass/fail verdict depends on its version: Ruff, mypy,
build, twine. Install with

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

matches CI. It is applied on CI only in the single-Python `lint` and
`build`/release jobs. It is **not** applied to the multi-version `test` matrix
(pytest 9 has no Python 3.9 build, and a coverage number does not depend on the
pytest patch version), and it does not pin NumPy/SciPy — the library supports a
range and each supported Python resolves a compatible build.

## `science.txt` (not present yet)

An exact runtime pin (Python, NumPy, SciPy and their transitive builds) for a
versioned simulation campaign. Add it alongside the campaign's data deposit
(`validation-vX.Y.Z`), naming the interpreter and platform it was resolved on.
