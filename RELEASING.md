# Releasing nsevt

This checklist keeps the Git tag, GitHub release, PyPI files, citation metadata,
and Zenodo archive tied to the same tested source revision. PyPI and published
GitHub releases are treated as immutable: a mistake is corrected in a new
version, never by replacing an existing distribution.

## One-time service configuration

1. Keep the GitHub--Zenodo integration enabled for `GaiskaSalomon/nsevt`.
2. Configure the PyPI Trusted Publisher for repository `GaiskaSalomon/nsevt`,
   workflow `publish.yml`, and GitHub environment `pypi`.
3. Configure the equivalent TestPyPI publisher for workflow `test-pypi.yml` and
   environment `testpypi`.
4. Require approval for the `pypi` environment. TestPyPI may remain unprotected.

Workflow actions are pinned to full commit hashes. Dependabot checks them
monthly; review the upstream release and keep the human-readable version comment
beside each updated hash.

No long-lived PyPI token belongs in GitHub secrets or a developer workstation.

## Prepare a version

Use semantic versioning. Patch releases correct compatible behavior or release
metadata, minor releases add backward-compatible functionality, and version 2
is reserved for changes to the stable public API or documented return schemas.

1. Start from a clean, current `main` branch.
2. Update `src/nsevt/_version.py`, which is the package version source.
3. Move the completed entries from `Unreleased` to a dated heading in
   `CHANGELOG.md`. Use the intended UTC publication date.
4. Update `version` and `date-released` in `CITATION.cff`, `version` in
   `.zenodo.json`, and the current-version text in `README.md`.
5. Run the local release gate in an environment built with the CI tooling pin
   (`pip install -e ".[dev,demo]" -c constraints/ci.txt`):

   ```bash
   python tools/release_check.py
   ruff check src tests demo tools
   mypy src/nsevt
   pytest -q --cov=nsevt --cov-report=term-missing --cov-fail-under=80
   python -m build
   python -m twine check dist/*
   python tools/write_checksums.py dist
   python tools/check_distributions.py dist
   ```

6. Exercise a release candidate through the manual `Publish to TestPyPI`
   workflow. TestPyPI versions are immutable too, so increment an `rcN` suffix
   before repeating a candidate.
7. Merge the reviewed release preparation. Do not tag an unmerged commit.

## Prepare the immutable artifacts

Create and push an annotated tag at the reviewed commit:

```bash
git tag -a vX.Y.Z -m "nsevt X.Y.Z"
git push origin vX.Y.Z
```

The `Prepare GitHub release` workflow reruns the quality gates, builds one wheel
and one source distribution, installs and smoke-tests both, writes
`SHA256SUMS`, verifies their embedded metadata, and creates a **draft** GitHub
release containing those three files. A failed job must be fixed through a new
commit and tag; do not move a published tag.

Before publishing the draft, verify:

- the tag and release title are `vX.Y.Z` and `nsevt X.Y.Z`;
- the workflow completed for every supported Python version;
- the release contains one `.whl`, one `.tar.gz`, and `SHA256SUMS`;
- the changelog describes every user-visible change and its claim boundary.

## Publish and verify

Publishing the GitHub draft triggers two independent consumers:

- `Publish to PyPI` downloads the attached artifacts, verifies their hashes and
  embedded versions, smoke-tests them again, and uploads those exact files via
  OIDC;
- Zenodo ingests the tagged software revision and creates an immutable
  version-specific DOI under concept DOI `10.5281/zenodo.21858232`.

After both finish:

1. Install `nsevt==X.Y.Z` from PyPI in a fresh environment and run
   `tests/smoke_release.py`.
2. Confirm the PyPI hashes equal `SHA256SUMS` from the GitHub release.
3. Confirm the Zenodo record reports version `X.Y.Z`, the correct tag, author,
   ORCID, MIT license, and release description.
4. Record the version-specific Zenodo DOI in the release notes or a subsequent
   documentation-only commit on `main`; keep the concept DOI in `CITATION.cff`
   so citations can resolve across versions.

## Validation data

Large simulation outputs are versioned separately from the Python package. A
validation deposit must identify the exact nsevt version and commit and include
its protocol, configurations, software environment, random-stream scheme,
requested/successful/failed replicate counts, stopping status, file hashes, and
reproduction commands. Name the first deposit for a release
`validation-vX.Y.Z`; a data-only correction becomes `validation-vX.Y.Z-r2`
without changing the software version.
