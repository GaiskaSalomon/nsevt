#!/usr/bin/env python3
"""Fail unless every release-facing version field is internally consistent."""
from __future__ import annotations

import argparse
import ast
import json
import re
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION_FILE = ROOT / "src" / "nsevt" / "_version.py"
VERSION_PATTERN = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+(?:rc[0-9]+)?$")
CONCEPT_DOI = "10.5281/zenodo.21858232"


def source_version() -> str:
    """Read ``__version__`` without importing nsevt or its dependencies."""
    tree = ast.parse(VERSION_FILE.read_text(encoding="utf-8"), filename=str(VERSION_FILE))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets):
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                return node.value.value
    raise RuntimeError(f"{VERSION_FILE} must assign a literal string to __version__")


def quoted_field(text: str, field: str) -> str | None:
    match = re.search(rf'^\s*{re.escape(field)}:\s*["\']([^"\']+)["\']\s*$', text, re.MULTILINE)
    return None if match is None else match.group(1)


def check_release(tag: str | None = None) -> list[str]:
    version = source_version()
    errors: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            errors.append(message)

    require(bool(VERSION_PATTERN.fullmatch(version)), f"unsupported version syntax: {version!r}")

    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    require(
        bool(re.search(r'^dynamic\s*=\s*\[\s*["\']version["\']\s*\]\s*$', pyproject, re.MULTILINE)),
        "pyproject.toml must declare project.version as dynamic",
    )
    require(
        bool(re.search(
            r'^version\s*=\s*\{\s*attr\s*=\s*["\']nsevt\._version\.__version__["\']\s*\}\s*$',
            pyproject,
            re.MULTILINE,
        )),
        "pyproject.toml must source the version from nsevt._version.__version__",
    )

    package_init = (ROOT / "src" / "nsevt" / "__init__.py").read_text(encoding="utf-8")
    require(
        "from ._version import __version__" in package_init,
        "nsevt.__init__ must re-export the canonical version",
    )
    require(
        not bool(re.search(r"^__version__\s*=", package_init, re.MULTILINE)),
        "nsevt.__init__ must not define a second version",
    )

    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    citation_version = quoted_field(citation, "version")
    release_date = quoted_field(citation, "date-released")
    require(citation_version == version, f"CITATION.cff version is {citation_version!r}, expected {version!r}")
    require(quoted_field(citation, "doi") == CONCEPT_DOI, "CITATION.cff must cite the Zenodo concept DOI")
    try:
        parsed_date = date.fromisoformat(release_date or "")
    except ValueError:
        parsed_date = None
    require(parsed_date is not None, "CITATION.cff date-released must be an ISO date")

    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    require(zenodo.get("version") == version, f".zenodo.json version is {zenodo.get('version')!r}, expected {version!r}")

    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    expected_heading = f"## [{version}] - {release_date}"
    require(expected_heading in changelog, f"CHANGELOG.md is missing {expected_heading!r}")
    require(changelog.count(expected_heading) == 1, "CHANGELOG.md must contain exactly one release heading")

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    require(f"current release is **nsevt {version}**" in readme, "README.md current-release text is stale")
    require(f"pip install nsevt=={version}" in readme, "README.md exact-install command is stale")
    require(CONCEPT_DOI in readme, "README.md must expose the Zenodo concept DOI")

    if tag is not None:
        require(tag == f"v{version}", f"release tag is {tag!r}, expected 'v{version}'")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", help="Git tag to compare with the canonical version")
    args = parser.parse_args()
    errors = check_release(args.tag)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"nsevt {source_version()} release metadata: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
