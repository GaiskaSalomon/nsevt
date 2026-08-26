#!/usr/bin/env python3
"""Verify release filenames, package metadata, and the SHA-256 manifest."""
from __future__ import annotations

import argparse
import hashlib
import tarfile
import zipfile
from email.parser import Parser
from pathlib import Path

from release_check import source_version


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata_from_wheel(path: Path) -> str:
    with zipfile.ZipFile(path) as archive:
        names = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(names) != 1:
            raise RuntimeError(f"{path.name}: expected one METADATA file")
        return archive.read(names[0]).decode("utf-8")


def metadata_from_sdist(path: Path) -> str:
    with tarfile.open(path, "r:gz") as archive:
        members = [
            member
            for member in archive.getmembers()
            if member.name.endswith("/PKG-INFO") and member.name.count("/") == 1
        ]
        if len(members) != 1:
            raise RuntimeError(f"{path.name}: expected one top-level PKG-INFO file")
        handle = archive.extractfile(members[0])
        if handle is None:
            raise RuntimeError(f"{path.name}: could not read PKG-INFO")
        return handle.read().decode("utf-8")


def expected_checksums(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        digest, separator, filename = line.partition("  ")
        if separator != "  " or len(digest) != 64 or not filename:
            raise RuntimeError(f"{path}: malformed checksum line {line!r}")
        entries[filename] = digest
    return entries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="dist", type=Path)
    parser.add_argument("--checksums", type=Path)
    args = parser.parse_args()
    version = source_version()
    wheels = sorted(args.directory.glob("*.whl"))
    sdists = sorted(args.directory.glob("*.tar.gz"))
    if len(wheels) != 1 or len(sdists) != 1:
        raise SystemExit("expected exactly one wheel and one source distribution")
    files = wheels + sdists
    expected_prefix = f"nsevt-{version}"
    if any(not path.name.startswith(expected_prefix) for path in files):
        raise SystemExit(f"distribution filenames must start with {expected_prefix!r}")

    for path, raw_metadata in (
        (wheels[0], metadata_from_wheel(wheels[0])),
        (sdists[0], metadata_from_sdist(sdists[0])),
    ):
        metadata = Parser().parsestr(raw_metadata)
        if metadata.get("Name") != "nsevt" or metadata.get("Version") != version:
            raise SystemExit(f"{path.name}: embedded name/version does not match nsevt {version}")

    checksums_path = args.checksums or args.directory / "SHA256SUMS"
    checksums = expected_checksums(checksums_path)
    if set(checksums) != {path.name for path in files}:
        raise SystemExit("SHA256SUMS must list exactly the wheel and source distribution")
    for path in files:
        if sha256(path) != checksums[path.name]:
            raise SystemExit(f"{path.name}: SHA-256 mismatch")
    print(f"nsevt {version} wheel, sdist, and checksums: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
