#!/usr/bin/env python3
"""Write a deterministic SHA-256 manifest for one wheel and one sdist."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="dist", type=Path)
    args = parser.parse_args()
    files = sorted((*args.directory.glob("*.whl"), *args.directory.glob("*.tar.gz")))
    if len(files) != 2 or sum(path.suffix == ".whl" for path in files) != 1:
        raise SystemExit("expected exactly one wheel and one source distribution")
    output = args.directory / "SHA256SUMS"
    output.write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in files),
        encoding="utf-8",
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
