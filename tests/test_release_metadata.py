"""Release metadata must agree with the installed distribution."""
from __future__ import annotations

import importlib.metadata as metadata
import subprocess
import sys
from pathlib import Path

import nsevt

ROOT = Path(__file__).resolve().parents[1]


def test_installed_and_runtime_versions_agree():
    assert metadata.version("nsevt") == nsevt.__version__


def test_release_metadata_agree():
    result = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "release_check.py")],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
