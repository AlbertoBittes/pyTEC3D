"""Custom build helpers for vendoring openpyxl alongside the API modules."""
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys

from setuptools.command.build_py import build_py as _build_py


class BuildPyWithVendoredOpenpyxl(_build_py):
    """Ensure openpyxl is available in ``api/local_packages`` before building."""

    description = "build Python modules and vendor openpyxl into api/local_packages"

    def run(self) -> None:
        repo_root = pathlib.Path(__file__).resolve().parent
        local_packages = repo_root / "api" / "local_packages"
        self._vendor_openpyxl(local_packages)
        super().run()
        self._copy_local_packages(local_packages)

    def _vendor_openpyxl(self, local_packages: pathlib.Path) -> None:
        local_packages.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        env.setdefault("PIP_DISABLE_PIP_VERSION_CHECK", "1")
        cmd = [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--upgrade",
            "openpyxl",
            f"--target={local_packages}",
        ]
        print("Vendoring openpyxl into", local_packages)
        subprocess.check_call(cmd, env=env)

    def _copy_local_packages(self, source: pathlib.Path) -> None:
        if not source.exists():
            return
        target = pathlib.Path(self.build_lib) / "local_packages"
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)
        print("Copied", source, "to", target)
