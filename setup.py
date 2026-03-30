from __future__ import annotations

import glob
import pathlib
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


class BuildPyWithSolverExtension(_build_py):
    """Copy a prebuilt solver extension into the Python build output."""

    def run(self) -> None:
        repo_root = pathlib.Path(__file__).resolve().parent
        package_dir = repo_root / "python" / "apec_m"
        super().run()
        self._copy_solver_extension(package_dir)

    def _copy_solver_extension(self, package_dir: pathlib.Path) -> None:
        patterns = [
            str(package_dir / "apec_m_solver*.so"),
            str(package_dir / "apec_m_solver*.pyd"),
            str(package_dir / "apec_m_solver*.dylib"),
        ]
        candidates = []
        for pattern in patterns:
            candidates.extend(glob.glob(pattern))

        if not candidates:
            print("No prebuilt apec_m_solver extension found under", package_dir)
            return

        build_lib = pathlib.Path(self.build_lib) / "apec_m"
        build_lib.mkdir(parents=True, exist_ok=True)
        for candidate in candidates:
            source = pathlib.Path(candidate)
            target = build_lib / source.name
            shutil.copy2(source, target)
            print("Copied", source, "to", target)


setup(
    cmdclass={"build_py": BuildPyWithSolverExtension},
)
