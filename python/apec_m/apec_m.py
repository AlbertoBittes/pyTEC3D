import os
import sys
import re
import gc
import importlib
import shutil
import atexit
from os import PathLike
from typing import Any, Optional, Union
import vtk

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOCAL_PACKAGES_PATH = os.path.join(CURRENT_DIR, "local_packages")

from .logging_utils import (
    get_logger,
    log_mesh_summary_info,
    log_python_stdio_redirect,
    log_solver_runtime_shutdown_failure,
    log_solver_shutdown_import_error,
)
from .apec_algorithm import ApecAlgorithm
from .apec_m_model import APEC_M_Model
from .model_edit import ModelEdit

LOGGER = get_logger("apec_m.APEC_M")

class APEC_M:
    """APEC solver/session manager that owns runtime/session and a Python model."""

    def __init__(self, num_threads: int = 1) -> None:
        if isinstance(num_threads, bool) or int(num_threads) < 1:
            raise ValueError("num_threads must be a positive integer.")

        self._num_threads = int(num_threads)
        self.model: APEC_M_Model = APEC_M_Model()
        self._mesh_filepath: Optional[str] = None
        self._solver_session: Optional[Any] = None
        self._null_output_stream = None
        self._closed = False

        self.algorithm = ApecAlgorithm(self)
        self._configure_python_output_root_only()
        atexit.register(self._cleanup_at_exit)

    @staticmethod
    def _world_rank() -> int:
        try:
            from mpi4py import MPI
        except ImportError:
            return 0
        return MPI.COMM_WORLD.Get_rank()

    def _configure_python_output_root_only(self) -> None:
        rank = self._world_rank()
        if rank == 0:
            return

        if self._null_output_stream is None or self._null_output_stream.closed:
            self._null_output_stream = open(os.devnull, "w", encoding="utf-8")
        sys.stdout = self._null_output_stream
        sys.stderr = self._null_output_stream
        log_python_stdio_redirect(LOGGER, rank)

    @property
    def num_threads(self) -> int:
        return self._num_threads

    @property
    def solver(self) -> Optional[Any]:
        return self._solver_session

    @property
    def mesh_filepath(self) -> Optional[str]:
        return self._mesh_filepath

    @mesh_filepath.setter
    def mesh_filepath(self, value: Optional[str]) -> None:
        if value is None:
            self._mesh_filepath = None
            return
        self._mesh_filepath = os.path.abspath(os.fspath(value))

    def __del__(self) -> None:
        try:
            self.close()
        except Exception as exc:
            LOGGER.debug("Exception while closing APEC_M in __del__: %s", exc)
        try:
            if self._null_output_stream is not None and not self._null_output_stream.closed:
                self._null_output_stream.close()
        except Exception as exc:
            LOGGER.debug("Exception while closing null output stream in __del__: %s", exc)

    def _cleanup_at_exit(self) -> None:
        self.close()

    def _close_solver_session(self) -> None:
        if self._solver_session is None:
            return

        try:
            close_method = getattr(self._solver_session, "close", None)
            if close_method is not None:
                close_method()
        finally:
            self._solver_session = None
            gc.collect()

    def _shutdown_solver_runtime(self) -> None:
        try:
            solver_module = self._require_solver_module()
        except Exception as exc:
            log_solver_shutdown_import_error(LOGGER, exc)
            return

        shutdown_runtime = getattr(solver_module, "shutdown_runtime", None)
        if shutdown_runtime is not None:
            try:
                shutdown_runtime()
            except Exception as exc:
                log_solver_runtime_shutdown_failure(LOGGER, exc)
                raise

    def close(self) -> None:
        if self._closed:
            return

        LOGGER.info("Closing APEC_M runtime resources.")
        self._close_solver_session()
        self._shutdown_solver_runtime()
        self._closed = True
        LOGGER.info("APEC_M runtime resources closed.")

    @classmethod
    def from_grid(
        cls,
        mesh: vtk.vtkUnstructuredGrid,
        num_threads: int = 1,
        mesh_name: Optional[str] = None,
    ) -> "APEC_M":
        if not isinstance(mesh, vtk.vtkUnstructuredGrid):
            raise TypeError("mesh must be a vtk.vtkUnstructuredGrid.")

        instance = cls(num_threads=num_threads)
        instance.model = APEC_M_Model.from_grid(mesh, mesh_name=mesh_name)
        instance.mesh_filepath = None
        instance._close_solver_session()
        log_mesh_summary_info(
            LOGGER,
            instance.model.mesh,
            instance.model._mesh_name or "APEC_Mesh",
            instance.mesh_filepath,
            os.getcwd(),
        )
        return instance

    def read_model(
        self,
        vtk_file_path: Union[str, PathLike[str]],
        reset_arrays: bool = True,
    ) -> None:
        if not isinstance(reset_arrays, bool):
            raise TypeError("reset_arrays must be a bool.")

        input_path = os.fspath(vtk_file_path)
        if os.path.isabs(input_path):
            mesh_path = input_path
        else:
            mesh_path = os.path.join(os.getcwd(), input_path)
        mesh_path = os.path.abspath(mesh_path)

        LOGGER.info("Loading model from '%s' (reset_arrays=%s).", mesh_path, reset_arrays)

        loaded_mesh = ModelEdit.read_mesh(mesh_path, reset_arrays=reset_arrays)
        resolved_name = os.path.splitext(os.path.basename(mesh_path))[0]
        loaded_model = APEC_M_Model.from_grid(loaded_mesh, mesh_name=resolved_name)
        base_cwd = os.getcwd()
        emat_dir = os.path.join(base_cwd, "emat")
        output_dir = os.path.join(base_cwd, "output")
        if self._world_rank() == 0:
            os.makedirs(emat_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)

            for entry_name in os.listdir(emat_dir):
                entry_path = os.path.join(emat_dir, entry_name)
                if os.path.isdir(entry_path):
                    shutil.rmtree(entry_path)
                elif os.path.exists(entry_path):
                    os.unlink(entry_path)

            for entry_name in os.listdir(output_dir):
                entry_path = os.path.join(output_dir, entry_name)
                if os.path.isfile(entry_path):
                    os.unlink(entry_path)

        self.model = loaded_model
        self.model._mesh_name = resolved_name
        self.mesh_filepath = mesh_path
        self._close_solver_session()
        log_mesh_summary_info(
            LOGGER,
            self.model.mesh,
            resolved_name,
            self.mesh_filepath,
            base_cwd,
        )

    def save_model(
        self,
        output_file: Optional[Union[str, PathLike[str]]] = None,
    ) -> None:
        if self._world_rank() != 0:
            return

        if output_file is None and self.model._mesh_name is None:
            raise ValueError("No mesh is loaded. Load a mesh before saving.")

        base_cwd = os.getcwd()
        if output_file is None:
            output_file = os.path.join(base_cwd, f"{self.model._mesh_name}.vtu")
        else:
            output_file = os.path.abspath(os.fspath(output_file))

        self.model._save_model(output_file, binary=True)
        LOGGER.info("Model saved to %s.", output_file)

    def write_output(self) -> None:
        if self.model._mesh_name is None:
            raise ValueError("No mesh is loaded. Load a mesh before writing output.")

        if self._world_rank() != 0:
            return

        mesh_name = self.model._mesh_name
        base_cwd = os.getcwd()
        output_dir = os.path.join(base_cwd, "output")
        os.makedirs(output_dir, exist_ok=True)

        pattern = re.compile(rf"^{re.escape(mesh_name)}_out_(\d+)\.vtu$")
        last_time_step = max(
            (
                int(match.group(1))
                for filename in os.listdir(output_dir)
                if (match := pattern.match(filename))
            ),
            default=0,
        )

        next_time_step = last_time_step + 1
        output_file = os.path.join(output_dir, f"{mesh_name}_out_{next_time_step}.vtu")
        self.save_model(output_file=output_file)

    @staticmethod
    def _require_solver_module() -> Any:
        try:
            return importlib.import_module("apec_m_core.apec_m_solver")
        except ImportError as exc:
            raise ImportError(
                "Could not import the C++ solver module 'apec_m_core.apec_m_solver'. "
                "Build/install the pybind11 extension before launching the solver."
            ) from exc

    def validate(self) -> None:
        from .validation_utils import validate_apec_model

        validate_apec_model(self.model)

    def load_model_into_solver(self) -> None:
        if self._closed:
            raise RuntimeError("This APEC_M instance is closed and cannot create a new solver session.")

        LOGGER.info("Validating model before creating solver session.")
        self.validate()

        LOGGER.info("Creating solver session from currently loaded model.")
        self._solver_session = self._load_model_into_solver()
        LOGGER.info("Solver session created.")

    def _load_model_into_solver(self) -> Any:
        solver_module = self._require_solver_module()

        if self.model._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        self._close_solver_session()

        solver = solver_module.ApecMSolver(self._num_threads)
        solver.load_model(self.model)
        return solver

    def solve(self) -> None:
        if self.model._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        if self._solver_session is None:
            raise RuntimeError("No solver session is loaded. Call load_model_into_solver() first.")

        if not hasattr(self._solver_session, "solve"):
            raise RuntimeError("Loaded solver session does not expose a solve() method.")

        LOGGER.info("Starting solver session solve().")
        self._solver_session.solve()
        LOGGER.info("Solver session solve() completed.")


if __name__ == "__main__":
    apec = APEC_M()
