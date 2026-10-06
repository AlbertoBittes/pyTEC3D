import copy
import os
import sys
import re
import gc
import importlib
import shutil
import atexit
from os import PathLike
from typing import Any, Dict, Mapping, Optional, Union
import vtk

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOCAL_PACKAGES_PATH = os.path.join(CURRENT_DIR, "local_packages")

from pytec3d.logging_utils import (
    get_logger,
    log_mesh_summary_info,
    log_python_stdio_redirect,
    log_solver_runtime_shutdown_failure,
    log_solver_shutdown_import_error,
)
from pytec3d.algorithm import pyTEC3DAlgorithm
from pytec3d.constants import DEFAULT_SOLVER_SETTINGS, DEFAULT_TIME_PARAMETERS
from pytec3d.model import pyTEC3D_Model
from pytec3d.model_edit import ModelEdit

LOGGER = get_logger("pytec3d.Solver")

class pyTEC3DSolver:
    """pyTEC3D session manager that owns runtime/session, a Python model and the run settings.

    Attributes:
        time_step: Dict with ``time`` (end time of the step), ``theta`` (time integration parameter in
            [0, 1]) and ``transient_effects`` (False for a static solve).
        solver_settings: Dict with the PETSc ``ksp_type``, ``pc_type`` and the linear tolerances and
            iteration limit. The solver reads it at every solve.
    """

    def __init__(self, num_threads: int = 1) -> None:
        """Initialize a session manager with an empty model.

        Args:
            num_threads: Number of threads used by the C++ solver backend.

        Raises:
            ValueError: If num_threads is not a positive integer.
        """
        if isinstance(num_threads, bool) or int(num_threads) < 1:
            raise ValueError("num_threads must be a positive integer.")

        self._num_threads = int(num_threads)
        self.model: pyTEC3D_Model = pyTEC3D_Model()
        self.time_step: Dict[str, Any] = copy.deepcopy(DEFAULT_TIME_PARAMETERS)
        self._solver_settings: Dict[str, Any] = copy.deepcopy(DEFAULT_SOLVER_SETTINGS)
        self._mesh_filepath: Optional[str] = None
        self._solver_session: Optional[Any] = None
        self._null_output_stream = None
        self._closed = False

        self.algorithm = pyTEC3DAlgorithm(self)
        self._configure_python_output_root_only()
        atexit.register(self._cleanup_at_exit)

    @staticmethod
    def _world_rank() -> int:
        """Return MPI rank for COMM_WORLD, or 0 when mpi4py is unavailable."""
        try:
            from mpi4py import MPI
        except ImportError:
            return 0
        return MPI.COMM_WORLD.Get_rank()

    def _configure_python_output_root_only(self) -> None:
        """
         configure python output root only.
        """
        rank = self._world_rank()
        if rank == 0:
            return

        if self._null_output_stream is None or self._null_output_stream.closed:
            self._null_output_stream = open(os.devnull, "w", encoding="utf-8")
        sys.stdout = self._null_output_stream
        sys.stderr = self._null_output_stream
        log_python_stdio_redirect(LOGGER, rank)

    @property
    def solver_settings(self) -> Dict[str, Any]:
        """PETSc settings dict; the same object is shared with the C++ solver session."""
        return self._solver_settings

    @solver_settings.setter
    def solver_settings(self, value: Mapping[str, Any]) -> None:
        # Update in place: a loaded solver session holds this dict.
        self._solver_settings.clear()
        self._solver_settings.update(value)

    @property
    def num_threads(self) -> int:
        """Thread count configured for solver session creation."""
        return self._num_threads

    @property
    def solver(self) -> Optional[Any]:
        """Current active solver session, if one has been created."""
        return self._solver_session

    @property
    def mesh_filepath(self) -> Optional[str]:
        """Absolute path of the currently loaded mesh file, if known."""
        return self._mesh_filepath

    @mesh_filepath.setter
    def mesh_filepath(self, value: Optional[str]) -> None:
        """Store mesh file path as an absolute path or clear it with None."""
        if value is None:
            self._mesh_filepath = None
            return
        self._mesh_filepath = os.path.abspath(os.fspath(value))

    def __del__(self) -> None:
        """Best-effort resource cleanup during garbage collection."""
        try:
            self.close()
        except Exception as exc:
            LOGGER.debug("Exception while closing pyTEC3D in __del__: %s", exc)
        try:
            if self._null_output_stream is not None and not self._null_output_stream.closed:
                self._null_output_stream.close()
        except Exception as exc:
            LOGGER.debug("Exception while closing null output stream in __del__: %s", exc)

    def _cleanup_at_exit(self) -> None:
        """
         cleanup at exit.
        """
        self.close()

    def _close_solver_session(self) -> None:
        """
         close solver session.
        """
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
        """
         shutdown solver runtime.
        """
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
        """
        Close.
        """
        if self._closed:
            return

        LOGGER.info("Closing pyTEC3D runtime resources.")
        self._close_solver_session()
        self._shutdown_solver_runtime()
        self._closed = True
        LOGGER.info("pyTEC3D runtime resources closed.")

    @classmethod
    def from_grid(
        cls,
        mesh: vtk.vtkUnstructuredGrid,
        num_threads: int = 1,
        mesh_name: Optional[str] = None,
    ) -> "pyTEC3DSolver":
        """Create an pyTEC3D instance directly from a VTK grid.

        Args:
            mesh: Input unstructured grid.
            num_threads: Solver thread count.
            mesh_name: Optional logical mesh name used in outputs.

        Returns:
            pyTEC3D: Initialized instance with model loaded from mesh.

        Raises:
            TypeError: If mesh is not a vtk.vtkUnstructuredGrid.
        """
        if not isinstance(mesh, vtk.vtkUnstructuredGrid):
            raise TypeError("mesh must be a vtk.vtkUnstructuredGrid.")

        instance = cls(num_threads=num_threads)
        instance.model = pyTEC3D_Model.from_grid(mesh, mesh_name=mesh_name)
        instance.mesh_filepath = None
        instance._close_solver_session()
        log_mesh_summary_info(
            LOGGER,
            instance.model.mesh,
            instance.model._mesh_name or "pyTEC3Desh",
            instance.mesh_filepath,
            os.getcwd(),
        )
        return instance

    def read_model(
        self,
        vtk_file_path: Union[str, PathLike[str]],
        reset_arrays: bool = True,
    ) -> None:
        """Load a VTU model file into this instance and reset runtime outputs.

        Args:
            vtk_file_path: Path to the mesh/model file to load.
            reset_arrays: Whether to reinitialize required VTK arrays on load.

        Raises:
            TypeError: If reset_arrays is not a bool.
        """
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
        loaded_model = pyTEC3D_Model.from_grid(loaded_mesh, mesh_name=resolved_name)
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
        """Persist the current model to a VTU file from rank 0 only.

        Args:
            output_file: Optional output path. Defaults to <mesh_name>.vtu in cwd.

        Raises:
            ValueError: If no mesh is loaded and no output path can be inferred.
        """
        if self._world_rank() != 0:
            return

        if output_file is None and self.model._mesh_name is None:
            raise ValueError("No mesh is loaded. Load a mesh before saving.")

        base_cwd = os.getcwd()
        if output_file is None:
            output_file = os.path.join(base_cwd, f"{self.model._mesh_name}.vtu")
        else:
            output_file = os.path.abspath(os.fspath(output_file))

        self.model.save_model(output_file, binary=True)
        LOGGER.info("Model saved to %s.", output_file)

    def write_output(self) -> None:
        """
        Write output.
        
        Raises:
            ValueError: If an error condition is met.
        """
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
        """Import and return the compiled pybind solver extension module.

        Raises:
            ImportError: If the extension cannot be imported.
        """
        try:
            return importlib.import_module("pytec3d_solver.pytec3d_solver")
        except ImportError as exc:
            raise ImportError(
                "Could not import the C++ solver module 'pytec3d_solver.pytec3d_solver'. "
                "Build/install the pybind11 extension before launching the solver."
            ) from exc

    def validate(self) -> None:
        """
        Validate.
        """
        from pytec3d.validation_utils import validate_pytec3dodel, validate_solver_settings

        validate_solver_settings(self._solver_settings)
        validate_pytec3dodel(self.model)

    def validate_pre_solve(self) -> None:
        """
        Validate pre solve.
        """
        from pytec3d.validation_utils import validate_pre_solve_inputs

        validate_pre_solve_inputs(self.model, self.time_step)

    def load_model_into_solver(self) -> None:
        """
        Load model into solver.
        
        Raises:
            RuntimeError: If an error condition is met.
        """
        if self._closed:
            raise RuntimeError("This pyTEC3D instance is closed and cannot create a new solver session.")

        LOGGER.info("Validating model before creating solver session.")
        self.validate()

        LOGGER.info("Creating solver session from currently loaded model.")
        self._solver_session = self._load_model_into_solver()
        LOGGER.info("Solver session created.")

    def _load_model_into_solver(self) -> Any:
        """Create a fresh backend solver and load the current model into it.

        Returns:
            Any: Backend solver object exposed by the pybind module.

        Raises:
            RuntimeError: If no mesh is loaded in the current model.
        """
        solver_module = self._require_solver_module()

        if self.model._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        self._close_solver_session()

        solver = solver_module.pyTEC3DSolver(self._num_threads)
        solver.load_model(self.model, self._solver_settings)
        return solver

    def solve(self) -> None:
        """
        Solve.
        
        Raises:
            RuntimeError: If an error condition is met.
        """
        if self.model._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        if self._solver_session is None:
            raise RuntimeError("No solver session is loaded. Call load_model_into_solver() first.")

        if not hasattr(self._solver_session, "solve"):
            raise RuntimeError("Loaded solver session does not expose a solve() method.")

        self.validate_pre_solve()
        LOGGER.info("Starting solver session solve().")
        self._solver_session.solve(self.time_step)
        LOGGER.info("Solver session solve() completed.")


if __name__ == "__main__":
    pytec3d = pyTEC3DSolver()
