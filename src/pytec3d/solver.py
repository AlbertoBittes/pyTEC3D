"""The solver: a model, the settings of its solution, and the C++ solver session that solves it."""

from __future__ import annotations

import atexit
import contextlib
import copy
import csv
import functools
import importlib
import importlib.metadata
import os
import re
import sys
import time
import weakref
from types import ModuleType, TracebackType
from typing import Any, TextIO

import numpy as np
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

from . import mesh as meshes
from ._logging import close_log_file, get_logger, message_counts, open_log_file, world_rank
from ._report import Table, print_table, provenance_table
from ._typing import Location, PathType
from .errors import SolverStateError
from .model import Model
from .validation import check_before_solve, check_model, check_solver_settings, load_time_arrays

LOGGER = get_logger("pytec3d.solver")

#: Name of the compiled solver module (pybind11 extension).
BACKEND_MODULE = "pytec3d_solver.pytec3d_solver"

#: Default solver settings. ``linear``: the PETSc KSP (a direct LU factorisation with MUMPS by default), used alone
#: or inside each Newton iteration. ``nonlinear``: the libMesh Newton solver of the models with nonlinear laws.
#: ``petsc_options``: PETSc options set before each linear solve ("" for a flag); None removes a default one.
DEFAULT_SOLVER_SETTINGS: dict[str, Any] = {
    "linear": {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "relative_tolerance": 1e-6,
        "absolute_tolerance": 1e-10,
        "divergence_tolerance": None,  # the PETSc default
        "max_iterations": 100,
        "initial_guess_nonzero": False,
        "factor_solver_type": "mumps",  # of lu / cholesky; None for the PETSc default
    },
    "nonlinear": {
        "max_iterations": 50,
        "absolute_residual_tolerance": 1e-35,
        "relative_residual_tolerance": 1e-8,
        "absolute_step_tolerance": 1e-8,
        "relative_step_tolerance": 1e-8,
        "require_residual_reduction": False,
        "require_finite_residual": True,
        "brent_line_search": True,  # False: plain backtracking
        "minimum_step_length": 1e-5,
        "linear_tolerance_multiplier": 1e-3,
        "continue_after_max_iterations": True,
        "continue_after_backtrack_failure": False,
        "verbose": True,
    },
    "petsc_options": {
        "-ksp_monitor": "",
        "-ksp_view": "",
        "-ksp_converged_reason": "",
        "-log_view": "",
        "-ksp_monitor_true_residual": "",
        "-ksp_monitor_solution": "",
        "-mat_mumps_icntl_4": "2",
    },
}

#: Default time step: theta for first-order formulations, Newmark beta and gamma (average acceleration) for
#: second-order ones (the "time_order" of the element types).
DEFAULT_TIME_STEP: dict[str, Any] = {
    "time": 1.0,
    "transient_effects": True,
    "theta": 1.0,
    "beta": 0.25,
    "gamma": 0.5,
}

#: Directory, in the working directory, that receives the output meshes.
OUTPUT_DIRECTORY = "output"


def _backend() -> ModuleType:
    """The compiled solver module.

    Raises:
        ImportError: If the extension has not been built.
    """
    try:
        return importlib.import_module(BACKEND_MODULE)
    except ImportError as exc:
        raise ImportError(
            f"Could not import the C++ solver module '{BACKEND_MODULE}': build it first (build.sh)."
        ) from exc


class Solver:
    """A pyTEC3D model, the settings of its solution and the C++ solver session that solves it.

    Use it as a context manager, or call :meth:`close` when done::

        from pytec3d import Solver

        with Solver(num_threads=4) as solver:
            solver.read_model("model.vtu")
            solver.model.read_json("model_data.json")
            solver.load_model_into_solver()
            solver.time_step["time"] = 1e-3
            solver.solve()
            solver.write_output()

    Attributes:
        time_step: ``time`` (end time of the step), ``theta`` (theta method, first-order formulations), ``beta`` and
            ``gamma`` (Newmark, second-order formulations) and ``transient_effects`` (False for a static solve).
            Read at every :meth:`solve`.
        solver_settings: ``{"linear": {...}, "nonlinear": {...}, "petsc_options": {...}}``, see
            :data:`DEFAULT_SOLVER_SETTINGS`; keys left out keep their default. Copied when the model is loaded:
            changes take effect at the next :meth:`load_model_into_solver`.
    """

    def __init__(self, num_threads: int = 1, root_only_output: bool = True, log_file: PathType | bool = True) -> None:
        """Create a solver with an empty model.

        Args:
            num_threads: Threads per MPI process of the C++ solver. The first solver session of the process sets it
                for the whole process.
            root_only_output: Silence ``sys.stdout`` and ``sys.stderr`` on the MPI ranks other than 0, until
                :meth:`close`.
            log_file: Also write the log of the run (Python and C++ messages and tables) to a file: True for
                ``<mesh_name>.log`` in the working directory, a path, or False for none. The file is started when
                the model is read or loaded, overwriting an older one. The raw output of PETSc and libMesh (e.g.
                ``-ksp_monitor``, the Newton iterations) only goes to the standard output.

        Raises:
            ValueError: If ``num_threads`` is not a positive integer.
            TypeError: If ``log_file`` is not a bool or a path.
        """
        if isinstance(num_threads, bool) or not isinstance(num_threads, int) or num_threads < 1:
            raise ValueError(f"num_threads must be a positive integer, got {num_threads!r}.")
        if not isinstance(log_file, (bool, str, os.PathLike)):
            raise TypeError(f"log_file must be a bool or a path, got {type(log_file).__name__}.")
        self._num_threads = num_threads
        self._log_file_option = log_file
        self._log_file: str | None = None
        self._steps: list[dict[str, Any]] = []
        self._cleared_outputs: set[str] = set()  # models whose previous outputs write_output() has deleted
        self._start = (time.perf_counter(), time.process_time())
        self._counts_at_start = _all_message_counts()
        self._model = Model()
        self.time_step: dict[str, Any] = copy.deepcopy(DEFAULT_TIME_STEP)
        self.solver_settings: dict[str, Any] = copy.deepcopy(DEFAULT_SOLVER_SETTINGS)
        self._mesh_filepath: str | None = None
        self._session: Any | None = None
        self._load_time_arrays: dict[tuple[Location, str], np.ndarray] = {}  # copies taken at the load
        self._saved_stdio: tuple[TextIO, TextIO] | None = None
        self._closed = False
        _register_shutdown()
        if root_only_output:
            self._silence_non_root_output()
        _open_solvers.add(self)

    def __repr__(self) -> str:
        state = "closed" if self._closed else "loaded" if self.is_loaded else "not loaded"
        return f"Solver({self._model!r}, {state})"

    # ------------------------------------------------------------------------------------------- lifecycle

    def __enter__(self) -> Solver:
        return self

    def __exit__(
        self, exc_type: type[BaseException] | None, exc: BaseException | None, traceback: TracebackType | None
    ) -> None:
        self.close()

    def __del__(self) -> None:
        with contextlib.suppress(Exception):
            self.close()

    def close(self) -> None:
        """Close the solver session and restore the standard streams. Idempotent.

        The process-wide C++ runtime stays up for other solvers; it is shut down at interpreter exit.
        """
        if self._closed:
            return
        try:
            if self._steps:
                self._report_run()
        finally:
            self.close_session()
            self._close_log_file()
            self._restore_output()
            self._closed = True
            _open_solvers.discard(self)

    def close_session(self) -> None:
        """Close the C++ solver session, if any, and unlock the model."""
        if self._session is not None:
            LOGGER.info("Closing the solver session of '%s'.", self._model.mesh_name)
            try:
                self._session.close()
            finally:
                self._session = None
        self._model.unlock()

    def _silence_non_root_output(self) -> None:
        rank = world_rank()
        if rank == 0:
            return
        self._saved_stdio = (sys.stdout, sys.stderr)
        null_stream = open(os.devnull, "w", encoding="utf-8")  # noqa: SIM115 - closed in _restore_output
        sys.stdout = sys.stderr = null_stream
        LOGGER.debug("stdout and stderr silenced on rank %d.", rank)

    def _restore_output(self) -> None:
        if self._saved_stdio is None:
            return
        null_stream = sys.stdout
        sys.stdout, sys.stderr = self._saved_stdio
        self._saved_stdio = None
        null_stream.close()

    # ------------------------------------------------------------------------------------------- state

    @property
    def rank(self) -> int:
        """MPI rank of this process in ``COMM_WORLD`` (0 without MPI)."""
        return world_rank()

    @property
    def num_threads(self) -> int:
        """Threads per MPI process requested for the solver session."""
        return self._num_threads

    @property
    def model(self) -> Model:
        """The model. Assigning another one closes the solver session.

        Raises:
            TypeError: When assigning something else than a :class:`~pytec3d.model.Model`.
        """
        return self._model

    @model.setter
    def model(self, value: Model) -> None:
        if not isinstance(value, Model):
            raise TypeError(f"model must be a pytec3d.Model, got {type(value).__name__}.")
        self.close_session()
        self._model = value

    @property
    def is_loaded(self) -> bool:
        """True when the model is loaded into a solver session."""
        return self._session is not None

    @property
    def mesh_filepath(self) -> str | None:
        """Absolute path of the file the model was read from, if any."""
        return self._mesh_filepath

    @property
    def log_file(self) -> str | None:
        """Absolute path of the log file of the run, once started (see ``log_file`` of :class:`Solver`)."""
        return self._log_file

    @property
    def steps(self) -> list[dict[str, Any]]:
        """One record per solved step: ``time``, ``index``, iterations, convergence, ``wall_seconds``, ``output``
        (the file written after it), and ``report``, the full report of the C++ solver."""
        return self._steps

    # ------------------------------------------------------------------------------------------- run log

    def _open_log_file(self) -> None:
        """Start the log file of the run, once: the Python log and the C++ solver log are appended to it."""
        if self._log_file is not None or self._log_file_option is False:
            return
        option = self._log_file_option
        path = os.path.abspath(f"{self._model.mesh_name}.log" if option is True else os.fspath(option))
        if self.rank == 0:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w", encoding="utf-8").close()  # noqa: SIM115 - truncate the log of a previous run
        open_log_file(path)
        self._log_file = path
        try:
            _backend().set_log_file(path)
        except ImportError as exc:
            LOGGER.warning("The C++ solver log will not be in the log file: %s", exc)
        LOGGER.info("Logging to '%s'.", path)

    def _close_log_file(self) -> None:
        if self._log_file is None:
            return
        close_log_file()
        backend = sys.modules.get(BACKEND_MODULE)
        if backend is not None:
            with contextlib.suppress(Exception):
                backend.set_log_file("")
        self._log_file = None

    def _report_run(self) -> None:
        """Print the summary of the run (a row per step) and write it to ``output/<mesh_name>_steps.csv``."""
        wall = time.perf_counter() - self._start[0]
        cpu = time.process_time() - self._start[1]
        warnings, errors = (now - start for now, start in zip(_all_message_counts(), self._counts_at_start))
        csv_path = None
        if self.rank == 0:
            directory = os.path.join(os.getcwd(), OUTPUT_DIRECTORY)
            os.makedirs(directory, exist_ok=True)
            csv_path = os.path.join(directory, f"{self._model.mesh_name}_steps.csv")
            columns = ["index", "time", "nonlinear", "newton_iterations", "linear_iterations", "convergence",
                       "true_residual", "max_law_input", "wall_seconds", "output"]
            with open(csv_path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(self._steps)

        table = Table(f"Run summary: {self._model.mesh_name}")
        table.header(["Step", "Time", "Newton its", "Linear its", "Convergence", "True residual", "Wall (s)", "Output"])
        for step in self._steps:
            table.row([
                step["index"], f"{step['time']:.4e}", step["newton_iterations"] if step["nonlinear"] else "-",
                step["linear_iterations"], step["convergence"],
                "-" if step["true_residual"] is None else f"{step['true_residual']:.2e}",
                f"{step['wall_seconds']:.2f}", os.path.basename(step["output"]) if step["output"] else "-",
            ])
        totals = Table(f"Run totals: {self._model.mesh_name}")
        totals.row(["Steps", len(self._steps)])
        totals.row(["Wall time (s)", f"{wall:.2f}"])
        totals.row(["CPU time, rank 0 (s)", f"{cpu:.2f}"])
        totals.row(["Warnings", warnings])
        totals.row(["Errors", errors])
        totals.row(["Log file", self._log_file or "-"])
        totals.row(["Steps (CSV)", csv_path or "-"])
        print_table(f"{table}\n{totals}")

    # ------------------------------------------------------------------------------------------- model files

    @classmethod
    def from_grid(cls, mesh: vtkUnstructuredGrid, num_threads: int = 1, mesh_name: str | None = None) -> Solver:
        """A solver whose model is built around a mesh (:meth:`pytec3d.Model.from_grid`).

        Raises:
            TypeError: If ``mesh`` is not a ``vtkUnstructuredGrid``.
        """
        solver = cls(num_threads=num_threads)
        solver.model = Model.from_grid(meshes.check_grid(mesh), mesh_name=mesh_name)
        LOGGER.info(meshes.summary(solver.model.mesh, solver.model.mesh_name))
        return solver

    def read_model(self, path: PathType, reset_arrays: bool = True) -> None:
        """Replace the model by the one stored in a ``.vtu`` / ``.pvtu`` file; closes the solver session.

        Args:
            path: File to read, absolute or relative to the working directory.
            reset_arrays: Set the arrays of the element types back to their default values, clear the selection and
                compute the midside-node mask (:meth:`pytec3d.Model.read_model`).

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file is not a valid mesh.
        """
        file_path = os.path.abspath(os.fspath(path))
        LOGGER.info("Reading the model from '%s' (reset_arrays=%s).", file_path, reset_arrays)
        self.model = Model.read_model(file_path, reset_arrays=reset_arrays)
        self._mesh_filepath = file_path
        self._open_log_file()
        LOGGER.info(meshes.summary(self.model.mesh, self.model.mesh_name, file_path))

    def save_model(self, path: PathType | None = None) -> None:
        """Save the model to a ``.vtu`` file, on rank 0 only.

        Args:
            path: File to write; ``<mesh_name>.vtu`` in the working directory by default.
        """
        if self.rank != 0:
            return
        file_path = os.path.abspath(os.fspath(path) if path is not None else f"{self.model.mesh_name}.vtu")
        self.model.save_model(file_path, binary=True)
        LOGGER.info("Model saved to '%s'.", file_path)

    def write_output(self) -> str | None:
        """Save the model as ``output/<mesh_name>_out_<n>.vtu`` in the working directory, on rank 0 only.

        ``n`` is the number of steps this solver has solved (0 before the first :meth:`solve`): writing again
        after the same step overwrites the file. The directory is created if needed, and on the first write the
        outputs of a previous run of the same model (``<mesh_name>_out_*.vtu``, ``<mesh_name>_steps.csv``) are
        deleted, so that the series only holds the steps of this run.

        Returns:
            The path written (None on the other ranks).

        Raises:
            SolverStateError: If the model has no mesh.
        """
        if not self.model.has_mesh:
            raise SolverStateError("The model has no mesh: nothing to write.")
        if self.rank != 0:
            return None
        name = self.model.mesh_name
        directory = os.path.join(os.getcwd(), OUTPUT_DIRECTORY)
        os.makedirs(directory, exist_ok=True)
        if name not in self._cleared_outputs:
            self._clear_previous_outputs(directory, name)
            self._cleared_outputs.add(name)
        path = os.path.join(directory, f"{name}_out_{len(self._steps)}.vtu")
        self.save_model(path)
        if self._steps and self._steps[-1]["output"] is None:
            self._steps[-1]["output"] = path
        return path

    @staticmethod
    def _clear_previous_outputs(directory: str, name: str) -> None:
        """Delete the output meshes and the steps CSV of a previous run of the model ``name``."""
        pattern = re.compile(rf"^{re.escape(name)}_(out_\d+\.vtu|steps\.csv)$")
        removed = [f for f in sorted(os.listdir(directory)) if pattern.match(f)]
        for file_name in removed:
            os.remove(os.path.join(directory, file_name))
        if removed:
            LOGGER.info("Deleted %d output file(s) of a previous run from '%s'.", len(removed), directory)

    # ------------------------------------------------------------------------------------------- solution

    def effective_solver_settings(self) -> dict[str, Any]:
        """:attr:`solver_settings` over the defaults: each group merged key by key, None PETSc options removed."""
        merged = copy.deepcopy(DEFAULT_SOLVER_SETTINGS)
        for group, values in self.solver_settings.items():
            if isinstance(merged.get(group), dict) and isinstance(values, dict):
                merged[group].update(values)
            else:
                merged[group] = values
        if isinstance(merged.get("petsc_options"), dict):
            merged["petsc_options"] = {k: v for k, v in merged["petsc_options"].items() if v is not None}
        return merged

    def validate(self) -> None:
        """Check the solver settings and the model (mesh arrays and model data).

        Raises:
            ValidationError: Listing every problem.
        """
        check_solver_settings(self.effective_solver_settings())
        check_model(self.model)

    def load_model_into_solver(self) -> None:
        """Validate the model and load it into a new solver session (closing the previous one).

        The session copies the model data and the solver settings, and keeps the mesh: until the session is closed
        the model is locked (only the values of its arrays may change).

        Raises:
            SolverStateError: If the solver is closed or the model has no mesh.
            ValidationError: If the model or the settings are invalid.
            RuntimeError: If the C++ solver rejects the model (in a parallel run it aborts instead).
        """
        if self._closed:
            raise SolverStateError("This solver is closed: create a new one.")
        if not self.model.has_mesh:
            raise SolverStateError("The model has no mesh: read one first (read_model).")

        self.close_session()
        self._open_log_file()
        self.model.update_grid()
        LOGGER.info("Validating the model.")
        self.validate()
        script = os.path.abspath(sys.argv[0]) if sys.argv and sys.argv[0] else None
        print_table(provenance_table(
            {"Mesh": self._mesh_filepath, "Model data": self.model.data_file,
             "Script": script if script and os.path.isfile(script) else None},
            _package_version(),
        ))

        settings = self.effective_solver_settings()
        changed = _changed_settings(settings, DEFAULT_SOLVER_SETTINGS)
        LOGGER.info(
            "Creating the solver session (%d thread(s) per rank); solver settings: %s.",
            self._num_threads, "defaults" if not changed else ", ".join(changed),
        )
        start = time.perf_counter()
        session = _backend().Session(self._num_threads)
        try:
            session.load_model(self.model, settings)
        except BaseException:
            session.close()
            raise
        self._session = session
        self._load_time_arrays = load_time_arrays(self.model.mesh)
        self.model.lock()
        LOGGER.info("Model '%s' loaded into the solver session in %.2f s.", self.model.mesh_name,
                    time.perf_counter() - start)

    def solve(self) -> None:
        """Solve the time step :attr:`time_step`; the solution is written into the mesh arrays.

        Raises:
            SolverStateError: If the model is not loaded into a solver session.
            ValidationError: If the time step is invalid or arrays are missing.
            RuntimeError: If the C++ solver fails (in a parallel run it aborts instead).
        """
        if self._session is None:
            raise SolverStateError("No solver session: call load_model_into_solver() first.")
        check_before_solve(self.model, self.time_step, self._load_time_arrays)
        step = self.time_step
        LOGGER.info(
            "Solving t = %g (%s; theta %g, beta %g, gamma %g).", step["time"],
            "transient" if step.get("transient_effects") else "static", step.get("theta", 1.0),
            step.get("beta", 0.25), step.get("gamma", 0.5),
        )
        start = time.perf_counter()
        self._session.solve(step)
        seconds = time.perf_counter() - start
        LOGGER.info("t = %g solved in %.2f s.", step["time"], seconds)
        self._record_step(self._session.last_step_report(), seconds)

    def _record_step(self, report: dict[str, Any], seconds: float) -> None:
        nonlinear = bool(report["nonlinear"])
        ranges = report["law_ranges"]
        self._steps.append({
            "index": report["index"],
            "time": report["time"],
            "nonlinear": nonlinear,
            "newton_iterations": report["newton_iterations"],
            "linear_iterations": report["linear_iterations"],
            "convergence": report["newton_result"] if nonlinear else report["linear_reason"],
            "true_residual": report["true_residual"],
            "max_law_input": "; ".join(
                f"{r['law']} {r['source']} {r['id']}: {r['max_input']:.4g}" for r in ranges
            ),
            "wall_seconds": seconds,
            "output": None,
            "report": report,
        })


def _package_version() -> str:
    try:
        return importlib.metadata.version("pytec3d")
    except importlib.metadata.PackageNotFoundError:
        return "0.0.0 (source tree)"


def _all_message_counts() -> tuple[int, int]:
    """Warnings and errors logged so far by the package and by the C++ solver (if loaded)."""
    warnings, errors = message_counts()
    backend = sys.modules.get(BACKEND_MODULE)
    if backend is not None:
        with contextlib.suppress(Exception):
            solver_warnings, solver_errors = backend.message_counts()
            warnings += solver_warnings
            errors += solver_errors
    return warnings, errors


def _changed_settings(settings: dict[str, Any], defaults: dict[str, Any], prefix: str = "") -> list[str]:
    """``group.key=value`` of every setting that differs from its default, for the log."""
    changed = []
    for key in sorted(set(settings) | set(defaults)):
        value, default = settings.get(key), defaults.get(key)
        if isinstance(value, dict) and isinstance(default, dict):
            changed += _changed_settings(value, default, f"{prefix}{key}.")
        elif value != default:
            changed.append(f"{prefix}{key}={value!r}")
    return changed


# Solvers not closed yet: closed at interpreter exit, before the C++ runtime is shut down.
_open_solvers: weakref.WeakSet[Solver] = weakref.WeakSet()


def _shutdown_at_exit() -> None:
    for solver in list(_open_solvers):
        solver.close()
    backend = sys.modules.get(BACKEND_MODULE)
    if backend is None:
        return
    try:
        backend.shutdown_runtime()
    except Exception:
        LOGGER.exception("The solver runtime failed to shut down.")


@functools.cache
def _register_shutdown() -> None:
    """Register the exit hook, once, when the first solver is created (importing the package starts no MPI)."""
    # MPI must still be up when the runtime shuts down: import mpi4py (which finalizes MPI at exit) before registering
    # the hook, so that the hook runs first.
    with contextlib.suppress(ImportError):
        from mpi4py import MPI  # noqa: F401
    atexit.register(_shutdown_at_exit)
