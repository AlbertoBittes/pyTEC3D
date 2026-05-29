import os
import sys
import re
import copy
import gc
import importlib
import shutil
import atexit
from os import PathLike
from typing import (
    Any,
    Dict,
    Optional,
    Sequence,
    Union,
)

import pandas as pd
import vtk

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
LOCAL_PACKAGES_PATH = os.path.join(CURRENT_DIR, "local_packages")

from .logging_utils import get_logger
from .model_edit import ModelEdit
from .excel import Excel
from .constants import *
from .apec_algorithm import ApecAlgorithm


LOGGER = get_logger("apec_m.APEC_M")


# Prefer environment openpyxl; fall back to vendored local_packages when needed.
try:
    import openpyxl
except ImportError:
    if os.path.isdir(LOCAL_PACKAGES_PATH) and LOCAL_PACKAGES_PATH not in sys.path:
        sys.path.insert(0, LOCAL_PACKAGES_PATH)
    try:
        import openpyxl
    except ImportError:
        LOGGER.warning(
            "Failed to import openpyxl. Ensure it is installed or vendored under local_packages."
        )


class APEC_M:
    """
    APEC class to manage VTK unstructured grids, metadata, and IO utilities.
    """

    def __init__(
        self,
        mpi_communicator: Optional[Any] = None,
        num_processors: Optional[int] = None,
    ) -> None:
        """
        Initialize the APEC object without requiring a file path.

        Args:
            mpi_communicator:
                Optional mpi4py communicator to use in solver sessions.
                If None, a communicator is built from COMM_WORLD according to
                num_processors.
            num_processors:
                Number of ranks to include in the internally-built solver
                communicator when mpi_communicator is None. If None, use all
                ranks in COMM_WORLD.

        Returns:
            None
        """
        self._mesh: Optional[vtk.vtkUnstructuredGrid] = None
        self._mesh_name: Optional[str] = None
        self.mesh_filepath: Optional[str] = None
        self._solver_session: Optional[Any] = None
        if mpi_communicator is not None:
            self._solver_mpi_communicator: Optional[Any] = mpi_communicator
            self._owns_solver_mpi_communicator: bool = False
        else:
            self._solver_mpi_communicator = self._build_solver_mpi_communicator(
                num_processors=num_processors
            )
            self._owns_solver_mpi_communicator = True
        self._null_output_stream = None
        self._closed: bool = False

        self.element_types: pd.DataFrame = pd.DataFrame(
            columns=list(ELEMENT_TYPES_COLUMNS_TYPES.keys())
        )
        self.materials: pd.DataFrame = pd.DataFrame(
            columns=list(MATERIALS_COLUMNS_TYPES.keys())
        )
        self.constants: pd.DataFrame = pd.DataFrame(
            columns=list(CONSTANTS_COLUMNS_TYPES.keys())
        )
        self.time_step: Dict[str, Any] = copy.deepcopy(DEFAULT_TIME_PARAMETERS)
        self.solver_settings: Dict[str, Any] = copy.deepcopy(DEFAULT_SOLVER_SETTINGS)
        self.dof_couplings: pd.DataFrame = pd.DataFrame(
            columns=list(DOF_COUPLINGS_COLUMNS_TYPES.keys())
        )
        self.algorithm = ApecAlgorithm(self)
        self._configure_python_output_root_only()
        atexit.register(self._cleanup_at_exit)

    @staticmethod
    def _world_rank() -> int:
        """Return MPI rank in COMM_WORLD when available.

        Returns:
            int: MPI rank, or 0 when mpi4py is unavailable.
        """
        try:
            from mpi4py import MPI
        except ImportError:
            return 0
        return MPI.COMM_WORLD.Get_rank()

    def _configure_python_output_root_only(self) -> None:
        """Enforce Python stdout/stderr emission only on MPI rank 0.

        Non-root ranks are permanently redirected to os.devnull for this process.

        Returns:
            None
        """
        rank = self._world_rank()
        if rank == 0:
            return

        if self._null_output_stream is None or self._null_output_stream.closed:
            self._null_output_stream = open(os.devnull, "w", encoding="utf-8")
        sys.stdout = self._null_output_stream
        sys.stderr = self._null_output_stream
        LOGGER.debug("Python stdout/stderr redirected to os.devnull for rank %d.", rank)

    @property
    def solver_mpi_communicator(self) -> Optional[Any]:
        """
        Solver communicator used when creating C++ solver sessions.

        Returns:
            Optional[Any]: MPI communicator used for solver creation.
        """
        return self._solver_mpi_communicator

    def __del__(self) -> None:
        """Best-effort destructor that releases runtime resources.

        Returns:
            None
        """
        try:
            self.close()
        except Exception:
            pass
        try:
            if self._null_output_stream is not None and not self._null_output_stream.closed:
                self._null_output_stream.close()
        except Exception:
            pass

    def _cleanup_at_exit(self) -> None:
        """Atexit hook that closes runtime resources.

        Returns:
            None
        """
        self.close()

    def _shutdown_solver_runtime(self) -> None:
        """Request C++ module runtime shutdown when available.

        Returns:
            None
        """
        try:
            solver_module = self._require_solver_module()
        except Exception:
            return

        shutdown_runtime = getattr(solver_module, "shutdown_runtime", None)
        if shutdown_runtime is not None:
            try:
                shutdown_runtime()
            except Exception:
                pass

    def _close_solver_session(self) -> None:
        """Close and clear the active solver session.

        Returns:
            None
        """
        if self._solver_session is None:
            return

        try:
            close_method = getattr(self._solver_session, "close", None)
            if close_method is not None:
                close_method()
        except Exception:
            pass
        finally:
            self._solver_session = None
            # Ensure pybind-wrapped C++ objects are finalized before runtime shutdown.
            gc.collect()

    def _release_runtime_references(self) -> None:
        """Release Python-held references that may keep C++ runtime objects alive.

        Returns:
            None
        """
        try:
            self.algorithm.equilibrium_mesh = None
            self.algorithm._b_target_axis_cache = None
        except Exception:
            pass

        self._mesh = None

    def _free_owned_solver_communicator(self) -> None:
        """Free internally-owned MPI communicator when applicable.

        Returns:
            None
        """
        if not self._owns_solver_mpi_communicator:
            return

        comm = self._solver_mpi_communicator
        if comm is None:
            return

        try:
            from mpi4py import MPI
        except ImportError:
            self._solver_mpi_communicator = None
            return

        try:
            if comm not in (MPI.COMM_NULL, MPI.COMM_WORLD, MPI.COMM_SELF):
                comm.Free()
        except Exception:
            pass
        finally:
            self._solver_mpi_communicator = None

    def close(self) -> None:
        """Close solver session and runtime resources.

        Returns:
            None
        """
        if self._closed:
            return

        LOGGER.info("Closing APEC_M runtime resources.")

        self._close_solver_session()
        self._release_runtime_references()
        self._shutdown_solver_runtime()
        self._free_owned_solver_communicator()

        self._closed = True

        LOGGER.info("APEC_M runtime resources closed.")

    @property
    def mesh(self) -> vtk.vtkUnstructuredGrid:
        """
        Mesh accessor.

        Returns:
            vtk.vtkUnstructuredGrid: Loaded mesh.
        """
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")
        return self._mesh

    @mesh.setter
    def mesh(self, value: vtk.vtkUnstructuredGrid) -> None:
        """
        Direct mesh assignment is disallowed.

        Args:
            value: Candidate mesh value.

        Returns:
            None
        """
        raise AttributeError("Direct mesh assignment is not allowed. Use read_mesh().")

    @property
    def solver(self) -> Optional[Any]:
        """
        Solver session accessor.

        Returns:
            Optional[Any]: Active solver session object.
        """
        return self._solver_session

    #
    #
    #
    #
    #
    #  Basic methods: read, initialize and save model
    #
    #
    #
    #
    #
    def read_mesh(
        self,
        vtk_file_path: Union[str, PathLike[str]],
        reset_arrays: bool = True,
    ) -> None:
        """
        Reads a VTK unstructured grid from a file in XML format.

        Args:
            vtk_file_path: Path to the VTK file in XML format.
            reset_arrays: When True, reset initialized arrays to defaults after loading.

        Returns:
            None
        """
        if not isinstance(reset_arrays, bool):
            raise TypeError("reset_arrays must be a bool.")

        LOGGER.info(
            "Loading mesh from '%s' (reset_arrays=%s).",
            os.fspath(vtk_file_path),
            reset_arrays,
        )

        loaded_mesh = ModelEdit.read_mesh(vtk_file_path, reset_arrays=reset_arrays)

        mesh_path = os.path.abspath(os.fspath(vtk_file_path))
        resolved_name = os.path.splitext(os.path.basename(mesh_path))[0]
        mesh_folder = os.path.dirname(mesh_path) or os.getcwd()
        os.chdir(mesh_folder)

        emat_dir = os.path.join(os.getcwd(), "emat")
        output_dir = os.path.join(os.getcwd(), "output")
        if self._world_rank() == 0:
            os.makedirs(emat_dir, exist_ok=True)
            os.makedirs(output_dir, exist_ok=True)

            # Always start from a clean emat folder.
            for entry_name in os.listdir(emat_dir):
                entry_path = os.path.join(emat_dir, entry_name)
                if os.path.isdir(entry_path):
                    shutil.rmtree(entry_path)
                elif os.path.exists(entry_path):
                    os.unlink(entry_path)

            # Always start from a clean output folder (files only).
            for entry_name in os.listdir(output_dir):
                entry_path = os.path.join(output_dir, entry_name)
                if os.path.isfile(entry_path):
                    os.unlink(entry_path)

        self._mesh = loaded_mesh
        self._mesh_name = resolved_name
        self.mesh_filepath = mesh_path
        self._close_solver_session()
        LOGGER.info(APEC_M._format_mesh_summary(loaded_mesh, resolved_name, mesh_path, os.getcwd()))

    def save_mesh(
        self, output_file: Optional[Union[str, PathLike[str]]] = None
    ) -> None:
        """
        Saves the unstructured grid in self.mesh to a .vtu file.
        Before writing, attaches the element_types DataFrame as field data to the grid.
        Validates that the first column of element_types has unique values, all values have a correspondence
        in the 'ansys_elem_num' array, and there are no undefined element types in the grid.

        Args:
            output_file: Optional output path. Defaults to '<mesh_name>.vtu' in current cwd.

        Raises:
            ValueError: If output_file is not provided and no mesh is loaded.

        Returns:
            None
        """
        if self._world_rank() != 0:
            return

        if output_file is None:
            if self._mesh_name is None:
                raise ValueError("No mesh is loaded. Load a mesh before saving.")

        # Define the output file name
        if output_file is None:
            output_file = os.path.join(os.getcwd(), f"{self._mesh_name}.vtu")

        # Write the unstructured grid to the .vtu file
        ModelEdit.save_mesh(self.mesh, output_file, binary=True)

        LOGGER.info("Mesh saved to '%s'.", output_file)

    def write_output(self) -> None:
        """
        Save the mesh into the output folder using the next timestep index.

        The target filename is:
            output/<mesh_name>_out_<time_step>.vtu

        The time step is computed at runtime by scanning existing files in the
        output directory and taking max(existing_step) + 1. If no output files
        are present, the first generated time step is 1.

        Returns:
            None
        """
        if self._mesh_name is None:
            raise ValueError("No mesh is loaded. Load a mesh before writing output.")

        if self._world_rank() != 0:
            return

        mesh_name = self._mesh_name
        output_dir = os.path.join(os.getcwd(), "output")
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
        self.save_mesh(output_file=output_file)

    #
    #
    #
    #
    #
    #  Excel I/O methods
    #
    #
    #
    #
    #
    def read_excel(
        self,
        excel_file_path: Union[str, PathLike[str]],
        data_to_read: Sequence[str] = tuple(DEFAULT_COLUMNS_ATTRIBUTES_TYPES.keys()),
    ) -> None:
        """
        Reads data from an Excel file and processes it based on the provided labels.

        Args:
            excel_file_path: Path to the Excel file.
            data_to_read: Iterable of sheet labels to read and process.

        Returns:
            None

        Raises:
            FileNotFoundError: If the Excel file is missing.
        """
        if not isinstance(data_to_read, Sequence):
            raise TypeError("data_to_read must be a sequence of sheet names.")

        LOGGER.info("Reading Excel workbook '%s'.", excel_file_path)

        # Loop through the data_to_read labels and process each one
        for label in data_to_read:
            try:
                self.read_excel_sheet(excel_file_path, label)
            except ValueError as e:
                LOGGER.warning("%s", e)

    def read_excel_sheet(
        self, excel_file: Union[str, PathLike[str]], sheet_name: str
    ) -> None:
        """
        Processes an Excel sheet by determining the appropriate method to call
        (read_excel_df, read_excel_dict, or update_vtk_grid_from_excel) based on a global dictionary.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Name of the sheet to process.

        Returns:
            None

        Raises:
            ValueError: If the sheet name is not recognized or mapped method is missing.
        """
        if sheet_name not in SHEET_PROCESSING_MAP:
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(SHEET_PROCESSING_MAP.keys())}."
            )
        processing_info = SHEET_PROCESSING_MAP[sheet_name]
        LOGGER.info("Reading sheet '%s' from '%s'.", sheet_name, excel_file)
        if isinstance(processing_info, str):
            func = getattr(Excel, processing_info)
            result = func(excel_file, sheet_name)
            if result is not None:
                setattr(self, sheet_name, result)
        elif isinstance(processing_info, tuple):
            if self._mesh is None:
                raise RuntimeError(
                    "No mesh is loaded. Load a mesh before reading mesh-dependent Excel sheets."
                )
            method_name, data_type = processing_info
            func = getattr(Excel, method_name)
            func(self.mesh, excel_file, sheet_name, data_type)

    #
    #
    #
    #
    #
    #  Method to plot the model in ParaView
    #
    #
    #
    #
    #
    def plot(self) -> None:
        """
        Replot the mesh in the active view in Paraview.

        Returns:
            None
        """
        from paraview.simple import TrivialProducer, RenameSource, Show, Render  # type: ignore[import-not-found]

        # Assuming you have a function to replot the mesh in your application
        # This is a placeholder for the actual replotting logic
        LOGGER.info("Plotting mesh in ParaView.")
        producer = TrivialProducer()
        producer.GetClientSideObject().SetOutput(self.mesh)
        RenameSource(self._mesh_name or "APEC_Mesh", producer)
        Show(producer)
        Render()

    #
    #
    #
    #
    #
    #  C++ solver bridge (pybind11 module)
    #
    #
    #
    #
    #
    @staticmethod
    def _format_mesh_summary(
        mesh: vtk.vtkUnstructuredGrid,
        mesh_name: str,
        mesh_path: str,
        cwd: str,
    ) -> str:
        """Build a formatted multi-line mesh summary for logging.

        Args:
            mesh: Loaded VTK unstructured grid.
            mesh_name: Job / mesh name (stem of the file).
            mesh_path: Absolute path to the mesh file.
            cwd: Working directory after mesh load.

        Returns:
            str: Multi-line summary string ready to pass to a logger.
        """
        W = 70
        sep = " " + "-" * W

        n_pts = mesh.GetNumberOfPoints()
        n_cells = mesh.GetNumberOfCells()
        bounds = mesh.GetBounds()   # xmin xmax ymin ymax zmin zmax
        z_span = bounds[5] - bounds[4]
        spatial_dim = "2D" if z_span < 1e-10 else "3D"

        def _bstr(lo: float, hi: float) -> str:
            return f"[{lo:.6g}, {hi:.6g}]"

        def _array_lines(data_obj: Any, indent: str = "    ") -> list[str]:
            rows: list[str] = []
            for i in range(data_obj.GetNumberOfArrays()):
                arr = data_obj.GetArray(i)
                if arr is None:
                    continue
                name = arr.GetName() or "(unnamed)"
                nc = arr.GetNumberOfComponents()
                dtype = arr.GetDataTypeAsString()
                comp_label = f"{nc} comp" + ("s" if nc != 1 else " ")
                rows.append(f"{indent}{name:<28}{comp_label:<10}{dtype}")
            return rows

        def _field_lines(fd: Any, indent: str = "    ") -> list[str]:
            rows: list[str] = []
            for i in range(fd.GetNumberOfArrays()):
                arr = fd.GetAbstractArray(i)
                if arr is None:
                    continue
                name = arr.GetName() or "(unnamed)"
                nc = arr.GetNumberOfComponents()
                try:
                    dtype = arr.GetDataTypeAsString()
                except Exception:
                    dtype = "unknown"
                comp_label = f"{nc} comp" + ("s" if nc != 1 else " ")
                rows.append(f"{indent}{name:<28}{comp_label:<10}{dtype}")
            return rows

        pd_obj = mesh.GetPointData()
        cd_obj = mesh.GetCellData()
        fd_obj = mesh.GetFieldData()
        n_pd = pd_obj.GetNumberOfArrays()
        n_cd = cd_obj.GetNumberOfArrays()
        n_fd = fd_obj.GetNumberOfArrays()

        L: list[str] = []
        L.append("")
        L.append(sep)
        L.append(f"  {'APEC-M  Mesh Summary':^{W}}")
        L.append(sep)
        L.append(f"  {'Job name':<18}: {mesh_name}")
        L.append(f"  {'Working dir':<18}: {cwd}")
        L.append(f"  {'Mesh file':<18}: {mesh_path}")
        L.append(sep)
        L.append("  Geometry")
        L.append(f"    {'Points':<20}: {n_pts:,}")
        L.append(f"    {'Cells':<20}: {n_cells:,}")
        L.append(f"    {'Spatial dimension':<20}: {spatial_dim}")
        L.append(f"    {'Bounds X':<20}: {_bstr(bounds[0], bounds[1])}")
        L.append(f"    {'Bounds Y':<20}: {_bstr(bounds[2], bounds[3])}")
        L.append(f"    {'Bounds Z':<20}: {_bstr(bounds[4], bounds[5])}")
        L.append(sep)
        L.append(f"  Point arrays ({n_pd})")
        L.extend(_array_lines(pd_obj))
        L.append(sep)
        L.append(f"  Cell arrays ({n_cd})")
        L.extend(_array_lines(cd_obj))
        L.append(sep)
        L.append(f"  Field arrays ({n_fd})")
        L.extend(_field_lines(fd_obj))
        L.append(sep)
        return "\n".join(L)

    @staticmethod
    def _require_solver_module() -> Any:
        """Import and return the pybind solver module.

        Returns:
            Any: Imported solver module.
        """
        package_name = __package__ or "apec_m"
        try:
            solver_module = importlib.import_module(f"{package_name}.apec_m_solver")
        except ImportError as exc:
            raise ImportError(
                "Could not import the C++ solver module 'apec_m.apec_m_solver'. "
                "Build/install the pybind11 extension before launching the solver."
            ) from exc
        return solver_module

    def validate(self) -> None:
        """Run a comprehensive consistency check on the loaded APEC model.

        Validates the ``element_types``, ``dof_couplings``, ``materials``, and
        ``constants`` DataFrames against each other and against the loaded VTK
        mesh.  All violations are reported together in a single ``ValueError``.
        Unused identifiers (TYPE_IDs, MATERIAL_IDs, CONSTANTS_IDs not referenced
        by any cell) are reported as warnings rather than errors.

        Returns:
            None

        Raises:
            ValueError: If the mesh is not loaded or any consistency check fails.
        """
        from .validation_utils import validate_apec_model
        validate_apec_model(self)

    def load_mesh_into_solver(self) -> None:
        """
        Build a C++ solver session from the loaded VTK grid and current APEC attributes.

        Runs a full model consistency check via :meth:`validate` before creating the
        session.  The created session is accessible via the `solver` property after
        this call.

        Returns:
            None
        """
        if self._closed:
            raise RuntimeError("This APEC_M instance is closed and cannot create a new solver session.")

        LOGGER.info("Validating model before creating solver session.")
        self.validate()

        LOGGER.info("Creating solver session from currently loaded mesh.")
        self._solver_session = self._load_mesh_into_solver()
        LOGGER.info("Solver session created.")

    def _load_mesh_into_solver(self) -> Any:
        """
        Build a C++ solver session from the loaded VTK grid and current APEC attributes.

        Returns:
            Any: Newly created solver session.
        """
        solver_module = self._require_solver_module()
        mpi_comm = self._solver_mpi_communicator

        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        self._close_solver_session()

        solver = solver_module.ApecMSolver(mpi_comm)
        solver.load_model(self)
        return solver

    def solve(self) -> None:
        """
        Solve using the loaded mesh solver session.

        Returns:
            None
        """
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        if self._solver_session is None:
            raise RuntimeError(
                "No solver session is loaded for the current mesh. "
                "Call load_mesh_into_solver() first."
            )

        if not hasattr(self._solver_session, "solve"):
            raise RuntimeError("Loaded solver session does not expose a solve() method.")

        LOGGER.info("Starting solver session solve().")
        self._solver_session.solve()
        LOGGER.info("Solver session solve() completed.")

    @staticmethod
    def _build_solver_mpi_communicator(num_processors: Optional[int]) -> Optional[Any]:
        """
        Create an mpi4py communicator for the solver.

        Notes:
            - If mpi4py is unavailable, only num_processors=None or 1 is supported.
            - If num_processors > COMM_WORLD size, an error is raised.
            - When COMM_WORLD has more ranks than requested, only the first
              num_processors ranks are included (others get COMM_NULL and must
              not invoke the solver).

        Args:
            num_processors: Requested number of participating ranks.

        Returns:
            Optional[Any]: Duplicated/split communicator, or None when mpi4py is unavailable.
        """
        try:
            from mpi4py import MPI
        except ImportError:
            if num_processors not in (None, 1):
                raise ImportError(
                    "mpi4py is required when num_processors > 1."
                ) from None
            return None

        world = MPI.COMM_WORLD
        world_size = world.Get_size()
        world_rank = world.Get_rank()

        if num_processors is None:
            num_processors = world_size
        elif isinstance(num_processors, bool):
            raise TypeError("num_processors must be an integer >= 1 or None.")
        elif num_processors < 1:
            raise ValueError("num_processors must be >= 1.")

        if num_processors > world_size:
            raise ValueError(
                f"Requested num_processors={num_processors}, but COMM_WORLD has "
                f"size={world_size}. Launch Python with at least {num_processors} "
                "MPI ranks."
            )

        if world_size == num_processors:
            return world.Dup()

        color = 0 if world_rank < num_processors else MPI.UNDEFINED
        comm = world.Split(color, world_rank)
        if comm == MPI.COMM_NULL:
            raise RuntimeError(
                f"Rank {world_rank} is outside the solver communicator (first "
                f"{num_processors} ranks). Invoke the solver only on participating ranks."
            )
        return comm


# Only execute this block if the script is run directly
if __name__ == "__main__":
    # Example usage
    apec = APEC_M()
    # print("APEC instance created.")
