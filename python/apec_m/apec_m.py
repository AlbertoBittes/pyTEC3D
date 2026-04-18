import os
import sys
import re
import copy
import importlib
import shutil
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

from .model_edit import ModelEdit
from .excel import Excel
from .constants import *
from .apec_algorithm import ApecAlgorithm


# Prefer environment openpyxl; fall back to vendored local_packages when needed.
try:
    import openpyxl
except ImportError:
    if os.path.isdir(LOCAL_PACKAGES_PATH) and LOCAL_PACKAGES_PATH not in sys.path:
        sys.path.insert(0, LOCAL_PACKAGES_PATH)
    try:
        import openpyxl
    except ImportError:
        print(
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
        """
        self._mesh: Optional[vtk.vtkUnstructuredGrid] = None
        self._mesh_name: Optional[str] = None
        self.mesh_filepath: Optional[str] = None
        self._solver_session: Optional[Any] = None
        if mpi_communicator is not None:
            self._solver_mpi_communicator: Optional[Any] = mpi_communicator
        else:
            self._solver_mpi_communicator = self._build_solver_mpi_communicator(
                num_processors=num_processors
            )
        self._python_output_root_only: bool = False
        self._saved_stdout = None
        self._saved_stderr = None
        self._null_output_stream = None

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
        self.set_python_output_root_only(True)

    @staticmethod
    def _world_rank() -> int:
        try:
            from mpi4py import MPI
        except ImportError:
            return 0
        return MPI.COMM_WORLD.Get_rank()

    def set_python_output_root_only(self, enabled: bool = True) -> None:
        """
        Enable/disable Python output only on MPI rank 0.

        When enabled, non-root ranks redirect both stdout and stderr to os.devnull.
        """
        enabled = bool(enabled)
        if enabled == self._python_output_root_only:
            return

        rank = self._world_rank()
        if enabled and rank != 0:
            if self._saved_stdout is None:
                self._saved_stdout = sys.stdout
            if self._saved_stderr is None:
                self._saved_stderr = sys.stderr
            if self._null_output_stream is None or self._null_output_stream.closed:
                self._null_output_stream = open(os.devnull, "w")
            sys.stdout = self._null_output_stream
            sys.stderr = self._null_output_stream
        else:
            if self._saved_stdout is not None:
                sys.stdout = self._saved_stdout
            if self._saved_stderr is not None:
                sys.stderr = self._saved_stderr

        self._python_output_root_only = enabled

    @property
    def solver_mpi_communicator(self) -> Optional[Any]:
        """
        Solver communicator used when creating C++ solver sessions.
        """
        return self._solver_mpi_communicator

    def __del__(self) -> None:
        try:
            self.set_python_output_root_only(False)
        except Exception:
            pass
        try:
            if self._null_output_stream is not None and not self._null_output_stream.closed:
                self._null_output_stream.close()
        except Exception:
            pass
    @property
    def mesh(self) -> vtk.vtkUnstructuredGrid:
        """
        Mesh accessor.
        """
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")
        return self._mesh

    @mesh.setter
    def mesh(self, value: vtk.vtkUnstructuredGrid) -> None:
        """
        Direct mesh assignment is disallowed.
        """
        raise AttributeError("Direct mesh assignment is not allowed. Use read_mesh().")

    @property
    def solver(self) -> Optional[Any]:
        """
        Solver session accessor.
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
        """
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
        self._solver_session = None

    def save_mesh(
        self, output_file: Optional[Union[str, PathLike[str]]] = None
    ) -> None:
        """
        Saves the unstructured grid in self.mesh to a .vtu file.
        Before writing, attaches the element_types DataFrame as field data to the grid.
        Validates that the first column of element_types has unique values, all values have a correspondence
        in the 'ansys_elem_num' array, and there are no undefined element types in the grid.

        Returns:
            None

        Raises:
            ValueError: If output_file is not provided and no mesh is loaded.
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

        print(f"Mesh saved to '{output_file}'.")

    def write_output(self) -> None:
        """
        Save the mesh into the output folder using the next timestep index.

        The target filename is:
            output/<mesh_name>_out_<time_step>.vtu

        The time step is computed at runtime by scanning existing files in the
        output directory and taking max(existing_step) + 1. If no output files
        are present, the first generated time step is 1.
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
        # Check if the Excel file exists
        if not os.path.exists(excel_file_path):
            raise FileNotFoundError(f"Excel file '{excel_file_path}' not found.")

        # Loop through the data_to_read labels and process each one
        for label in data_to_read:
            try:
                self.read_excel_sheet(excel_file_path, label)
            except ValueError as e:
                print(f"Warning: {e}")

    def read_excel_sheet(
        self, excel_file: Union[str, PathLike[str]], sheet_name: str
    ) -> None:
        """
        Processes an Excel sheet by determining the appropriate method to call
        (read_excel_df, read_excel_dict, or update_vtk_grid_from_excel) based on a global dictionary.

        Args:
            excel_file: Path to the Excel file.
            sheet_name: Name of the sheet to process.

        Raises:
            ValueError: If the sheet name is not recognized or mapped method is missing.
        """
        if sheet_name not in SHEET_PROCESSING_MAP:
            raise ValueError(
                f"Invalid sheet name '{sheet_name}'. Must be one of {list(SHEET_PROCESSING_MAP.keys())}."
            )
        processing_info = SHEET_PROCESSING_MAP[sheet_name]
        print(f"Reading sheet '{sheet_name}' from '{excel_file}'...")
        if isinstance(processing_info, str):
            func = getattr(Excel, processing_info)
            result = func(excel_file, sheet_name)
            if result is not None:
                setattr(self, sheet_name, result)
        elif isinstance(processing_info, tuple):
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
        """
        from paraview.simple import TrivialProducer, RenameSource, Show, Render

        # Assuming you have a function to replot the mesh in your application
        # This is a placeholder for the actual replotting logic
        print("Plotting the mesh...")
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
    def _require_solver_module() -> Any:
        package_name = __package__ or "apec_m"
        try:
            solver_module = importlib.import_module(f"{package_name}.apec_m_solver")
        except ImportError as exc:
            raise ImportError(
                "Could not import the C++ solver module 'apec_m.apec_m_solver'. "
                "Build/install the pybind11 extension before launching the solver."
            ) from exc
        return solver_module

    def load_mesh_into_solver(self) -> Any:
        """
        Build a C++ solver session from the loaded VTK grid and current APEC attributes.

        """
        self._solver_session = self._load_mesh_into_solver()

    def _load_mesh_into_solver(self) -> Any:
        """
        Build a C++ solver session from the loaded VTK grid and current APEC attributes.
        """
        solver_module = self._require_solver_module()
        mpi_comm = self._solver_mpi_communicator

        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        solver = solver_module.ApecMSolver(mpi_comm)
        solver.load_model(self)
        return solver

    def solve(self, nodal_solution: bool = True) -> None:
        """
        Solve using the loaded mesh solver session.

        Args:
            nodal_solution:
                Whether nodal solution output should be computed.
        """
        if self._mesh is None:
            raise RuntimeError("No mesh is loaded.")

        if self._solver_session is None:
            raise RuntimeError(
                "No solver session is loaded for the current mesh. "
                "Call load_mesh_into_solver() first."
            )

        self._solver_session.solve(nodal_solution)

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
