import os
import sys
import re
import copy
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

from model_edit import ModelEdit
from excel import Excel
from constants import *


# Add the current directory to the Python path to ensure imports work
try:
    import openpyxl
except ImportError:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    local_packages_path = os.path.join(current_dir, "local_packages")
    sys.path.append(local_packages_path)
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    try:
        import openpyxl

        # print("openpyxl imported successfully.")
    except ImportError:
        print(
            "Failed to import openpyxl. Ensure the package is downloaded and the path is correct."
        )


class APEC_M:
    """
    APEC class to manage VTK unstructured grids, metadata, and IO utilities.
    """

    def __init__(self) -> None:
        """
        Initialize the APEC object without requiring a file path.
        """
        self.meshes: Dict[str, vtk.vtkUnstructuredGrid] = {}
        self.mesh_filepaths: Dict[str, str] = {}
        self._active_mesh_name: Optional[str] = None

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
        self._solver_session: Optional[Any] = None

    @property
    def mesh(self) -> vtk.vtkUnstructuredGrid:
        """
        Active mesh accessor.
        """
        if not self._active_mesh_name:
            raise RuntimeError("No active mesh is set.")
        if self._active_mesh_name not in self.meshes:
            raise RuntimeError(
                f"Active mesh '{self._active_mesh_name}' is not available."
            )
        return self.meshes[self._active_mesh_name]

    @mesh.setter
    def mesh(self, value: vtk.vtkUnstructuredGrid) -> None:
        """
        Direct mesh assignment is disallowed.
        """
        raise AttributeError("Direct mesh assignment is not allowed. Use read_mesh().")

    @property
    def active_mesh_name(self) -> Optional[str]:
        """
        Name of the active mesh.
        """
        return self._active_mesh_name

    def set_active_mesh(self, mesh_name: str) -> bool:
        """
        Set the active mesh by name.

        Returns:
            True when the mesh exists and is activated, otherwise False.
        """
        if mesh_name not in self.meshes:
            return False

        self._active_mesh_name = mesh_name
        return True

    def get_mesh(self, mesh_name: str) -> Optional[vtk.vtkUnstructuredGrid]:
        """
        Fetch a mesh by name.
        """
        return self.meshes.get(mesh_name)

    def clear_mesh(self, mesh_name: str) -> None:
        """
        Remove a mesh by name and update the active mesh if needed.
        """
        self.meshes.pop(mesh_name, None)
        self.mesh_filepaths.pop(mesh_name, None)

        if self._active_mesh_name == mesh_name:
            if self.meshes:
                next_mesh = next(iter(self.meshes))
                self._active_mesh_name = next_mesh
            else:
                self._active_mesh_name = None

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
        first_mesh = len(self.meshes) == 0
        loaded_mesh = ModelEdit.read_mesh(vtk_file_path, reset_arrays=reset_arrays)

        mesh_path = os.path.abspath(os.fspath(vtk_file_path))
        resolved_name = os.path.splitext(os.path.basename(mesh_path))[0]
        mesh_folder = os.path.dirname(mesh_path) or os.getcwd()

        if first_mesh:
            os.chdir(mesh_folder)

        model_folder = os.path.join(os.getcwd(), resolved_name)
        os.makedirs(model_folder, exist_ok=True)
        os.makedirs(os.path.join(model_folder, "emat"), exist_ok=True)
        os.makedirs(os.path.join(model_folder, "output"), exist_ok=True)

        self.meshes[resolved_name] = loaded_mesh
        self.mesh_filepaths[resolved_name] = mesh_path
        self.set_active_mesh(resolved_name)

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
            ValueError: If output_file is not provided and no active mesh is set.
        """
        if output_file is None:
            if not self._active_mesh_name:
                raise ValueError(
                    "No active mesh is set. Activate a mesh before saving."
                )

        # Define the output file name
        if output_file is None:
            output_file = os.path.join(os.getcwd(), f"{self._active_mesh_name}.vtu")

        # Write the unstructured grid to the .vtu file
        ModelEdit.save_mesh(self.mesh, output_file, binary=True)

        print(f"Mesh saved to '{output_file}'.")

    def write_output(self) -> None:
        """
        Save the active mesh into its output folder using the next timestep index.

        The target filename is:
            <mesh_name>/output/<mesh_name>_out_<time_step>.vtu

        The time step is computed at runtime by scanning existing files in the
        output directory and taking max(existing_step) + 1. If no output files
        are present, the first generated time step is 1.
        """
        if not self._active_mesh_name:
            raise ValueError(
                "No active mesh is set. Activate a mesh before writing output."
            )

        mesh_name = self._active_mesh_name
        output_dir = os.path.join(os.getcwd(), mesh_name, "output")
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
        RenameSource(self.active_mesh_name or "APEC_Mesh", producer)
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
        try:
            import apec_m_solver as solver_module
        except ImportError as exc:
            raise ImportError(
                "Could not import the C++ solver module 'apec_m_solver'. "
                "Build/install the pybind11 extension before launching the solver."
            ) from exc
        return solver_module

    def load_mesh_into_solver(
        self,
        num_processors: int = 1,
        mesh_name: Optional[str] = None,
    ) -> Any:
        """
        Build a C++ solver session from the active VTK grid and current APEC attributes.

        Args:
            num_processors:
                Number of MPI ranks to use for the solver communicator (default: 1).
            mesh_name:
                Optional mesh to activate before launching.
        """
        solver_module = self._require_solver_module()
        mpi_comm = self._build_solver_mpi_communicator(num_processors)

        if mesh_name is not None:
            if not self.set_active_mesh(mesh_name):
                raise ValueError(f"Mesh '{mesh_name}' does not exist.")

        solver = solver_module.ApecMSolver()
        solver.load_model(self, mpi_comm)

        self._solver_session = solver
        return solver

    @staticmethod
    def _build_solver_mpi_communicator(num_processors: int) -> Optional[Any]:
        """
        Create an mpi4py communicator for the solver.

        Notes:
            - If mpi4py is unavailable, only num_processors=1 is supported.
            - If num_processors > COMM_WORLD size, an error is raised.
            - When COMM_WORLD has more ranks than requested, only the first
              num_processors ranks are included (others get COMM_NULL and must
              not invoke the solver).
        """
        if num_processors < 1:
            raise ValueError("num_processors must be >= 1.")

        try:
            from mpi4py import MPI
        except ImportError:
            if num_processors != 1:
                raise ImportError(
                    "mpi4py is required when num_processors > 1."
                ) from None
            return None

        world = MPI.COMM_WORLD
        world_size = world.Get_size()
        world_rank = world.Get_rank()

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
