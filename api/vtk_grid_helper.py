import vtk
import vtk.util.numpy_support as vtknp
import numpy as np
from constants import VTK_ARRAYS_TO_INIT

class VTKGridHelper:
    @staticmethod
    def initialize_vtk_arrays(mesh, arrays_to_init=None):
        """
        Initializes required VTK arrays in cell/point data as specified in arrays_to_init.

        Parameters:
            mesh: The VTK unstructured grid to initialize arrays for.
            arrays_to_init (dict, optional): Dictionary specifying arrays to initialize.
                If None, uses VTK_ARRAYS_TO_INIT from constants.
        """
        if arrays_to_init is None:
            from constants import VTK_ARRAYS_TO_INIT
            arrays_to_init = VTK_ARRAYS_TO_INIT

        cell_data = mesh.GetCellData()
        point_data = mesh.GetPointData()

        for name, info in arrays_to_init.items():
            data_obj = cell_data if info["location"] == "cell" else point_data
            if data_obj.GetArray(name) is None:
                num = (
                    mesh.GetNumberOfCells()
                    if info["location"] == "cell"
                    else mesh.GetNumberOfPoints()
                )
                shape = (num, info["components"]) if info["components"] > 1 else num
                arr = np.full(shape, info.get("fill", 0))
                vtk_arr = vtknp.numpy_to_vtk(arr, deep=True, array_type=info.get("dtype", vtk.VTK_DOUBLE))
                vtk_arr.SetName(name)
                if info["components"] > 1:
                    vtk_arr.SetNumberOfComponents(info["components"])
                else:
                    vtk_arr.SetNumberOfComponents(1)
                data_obj.AddArray(vtk_arr)

    @staticmethod
    def initialize_vtk_array(mesh, name, location="point", components=1, fill=0, dtype=vtk.VTK_DOUBLE):
        """
        Initializes a single VTK array in cell or point data.

        Parameters:
            mesh: The VTK unstructured grid to initialize the array for.
            name (str): Name of the array.
            location (str): "cell" or "point" (default "point").
            components (int): Number of components (default 1).
            fill (numeric): Default value to fill the array (default 0).
            dtype: VTK data type (default vtk.VTK_DOUBLE).
        """
        arrays_to_init = {
            name: {
                "location": location,
                "components": components,
                "fill": fill,
                "dtype": dtype,
            }
        }
        VTKGridHelper.initialize_vtk_arrays(mesh, arrays_to_init)

    @staticmethod
    def safe_remove_array(mesh, array_name, location=None):
        """
        Safely removes an array from CellData or PointData of the mesh if it exists.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh from which to remove the array.
            array_name (str): The name of the array to remove.
            location (str, optional): "cell", "point", or None. If None, tries to guess location.

        Side Effects:
            Removes the array from the detected or specified data location if present.
        """
        # Try to guess location if not provided
        if location is None:
            cell_data = mesh.GetCellData()
            point_data = mesh.GetPointData()
            if cell_data.HasArray(array_name):
                cell_data.RemoveArray(array_name)
                return
            elif point_data.HasArray(array_name):
                point_data.RemoveArray(array_name)
                return
            else:
                print(f"Array '{array_name}' not found in CellData or PointData.")
                return
        elif location == "cell":
            data = mesh.GetCellData()
        elif location == "point":
            data = mesh.GetPointData()
        else:
            raise ValueError(f"Unknown location '{location}'. Use 'cell', 'point', or None.")
        if data.HasArray(array_name):
            data.RemoveArray(array_name)


    @staticmethod
    def validate_dataframe_with_vtk_data(vtk_data, df, array_name):
        """
        Validates that the first column of a DataFrame has unique values,
        all values have a correspondence in the specified VTK data array,
        and there are no undefined values in the VTK data array. Useful
        to check 1-to-1 correspondence between IDs.

        Parameters:
            vtk_data (vtk.vtkDataSetAttributes): The VTK cell or point data to validate against.
            df (pd.DataFrame): The DataFrame to validate.
            array_name (str): The name of the VTK data array to check.

        Raises:
            ValueError: If the DataFrame is empty, the first column is not unique,
                        the VTK data array is missing, or the checks fail.
        """
        if df.empty:
            print("The DataFrame is empty.")
            return
        first_column = df.iloc[:, 0]
        if not first_column.is_unique:
            raise ValueError(
                "The first column of the DataFrame must have unique values."
            )
        vtk_array = vtk_data.GetArray(array_name)
        if vtk_array is None:
            raise ValueError(f"The '{array_name}' array is missing in the VTK data.")
        vtk_array_np = vtknp.vtk_to_numpy(vtk_array)
        if not set(first_column).issubset(set(vtk_array_np)):
            raise ValueError(
                f"Some values in the first column of the DataFrame do not have a correspondence in the '{array_name}' array."
            )
        undefined_values = set(vtk_array_np) - set(first_column)
        if undefined_values:
            raise ValueError(
                f"The VTK data contains undefined values: {undefined_values}."
            )
        print(f"Validation passed for '{array_name}' array and DataFrame.")

    @staticmethod
    def update_vtk_grid(data, df):
        """
        Updates the VTK grid's cell/point data arrays using the DataFrame resulting from `replace_ansys_with_vtk_ids`.
        The function assigns values to the appropriate components at the corresponding VTK IDs.

        Args:
            data: VTK DataSetAttributes (CellData or PointData)
            df (pd.DataFrame): The DataFrame containing columns 'vtk_id', 'label', 'component', and 'value'.
        """

        loads_bc_labels = df["label"].unique()
        for load_bc_label in loads_bc_labels:
            load_bc_df = df[df["label"] == load_bc_label]
            vtk_array = data.GetArray(load_bc_label)
            if vtk_array is None:
                raise ValueError(f"Array '{load_bc_label}' not found in the VTK file.")
            vtk_array_np = vtknp.vtk_to_numpy(vtk_array)
            if vtk_array_np.ndim == 1:
                vtk_array_np[:] = np.nan
            else:
                vtk_array_np[:] = np.full(vtk_array_np.shape, np.nan)
            if vtk_array_np.ndim == 1:
                for component in load_bc_df["component"].unique():
                    component_mask = load_bc_df["component"] == component
                    vtk_ids = load_bc_df.loc[component_mask, "vtk_id"].to_numpy()
                    values = load_bc_df.loc[component_mask, "value"].to_numpy()
                    vtk_array_np[vtk_ids] = values
            else:
                for component in load_bc_df["component"].unique():
                    component_mask = load_bc_df["component"] == component
                    vtk_ids = load_bc_df.loc[component_mask, "vtk_id"].to_numpy()
                    values = load_bc_df.loc[component_mask, "value"].to_numpy()
                    vtk_array_np[vtk_ids, component] = values
            vtk_array.Modified()

    @staticmethod
    def compute_cell_centroids(mesh):
        """
        Computes the geometric centroids for all cells in the given VTK unstructured grid.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh whose cell centroids are to be computed.

        Returns:
            np.ndarray: Array of shape (num_cells, 3) containing the centroid coordinates (x, y, z) for each cell.

        Notes:
            Uses vtkCellCenters to compute centroids. The output array is ordered by cell index.
        """
        cellCenters = vtk.vtkCellCenters()
        cellCenters.SetInputData(mesh)
        cellCenters.Update()
        vtk_points = cellCenters.GetOutput().GetPoints().GetData()
        np_centers = vtknp.vtk_to_numpy(vtk_points)
        return np_centers

    @staticmethod
    def store_cell_centroids(mesh, array_name="cell_centroids"):
        """
        Computes cell centroids and stores them in a VTK array in CellData.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh whose cell centroids are to be computed and stored.
            array_name (str): Name of the array to store centroids in CellData. Default is "cell_centroids".

        The resulting array will have 3 components per cell, representing the centroid coordinates (x, y, z).
        The array is added to mesh.GetCellData().
        """
        np_centers = VTKGridHelper.compute_cell_centroids(mesh)
        vtk_centers = vtknp.numpy_to_vtk(np_centers, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_centers.SetName(array_name)
        vtk_centers.SetNumberOfComponents(3)
        mesh.GetCellData().AddArray(vtk_centers)

    @staticmethod
    def compute_cell_volumes(mesh):
        """
        Computes the volume (for 3D cells) or area (for 2D cells) for each cell in the mesh.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh whose cell volumes/areas are to be computed.

        Returns:
            np.ndarray: Array of shape (num_cells,) containing the volume (or area) for each cell.

        Notes:
            Uses ComputeVolume or ComputeArea if available on the cell object. If not available, sets value to 0.0.
            The output array is ordered by cell index.
        """
        num_cells = mesh.GetNumberOfCells()
        volumes = np.zeros(num_cells)
        for i in range(num_cells):
            cell = mesh.GetCell(i)
            if hasattr(cell, "ComputeVolume"):
                try:
                    volumes[i] = cell.ComputeVolume()
                except Exception:
                    volumes[i] = 0.0
            elif hasattr(cell, "ComputeArea"):
                try:
                    volumes[i] = cell.ComputeArea()
                except Exception:
                    volumes[i] = 0.0
            else:
                volumes[i] = 0.0
        return volumes

    @staticmethod
    def store_cell_volumes(mesh, array_name="cell_volume"):
        """
        Computes cell volumes/areas and stores them in a VTK array in CellData.

        Parameters:
            mesh (vtk.vtkUnstructuredGrid): The mesh whose cell volumes/areas are to be computed and stored.
            array_name (str): Name of the array to store volumes/areas in CellData. Default is "cell_volume".

        The resulting array will have one component per cell, representing the volume or area.
        The array is added to mesh.GetCellData().
        """
        volumes = VTKGridHelper.compute_cell_volumes(mesh)
        vtk_vol_arr = vtknp.numpy_to_vtk(volumes, deep=True, array_type=vtk.VTK_DOUBLE)
        vtk_vol_arr.SetName(array_name)
        vtk_vol_arr.SetNumberOfComponents(1)
        mesh.GetCellData().AddArray(vtk_vol_arr)
