import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkDataArraySelection
from paraview.util.vtkAlgorithm import smproxy, smproperty, smdomain
from numbers import Real

# same imports as earlier.
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase

from .logging_utils import get_logger


LOGGER = get_logger("apec_m.paraview_plugin")


# ------------------------------------------------------------------------------
# A reader example.
# ------------------------------------------------------------------------------
def _create_modified_callback(anobject):
    """Return a callback that marks the owning VTK algorithm as modified.

    Args:
        anobject: VTK object with a Modified() method.

    Returns:
        callable: Callback function.
    """
    import weakref

    weakref_obj = weakref.ref(anobject)
    anobject = None

    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()

    return _markmodified


@smproxy.filter(name="SetValueOnSelection", label="Set Value On Selection")
@smproperty.input(name="Input Grid", port_index=0)
@smproperty.input(name="Extracted Selection", port_index=1)
class SetValueOnSelection(VTKPythonAlgorithmBase):
    def __init__(self):
        """Initialize plugin state and array selection controls.

        Returns:
            None
        """
        self.__component = 0
        self.__value = 0.0
        self.__array_names = []
        self.__arrayselection = vtkDataArraySelection()
        self.__arrayselection.AddObserver("ModifiedEvent", _create_modified_callback(self))
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1)

    def RequestInformation(self, request, inInfo, outInfo):
        """Populate array-selection metadata from the input dataset.

        Args:
            request: VTK pipeline request.
            inInfo: Input info vectors.
            outInfo: Output info vectors.

        Returns:
            int: 1 on success, 0 when required input is missing.
        """
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        if input_data is None:
            LOGGER.error("SetValueOnSelection received no input grid on port 0.")
            return 0

        point_data_arrays = [
            input_data.GetPointData().GetArrayName(i)
            for i in range(input_data.GetPointData().GetNumberOfArrays())
        ]
        cell_data_arrays = [
            input_data.GetCellData().GetArrayName(i)
            for i in range(input_data.GetCellData().GetNumberOfArrays())
        ]
        self.__array_names = point_data_arrays + cell_data_arrays

        for array_name in self.__array_names:
            self.__arrayselection.AddArray(array_name)
            # self.__arrayselection.EnableArray(array_name)  # Enable the array if needed

        self.__arrayselection.Modified()
        LOGGER.debug("SetValueOnSelection discovered %d candidate arrays.", len(self.__array_names))

        return 1

    def RequestData(self, request, inInfo, outInfo):
        """Apply value updates to selected point/cell arrays and forward output.

        Args:
            request: VTK pipeline request.
            inInfo: Input info vectors.
            outInfo: Output info vectors.

        Returns:
            int: 1 on success.
        """
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        selection_data = vtk.vtkUnstructuredGrid.GetData(inInfo[1])
        if input_data is None:
            raise ValueError("Input grid (port 0) is required.")
        if selection_data is None:
            raise ValueError("Extracted selection grid (port 1) is required.")

        self.__arrayselection.Modified()
        for i in range(self.__arrayselection.GetNumberOfArrays()):
            array_name = self.__arrayselection.GetArrayName(i)
            if not self.__arrayselection.ArrayIsEnabled(array_name):
                continue

            # Try to get the array from point data first
            vtk_data_array = input_data.GetPointData().GetArray(array_name)
            if vtk_data_array is not None:
                data_array = vtk_to_numpy(vtk_data_array)
                ids_array = selection_data.GetPointData().GetArray("vtkOriginalPointIds")
                if ids_array is None:
                    raise ValueError(
                        "Selection is missing 'vtkOriginalPointIds' for point-data update."
                    )
                original_ids = vtk_to_numpy(ids_array)
                is_point_data = True
            else:
                # If not found in point data, try to get it from cell data
                vtk_data_array = input_data.GetCellData().GetArray(array_name)
                if vtk_data_array is not None:
                    data_array = vtk_to_numpy(vtk_data_array)
                    ids_array = selection_data.GetCellData().GetArray("vtkOriginalCellIds")
                    if ids_array is None:
                        raise ValueError(
                            "Selection is missing 'vtkOriginalCellIds' for cell-data update."
                        )
                    original_ids = vtk_to_numpy(ids_array)
                    is_point_data = False
                else:
                    raise ValueError(
                        f"The data array '{array_name}' is not available in either point data or cell data."
                    )

            if original_ids is None or data_array is None:
                raise ValueError(
                    f"The selection or data array '{array_name}' is not available."
                )

            original_ids = np.asarray(original_ids, dtype=np.int64).reshape(-1)
            if original_ids.size == 0:
                continue
            if np.any(original_ids < 0) or np.any(original_ids >= data_array.shape[0]):
                raise ValueError(
                    f"Selection IDs for '{array_name}' are out of bounds for target array length {data_array.shape[0]}."
                )

            # Check if the array is a vector
            if vtk_data_array.GetNumberOfComponents() > 1:
                data_array[original_ids, self.__component] = self.__value
            else:
                data_array[original_ids] = self.__value

            vtk_array = numpy_to_vtk(data_array)
            vtk_array.SetName(array_name)
            vtk_array.Modified()
            if is_point_data:
                input_data.GetPointData().AddArray(vtk_array)
            else:
                input_data.GetCellData().AddArray(vtk_array)

        input_data.Modified()

        return 1

    def _get_array_selection(self):
        """Return the internal VTK data-array selection object.

        Returns:
            vtkDataArraySelection: Array-selection object.
        """
        return self.__arrayselection

    # Array selection API is typical with readers in VTK
    # This is intended to allow ability for users to choose which arrays to
    # load. To expose that in ParaView, simply use the
    # smproperty.dataarrayselection().
    # This method **must** return a `vtkDataArraySelection` instance.
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        """Expose array-selection API to ParaView properties panel.

        Returns:
            vtkDataArraySelection: Array-selection object.
        """
        return self._get_array_selection()

    @smproperty.intvector(name="Component: ", default_values=0)
    @smdomain.intrange(min=0, max=2)
    def SetComponent(self, component):
        """Set vector component index used for updates on multi-component arrays.

        Args:
            component: Component index in [0, 2].

        Returns:
            None
        """
        if not isinstance(component, int):
            raise TypeError("component must be an integer.")
        if component < 0 or component > 2:
            raise ValueError("component must be in [0, 2].")
        self.__component = component
        self.Modified()

    @smproperty.doublevector(name="Value: ", default_values=0.0)
    @smdomain.doublerange()
    def SetValue(self, value):
        """Set scalar value written to the selected IDs/components.

        Args:
            value: Scalar numeric value.

        Returns:
            None
        """
        if not isinstance(value, Real):
            raise TypeError("value must be a real scalar.")
        if not np.isfinite(float(value)):
            raise ValueError("value must be finite.")
        self.__value = value
        self.Modified()
