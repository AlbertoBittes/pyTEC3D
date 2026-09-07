# pyright: reportMissingImports=false, reportPossiblyUnboundVariable=false, reportIncompatibleMethodOverride=false, reportArgumentType=false
"""ParaView integration entrypoint for APEC-M.

This module serves two use-cases:
1) Console bootstrap: import it in ParaView Python console to make apec_m importable.
2) Plugin file: load this file in ParaView's plugin manager to register filters.
"""

from __future__ import annotations

import os
import sys
from numbers import Real
from typing import Any

# When loaded directly as a plugin file, ensure the top-level python/ root is importable.
_THIS_FILE = os.path.abspath(__file__)
_PACKAGE_DIR = os.path.dirname(_THIS_FILE)
_PACKAGE_PARENT = os.path.dirname(_PACKAGE_DIR)
if _PACKAGE_PARENT not in sys.path:
    sys.path.insert(0, _PACKAGE_PARENT)

from apec_m_solver.apec_m import APEC_M  # type: ignore[import-not-found]
from apec_m.logging_utils import (  # type: ignore[import-not-found]
    get_logger,
    log_paraview_console_package_root,
    log_paraview_init,
    log_paraview_working_directory,
)
from apec_m.selector import Selector  # type: ignore[import-not-found]

import numpy as np
import vtk
import vtk.util.numpy_support as vtknp  # type: ignore[import-not-found]

LOGGER = get_logger("apec_m.paraview")


def setup_paraview_console() -> str:
    """Ensure the package parent is importable in console/plugin contexts."""
    this_file = os.path.abspath(__file__)
    package_dir = os.path.dirname(this_file)
    package_parent = os.path.dirname(package_dir)
    if package_parent not in sys.path:
        sys.path.insert(0, package_parent)
    log_paraview_console_package_root(LOGGER, package_parent)
    return package_parent


def start_apec_for_paraview(num_threads: int = 1) -> APEC_M:
    """Create and return an APEC_M instance for ParaView console workflows."""
    setup_paraview_console()
    log_paraview_init(LOGGER)
    apec = APEC_M(num_threads=num_threads)
    log_paraview_working_directory(LOGGER, os.getcwd())
    return apec


def plot(mesh: vtk.vtkUnstructuredGrid) -> None:
    """Plot an unstructured grid in ParaView. Usage: plot(apec.model.mesh)."""
    if not isinstance(mesh, vtk.vtkUnstructuredGrid):
        raise TypeError("mesh must be a vtk.vtkUnstructuredGrid.")

    from paraview.simple import Render, RenameSource, Show, TrivialProducer  # type: ignore[import-not-found]

    producer = TrivialProducer()
    producer.GetClientSideObject().SetOutput(mesh)
    RenameSource("APEC_Mesh", producer)
    Show(producer)
    Render()


vtkDataArraySelection: Any
VTKPythonAlgorithmBase: Any
smdomain: Any
smproperty: Any
smproxy: Any

try:
    from vtkmodules.vtkCommonCore import vtkDataArraySelection as _vtkDataArraySelection
    from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase as _VTKPythonAlgorithmBase
    from paraview.util.vtkAlgorithm import smdomain as _smdomain, smproperty as _smproperty, smproxy as _smproxy  # type: ignore[import-not-found]
    vtkDataArraySelection = _vtkDataArraySelection
    VTKPythonAlgorithmBase = _VTKPythonAlgorithmBase
    smdomain = _smdomain
    smproperty = _smproperty
    smproxy = _smproxy
    _PARAVIEW_PLUGIN_AVAILABLE = True
except Exception:  # pragma: no cover
    _PARAVIEW_PLUGIN_AVAILABLE = False

    class _NoOpDecoratorNamespace:
        @staticmethod
        def _identity(*args, **kwargs):
            def _wrap(obj):
                return obj

            return _wrap

        filter = _identity
        input = _identity
        datatype = _identity
        xml = _identity
        intvector = _identity
        stringvector = _identity
        dataarrayselection = _identity
        doublevector = _identity
        intrange = _identity
        doublerange = _identity

    class _FallbackVTKPythonAlgorithmBase:
        def __init__(self, *args, **kwargs):
            pass

        def Modified(self):
            pass

        @staticmethod
        def INPUT_REQUIRED_DATA_TYPE():
            return None

        @staticmethod
        def INPUT_IS_OPTIONAL():
            return None

    class _FallbackDataArraySelection:
        def AddObserver(self, *args, **kwargs):
            return None

        def AddArray(self, *args, **kwargs):
            return None

        def Modified(self):
            return None

        def ArrayIsEnabled(self, *args, **kwargs):
            return True

    VTKPythonAlgorithmBase = _FallbackVTKPythonAlgorithmBase
    vtkDataArraySelection = _FallbackDataArraySelection
    smdomain = _NoOpDecoratorNamespace()
    smproperty = _NoOpDecoratorNamespace()
    smproxy = _NoOpDecoratorNamespace()


def _create_modified_callback(anobject):
    """
     create modified callback.
    
    Args:
        anobject (Any): Description.
    
    Returns:
        Any: Description.
    """
    import weakref

    weakref_obj = weakref.ref(anobject)
    anobject = None

    def _markmodified(*args, **kwargs):
        o = weakref_obj()
        if o is not None:
            o.Modified()

    return _markmodified


@smproxy.filter(name="SelectExtractedData", label="Select Extracted Data")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
@smproperty.input(name="Extracted Selection", port_index=1)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class SelectExtractedData(VTKPythonAlgorithmBase):
    """SelectExtractedData class."""
    def __init__(self):
        """
        Special method `__init__`.
        
        Returns:
            Any: Description.
        """
        self.__mode = "s"
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1)

    def FillInputPortInformation(self, port, info):
        """
        Fillinputportinformation.
        
        Args:
            port (Any): Description.
            info (Any): Description.
        
        Returns:
            Any: Description.
        """
        if port == 0:
            info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        elif port == 1:
            info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet")
            info.Set(self.INPUT_IS_OPTIONAL(), 1)
        return 1

    def RequestData(self, request, inInfo, outInfo):
        """
        Requestdata.
        
        Args:
            request (Any): Description.
            inInfo (Any): Description.
            outInfo (Any): Description.
        
        Returns:
            Any: Description.
        
        Raises:
            ValueError: If an error condition is met.
        """
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        if input_data is None:
            raise ValueError("Input grid (port 0) is required.")

        output = vtk.vtkUnstructuredGrid.GetData(outInfo)
        output.ShallowCopy(input_data)
        Selector.ensure_selection_arrays(output)

        extracted = vtk.vtkDataSet.GetData(inInfo[1])
        if extracted is None:
            LOGGER.warning("SelectExtractedData: no dataset connected on port 1; selection unchanged.")
            return 1

        point_mask = np.zeros(output.GetNumberOfPoints(), dtype=bool)
        point_ids_array = extracted.GetPointData().GetArray("vtkOriginalPointIds")
        if point_ids_array is not None:
            point_ids = vtknp.vtk_to_numpy(point_ids_array).astype(np.int64, copy=False)
            point_ids = point_ids[(point_ids >= 0) & (point_ids < point_mask.shape[0])]
            point_mask[point_ids] = True
        Selector.apply_boolean_selection_mask(point_mask, output.GetPointData(), "node", self.__mode)

        cell_mask = np.zeros(output.GetNumberOfCells(), dtype=bool)
        cell_ids_array = extracted.GetCellData().GetArray("vtkOriginalCellIds")
        if cell_ids_array is not None:
            cell_ids = vtknp.vtk_to_numpy(cell_ids_array).astype(np.int64, copy=False)
            cell_ids = cell_ids[(cell_ids >= 0) & (cell_ids < cell_mask.shape[0])]
            cell_mask[cell_ids] = True
        Selector.apply_boolean_selection_mask(cell_mask, output.GetCellData(), "elem", self.__mode)
        return 1

    @smproperty.xml("""
        <IntVectorProperty name="Mode"
                            command="SetMode"
                            number_of_elements="1"
                            default_values="0">
          <EnumerationDomain name="enum">
            <Entry value="0" text="Set"/>
            <Entry value="1" text="Add"/>
            <Entry value="2" text="Reselect"/>
            <Entry value="3" text="Unselect"/>
          </EnumerationDomain>
        </IntVectorProperty>
    """)
    def SetMode(self, mode):
        """
        Setmode.
        
        Args:
            mode (Any): Description.
        
        Returns:
            Any: Description.
        """
        self.__mode = Selector.mode_label(int(mode))
        self.Modified()


@smproxy.filter(name="SelectByCriteria", label="Select")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class Select(VTKPythonAlgorithmBase):
    """Select class."""
    def __init__(self):
        """
        Special method `__init__`.
        
        Returns:
            Any: Description.
        """
        self.__array_name = ""
        self.__type = "node"
        self.__mode = "s"
        self.__component = 0
        self.__values_text = ""
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1)

    def FillInputPortInformation(self, port, info):
        """
        Fillinputportinformation.
        
        Args:
            port (Any): Description.
            info (Any): Description.
        
        Returns:
            Any: Description.
        """
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    @staticmethod
    def _parse_values(text: str):
        """
         parse values.
        
        Args:
            text (str): Description.
        
        Returns:
            Any: Description.
        
        Raises:
            ValueError: If an error condition is met.
        """
        text = text.strip()
        if not text:
            raise ValueError("The 'Values' field is empty.")
        if ":" in text:
            parts = [p.strip() for p in text.split(":")]
            if len(parts) != 2:
                raise ValueError("A range must be written as 'min:max'.")
            return (float(parts[0]), float(parts[1]))
        if "," in text:
            parts = [p.strip() for p in text.split(",") if p.strip()]
            if not parts:
                raise ValueError("No values found in the comma-separated list.")
            return [float(p) for p in parts]
        return float(text)

    def RequestData(self, request, inInfo, outInfo):
        """
        Requestdata.
        
        Args:
            request (Any): Description.
            inInfo (Any): Description.
            outInfo (Any): Description.
        
        Returns:
            Any: Description.
        
        Raises:
            ValueError: If an error condition is met.
        """
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        if input_data is None:
            raise ValueError("Input grid (port 0) is required.")

        output = vtk.vtkUnstructuredGrid.GetData(outInfo)
        output.ShallowCopy(input_data)
        Selector.ensure_selection_arrays(output)

        if not self.__array_name:
            LOGGER.warning("Select: no array chosen; selection masks left unchanged.")
            return 1

        value = self._parse_values(self.__values_text)
        component = "magnitude" if self.__component < 0 else self.__component
        Selector.select(
            output,
            self.__array_name,
            value,
            component=component,
            type=self.__type,
            mode=self.__mode,
        )
        return 1

    @smproperty.xml("""
        <StringVectorProperty name="ArrayName"
                               command="SetArrayName"
                               number_of_elements="1"
                               animateable="0">
          <ArrayListDomain name="array_list"
                            attribute_type="Scalars"
                            input_domain_name="input_array">
            <RequiredProperties>
              <Property name="Input" function="Input"/>
            </RequiredProperties>
          </ArrayListDomain>
        </StringVectorProperty>
    """)
    def SetArrayName(self, name):
        """
        Setarrayname.
        
        Args:
            name (Any): Description.
        
        Returns:
            Any: Description.
        """
        self.__array_name = name
        self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="Type"
                            command="SetSelectionType"
                            number_of_elements="1"
                            default_values="0">
          <EnumerationDomain name="enum">
            <Entry value="0" text="Node"/>
            <Entry value="1" text="Element"/>
          </EnumerationDomain>
        </IntVectorProperty>
    """)
    def SetSelectionType(self, value):
        """
        Setselectiontype.
        
        Args:
            value (Any): Description.
        
        Returns:
            Any: Description.
        """
        self.__type = Selector.type_label(int(value))
        self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="Mode"
                            command="SetMode"
                            number_of_elements="1"
                            default_values="0">
          <EnumerationDomain name="enum">
            <Entry value="0" text="Set"/>
            <Entry value="1" text="Add"/>
            <Entry value="2" text="Reselect"/>
            <Entry value="3" text="Unselect"/>
          </EnumerationDomain>
        </IntVectorProperty>
    """)
    def SetMode(self, value):
        """
        Setmode.
        
        Args:
            value (Any): Description.
        
        Returns:
            Any: Description.
        """
        self.__mode = Selector.mode_label(int(value))
        self.Modified()

    @smproperty.intvector(name="Component", default_values=0)
    @smdomain.intrange(min=-1, max=2)
    def SetComponent(self, component):
        """
        Setcomponent.
        
        Args:
            component (Any): Description.
        
        Returns:
            Any: Description.
        
        Raises:
            TypeError: If an error condition is met.
            ValueError: If an error condition is met.
        """
        if not isinstance(component, int):
            raise TypeError("component must be an integer.")
        if component < -1 or component > 2:
            raise ValueError("component must be in [-1, 2] (-1 = magnitude).")
        self.__component = component
        self.Modified()

    @smproperty.stringvector(name="Values", default_values="0")
    def SetValues(self, text):
        """
        Setvalues.
        
        Args:
            text (Any): Description.
        
        Returns:
            Any: Description.
        """
        self.__values_text = text
        self.Modified()


@smproxy.filter(name="SetValueOnSelection", label="Set Value On Selection")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class SetValueOnSelection(VTKPythonAlgorithmBase):
    """SetValueOnSelection class."""
    def __init__(self):
        """
        Special method `__init__`.
        
        Returns:
            Any: Description.
        """
        self.__component = 0
        self.__value = 0.0
        self.__array_names = []
        self.__arrayselection = vtkDataArraySelection()
        self.__arrayselection.AddObserver("ModifiedEvent", _create_modified_callback(self))
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1)

    def FillInputPortInformation(self, port, info):
        """
        Fillinputportinformation.
        
        Args:
            port (Any): Description.
            info (Any): Description.
        
        Returns:
            Any: Description.
        """
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    def RequestInformation(self, request, inInfo, outInfo):  # pyright: ignore[reportIncompatibleMethodOverride]
        """
        Requestinformation.
        
        Args:
            request (Any): Description.
            inInfo (Any): Description.
            outInfo (Any): Description.
        
        Returns:
            Any: Description.
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
        self.__array_names = [a for a in (point_data_arrays + cell_data_arrays) if a]

        for array_name in self.__array_names:
            self.__arrayselection.AddArray(array_name)

        self.__arrayselection.Modified()
        return 1

    def RequestData(self, request, inInfo, outInfo):
        """
        Requestdata.
        
        Args:
            request (Any): Description.
            inInfo (Any): Description.
            outInfo (Any): Description.
        
        Returns:
            Any: Description.
        
        Raises:
            ValueError: If an error condition is met.
        """
        input_data = vtk.vtkUnstructuredGrid.GetData(inInfo[0])
        if input_data is None:
            raise ValueError("Input grid (port 0) is required.")

        output = vtk.vtkUnstructuredGrid.GetData(outInfo)
        output.ShallowCopy(input_data)
        Selector.ensure_selection_arrays(output)

        self.__arrayselection.Modified()
        for array_name in self.__array_names:
            if not self.__arrayselection.ArrayIsEnabled(array_name):
                continue

            if output.GetPointData().GetArray(array_name) is not None:
                Selector.set_value_on_selection(
                    output,
                    array_name,
                    self.__value,
                    component=self.__component,
                    type="node",
                )

            if output.GetCellData().GetArray(array_name) is not None:
                Selector.set_value_on_selection(
                    output,
                    array_name,
                    self.__value,
                    component=self.__component,
                    type="elem",
                )

        return 1

    def _get_array_selection(self):
        """
         get array selection.
        
        Returns:
            Any: Description.
        """
        return self.__arrayselection

    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        """
        Getdataarrayselection.
        
        Returns:
            Any: Description.
        """
        return self._get_array_selection()

    @smproperty.intvector(name="Component", default_values=0)
    @smdomain.intrange(min=0, max=2)
    def SetComponent(self, component):
        """
        Setcomponent.
        
        Args:
            component (Any): Description.
        
        Returns:
            Any: Description.
        
        Raises:
            TypeError: If an error condition is met.
            ValueError: If an error condition is met.
        """
        if not isinstance(component, int):
            raise TypeError("component must be an integer.")
        if component < 0 or component > 2:
            raise ValueError("component must be in [0, 2].")
        self.__component = component
        self.Modified()

    @smproperty.doublevector(name="Value", default_values=0.0)
    @smdomain.doublerange()
    def SetValue(self, value):
        """
        Setvalue.
        
        Args:
            value (Any): Description.
        
        Returns:
            Any: Description.
        
        Raises:
            TypeError: If an error condition is met.
            ValueError: If an error condition is met.
        """
        if not isinstance(value, Real):
            raise TypeError("value must be a real scalar.")
        if not np.isfinite(float(value)):
            raise ValueError("value must be finite.")
        self.__value = float(value)
        self.Modified()


if __name__ == "__main__":
    setup_paraview_console()
    apec = start_apec_for_paraview()
    LOGGER.info("apec_m_paraview console bootstrap complete.")
