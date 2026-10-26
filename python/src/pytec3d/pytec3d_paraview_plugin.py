"""ParaView integration of pytec3d.

Two uses, which can be combined:

* Console: *Run Script* of this file in the ParaView Python shell. The shell then has pytec3d and its tools, with no
  import to type: ``selection``, ``meshes`` (:mod:`pytec3d.mesh`), ``model_edit``, ``Model``, the coordinate systems,
  ``np``, and :func:`plot` to show a mesh (:data:`CONSOLE_NAMES`).
* Plugin: load this file in ParaView (Tools > Manage Plugins > Load New) to get the filters *Select*,
  *Select Extracted Data* and *Set Value On Selection*, which edit the selection and the arrays of a model mesh as
  :mod:`pytec3d.selection` does. Loaded this way it also puts the names of the console in the Python shell, when the
  shell already exists.
"""

# No ``from __future__`` import: *Run Script* runs this file after a first line of its own (``__file__ = ...``).

import os
import sys
import weakref
from typing import Any

# Loaded as a plugin file, the package directory is not on the path yet. Installed, the element library
# (pytec3d_element_library) is next to pytec3d; in the repository (python/src/pytec3d), it is element_library/.
_PACKAGE_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_REPOSITORY_ELEMENT_LIBRARY = os.path.join(_PACKAGE_PARENT, "..", "..", "element_library")
for _directory in (_PACKAGE_PARENT, os.path.normpath(_REPOSITORY_ELEMENT_LIBRARY)):
    if os.path.isdir(_directory) and _directory not in sys.path:
        sys.path.insert(0, _directory)

import numpy as np
from vtkmodules.util import numpy_support as vtknp
from vtkmodules.vtkCommonDataModel import vtkDataSet, vtkUnstructuredGrid

# The console names (CONSOLE_NAMES): with *Run Script* every top-level name of this file is a name of the shell.
import pytec3d  # noqa: F401
from pytec3d import (  # noqa: F401
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    Model,
    SphericalCoordinateSystem,
    XAxisCylindricalCoordinateSystem,
    YAxisCylindricalCoordinateSystem,
    ZAxisCylindricalCoordinateSystem,
    element_library,
    model_edit,
    selection,
    tree_gauge,
)
from pytec3d import mesh as meshes  # noqa: F401
from pytec3d._logging import get_logger
from pytec3d._typing import Location, SelectionMode

LOGGER = get_logger("pytec3d.paraview")

#: Selection modes and locations as indexed by the enumerations of the filter properties.
_MODES: tuple[SelectionMode, ...] = ("s", "a", "r", "u")
_LOCATIONS: tuple[Location, ...] = ("node", "elem")

try:
    from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
    from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
    from vtkmodules.vtkCommonCore import vtkDataArraySelection

    PARAVIEW_AVAILABLE = True
except ImportError:  # Outside ParaView: the decorators do nothing and the filters cannot run.
    PARAVIEW_AVAILABLE = False

    class _Decorators:
        def __getattr__(self, _name: str) -> Any:
            return lambda *args, **kwargs: (lambda obj: obj)

    smdomain = smproperty = smproxy = _Decorators()

    class VTKPythonAlgorithmBase:  # type: ignore[no-redef]
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            pass

        def Modified(self) -> None:
            pass

    class vtkDataArraySelection:  # type: ignore[no-redef]  # noqa: N801
        def AddObserver(self, *args: Any) -> None:
            pass


_MODE_XML = """
    <IntVectorProperty name="Mode" command="SetMode" number_of_elements="1" default_values="0">
      <EnumerationDomain name="enum">
        <Entry value="0" text="Set"/>
        <Entry value="1" text="Add"/>
        <Entry value="2" text="Reselect"/>
        <Entry value="3" text="Unselect"/>
      </EnumerationDomain>
    </IntVectorProperty>
"""


#: The names this file puts in the ParaView Python shell.
CONSOLE_NAMES: tuple[str, ...] = (
    "np",
    "pytec3d",
    "Model",
    "selection",
    "meshes",
    "model_edit",
    "element_library",
    "tree_gauge",
    "CartesianCoordinateSystem",
    "CylindricalCoordinateSystem",
    "SphericalCoordinateSystem",
    "XAxisCylindricalCoordinateSystem",
    "YAxisCylindricalCoordinateSystem",
    "ZAxisCylindricalCoordinateSystem",
    "plot",
)


def setup_paraview_console(namespace: dict[str, Any] | None = None) -> dict[str, Any] | None:
    """Put the names of :data:`CONSOLE_NAMES` in the ParaView Python shell.

    Args:
        namespace: The namespace of the shell. By default the one of ParaView's Python shell (``__vtkConsoleLocals``
            of ``__main__``), when the shell exists.

    Returns:
        The namespace, or None when there is no shell (outside ParaView, or a plugin loaded before the shell was
        created: then *Run Script* of this file gives the names).
    """
    if namespace is None:
        import __main__

        namespace = getattr(__main__, "__vtkConsoleLocals", None)
        if namespace is None:
            return None
    names = globals()
    namespace.update({name: names[name] for name in CONSOLE_NAMES})
    LOGGER.info(
        "pytec3d %s in the ParaView Python shell (from %s): %s.",
        pytec3d.__version__,
        _PACKAGE_PARENT,
        ", ".join(CONSOLE_NAMES),
    )
    return namespace


def plot(mesh: vtkUnstructuredGrid, name: str = "Mesh") -> None:
    """Show a mesh in the ParaView view, e.g. ``plot(solver.model.mesh)``.

    Raises:
        TypeError: If ``mesh`` is not a ``vtkUnstructuredGrid``.
        ImportError: Outside ParaView.
    """
    if not isinstance(mesh, vtkUnstructuredGrid):
        raise TypeError(f"mesh must be a vtkUnstructuredGrid, got {type(mesh).__name__}.")
    from paraview.simple import RenameSource, Render, Show, TrivialProducer

    producer = TrivialProducer()
    producer.GetClientSideObject().SetOutput(mesh)
    RenameSource(name, producer)
    Show(producer)
    Render()


def _input_grid(in_info: Any) -> vtkUnstructuredGrid:
    grid = vtkUnstructuredGrid.GetData(in_info[0])
    if grid is None:
        raise ValueError("The filter needs an unstructured grid on its input.")
    selection.ensure_selection_arrays(grid)
    return grid


def _modified_callback(algorithm: Any) -> Any:
    """An observer that marks ``algorithm`` modified, without keeping it alive."""
    reference = weakref.ref(algorithm)

    def mark_modified(*_args: Any) -> None:
        target = reference()
        if target is not None:
            target.Modified()

    return mark_modified


@smproxy.filter(name="SelectExtractedData", label="Select Extracted Data")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
@smproperty.input(name="Extracted Selection", port_index=1)
@smdomain.datatype(dataTypes=["vtkDataSet"], composite_data_supported=False)
class SelectExtractedData(VTKPythonAlgorithmBase):
    """Select on the model mesh (port 0) the points and cells of an *Extract Selection* output (port 1).

    The extracted data identifies them by its ``vtkOriginalPointIds`` / ``vtkOriginalCellIds`` arrays.
    """

    def __init__(self) -> None:
        self._mode: SelectionMode = "s"
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=2, nOutputPorts=1)

    def FillInputPortInformation(self, port: int, info: Any) -> int:
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" if port == 0 else "vtkDataSet")
        if port == 1:
            info.Set(self.INPUT_IS_OPTIONAL(), 1)
        return 1

    def RequestData(self, request: Any, in_info: Any, out_info: Any) -> int:
        grid = _input_grid(in_info)
        extracted = vtkDataSet.GetData(in_info[1])
        if extracted is None:
            LOGGER.warning("Select Extracted Data: nothing connected on port 1; selection unchanged.")
            return 1
        for location, data, ids_name in (
            ("node", extracted.GetPointData(), "vtkOriginalPointIds"),
            ("elem", extracted.GetCellData(), "vtkOriginalCellIds"),
        ):
            mask = np.zeros(grid.GetNumberOfPoints() if location == "node" else grid.GetNumberOfCells(), dtype=bool)
            ids_array = data.GetArray(ids_name)
            if ids_array is not None:
                ids = vtknp.vtk_to_numpy(ids_array).astype(np.int64, copy=False)
                mask[ids[(ids >= 0) & (ids < mask.size)]] = True
            selection.apply_mask(grid, mask, location, self._mode)  # type: ignore[arg-type]
        grid.Modified()
        return 1

    @smproperty.xml(_MODE_XML)
    def SetMode(self, mode: int) -> None:
        self._mode = _MODES[int(mode)]
        self.Modified()


@smproxy.filter(name="SelectByCriteria", label="Select")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class Select(VTKPythonAlgorithmBase):
    """Select points or cells by the value of an array (:func:`pytec3d.selection.select`).

    *Values* is one number, a range ``min:max`` or a list ``a, b, c``. *Component* -1 selects on the magnitude.
    """

    def __init__(self) -> None:
        self._array_name = ""
        self._location: Location = "node"
        self._mode: SelectionMode = "s"
        self._component = 0
        self._values_text = "0"
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1)

    def FillInputPortInformation(self, port: int, info: Any) -> int:
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    @staticmethod
    def parse_values(text: str) -> float | tuple[float, float] | list[float]:
        """``"1"`` -> 1.0, ``"0:2"`` -> (0.0, 2.0), ``"1, 3"`` -> [1.0, 3.0].

        Raises:
            ValueError: If the text is empty or not numbers.
        """
        text = text.strip()
        if not text:
            raise ValueError("The 'Values' field is empty.")
        if ":" in text:
            low, _, high = text.partition(":")
            return float(low), float(high)
        if "," in text:
            return [float(part) for part in text.split(",") if part.strip()]
        return float(text)

    def RequestData(self, request: Any, in_info: Any, out_info: Any) -> int:
        grid = _input_grid(in_info)
        if not self._array_name:
            LOGGER.warning("Select: no array chosen; selection unchanged.")
            return 1
        selection.select(
            grid,
            self._array_name,
            self.parse_values(self._values_text),
            component="magnitude" if self._component < 0 else self._component,
            location=self._location,
            mode=self._mode,
        )
        grid.Modified()
        return 1

    @smproperty.xml("""
        <StringVectorProperty name="ArrayName" command="SetArrayName" number_of_elements="1" animateable="0">
          <ArrayListDomain name="array_list" attribute_type="Scalars" input_domain_name="input_array">
            <RequiredProperties><Property name="Input" function="Input"/></RequiredProperties>
          </ArrayListDomain>
        </StringVectorProperty>
    """)
    def SetArrayName(self, name: str) -> None:
        self._array_name = name
        self.Modified()

    @smproperty.xml("""
        <IntVectorProperty name="Type" command="SetSelectionType" number_of_elements="1" default_values="0">
          <EnumerationDomain name="enum">
            <Entry value="0" text="Node"/>
            <Entry value="1" text="Element"/>
          </EnumerationDomain>
        </IntVectorProperty>
    """)
    def SetSelectionType(self, value: int) -> None:
        self._location = _LOCATIONS[int(value)]
        self.Modified()

    @smproperty.xml(_MODE_XML)
    def SetMode(self, value: int) -> None:
        self._mode = _MODES[int(value)]
        self.Modified()

    @smproperty.intvector(name="Component", default_values=0)
    @smdomain.intrange(min=-1, max=2)
    def SetComponent(self, component: int) -> None:
        if not -1 <= int(component) <= 2:
            raise ValueError("Component must be in [-1, 2] (-1 = magnitude).")
        self._component = int(component)
        self.Modified()

    @smproperty.stringvector(name="Values", default_values="0")
    def SetValues(self, text: str) -> None:
        self._values_text = text
        self.Modified()


@smproxy.filter(name="SetValueOnSelection", label="Set Value On Selection")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class SetValueOnSelection(VTKPythonAlgorithmBase):
    """Set a value in the checked arrays at the selected points / cells.

    As :func:`pytec3d.selection.set_value_on_selection`, for every checked array, on points and cells.
    """

    def __init__(self) -> None:
        self._component = 0
        self._value = 0.0
        self._array_names: list[str] = []
        self._arrays = vtkDataArraySelection()
        self._arrays.AddObserver("ModifiedEvent", _modified_callback(self))  # type: ignore[arg-type]
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1)

    def FillInputPortInformation(self, port: int, info: Any) -> int:
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    def RequestInformation(self, request: Any, in_info: Any, out_info: Any) -> int:
        grid = vtkUnstructuredGrid.GetData(in_info[0])
        if grid is None:
            LOGGER.error("Set Value On Selection: no unstructured grid on the input.")
            return 0
        names = []
        for data in (grid.GetPointData(), grid.GetCellData()):
            names += [data.GetArrayName(i) for i in range(data.GetNumberOfArrays())]
        self._array_names = [name for name in names if name]
        for name in self._array_names:
            self._arrays.AddArray(name)
        self._arrays.Modified()
        return 1

    def RequestData(self, request: Any, in_info: Any, out_info: Any) -> int:
        grid = _input_grid(in_info)
        for name in self._array_names:
            if not self._arrays.ArrayIsEnabled(name):
                continue
            for location, data in (("node", grid.GetPointData()), ("elem", grid.GetCellData())):
                if data.GetArray(name) is not None:
                    selection.set_value_on_selection(grid, name, self._value, self._component, location)  # type: ignore[arg-type]
        grid.Modified()
        return 1

    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self) -> Any:
        return self._arrays

    @smproperty.intvector(name="Component", default_values=0)
    @smdomain.intrange(min=0, max=2)
    def SetComponent(self, component: int) -> None:
        if not 0 <= int(component) <= 2:
            raise ValueError("Component must be in [0, 2].")
        self._component = int(component)
        self.Modified()

    @smproperty.doublevector(name="Value", default_values=0.0)
    @smdomain.doublerange()
    def SetValue(self, value: float) -> None:
        if not np.isfinite(float(value)):
            raise ValueError("Value must be finite.")
        self._value = float(value)
        self.Modified()


if PARAVIEW_AVAILABLE:
    setup_paraview_console()
