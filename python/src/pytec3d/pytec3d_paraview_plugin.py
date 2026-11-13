"""ParaView integration of pytec3d: one model mesh shown as the source ``Mesh``, edited from the console or the plugin.

Two uses, which can be combined:

* Console: load this file as a plugin (below), or *Run Script* of it in the ParaView Python shell. The shell then has
  all of pytec3d with nothing to import - ``Solver``, ``Model``, ``selection``, ``meshes`` (:mod:`pytec3d.mesh`),
  ``model_edit``, ``element_library``, ``tree_gauge``, the coordinate systems, ``np`` - and the helpers of this file
  (:data:`CONSOLE_NAMES`):

  - :func:`from_paraview` - a model whose mesh is the grid of a source of the pipeline, e.g. a ``.vtu`` dragged
    into ParaView: ``model = from_paraview()``;
  - :func:`plot` - show a model read in the console: ``model = Model.read_model(...); plot(model)``, or the model of
    a solver: ``plot(solver)``;
  - :func:`refresh` - show the current state of the model mesh after editing it in the console;
  - :func:`get_model`, :func:`get_solver` - the shown model, the shown solver (after ``plot(solver)``).

  A model is solved by attaching it to a solver: ``solver = Solver(); solver.model = model``.

  Either way the model mesh is the source :data:`MESH_NAME` (``"Mesh"``): its output *is* the model's grid, not a
  copy, so what the console changes is what ParaView shows after :func:`refresh`.
* Plugin: load this file in ParaView (Tools > Manage Plugins > Load New) to get the source *Set Value On Selection*.
  It writes a value into the checked arrays of the model mesh, on ParaView's selection of the ``Mesh`` source
  (interactive, or *Find Data*), or on the pyTEC3D selection arrays ``selection_point`` / ``selection_cell`` when
  ParaView has no selection, then refreshes the view. Loaded this way it also puts the console names in the Python
  shell, when the shell already exists.

The shown solver and mesh are kept in one place for the whole ParaView session (:data:`STATE_MODULE`), so the
plugin and a *Run Script* of this file share them.
"""

# No ``from __future__`` import: *Run Script* runs this file after a first line of its own (``__file__ = ...``).

import os
import sys
import types
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
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

# The console names (CONSOLE_NAMES): with *Run Script* every top-level name of this file is a name of the shell.
# ParaView's Python has no environment to import from, so the console gets all of pytec3d.
import pytec3d
from pytec3d import (  # noqa: F401
    CartesianCoordinateSystem,
    CylindricalCoordinateSystem,
    Model,
    Solver,
    SphericalCoordinateSystem,
    XAxisCylindricalCoordinateSystem,
    YAxisCylindricalCoordinateSystem,
    ZAxisCylindricalCoordinateSystem,
    element_library,
    model_edit,
    selection,
    tree_gauge,
)
from pytec3d import mesh as meshes
from pytec3d._logging import get_logger
from pytec3d._typing import Location

LOGGER = get_logger("pytec3d.paraview")

#: Name of the pipeline source that shows the model mesh.
MESH_NAME = "Mesh"

#: Name of the module (in ``sys.modules``) that keeps the shown solver and mesh for the whole session.
STATE_MODULE = "_pytec3d_paraview_state"

#: The names this file puts in the ParaView Python shell.
CONSOLE_NAMES: tuple[str, ...] = (
    "np",
    "pytec3d",
    "Solver",
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
    "from_paraview",
    "plot",
    "refresh",
    "get_model",
    "get_solver",
    "MESH_NAME",
)

try:
    from paraview.util.vtkAlgorithm import smdomain, smproperty, smproxy
    from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
    from vtkmodules.vtkCommonCore import vtkDataArraySelection

    # Not a mocked paraview (the documentation build): there the console setup must not run
    PARAVIEW_AVAILABLE = not getattr(smproxy, "__sphinx_mock__", False)
except ImportError:  # Outside ParaView: the decorators do nothing and the plugin cannot run.
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

    from vtkmodules.vtkCommonCore import vtkDataArraySelection


# ------------------------------------------------------------------------------------------------ shared state


def _state() -> types.SimpleNamespace:
    """The session state: ``shown``, what :data:`MESH_NAME` shows - a :class:`~pytec3d.Model`, a
    :class:`~pytec3d.Solver` (its model) or a ``vtkUnstructuredGrid`` - or None.

    Held by a module of its own in ``sys.modules``, so that every copy of this file (the plugin, *Run Script*)
    shares it.
    """
    module = sys.modules.get(STATE_MODULE)
    if module is None:
        module = types.ModuleType(STATE_MODULE, "pytec3d ParaView session state")
        sys.modules[STATE_MODULE] = module
    return module.__dict__.setdefault("session", types.SimpleNamespace(shown=None))


def get_model() -> Model:
    """The shown model (:func:`from_paraview`, or :func:`plot` of a model or a solver).

    Raises:
        RuntimeError: If no model is shown (nothing, or a bare grid).
    """
    shown = _state().shown
    if isinstance(shown, Solver):
        return shown.model
    if isinstance(shown, Model):
        return shown
    raise RuntimeError("No pyTEC3D model is shown: use model = from_paraview() or plot(model).")


def get_solver() -> Solver:
    """The shown solver, when :func:`plot` was given a solver.

    Raises:
        RuntimeError: If what is shown is not a solver.
    """
    shown = _state().shown
    if isinstance(shown, Solver):
        return shown
    raise RuntimeError("No pyTEC3D solver is shown: plot(solver) shows the model of a solver "
                       "(e.g. solver = Solver(); solver.model = get_model(); plot(solver)).")


def _grid_of(target: Any) -> vtkUnstructuredGrid:
    """The grid of a solver (of its current model), a model or a grid.

    Raises:
        TypeError: If ``target`` is none of them.
        ValueError: If the model has no mesh.
    """
    if isinstance(target, Solver):
        grid = target.model.mesh
    elif isinstance(target, Model):
        grid = target.mesh
    elif isinstance(target, vtkUnstructuredGrid):
        grid = target
    else:
        raise TypeError(f"Expected a pytec3d Model, Solver or vtkUnstructuredGrid, got {type(target).__name__}.")
    if grid is None:
        raise ValueError("The model has no mesh: read a model first (Model.read_model(...)).")
    return grid


def _current_grid() -> vtkUnstructuredGrid:
    """The grid to show: that of what is shown now (a solver may have read another model since :func:`plot`).

    Raises:
        RuntimeError: If nothing is shown.
    """
    shown = _state().shown
    if shown is None:
        raise RuntimeError("No model mesh is shown: use model = from_paraview() or plot(model).")
    return _grid_of(shown)


# ------------------------------------------------------------------------------------------------ console


def working_directory() -> str:
    """The working folder of pytec3d in ParaView: ``$PYTEC3D_PARAVIEW_WORKDIR``, else the home folder.

    Not ParaView's own folder (where it often starts, read-only on Windows): relative paths, ``solver.save_model()``
    and the ``<mesh_name>.log`` of :meth:`pytec3d.Solver.read_model` use the working folder.
    """
    return os.path.abspath(os.environ.get("PYTEC3D_PARAVIEW_WORKDIR") or os.path.expanduser("~"))


def setup_paraview_console(namespace: dict[str, Any] | None = None, change_directory: bool = True) -> dict[str, Any]:
    """Make :data:`CONSOLE_NAMES` available in the ParaView Python shell and move to :func:`working_directory`.

    The names go into the namespace of the shell when it exists, and into ``builtins``, so that a shell opened later
    (e.g. after a plugin *Auto Load* at the start of ParaView) has them too.

    Args:
        namespace: The namespace of the shell. By default the one of ParaView's Python shell (``__vtkConsoleLocals``
            of ``__main__``), when the shell exists.
        change_directory: Change the working folder to :func:`working_directory` (created when missing).

    Returns:
        The names and their values.
    """
    import builtins

    names = globals()
    values = {name: names[name] for name in CONSOLE_NAMES}
    if namespace is None:
        import __main__

        namespace = getattr(__main__, "__vtkConsoleLocals", None)
    if namespace is not None:
        namespace.update(values)
    for name, value in values.items():
        setattr(builtins, name, value)
    if change_directory:
        folder = working_directory()
        try:
            os.makedirs(folder, exist_ok=True)
            os.chdir(folder)
        except OSError as error:
            LOGGER.warning("Could not move to the working folder '%s' (%s); staying in '%s'.", folder, error,
                           os.getcwd())
    LOGGER.info(
        "pytec3d %s in the ParaView Python shell (from %s, working folder %s): %s.",
        pytec3d.__version__,
        _PACKAGE_PARENT,
        os.getcwd(),
        ", ".join(CONSOLE_NAMES),
    )
    return values


def _mesh_proxy(create: bool = False) -> Any:
    """The source :data:`MESH_NAME` (a TrivialProducer), or None; with ``create``, made (and shown) when missing.

    Another kind of source with that name (e.g. a reader the user renamed) is renamed out of the way.
    """
    from paraview.simple import FindSource, RenameSource, Show, TrivialProducer

    proxy = FindSource(MESH_NAME)
    if proxy is not None and proxy.GetXMLName() != "TrivialProducer":
        RenameSource(f"{MESH_NAME} (original)", proxy)
        proxy = None
    if proxy is None and create:
        proxy = TrivialProducer()
        proxy.GetClientSideObject().SetOutput(_current_grid())
        RenameSource(MESH_NAME, proxy)
        Show(proxy)
    return proxy


def refresh() -> None:
    """Show the current state of the model mesh in ParaView, after it changed (arrays, selection, another model).

    Raises:
        RuntimeError: If no model mesh is shown.
    """
    from paraview.simple import GetDisplayProperties, Render, UpdatePipeline

    grid = _current_grid()
    proxy = _mesh_proxy(create=True)
    display = GetDisplayProperties(proxy)
    proxy.GetClientSideObject().SetOutput(grid)
    proxy.MarkModified(proxy)
    UpdatePipeline(proxy=proxy)
    # Quadratic cells as they are, not split into linear pieces
    properties = display.ListProperties()
    if "Triangulate" in properties:
        display.Triangulate = 0
    if "NonlinearSubdivisionLevel" in properties:
        display.NonlinearSubdivisionLevel = 0
    Render()


def plot(target: Any = None) -> Any:
    """Show a model mesh as the source :data:`MESH_NAME`, e.g. ``plot(model)`` after ``model = Model.read_model(...)``.

    The source shows the grid itself, so later changes need only :func:`refresh`.

    Args:
        target: A :class:`~pytec3d.Model`, a :class:`~pytec3d.Solver` (its current model, also after it reads another
            one) or a ``vtkUnstructuredGrid``; None shows again the current one.

    Returns:
        The source :data:`MESH_NAME`.

    Raises:
        TypeError: If ``target`` is not a model, a solver or a grid.
        ImportError: Outside ParaView.
    """
    if target is not None:
        _grid_of(target)  # checked before it replaces what is shown
        _state().shown = target
    refresh()
    return _mesh_proxy()


def from_paraview(source: Any = None, mesh_name: str | None = None) -> Model:
    """A model whose mesh is the grid of a pipeline source, e.g. a ``.vtu`` dragged into ParaView, shown as
    :data:`MESH_NAME` (the source itself is hidden): ``model = from_paraview()``.

    The grid is copied from the source: the model is the copy, which :func:`refresh` shows. To solve it, attach it
    to a solver: ``solver = Solver(); solver.model = model``.

    Args:
        source: The source, or its name in the Pipeline Browser; the active source by default.
        mesh_name: Name of the model; the file name of a reader by default.

    Returns:
        The model (also :func:`get_model`).

    Raises:
        ValueError: If there is no such source.
        TypeError: If its output is not an unstructured grid.
    """
    from paraview import servermanager
    from paraview.simple import FindSource, GetActiveSource, Hide

    if isinstance(source, str):
        name = source
        source = FindSource(name)
        if source is None:
            raise ValueError(f"No source '{name}' in the pipeline.")
    elif source is None:
        source = GetActiveSource()
        if source is None:
            raise ValueError("No active source: select the dragged .vtu in the Pipeline Browser, or pass its name.")

    state = _state()
    if (state.shown is not None and source.GetXMLName() == "TrivialProducer"
            and source.GetClientSideObject().GetOutputDataObject(0) is _current_grid()):
        if isinstance(state.shown, vtkUnstructuredGrid):  # a bare grid shown: make it a model
            state.shown = Model.from_grid(state.shown, mesh_name=mesh_name)
        return get_model()

    source.UpdatePipeline()
    grid = servermanager.Fetch(source)
    if not isinstance(grid, vtkUnstructuredGrid):
        raise TypeError(f"The source gives a {type(grid).__name__}, not a vtkUnstructuredGrid (a .vtu mesh).")
    if mesh_name is None and "FileName" in source.ListProperties():
        file_name = source.FileName
        file_name = file_name[0] if isinstance(file_name, (list, tuple)) and file_name else file_name
        if isinstance(file_name, str) and file_name:
            mesh_name = os.path.splitext(os.path.basename(file_name))[0]
    model = Model.from_grid(grid, mesh_name=mesh_name)
    Hide(source)
    plot(model)
    return model


# ------------------------------------------------------------------------------------------------ selection


def _element_type(selection_proxy: Any) -> Location | None:
    """``"node"`` or ``"elem"``: what a ParaView selection selects, from its source properties; None when unknown.

    ``ElementType`` (ParaView >= 5.11) counts as vtkDataObject attribute types (POINT 0, CELL 1), the older
    ``FieldType`` as vtkSelectionNode field types (CELL 0, POINT 1); an enumeration may give its text.
    """
    properties = selection_proxy.ListProperties()
    for name, point_value in (("ElementType", 0), ("FieldType", 1)):
        if name not in properties:
            continue
        value = selection_proxy.GetPropertyValue(name)
        if isinstance(value, (list, tuple)):
            value = value[0] if value else None
        if isinstance(value, str):
            text = value.upper()
            if "POINT" in text:
                return "node"
            if "CELL" in text:
                return "elem"
            continue
        if value is not None:
            return "node" if int(value) == point_value else "elem"
    return None


def paraview_selection(proxy: Any) -> tuple[np.ndarray, np.ndarray | None] | None:
    """The points and cells of ParaView's current selection of the source ``proxy`` (interactive or *Find Data*).

    The selection is evaluated on the source by an *Extract Selection* that is not added to the pipeline.

    Returns:
        ``(point ids, cell ids)``: for a cell selection the selected cells and all their points, for a point selection
        the points and None; None when the source has no selection.
    """
    from paraview import servermanager

    selection_input = proxy.GetSelectionInput(proxy.Port)
    if selection_input is None:
        return None
    selection_proxy = servermanager._getPyProxy(selection_input)
    extract = servermanager.filters.ExtractSelection(Input=proxy, Selection=selection_proxy)
    try:
        extract.UpdatePipeline()
        extracted = servermanager.Fetch(extract)
    finally:
        del extract
    if extracted is None or not hasattr(extracted, "GetPointData"):
        LOGGER.warning("Set Value On Selection: the selection gave a %s; ignored.", type(extracted).__name__)
        return None

    def ids(data: Any, name: str) -> np.ndarray:
        array = data.GetArray(name)
        if array is None:
            return np.zeros(0, dtype=np.int64)
        return np.unique(vtknp.vtk_to_numpy(array).astype(np.int64, copy=False))

    point_ids = ids(extracted.GetPointData(), "vtkOriginalPointIds")
    cell_ids = ids(extracted.GetCellData(), "vtkOriginalCellIds")
    kind = _element_type(selection_proxy)
    if kind is None:  # e.g. a combination of selections: a cell selection extracts cells that are not vertices
        types = (vtknp.vtk_to_numpy(extracted.GetDistinctCellTypesArray())
                 if extracted.GetNumberOfCells() else np.zeros(0))
        kind = "elem" if cell_ids.size and np.any(types != 1) else "node"  # 1 = VTK_VERTEX
    return (point_ids, cell_ids) if kind == "elem" else (point_ids, None)


def set_value_at(grid: vtkUnstructuredGrid, name: str, location: Location, ids: np.ndarray, value: float,
                 component: int) -> None:
    """Write ``value`` in component ``component`` (ignored for a scalar array) of the array ``name`` at ``ids``.

    Raises:
        ValueError: If the array has no such component.
    """
    array = meshes.get_array(grid, name, location)
    values = vtknp.vtk_to_numpy(array)
    if values.ndim == 1:
        values[ids] = value
    else:
        if component >= values.shape[1]:
            raise ValueError(f"The array '{name}' has {values.shape[1]} components: no component {component}.")
        values[ids, component] = value
    array.Modified()


def set_value_on_selection(names: list[str], value: float, component: int = 0) -> dict[str, int]:
    """Write ``value`` in the arrays ``names`` of the shown mesh, on ParaView's selection of the source
    :data:`MESH_NAME`, or on the pyTEC3D selection arrays when ParaView has none, then :func:`refresh`.

    An array that exists on both points and cells is written on both. On a point selection, cell arrays are left
    unchanged (warning).

    Returns:
        The number of values written per array (``"<name> (point)"`` / ``"<name> (cell)"``).
    """
    grid = _current_grid()
    selected = paraview_selection(_mesh_proxy(create=True))
    if selected is None:
        source = "the pyTEC3D selection arrays"
        point_ids, cell_ids = selection.selected_ids(grid, "node"), selection.selected_ids(grid, "elem")
    else:
        source = "the ParaView selection"
        point_ids, cell_ids = selected
    written: dict[str, int] = {}
    for name in names:
        targets: tuple[tuple[Location, Any, np.ndarray | None], ...] = (
            ("node", grid.GetPointData(), point_ids), ("elem", grid.GetCellData(), cell_ids))
        for location, data, ids in targets:
            if data.GetArray(name) is None:
                continue
            if ids is None:
                LOGGER.warning("Set Value On Selection: '%s' is a cell array and the ParaView selection has points "
                               "only; not written.", name)
                continue
            set_value_at(grid, name, location, ids, value, component)
            written[f"{name} ({'point' if location == 'node' else 'cell'})"] = int(ids.size)
    grid.Modified()
    LOGGER.info("Set Value On Selection: %g (component %d) on %s: %s.", value, component, source,
                ", ".join(f"{key} {count}" for key, count in written.items()) or "no array written")
    refresh()
    return written


# ------------------------------------------------------------------------------------------------ plugin


def _modified_callback(algorithm: Any) -> Any:
    """An observer that marks ``algorithm`` modified, without keeping it alive."""
    reference = weakref.ref(algorithm)

    def mark_modified(*_args: Any) -> None:
        target = reference()
        if target is not None:
            target.Modified()

    return mark_modified


@smproxy.filter(name="SetValueOnSelection", label="Set Value On Selection")
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkUnstructuredGrid"], composite_data_supported=False)
class SetValueOnSelection(VTKPythonAlgorithmBase):
    """Set a value in the checked arrays of the model mesh (the source ``Mesh``, its input) on the current selection.

    The selection is ParaView's selection of ``Mesh`` (interactive, or *Find Data*), or the pyTEC3D selection arrays
    when there is none. Set the properties, *Apply*, then press *Set Value* (once per selection): the view is
    refreshed. The work is done by the button, not by the pipeline, so the output is empty and the view can be
    refreshed from it.
    """

    def __init__(self) -> None:
        self._component = 0
        self._value = 0.0
        self._arrays = vtkDataArraySelection()
        self._arrays.AddObserver("ModifiedEvent", _modified_callback(self))  # type: ignore[arg-type]
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkUnstructuredGrid")

    def FillInputPortInformation(self, port: int, info: Any) -> int:
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid")
        return 1

    def _input_grid(self) -> Any:
        return self.GetInputDataObject(0, 0) if self.GetNumberOfInputConnections(0) else None

    def _update_array_list(self) -> None:
        grid = self._input_grid()
        if grid is None:
            return
        for data in (grid.GetPointData(), grid.GetCellData()):
            for i in range(data.GetNumberOfArrays()):
                name = data.GetArrayName(i)
                if name and not self._arrays.ArrayExists(name):
                    self._arrays.AddArray(name, False)

    def RequestInformation(self, request: Any, in_info: Any, out_info: Any) -> int:
        self._update_array_list()
        return 1

    def RequestData(self, request: Any, in_info: Any, out_info: Any) -> int:
        return 1  # empty output: the work is done by the Set Value button

    def checked_arrays(self) -> list[str]:
        return [self._arrays.GetArrayName(i) for i in range(self._arrays.GetNumberOfArrays())
                if self._arrays.GetArraySetting(i)]

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

    @smproperty.xml("""
        <Property name="SetValueButton" label="Set Value" command="Run" panel_widget="command_button">
          <Documentation>Write the value on the current selection of Mesh (Apply the properties first).</Documentation>
        </Property>
    """)
    def Run(self) -> None:
        """The *Set Value* button: write the value in the checked arrays on the current selection, refresh."""
        try:
            shown = _current_grid()
            if self._input_grid() is not shown:
                raise RuntimeError(f"apply the filter on the source '{MESH_NAME}' (the model mesh): its input is "
                                   "another source. Show the model with from_paraview() or plot(solver).")
            names = self.checked_arrays()
            if not names:
                LOGGER.warning("Set Value On Selection: no array checked; nothing written.")
                return
            set_value_on_selection(names, self._value, self._component)
        except Exception as error:  # a button: report, do not let ParaView swallow it silently
            LOGGER.error("Set Value On Selection: %s", error)
            raise
        finally:
            self._update_array_list()


if PARAVIEW_AVAILABLE:
    try:
        setup_paraview_console()
    except Exception as _error:  # the plugin and its filter must load even so
        LOGGER.warning("pytec3d console setup failed: %s", _error)
