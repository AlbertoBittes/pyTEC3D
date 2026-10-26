ParaView Integration
====================

`ParaView <https://www.paraview.org/>`_ [Ahrens2005]_ (`documentation <https://docs.paraview.org/en/latest/>`_) is the
graphical front end of pyTEC3D: the model is a VTK grid, so ParaView shows it directly, and with
``python/src/pytec3d/pytec3d_paraview_plugin.py`` the mesh can also be **edited interactively** with the pyTEC3D tools, in
two complementary ways:

* the **Python Shell** of ParaView (:paraview:`ParaView and Python
  <Tutorials/SelfDirectedTutorial/batchPythonScripting.html>`): the whole :mod:`pytec3d.selection`,
  :mod:`pytec3d.mesh` and :mod:`pytec3d.model_edit` API, run on the grids loaded in the GUI, with the results shown at
  once;
* three **filters** (*Select*, *Select Extracted Data*, *Set Value On Selection*) that do the same operations from
  the pipeline browser, without code, and combine with the graphical selection tools of ParaView
  (:paraview:`selecting data <UsersGuide/selectingData.html>`).

The idea is to prepare the model where you see it: pick regions and boundaries with the mouse or by criteria, check
them visually, write values (Dirichlet arrays, loads, ``TYPE`` / ``MATERIAL`` ids, coupling ids), save the mesh, and
then run the analysis script in the ``pytec3d`` environment. The solution itself is computed outside ParaView (the
compiled solver module is not available in ParaView's Python).

Making pytec3d available in ParaView
------------------------------------

ParaView runs its own Python interpreter, which does **not** see the packages of the ``pytec3d`` conda environment.
The plugin file solves this: when it is executed it adds the folder that contains the ``pytec3d`` package to
``sys.path`` of ParaView's Python, so ``import pytec3d`` works from then on (ParaView's Python provides VTK and NumPy).

* **For the console**: open *View > Python Shell*, click **Run Script** and choose
  ``python/src/pytec3d/pytec3d_paraview_plugin.py`` (from the repository, or from the installed package folder), once per
  ParaView session. The shell then has pytec3d and its tools, **with nothing to import**: ``selection``, ``meshes``
  (:mod:`pytec3d.mesh`), ``model_edit``, ``element_library``, ``tree_gauge``, ``Model``, the coordinate systems
  (``CylindricalCoordinateSystem``...), ``np``, ``pytec3d`` itself and ``plot`` (the list is
  :data:`pytec3d.pytec3d_paraview_plugin.CONSOLE_NAMES`). ``selection`` is :mod:`pytec3d.selection`, not the selection
  module of ParaView.
* **For the filters**: *Tools > Manage Plugins > Load New...*, select the same file (tick *Auto Load* to have the
  filters at every start). Loaded this way, the plugin also puts the names of the console in the Python Shell if the
  shell is already open; otherwise (e.g. *Auto Load* at the start of ParaView) use *Run Script* as above.

The two can be combined: the filters from the Plugin Manager, the console from *Run Script*.

Use a ParaView version whose Python has the same major version as the environment (3.12) and NumPy available
(ParaView 5.12 or later binaries).

Working in the Python Shell
---------------------------

After *Run Script* of the plugin file, with no import:

.. code-block:: python

   # 1. Read a model with pytec3d and show it in the view
   model = Model.read_model("/path/to/coil_2D_avolt.vtu", reset_arrays=True)
   grid = model.mesh
   plot(grid, "coil")                         # a new source "coil" in the pipeline browser

   # 2. Select and inspect: the selection arrays are ordinary arrays, colour by them in the GUI
   selection.select(grid, "MATERIAL", 23, location="elem")
   selection.nsle(grid)
   print(selection.selected_ids(grid, "node").size, "nodes selected")
   grid.Modified(); Render()                  # refresh the view; colour by "selection_point"

   # 3. Write values and save the prepared mesh for the analysis script
   selection.set_value_on_selection(grid, "CP_0", 1)
   model.save_model("/path/to/coil_2D_prepared.vtu")

A grid already opened in the GUI (*File > Open*) can be fetched into the console as a VTK object:

.. code-block:: python

   from paraview import servermanager                # ParaView's own module
   grid = servermanager.Fetch(GetActiveSource())     # a copy of the data of the selected pipeline item
   selection.select_all(grid); selection.nsel_ext(grid)
   plot(grid, "boundary nodes")                     # show the edited copy

Every operation of :doc:`selection`, :doc:`boundary_conditions` and :doc:`coordinate_systems` is available, so the
console is a good place to try selections (by position in a cylindrical system, by array value, by topology) before
writing them in a macro.

The filters
-----------

After *Load New* of the plugin, three filters appear in the *Filters* menu (search with Ctrl+Space). They work on an
unstructured grid input and edit its selection arrays and values, as :mod:`pytec3d.selection` does:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Filter
     - Effect
   * - **Select**
     - Select points or cells by the value of an array. Properties: *ArrayName*, *Type* (Node / Element), *Mode*
       (Set, Add, Reselect, Unselect), *Component* (-1 for the magnitude), *Values* (``"1"`` one value, ``"0:2"`` an
       inclusive range, ``"1, 3"`` a list).
   * - **Select Extracted Data**
     - Turn an interactive ParaView selection into a pyTEC3D selection. Select cells or points in the view (rubber
       band, *Find Data*...), apply *Extract Selection*, then apply *Select Extracted Data* with the model as input 0
       and the extraction as input 1, with a *Mode*.
   * - **Set Value On Selection**
     - Write *Value* (and *Component*) in the checked *Arrays* at the selected points / cells.

Filters can be chained (e.g. *Select* the cells of a material, then *Select Extracted Data* to add a hand-picked
region, then *Set Value On Selection*). Save the last output with *File > Save Data* (``.vtu``): the analysis script
reads it with ``Solver.read_model(path, reset_arrays=False)`` so that the values written in ParaView are kept.

Viewing the results
-------------------

Open the output series ``output/<name>_out_*.vtu``; colour by ``B``, ``Jt``, ``A``; use *Glyph* for vector fields,
*Threshold* on ``TYPE`` or ``MATERIAL`` to isolate a region, *Plot Over Line*, *Integrate Variables* (see
:doc:`outputs`).
