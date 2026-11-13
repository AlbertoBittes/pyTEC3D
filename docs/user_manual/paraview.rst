ParaView Integration
====================

`ParaView <https://www.paraview.org/>`_ [Ahrens2005]_ (`documentation <https://docs.paraview.org/en/latest/>`_) is the
graphical front end of pyTEC3D: the model is a VTK grid, so ParaView shows it directly, and with
``python/src/pytec3d/pytec3d_paraview_plugin.py`` the mesh can also be **edited interactively** with the pyTEC3D tools, in
two complementary ways:

* the **Python Shell** of ParaView (:paraview:`ParaView and Python
  <Tutorials/SelfDirectedTutorial/batchPythonScripting.html>`): a :class:`~pytec3d.Model` whose mesh is the one
  shown in the view, edited with the whole pytec3d API, and :func:`~pytec3d.pytec3d_paraview_plugin.refresh` to see
  the result;
* the filter **Set Value On Selection**, which writes values on the selection made with the graphical tools of
  ParaView (:paraview:`selecting data <UsersGuide/selectingData.html>`: rubber band, *Find Data*...), without code.

Both work on **one model mesh**, shown as the pipeline source ``Mesh``
(:data:`~pytec3d.pytec3d_paraview_plugin.MESH_NAME`). Its output *is* the model's grid, not a copy: what the console
or the plugin writes is what the view shows after a refresh.

The idea is to prepare the model where you see it: pick regions and boundaries with the mouse or by criteria, check
them visually, write values (Dirichlet arrays, loads, ``TYPE`` / ``MATERIAL`` ids, coupling ids), save the mesh, and
then run the analysis script in the ``pytec3d`` environment. The solution itself is computed outside ParaView (the
compiled solver module is not available in ParaView's Python).

Making pytec3d available in ParaView
------------------------------------

ParaView runs its own Python interpreter, which does **not** see the packages of the ``pytec3d`` conda environment.
The plugin file solves this: when it is executed it adds the folder that contains the ``pytec3d`` package to
``sys.path`` of ParaView's Python, so ``import pytec3d`` works from then on (ParaView's Python provides VTK and NumPy).

**Load the plugin**: *Tools > Manage Plugins > Load New...*, select ``python/src/pytec3d/pytec3d_paraview_plugin.py``
(from the repository, or from the installed package folder); tick *Auto Load* to have it at every start. That is all:

* the filter *Set Value On Selection* is in the *Filters* menu;
* the Python Shell (*View > Python Shell*, also one opened later) has all of pytec3d, **with nothing to import**
  (ParaView's Python has no environment to import from): ``Solver``, ``Model``, ``selection``, ``meshes``
  (:mod:`pytec3d.mesh`), ``model_edit``, ``element_library``, ``tree_gauge``, the coordinate systems
  (``CylindricalCoordinateSystem``...), ``np``, ``pytec3d`` itself, and the helpers ``from_paraview``, ``plot``,
  ``refresh``, ``get_model``, ``get_solver`` and ``MESH_NAME`` (the list is :data:`pytec3d.pytec3d_paraview_plugin.CONSOLE_NAMES`).
  ``selection`` is :mod:`pytec3d.selection`, not the selection module of ParaView;
* the working folder is your **home folder** (or ``$PYTEC3D_PARAVIEW_WORKDIR``), not ParaView's own folder: relative
  paths, ``model.save_model()`` and the ``<mesh_name>.log`` of ``solver.read_model()`` go there.

Without loading the plugin, *Run Script* of the same file in the Python Shell gives the console part only. The
plugin and *Run Script* share the shown model for the whole session.

Use a ParaView version whose Python has the same major version as the environment (3.12) and NumPy available
(ParaView 5.12 or later binaries).

Working in the Python Shell
---------------------------

There are two ways to get the model into ParaView; either way it becomes the source ``Mesh``.

**A .vtu dragged into ParaView** (or opened with *File > Open*): select it in the Pipeline Browser, then

.. code-block:: python

   model = from_paraview()             # or from_paraview("coil_2D_avolt.vtu"): the source by name

The model is a copy of the reader's grid, named after the file; the reader is hidden and ``Mesh`` shows the model.

**A model read in the console**:

.. code-block:: python

   model = Model.read_model("/path/to/coil_2D_avolt.vtu")
   plot(model)

Then edit the model with the pytec3d API and call ``refresh()`` to update the view:

.. code-block:: python

   grid = model.mesh
   selection.select(grid, "MATERIAL", 23, location="elem")
   selection.nsle(grid)
   selection.set_value_on_selection(grid, "CP_0", 1)
   refresh()                           # colour by "selection_point" or "CP_0" in the GUI
   model.save_model("/path/to/coil_2D_prepared.vtu")

``refresh()`` sets the model's current grid as the output of ``Mesh``, marks it modified, updates the pipeline, shows
the quadratic cells as they are (*Triangulate* and *Nonlinear Subdivision Level* 0) and renders. ``get_model()``
returns the shown model, e.g. after a restart of the console namespace; ``plot(other)`` shows another model (or a
grid) in the same ``Mesh``.

To solve the model, attach it to a solver; ``plot(solver)`` then shows the solver's current model (also after
``solver.read_model`` of another file), and ``get_solver()`` returns it:

.. code-block:: python

   solver = Solver()
   solver.model = model
   plot(solver)

Every operation of :doc:`selection`, :doc:`boundary_conditions` and :doc:`coordinate_systems` is available, so the
console is a good place to try selections (by position in a cylindrical system, by array value, by topology) before
writing them in a macro.

Set Value On Selection
----------------------

After loading the plugin, *Set Value On Selection* is in the *Filters* menu (search with Ctrl+Space). It is applied on
``Mesh`` and acts on it; its output is empty (the work is done by its button, not by the pipeline). Applied on another
source, its button refuses.

1. Show the model (``from_paraview()`` or ``plot(model)``), select ``Mesh`` in the Pipeline Browser, add the filter
   *Set Value On Selection* and *Apply*: its *Arrays* list has the arrays of the model mesh.
2. Check the *Arrays* to write, set *Value* and *Component* (of a vector array; ignored for a scalar), *Apply*.
3. Select on ``Mesh`` with ParaView's tools (rubber band, *Find Data* with its query: no *Extract Selection*
   needed), then press **Set Value**. The view is refreshed. Select again and press *Set Value* again for the next
   region; *Apply* only when the properties change.

What is written:

* on a **cell selection**: the cell arrays on the selected cells, the point arrays on all the points of those cells;
* on a **point selection**: the point arrays on the selected points; cell arrays are left unchanged (a warning);
* with **no ParaView selection** on ``Mesh``: the pyTEC3D selection (arrays ``selection_point`` / ``selection_cell``,
  as :func:`pytec3d.selection.set_value_on_selection`), e.g. a selection made in the console.

An array that exists on both points and cells is written on both. Save the model from the console
(``model.save_model(path)``) or with *File > Save Data* on ``Mesh`` (``.vtu``): the analysis script reads it with
``Solver.read_model(path, reset_arrays=False)`` so that the values written in ParaView are kept.

Viewing the results
-------------------

Open the output series ``output/<name>_out_*.vtu``; colour by ``B``, ``Jt``, ``A``; use *Glyph* for vector fields,
*Threshold* on ``TYPE`` or ``MATERIAL`` to isolate a region, *Plot Over Line*, *Integrate Variables* (see
:doc:`outputs`).
