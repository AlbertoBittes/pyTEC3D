import os
import sys

project = "pyTEC3D"
author = "pyTEC3D Team"
release = "0.1.0"
version = release

# Add project root and python package to path so autodoc can find modules
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../python/src"))
sys.path.insert(0, os.path.abspath("../element_library"))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    # Jupyter notebooks of the user manual, shown with the outputs stored in them (not executed at build time: the
    # C++ solver is not available on the documentation builders).
    "myst_nb",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "refactor_open_points.md", "**.ipynb_checkpoints"]

nb_execution_mode = "off"
# $...$ and $$...$$ math in the markdown cells, as Jupyter renders it.
myst_enable_extensions = ["dollarmath", "amsmath"]
numfig = True
math_number_all = False

autosummary_generate = True
autosummary_generate_overwrite = True
# Each object is documented once, in its module (the names re-exported by pytec3d point there).
autosummary_imported_members = False
suppress_warnings = ["autodoc.mocked_object"]
# "Attributes:" sections of the docstrings as instance variable fields (no duplicate attribute entries).
napoleon_use_ivar = True

# Mock heavy optional dependencies to keep RTD builds light (numpy is installed: module-level arithmetic needs it)
autodoc_mock_imports = [
    "mpi4py",
    "vtk",
    "vtkmodules",
    "vtk.util",
    "pytec3d_solver.pytec3d_solver",
    "paraview",
    "ansys",
    "ansys.mapdl",
    "ansys.mapdl.reader",
]

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]        # scrollable tables

# Links to the documentation of the libraries: :petsc:`KSPGMRES <KSP/KSPGMRES>` (PETSc manual pages),
# :vtk:`vtkUnstructuredGrid` (VTK class reference), :paraview:`text <UsersGuide/selectingData.html>` (ParaView guide),
# :libmesh:`NewtonSolver` (libMesh class reference).
extlinks = {
    "petsc": ("https://petsc.org/release/manualpages/%s/", "%s"),
    "vtk": ("https://vtk.org/doc/nightly/html/class%s.html", "%s"),
    "paraview": ("https://docs.paraview.org/en/latest/%s", "%s"),
    "libmesh": ("https://libmesh.github.io/doxygen/classlibMesh_1_1%s.html", "%s"),
}
# Python, NumPy and pyvista objects named in the docstrings link to their documentation.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pyvista": ("https://docs.pyvista.org", None),
}
html_last_updated_fmt = "%Y-%m-%d"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False}
