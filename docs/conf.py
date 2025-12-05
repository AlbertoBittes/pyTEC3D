import os
import sys
from datetime import datetime

project = "APEC-M"
author = "APEC-M Team"
release = "0.1.0"
version = release

# Add project root and api package to path so autodoc can find modules
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../api"))

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autosummary_generate = True
autosummary_imported_members = True

# Mock heavy optional dependencies to keep RTD builds light
autodoc_mock_imports = [
    "numpy",
    "pandas",
    "vtk",
    "vtkmodules",
    "vtk.util",
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
html_last_updated_fmt = "%Y-%m-%d"
