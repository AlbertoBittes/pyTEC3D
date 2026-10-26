# pytec3d-element-library

The element library of pyTEC3D: the JSON definitions of the domains, fields, element types and nonlinear laws
(`pytec3d_element_library/data/`). It is the single source of the library, shared by:

- the C++ solver (`../solver`, distribution `pytec3d-solver`), which compiles it into its module when it is built;
- the Python API (`../python`, distribution `pytec3d`), which reads it at run time.

```python
import pytec3d_element_library as library
library.data_dir()                    # the directory of the JSON files
library.load("element_types.json")    # the content of one of them
```

`python -m pytec3d_element_library` prints the directory (the solver build uses it). See `docs/user_manual/element_library.rst`
of the repository for the format.
