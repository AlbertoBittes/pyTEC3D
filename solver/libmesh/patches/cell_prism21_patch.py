"""Prism21::connectivity: the VTK node order (libMesh has none: libmesh_not_implemented).

Usage: python cell_prism21_patch.py <libmesh>/src/geom/cell_prism21.C

VTK has no fixed-order 21-point wedge; its VTK_LAGRANGE_WEDGE of order 2 has 21 points: the PRISM15 points in the
VTK_QUADRATIC_WEDGE order, the centres of the triangles (bottom, top), of the quadrilaterals (0-1-4-3, 1-2-5-4, 2-0-3-5)
and of the cell. From libMesh's order it differs in the rings of edge nodes (VTK: bottom, top, vertical; libMesh:
bottom, vertical, top) and in the face nodes (libMesh: quadrilaterals, then triangles). Tecplot stays not implemented.
"""

import re
import sys

file_path = sys.argv[1]

with open(file_path, 'r') as f:
    content = f.read()

new_func = '''void Prism21::connectivity(const unsigned int libmesh_dbg_var(sc),
                           const IOPackage iop,
                           std::vector<dof_id_type> & conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case VTK:
      {
        // VTK_LAGRANGE_WEDGE of order 2 (21 points): conn[k] is the libMesh node of VTK point k.
        static const unsigned int vtk_to_libmesh[21] =
          {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 18, 19, 15, 16, 17, 20};
        conn.resize(21);
        for (unsigned int k = 0; k != 21; ++k)
          conn[k] = this->node_id(vtk_to_libmesh[k]);
        return;
      }

    default:
      libmesh_not_implemented(); // FIXME RHS: Tecplot
    }
}'''

# The definition, up to the closing brace at the start of a line (the old body is commented out with /* */ inside).
pattern = r'void\s+Prism21::connectivity\s*\([\s\S]*?^\}'

if not re.search(pattern, content, re.MULTILINE):
    print('Error: Could not locate Prism21::connectivity function in target file.', file=sys.stderr)
    sys.exit(1)
if 'vtk_to_libmesh' in content:
    print('Error: Prism21::connectivity is already patched.', file=sys.stderr)
    sys.exit(1)

new_content = re.sub(pattern, lambda m: new_func, content, count=1, flags=re.MULTILINE)

with open(file_path, 'w') as f:
    f.write(new_content)

print(f'Successfully patched Prism21::connectivity in {file_path}')
