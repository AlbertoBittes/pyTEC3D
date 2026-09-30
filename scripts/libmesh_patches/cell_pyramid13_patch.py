import sys, re

file_path = sys.argv[1]

with open(file_path, 'r') as f:
    content = f.read()

# Define the updated connectivity function
new_func = '''void Pyramid13::connectivity(const unsigned int sc,
                           const IOPackage iop,
                           std::vector< dof_id_type > & conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        // TODO
        libmesh_not_implemented();
      }

    case VTK:
      {
        switch (sc)
          {
          case 0:
            // libMesh and VTK node orderings match perfectly for Pyramid13.
            conn.resize(13);
            conn[0]  = this->node_id(0);
            conn[1]  = this->node_id(1);
            conn[2]  = this->node_id(2);
            conn[3]  = this->node_id(3);
            conn[4]  = this->node_id(4);
            conn[5]  = this->node_id(5);
            conn[6]  = this->node_id(6);
            conn[7]  = this->node_id(7);
            conn[8]  = this->node_id(8);
            conn[9]  = this->node_id(9);
            conn[10] = this->node_id(10);
            conn[11] = this->node_id(11);
            conn[12] = this->node_id(12);
            return;

          default:
            libmesh_error_msg(\"Unknown sc = \" << sc);
          }
      }

    default:
      libmesh_error_msg(\"Unsupported IO package \" << iop);
    }
}'''

# Regex to find and replace the void Pyramid13::connectivity method definition
pattern = r'void\s+Pyramid13::connectivity\s*\([\s\S]*?^\}'

if not re.search(pattern, content, re.MULTILINE):
    print('Error: Could not locate Pyramid13::connectivity function in target file.', file=sys.stderr)
    sys.exit(1)

new_content = re.sub(pattern, new_func, content, flags=re.MULTILINE)

with open(file_path, 'w') as f:
    f.write(new_content)

print(f'Successfully patched Pyramid13::connectivity in {file_path}')
