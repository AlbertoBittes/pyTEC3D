import sys, re

file_path = sys.argv[1]

with open(file_path, 'r') as f:
    content = f.read()

# Target specifically the second-order pyramid block using its unique 13-node assertion anchor
pattern = r'(libmesh_assert_less\s*\(\s*i\s*,\s*13\s*\)\s*;\s*const\s+Real\s+xi\s*=\s*p\(0\)\s*;\s*const\s+Real\s+eta\s*=\s*p\(1\)\s*;\s*)const\s+Real\s+zeta\s*=\s*p\(2\)\s*;'

replacement = r'\1const Real zeta = (p(2) > 1.0 - 1.e-6) ? (1.0 - 1.e-6) : p(2);'

if not re.search(pattern, content):
    print('Error: Could not locate second-order pyramid (PYRAMID13) block in target file.', file=sys.stderr)
    sys.exit(1)

new_content, count = re.subn(pattern, replacement, content)

with open(file_path, 'w') as f:
    f.write(new_content)

print(f'Successfully patched {count} instance(s) in second-order pyramid block of {file_path}')
