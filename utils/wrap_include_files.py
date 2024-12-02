import os
import re
import shutil
import sys

def find_includes(source_file):
    """Find all #include directives in a source file."""
    includes = []
    with open(source_file, 'r') as f:
        for line in f:
            match = re.match(r'#include\s*["<](.*)[">]', line)
            if match:
                includes.append(match.group(1))
    return includes

def copy_headers(source_dir, include_dir, include_paths):
    """Copy header files to the include directory."""
    if not os.path.exists(include_dir):
        os.makedirs(include_dir)

    for root, _, files in os.walk(source_dir):
        for file in files:
            if file.endswith('.cpp') or file.endswith('.h'):
                source_file = os.path.join(root, file)
                includes = find_includes(source_file)
                for include in includes:
                    for path in include_paths:
                        header_path = os.path.join(path, include)
                        if os.path.exists(header_path):
                            dest_path = os.path.join(include_dir, include)
                            dest_dir = os.path.dirname(dest_path)
                            if not os.path.exists(dest_dir):
                                os.makedirs(dest_dir)
                            shutil.copy(header_path, dest_path)
                            print(f"Copied {header_path} to {dest_path}")
                            break

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <source_dir> <include_dir> <include_paths>", file=sys.stderr)
        sys.exit(1)

    source_dir = sys.argv[1]
    include_dir = sys.argv[2]
    include_paths = sys.argv[3].split(':')

    copy_headers(source_dir, include_dir, include_paths)

if __name__ == "__main__":
    main()