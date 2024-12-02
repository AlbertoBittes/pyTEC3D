import os
import shutil
import subprocess
import sys

def get_shared_libraries(executable):
    """Get the list of shared libraries for the given executable using ldd."""
    try:
        result = subprocess.run(['ldd', executable], capture_output=True, text=True, check=True)
        lines = result.stdout.splitlines()
        libraries = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 3 and '=>' in parts:
                lib_path = parts[2]
                if os.path.isfile(lib_path):
                    libraries.append(lib_path)
        return libraries
    except subprocess.CalledProcessError as e:
        print(f"Error running ldd: {e}", file=sys.stderr)
        return []

def copy_libraries(libraries, dest_folder):
    """Copy the list of libraries to the destination folder."""
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
    
    for lib in libraries:
        try:
            shutil.copy(lib, dest_folder)
            print(f"Copied {lib} to {dest_folder}")
        except Exception as e:
            print(f"Error copying {lib}: {e}", file=sys.stderr)

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <executable> <destination_folder>", file=sys.stderr)
        sys.exit(1)
    
    executable = sys.argv[1]
    dest_folder = sys.argv[2]
    
    libraries = get_shared_libraries(executable)
    if libraries:
        copy_libraries(libraries, dest_folder)
    else:
        print("No shared libraries found or error occurred.", file=sys.stderr)

if __name__ == "__main__":
    main()