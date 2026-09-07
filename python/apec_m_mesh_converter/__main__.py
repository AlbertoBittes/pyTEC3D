"""CLI entry point for ``python -m apec_m_mesh_converter``."""

from __future__ import annotations

import sys

from .mesh_converter import convert_mesh_from_cdb


def main() -> int:
    """Run the mesh converter CLI."""
    if len(sys.argv) == 2:
        input_file = sys.argv[1]
    elif len(sys.argv) == 1:
        input_file = input("Enter the path of the file to be converted: ").strip()
    else:
        print("Usage: python -m apec_m_mesh_converter <path_to_mesh_file>")
        return 1

    try:
        convert_mesh_from_cdb(input_file)
    except Exception as exc:
        print(f"Mesh conversion failed: {exc}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
