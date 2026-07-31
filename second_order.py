import sys
import numpy as np
import vtk

# ==============================================================================
# Canonical VTK 1st-to-2nd Order Topology Configurations
# ==============================================================================
# Defines target quadratic cell type, number of corners, and the exact local
# edge node pairs (u, v) matching VTK quadratic node numbering standards.
CELL_CONFIGS = {
    vtk.VTK_LINE: {
        "target_type": vtk.VTK_QUADRATIC_EDGE,
        "num_corners": 2,
        "edges": [(0, 1)],
    },
    vtk.VTK_TRIANGLE: {
        "target_type": vtk.VTK_QUADRATIC_TRIANGLE,
        "num_corners": 3,
        "edges": [(0, 1), (1, 2), (2, 0)],
    },
    vtk.VTK_QUAD: {
        "target_type": vtk.VTK_QUADRATIC_QUAD,
        "num_corners": 4,
        "edges": [(0, 1), (1, 2), (2, 3), (3, 0)],
    },
    vtk.VTK_TETRA: {
        "target_type": vtk.VTK_QUADRATIC_TETRA,
        "num_corners": 4,
        "edges": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    },
    vtk.VTK_HEXAHEDRON: {
        "target_type": vtk.VTK_QUADRATIC_HEXAHEDRON,
        "num_corners": 8,
        "edges": [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),  # Bottom face edges (nodes 8..11)
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),  # Top face edges    (nodes 12..15)
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),  # Vertical edges    (nodes 16..19)
        ],
    },
    vtk.VTK_WEDGE: {
        "target_type": vtk.VTK_QUADRATIC_WEDGE,
        "num_corners": 6,
        "edges": [
            (0, 1),
            (1, 2),
            (2, 0),  # Bottom face edges (nodes 6..8)
            (3, 4),
            (4, 5),
            (5, 3),  # Top face edges    (nodes 9..11)
            (0, 3),
            (1, 4),
            (2, 5),  # Vertical edges    (nodes 12..14)
        ],
    },
    vtk.VTK_PYRAMID: {
        "target_type": vtk.VTK_QUADRATIC_PYRAMID,
        "num_corners": 5,
        "edges": [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),  # Base edges        (nodes 5..8)
            (0, 4),
            (1, 4),
            (2, 4),
            (3, 4),  # Slanted edges     (nodes 9..12)
        ],
    },
}


def fix_hex_node_orientation(
    corner_ids: list, points: vtk.vtkPoints
) -> list:
    """Checks and fixes inverted/twisted node ordering for a 8-node linear

    hexahedron to guarantee positive Jacobian volume in VTK.
    """
    p = [np.array(points.GetPoint(corner_ids[i])) for i in range(8)]

    # Compute scalar triple product of bottom face vectors and vertical vector
    v01 = p[1] - p[0]
    v03 = p[3] - p[0]
    v04 = p[4] - p[0]

    det = np.dot(np.cross(v01, v03), v04)

    # If determinant is negative, the winding order is inverted (left-handed)
    if det < 0:
        # Swap 1 <-> 3 and 5 <-> 7 to restore right-handed orientation
        return [
            corner_ids[0],
            corner_ids[3],
            corner_ids[2],
            corner_ids[1],
            corner_ids[4],
            corner_ids[7],
            corner_ids[6],
            corner_ids[5],
        ]
    return corner_ids


def promote_to_quadratic_grid(
    input_grid: vtk.vtkUnstructuredGrid,
) -> vtk.vtkUnstructuredGrid:
    """Promotes all linear cells in a vtkUnstructuredGrid to their corresponding

    quadratic cell types with correct VTK node ordering.
    """
    num_orig_points = input_grid.GetNumberOfPoints()
    orig_points = input_grid.GetPoints()

    new_points = vtk.vtkPoints()
    new_points.SetDataTypeToDouble()

    for i in range(num_orig_points):
        new_points.InsertNextPoint(orig_points.GetPoint(i))

    # Setup PointData interpolation
    in_pd = input_grid.GetPointData()
    out_grid = vtk.vtkUnstructuredGrid()
    out_pd = out_grid.GetPointData()

    estimated_pts = num_orig_points * 4
    out_pd.InterpolateAllocate(in_pd, estimated_pts)
    out_pd.CopyAllocate(in_pd, estimated_pts)

    for i in range(num_orig_points):
        out_pd.CopyData(in_pd, i, i)

    # Global shared edge cache: (min_id, max_id) -> new_point_id
    edge_to_midpoint = {}

    output_cell_array = vtk.vtkCellArray()
    output_cell_types = vtk.vtkUnsignedCharArray()

    num_cells = input_grid.GetNumberOfCells()

    for cell_idx in range(num_cells):
        cell = input_grid.GetCell(cell_idx)
        cell_type = cell.GetCellType()

        # If cell type is not registered for promotion, keep original
        if cell_type not in CELL_CONFIGS:
            output_cell_types.InsertNextValue(cell_type)
            pt_ids = vtk.vtkIdList()
            cell.GetPointIds(pt_ids)
            output_cell_array.InsertNextCell(pt_ids)
            continue

        config = CELL_CONFIGS[cell_type]
        quad_cell_type = config["target_type"]
        output_cell_types.InsertNextValue(quad_cell_type)

        num_corners = config["num_corners"]
        corner_ids = [cell.GetPointId(i) for i in range(num_corners)]

        # Apply orientation fix for hexahedra if inverted
        if cell_type == vtk.VTK_HEXAHEDRON:
            corner_ids = fix_hex_node_orientation(corner_ids, orig_points)

        # Build quadratic cell node array: [corners..., midsides...]
        quad_cell_ids = list(corner_ids)

        # Generate/Fetch midside nodes in canonical VTK edge order
        for u_local, v_local in config["edges"]:
            u_global = corner_ids[u_local]
            v_global = corner_ids[v_local]

            edge_key = (min(u_global, v_global), max(u_global, v_global))

            if edge_key not in edge_to_midpoint:
                # Calculate physical midpoint
                p_u = np.array(new_points.GetPoint(u_global))
                p_v = np.array(new_points.GetPoint(v_global))
                p_mid = 0.5 * (p_u + p_v)

                # Insert coordinate
                mid_pt_id = new_points.InsertNextPoint(p_mid)
                edge_to_midpoint[edge_key] = mid_pt_id

                # Interpolate field variables
                out_pd.InterpolateEdge(in_pd, mid_pt_id, u_global, v_global, 0.5)

            quad_cell_ids.append(edge_to_midpoint[edge_key])

        output_cell_array.InsertNextCell(len(quad_cell_ids), quad_cell_ids)

    # Set mesh topology and geometry
    out_grid.SetPoints(new_points)
    out_grid.SetCells(output_cell_types, output_cell_array)

    # Pass all CellData arrays (materials, tags, stresses) unchanged
    out_grid.GetCellData().PassData(input_grid.GetCellData())

    return out_grid


# ==============================================================================
# Main IO Execution
# ==============================================================================
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python convert_quadratic.py <input_mesh.vtu> <output_mesh.vtu>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    # Load input mesh
    if input_path.endswith(".vtu"):
        reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        reader = vtk.vtkGenericDataObjectReader()

    reader.SetFileName(input_path)
    reader.Update()
    mesh = reader.GetOutput()

    print(
        f"Input Mesh:  {mesh.GetNumberOfPoints()} points, "
        f"{mesh.GetNumberOfCells()} cells"
    )

    # Execute conversion
    quad_mesh = promote_to_quadratic_grid(mesh)

    print(
        f"Output Mesh: {quad_mesh.GetNumberOfPoints()} points, "
        f"{quad_mesh.GetNumberOfCells()} cells"
    )

    # Write output mesh
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(output_path)
    writer.SetInputData(quad_mesh)
    writer.Write()

    print(f"Successfully converted and saved quadratic mesh to {output_path}")