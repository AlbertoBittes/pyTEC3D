import vtk
import re


def parse_node_data(file_path):
    # Dictionary to store node data
    node_data = {}

    # Regular expression patterns
    node_pattern = re.compile(r"^\s*(\d+)\s+(\S+)\s+([\d\.\-E]+)\s+([\d\.\-E]+)")

    with open(file_path, "r") as file:
        for line in file:
            # Skip header lines
            if line.startswith("NODE") or line.startswith("//"):
                continue

            # Check for node data line
            node_match = node_pattern.match(line)
            if node_match:
                node_num = int(node_match.group(1))
                label = node_match.group(2)
                real = float(node_match.group(3))
                imag = float(node_match.group(4))

                if node_num not in node_data:
                    node_data[node_num] = {}

                node_data[node_num][label] = (real, imag)

    return node_data


def update_vtk_grid(vtk_file_path, node_data, time_step):
    # Read the VTK file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtk_file_path)
    reader.Update()
    grid = reader.GetOutput()

    # Get the dimension of the unstructured grid
    bounds = grid.GetBounds()
    if bounds[4] == bounds[5]:  # Check if Z bounds are the same
        spatial_dimension = 2
    else:
        spatial_dimension = 3

    # Get the "A" and "ansys_node_num" arrays
    a_array = grid.GetPointData().GetArray("A_" + str(time_step))
    ansys_node_num_array = grid.GetPointData().GetArray("ansys_node_num")

    # Update the "A" array based on the node data
    for i in range(grid.GetNumberOfPoints()):
        node_num = int(ansys_node_num_array.GetValue(i))
        if node_num in node_data:
            a_values = node_data[node_num]
            if spatial_dimension == 2:
                a_tuple = (a_values["AZ"][0],)  # Only Az exists in 2D
            else:
                # Initialize with current values
                current_tuple = a_array.GetTuple(i)
                a_tuple = [
                    a_values.get("AX", [current_tuple[0]])[0],
                    a_values.get("AY", [current_tuple[1]])[0],
                    a_values.get("AZ", [current_tuple[2]])[0],
                ]
            a_array.SetTuple(i, a_tuple)

    # Write the updated grid to a new file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(vtk_file_path)
    writer.SetInputData(grid)
    writer.Write()


if __name__ == "__main__":
    filename = "3D_loop_cond_bc"
    vtk_grid = "3D_loop_cond.vtu"
    filetype = "out"
    time_steps = 5
    for i in range(1, time_steps + 1):
        node_data = parse_node_data(f"{filename}_{i}." + filetype)
        update_vtk_grid(vtk_grid, node_data, i)
