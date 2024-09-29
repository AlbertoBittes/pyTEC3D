import vtk
import numpy as np


def parse_node_constraints(file_path):
    node_data = {}
    with open(file_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            parts = line.split()
            if len(parts) == 5 and parts[0].isdigit():
                node = int(parts[0])
                label = parts[1]
                real_value = float(parts[2])
                if node not in node_data:
                    node_data[node] = [np.nan, np.nan, np.nan]
                if label == "AX":
                    node_data[node][0] = real_value
                elif label == "AY":
                    node_data[node][1] = real_value
                elif label == "AZ":
                    node_data[node][2] = real_value
    return node_data


def update_vtk_grid(vtk_file_path, node_data):
    # Read the VTK file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtk_file_path)
    reader.Update()
    grid = reader.GetOutput()

    # Get the "A" and "ansys_node_num" arrays
    a_array = grid.GetPointData().GetArray("A")
    ansys_node_num_array = grid.GetPointData().GetArray("ansys_node_num")

    # Update the "A" array based on the node data
    for i in range(grid.GetNumberOfPoints()):
        node_num = int(ansys_node_num_array.GetValue(i))
        if node_num in node_data:
            a_values = node_data[node_num]
            current_values = a_array.GetTuple(i)
            updated_values = [
                a_values[j] if not np.isnan(a_values[j]) else current_values[j]
                for j in range(3)
            ]
            a_array.SetTuple(i, updated_values)

    # Write the updated grid to a new file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(vtk_file_path)
    writer.SetInputData(grid)
    writer.SetDataModeToBinary()
    writer.Write()


if __name__ == "__main__":
    node_data = parse_node_constraints(
        r"D:\PlasmaArc\simulations\3D_loop\3D_loop_bcs_1.txt"
    )
    update_vtk_grid(r"D:\PlasmaArc\simulations\3D_loop\3D_loop.vtu", node_data)
