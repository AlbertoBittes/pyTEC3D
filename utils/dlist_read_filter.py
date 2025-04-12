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
                if label.endswith("X"):
                    label = label[:-1] + "_x"
                elif label.endswith("Y"):
                    label = label[:-1] + "_y"
                elif label.endswith("Z"):
                    label = label[:-1] + "_z"
                real = float(node_match.group(3))
                imag = float(node_match.group(4))

                if node_num not in node_data:
                    node_data[node_num] = {}

                node_data[node_num][label] = (real, imag)

    return node_data


def update_vtk_grid(self, node_data):
    # Read the VTK file

    # Get the "A" and "ansys_node_num" arrays
    a_array = self.mesh.GetPointData().GetArray("A")
    ansys_node_num_array = self.mesh.GetPointData().GetArray("ansys_node_num")
    v_array = self.mesh.GetPointData().GetArray("V")

    # Update the "A" array based on the node data

    for i in range(self.mesh.GetNumberOfPoints()):
        node_num = int(ansys_node_num_array.GetValue(i))
        if node_num in node_data.keys():

            a_values = node_data[node_num]
            # Initialize with current values
            current_tuple = a_array.GetTuple(i)
            a_tuple = [
                a_values.get("A_x", [current_tuple[0]])[0],
                a_values.get("A_y", [current_tuple[1]])[0],
                a_values.get("A_z", [current_tuple[2]])[0],
            ]
            a_array.SetTuple(i, a_tuple)
            if "VOLT" in node_data[node_num]:
                volt_value = node_data[node_num]["VOLT"][0]
                v_array.SetValue(i, volt_value)


if __name__ == "__main__":
    filename = "3D_loop_cond_bc"
    vtk_grid = "3D_loop_cond.vtu"
    filetype = "out"
    time_steps = 5
    for i in range(1, time_steps + 1):
        node_data = parse_node_data(f"{filename}_{i}." + filetype)
        update_vtk_grid(vtk_grid, node_data, i)
