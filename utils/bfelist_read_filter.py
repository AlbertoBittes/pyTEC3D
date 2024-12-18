import vtk
import re


def parse_element_data(file_path):
    # Dictionary to store current densities for each element
    current_densities = {}

    # Regular expression patterns
    element_pattern = re.compile(r"ELEMENT=\s+(\d+)")
    density_pattern = re.compile(
        r"\s*([\d\.\-E]+)\s+([\d\.\-E]+)\s+([\d\.\-E]+)\s+([\d\.\-E]+)"
    )

    with open(file_path, "r") as file:
        current_element = None
        for line in file:
            # Check for element line
            element_match = element_pattern.search(line)
            if element_match:
                current_element = int(element_match.group(1))
                current_densities[current_element] = {
                    "JSX": None,
                    "JSY": None,
                    "JSZ": None,
                }
                continue

            # Check for density line (only the first one after the element line)
            if (
                current_element is not None
                and current_densities[current_element]["JSX"] is None
            ):
                density_match = density_pattern.match(line)
                if density_match:
                    current_densities[current_element]["JSX"] = float(
                        density_match.group(1)
                    )
                    current_densities[current_element]["JSY"] = float(
                        density_match.group(2)
                    )
                    current_densities[current_element]["JSZ"] = float(
                        density_match.group(3)
                    )

    return current_densities


def update_vtk_grid(vtk_file_path, element_data, time_step):
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

    # Get the "Js" and "ansys_elem_num" arrays
    js_array = grid.GetCellData().GetArray("Js_" + str(time_step))
    ansys_elem_num_array = grid.GetCellData().GetArray("ansys_elem_num")

    # Update the "Js" array based on the element data
    for i in range(grid.GetNumberOfCells()):
        elem_num = int(ansys_elem_num_array.GetValue(i))
        if elem_num in element_data:
            js_values = element_data[elem_num]
            if spatial_dimension == 2:
                js_tuple = (js_values["JSZ"],)  # Only JSZ exists in 2D
            else:
                js_tuple = (js_values["JSX"], js_values["JSY"], js_values["JSZ"])
            js_array.SetTuple(i, js_tuple)

    # Write the updated grid to a new file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(vtk_file_path)
    writer.SetInputData(grid)
    writer.Write()


if __name__ == "__main__":
    filename = "3D_loop_cond_loads"
    vtk_grid = "3D_loop_cond.vtu"
    filetype = "out"
    time_steps = 5
    for i in range(1, time_steps + 1):
        element_data = parse_element_data(f"{filename}_{i}." + filetype)
        update_vtk_grid(vtk_grid, element_data, i)
