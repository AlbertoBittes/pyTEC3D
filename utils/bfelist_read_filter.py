import vtk


def parse_element_data(file_path):
    element_data = {}
    with open(file_path, "r") as file:
        lines = file.readlines()
        element_number = None
        for line in lines:
            if line.startswith(" ELEMENT="):
                element_number = int(line.split("=")[1].strip())
            elif element_number is not None:
                values = list(map(float, line.split()))
                if len(values) == 4:
                    element_data[element_number] = values
                    element_number = None
    return element_data


def update_vtk_grid(vtk_file_path, element_data):
    # Read the VTK file
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtk_file_path)
    reader.Update()
    grid = reader.GetOutput()

    # Get the "Js" and "ansys_elem_num" arrays
    js_array = grid.GetCellData().GetArray("Js")
    ansys_elem_num_array = grid.GetCellData().GetArray("ansys_elem_num")

    # Update the "Js" array based on the element data
    for i in range(grid.GetNumberOfCells()):
        elem_num = int(ansys_elem_num_array.GetValue(i))
        if elem_num in element_data:
            js_values = element_data[elem_num]
            js_array.SetTuple(i, js_values)

    # Write the updated grid to a new file
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName("updated_" + vtk_file_path)
    writer.SetInputData(grid)
    writer.Write()


if __name__ == "__main__":
    element_data = parse_element_data("element_data.txt")
    update_vtk_grid("unstructured_grid.vtk", element_data)
