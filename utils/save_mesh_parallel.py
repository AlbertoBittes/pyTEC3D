import vtk
import os


def save_mesh_parallel():
    # Ask the user for the file path
    input_filepath = input("Enter the file path of the .vtu file: \n")

    # Check if the file exists
    if not os.path.isfile(input_filepath):
        print("Error: The specified file does not exist.")
        return

    # Ask the user for the number of partitions
    try:
        n = int(input("Enter the number of partitions: \n"))
        if n <= 0:
            raise ValueError("Number of partitions must be greater than 0.")
    except ValueError as e:
        print(f"Invalid input: {e}")
        return

    # Read the unstructured grid from the file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_filepath)
    reader.Update()

    # Check if the input file is valid
    if not reader.GetOutput():
        print("Error: Unable to read the input file.")
        return

    # Get the unstructured grid data
    unstructured_grid = reader.GetOutput()

    # Set the output directory to be a new folder with the same basename as the input file
    input_basename = os.path.splitext(os.path.basename(input_filepath))[0]
    output_dir = os.path.join(os.path.dirname(input_filepath), input_basename)

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Set the output file path with the same basename but with a .pvtu extension
    output_filepath = os.path.join(output_dir, f"{input_basename}.pvtu")

    # Write the partitioned data to .pvtu format
    writer = vtk.vtkXMLPUnstructuredGridWriter()
    writer.SetInputData(unstructured_grid)
    writer.SetNumberOfPieces(n)
    writer.SetFileName(output_filepath)
    writer.SetStartPiece(0)
    writer.SetEndPiece(n - 1)
    writer.SetDataModeToBinary()  # Set the data mode to binary
    writer.Write()

    print(f"Partitioned mesh saved to {output_dir}")


# Run the function
save_mesh_parallel()
