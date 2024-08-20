import vtk
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid
from vtkmodules.vtkCommonExecutionModel import vtkInformation, vtkInformationVector


class TimeDependentAlgorithm(VTKPythonAlgorithmBase):
    def __init__(self, time_steps):
        super().__init__(
            nInputPorts=1,
            inputType="vtkUnstructuredGrid",
            nOutputPorts=1,
            outputType="vtkUnstructuredGrid",
        )
        self.time_steps = time_steps

    def RequestInformation(self, request, inInfo, outInfo):
        # Set up time steps
        outInfo.GetInformationObject(0).Set(
            vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(),
            self.time_steps,
            len(self.time_steps),
        )
        outInfo.GetInformationObject(0).Set(
            vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(),
            [self.time_steps[0], self.time_steps[-1]],
            2,
        )
        return 1

    def RequestData(self, request, inInfo, outInfo):
        input_data = vtkUnstructuredGrid.GetData(inInfo[0])
        output_data = vtkUnstructuredGrid.GetData(outInfo)

        # Get the current time step
        current_time = request.Get(
            vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP()
        )

        # Process the input data to produce time-dependent output
        # For simplicity, we just copy the input to the output
        output_data.ShallowCopy(input_data)

        # Here you can modify the output_data based on the current_time if needed

        return 1


# Example usage
time_steps = [0.0, 1.0, 2.0, 3.0]  # Define your time steps here
algorithm = TimeDependentAlgorithm(time_steps)
