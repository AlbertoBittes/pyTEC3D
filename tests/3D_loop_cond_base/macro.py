# This example shows how to setup a 3D transient analysis for a passive
# conductive 3D current loop on top an identical active current source.
# Only one-fourth of the coils and surrounded air is modeled.
# The simulation includes far-field, symmetry, and voltage boundary conditions,
# as well as transient analysis with imposed source currents (Js) on the active coil elements.
#
# element types:
# - 21: 3D solid element for the surrounding air (ELEMENT TYPE = 1 (3D nodal), OPTION1 = 0 (Air), OPTION2 = 1 (1st Order))
# - 23: 3D solid element for the active coil elements (ELEMENT TYPE = 1 (3D nodal), OPTION1 = 3 (Stranded Conductor), OPTION2 = 1 (1st Order))
# - 24: 3D solid element for the passive coil elements (ELEMENT TYPE = 1 (3D nodal), OPTION1 = 4 (Solid Conductor), OPTION2 = 1 (1st Order))
#
# materials:
# - 21: mu_r (x,y,z) = 1, rho (x,y,z) = - for the surrounding air (not used)
# - 23: mu_r (x,y,z) = 1, rho (x,y,z) = - for the active coil elements (not used)
# - 24: mu_r (x,y,z) = 1, rho (x,y,z) = 1e-9 [Ohm m] for the passive coil elements


# 0. Import necessary modules
from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
import os

# 1. Initialize the program instance
pytec3d = pyTEC3DSolver()

# 2. Define the file name and path for the model and read it into the program
filename = "3D_loop_cond_base"
filepath = os.path.join(os.getcwd(), filename + ".vtu")
pytec3d.read_model(filepath)

# 3. Read the Excel file for new element types, material properties, real constants, solver settings etc.
pytec3d.model.read_json("model_data.json")

# 4. Save the mesh and the model data read from the Excel file to .vtu
pytec3d.save_model()

# 5. Load the model into the solver module
pytec3d.load_model_into_solver()

# 6. Boundary conditions
# Far-field boundary conditions: set magnetic vector potential A to zero
Selector.select(pytec3d.model.mesh, "Z_BOUNDS_1", 1)
Selector.select(pytec3d.model.mesh, "Z_BOUNDS_2", 1, mode="a")
Selector.select(pytec3d.model.mesh, "Y_BOUNDS_2", 1, mode="a")
Selector.select(pytec3d.model.mesh, "X_BOUNDS_2", 1, mode="a")
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# Symmetry boundary conditions: set parallel magnetic vector potential A to zero on the symmetry planes
Selector.select(pytec3d.model.mesh, "X_BOUNDS_1", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select(pytec3d.model.mesh, "Y_BOUNDS_1", 1)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 0)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# Voltage boundary conditions: set the electric potential V on the voltage boundaries
Selector.select(pytec3d.model.mesh, "VOLT_0_1", 1)
Selector.select(pytec3d.model.mesh, "VOLT_0_2", 1, mode="a")
Selector.set_value_on_selection(pytec3d.model.mesh, "_V", 0)

# sqrt(2)/2 factor for the corner element of the current loop
b = 2**0.5 / 2

# 6. Transient analysis: loop over the number of desired time-steps
for i in range(1, 3):
    # Set the current time for this time-step
    pytec3d.time_step["time"] = i * 0.0001
    # Set the transient effects flag
    # initial time-step is static to establish the initial conditions
    if i == 1:
        pytec3d.time_step["transient_effects"] = 0
    # subsequent time-steps have transient effects
    else:
        pytec3d.time_step["transient_effects"] = 1

    # 7. Loads: impose a fixed source current Js on the active coil elements
    Selector.select(pytec3d.model.mesh, "COIL_NOCOND_1", 1, type="elem")
    Selector.set_value_on_selection(pytec3d.model.mesh, "Js", 1 / i, 1, type="elem")
    Selector.select(pytec3d.model.mesh, "COIL_NOCOND_2", 1, type="elem")
    Selector.set_value_on_selection(pytec3d.model.mesh, "Js", -1 / i, 0, type="elem")
    Selector.select(pytec3d.model.mesh, "COIL_NOCOND_3", 1, type="elem")
    Selector.set_value_on_selection(pytec3d.model.mesh, "Js", -b / i, 0, type="elem")
    Selector.set_value_on_selection(pytec3d.model.mesh, "Js", b / i, 1, type="elem")

    # 9. Solve the system for this time-step
    pytec3d.solve()

    # 10. Write the output for this time-step
    pytec3d.write_output()
