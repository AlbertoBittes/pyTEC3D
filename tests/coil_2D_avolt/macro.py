# This example shows how to setup a 2D axisymmetric transient analysis for a passive
# conductive 3D current loop on top an identical active current source.
# The cross-section of the coils is modelled and the axisymmetric assumption is applied.
# The simulation includes far-field, symmetry, and voltage boundary conditions,
# as well as transient analysis with imposed current on a circuit coupled solid conductor
# (AMPS) on the coil master node.
#
# element types:
# - 21: 2D axisymmetric plane element for the surrounding air (ELEMENT TYPE = 1 (2D axisymmetric nodal), OPTION1 = 0 (Air), OPTION2 = 2 (2nd Order))
# - 23: 2D axisymmetric plane element for the active coil elements (ELEMENT TYPE = 1 (2D axisymmetric nodal), OPTION1 = 5 (Circuit Coupled Solid Conductor), OPTION2 = 2 (2nd Order))
# - 24: 2D axisymmetric plane element for the passive coil elements (ELEMENT TYPE = 1 (2D axisymmetric nodal), OPTION1 = 4 (Passive Solid Conductor), OPTION2 = 2 (2nd Order))
#
# materials:
# - 21: mu_r (x,y,z) = 1, rho (x,y,z) = - for the surrounding air (not used)
# - 23: mu_r (x,y,z) = 1, rho (x,y,z) = - for the active coil elements (not used)
# - 24: mu_r (x,y,z) = 1, rho (x,y,z) = 1e-9 [Ohm m] for the passive coil elements



from pytec3d_solver.solver import pyTEC3DSolver
from pytec3d.selector import Selector
import os
import vtk.util.numpy_support as vtknp

pytec3d = pyTEC3DSolver(4)

filename = "coil_2D_avolt"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_excel("input_template.xlsx")

Selector.select_all(pytec3d.model.mesh)
Selector.nsel_ext(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

Selector.select(pytec3d.model.mesh, "MATERIAL", 23, type="elem")
Selector.nsle(pytec3d.model.mesh)
Selector.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

pytec3d.load_model_into_solver()

Selector.select_all(pytec3d.model.mesh)
pytec3d.model.time_step["theta"] = 0.5

for i in range(1, 3):
    pytec3d.model.time_step["time"] = i * 0.00001
    if i == 1:
        Selector.select(pytec3d.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(pytec3d.model.mesh, "CURR", 100)
    else:
        Selector.select(pytec3d.model.mesh, "ID", master_node)
        Selector.set_value_on_selection(pytec3d.model.mesh, "CURR", 0)
        Selector.set_value_on_selection(pytec3d.model.mesh, "_VOLT", 10000)
    pytec3d.solve()
    pytec3d.write_output()

pytec3d.save_model()
