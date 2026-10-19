# =====================================================================================================================
# coil_2D_avolt - 2D axisymmetric massive coil with a terminal, driven first by its current then by its voltage, with a passive ring
# =====================================================================================================================
#
# What this example shows
# -----------------------
# A 2D axisymmetric transient eddy-current analysis of two identical rings of square cross-section. The lower ring is
# an active CIRCUIT_COUPLED_CONDUCTOR: a massive conductor with a terminal. Its current is not imposed point by point
# but follows from the electric field, and the ring is driven through its terminal, i.e. through one electric DOF
# shared by the whole cross-section (V, the voltage across its out-of-plane length).
# The upper ring is a passive SOLID_CONDUCTOR (massive, no terminal) in which only eddy currents flow. It shows:
#   * how to read a mesh and its model data (element types, materials, real constants, DOF couplings);
#   * how to build a boundary condition from a selection (A = 0 on the external boundary);
#   * how to merge the V DOFs of a conductor into a single terminal DOF with a DOF coupling (CP_0 array);
#   * how to switch the excitation of a conductor between an imposed total current (nodal load AMPS on the master
#     node) and an imposed voltage (Dirichlet value _V) from one time step to the next;
#   * the transient solution with the theta method (Crank-Nicolson, theta = 0.5) and the output files.
#
# Geometry and mesh (coil_2D_avolt.vtu)
# -------------------------------------
#   2D axisymmetric about the y axis: x is the radius r, y the axial coordinate z. Domain 0 <= r <= 1 m,
#   -1 <= z <= 1 m, 800 QUAD8 (second-order) cells of 0.05 m.
#   * Lower ring (TYPE / MATERIAL 23): 0.2 m x 0.2 m cross-section centred at r = 0.5 m, z = -0.5 m (16 cells).
#   * Upper ring (TYPE / MATERIAL 24): the same cross-section centred at r = 0.5 m, z = +0.5 m (16 cells).
#   * Air (TYPE / MATERIAL 21): the rest of the domain.
#
# Model data (model_data.json)
# ----------------------------
#   Element types: TYPE id -> element_type / option1 (domain) / option2 (FE order) / option3 (axisymmetric)
#   * 21: 0 = AV2D_nodal / 0 = NON_CONDUCTING            / 2 / 1   air: only A_z, no current
#   * 23: 0 = AV2D_nodal / 4 = CIRCUIT_COUPLED_CONDUCTOR / 2 / 1   active ring: A_z + V (terminal voltage)
#   * 24: 0 = AV2D_nodal / 5 = SOLID_CONDUCTOR           / 2 / 1   passive ring: A_z only; in 2D a closed, shorted
#                                                                  massive conductor carrying eddy currents
#   Materials (MATERIAL id): mur relative permeability, rho resistivity [Ohm m] (3x3 tensors, here isotropic)
#   * 21: mur = 1, rho = 1      (rho is not used: NON_CONDUCTING elements have no conductivity term)
#   * 23: mur = 1, rho = 1e-7
#   * 24: mur = 1, rho = 1e-8
#   Real constants (CONSTANT id)
#   * 1: thickness = 1: out-of-plane length L = thickness * 2 pi r (axisymmetric), i.e. the full 360 degree ring.
#   DOF couplings (CP id): 1 -> ["V"]: all the nodes with CP_0 = 1 share one V value.
#
# Boundary conditions
# -------------------
#   * A_z = 0 on the external boundary: on the axis (r = 0, where A_theta vanishes by symmetry) and on the outer
#     boundary (truncation of the open domain: flux parallel to it).
#   * Step 2 only: V = 10000 V on the coupled terminal of the lower ring (Dirichlet value _V).
#
# Loads and time history (theta = 0.5, dt = 1e-5 s, starting from zero fields)
# ----------------------------------------------------------------------------
#   * Step 1, t = 1e-5 s: total current AMPS = 100 A in the lower ring (current driven). The passive upper ring
#     reacts with eddy currents opposing the change of flux.
#   * Step 2, t = 2e-5 s: AMPS = 0 and V = 10 kV (voltage driven): the current of the lower ring now follows from
#     its resistance and from the induced voltage.
#
# Results (written in ./output)
# -----------------------------
#   coil_2D_avolt_out_1.vtu, coil_2D_avolt_out_2.vtu: A (A_z in component 2), V, B and H (nodal and element), Jt
#   (total current density, element, component 2 = toroidal); coil_2D_avolt_steps.csv and coil_2D_avolt.log: the
#   run report (model tables, constraints, solver statistics, timings).
#
# How to run:  cd tests/coil_2D_avolt && python macro.py
# Documentation: docs, "Examples" (2D axisymmetric coils), "Element library" (AV2D_nodal) and "Theory"
# (2D circuit-coupled conductor, theta method).
# =====================================================================================================================


# ---------------------------------------------------------------------------------------------------------------------
# 0. Imports
#    Solver: the analysis (model + solver session); selection: APDL-like selection of nodes and elements.
# ---------------------------------------------------------------------------------------------------------------------
from pytec3d import selection, Solver
import os
import vtk.util.numpy_support as vtknp

# ---------------------------------------------------------------------------------------------------------------------
# 1. Create the solver
#    One Solver per analysis. It owns the model (mesh + model data), the time step and the solver settings, and starts
#    the job log <mesh_name>.log in the working directory when the model is read.
# ---------------------------------------------------------------------------------------------------------------------
pytec3d = Solver()

# ---------------------------------------------------------------------------------------------------------------------
# 2. Read the model
#    The mesh is a VTK unstructured grid (.vtu). read_model() names the model after the file, adds the arrays the
#    element types need (A, _A, V, _V, AMPS, Js, ...) and clears the selection.
#    read_json() then replaces the model data (element types, materials, real constants, DOF couplings, nonlinear
#    laws, coordinate systems) by the content of model_data.json. The cells pick their rows through the cell arrays
#    TYPE, MATERIAL, CONSTANT and COORDINATE_SYSTEM already stored in the mesh.
# ---------------------------------------------------------------------------------------------------------------------
filename = "coil_2D_avolt"

filepath = os.path.join(os.getcwd(), filename + ".vtu")

pytec3d.read_model(filepath)

pytec3d.model.read_json("model_data.json")

# ---------------------------------------------------------------------------------------------------------------------
# 3. Boundary condition: A_z = 0 on the external boundary
#    Dirichlet values are written in the point array "_<field>": NaN means free, a number fixes that component.
#    select_all() selects every node and cell; nsel_ext() replaces the node selection by the nodes on the external
#    boundary of the selected cells. A is a 3-component array; in 2D only its component 2 (A_z, out of plane) is a DOF.
# ---------------------------------------------------------------------------------------------------------------------
selection.select_all(pytec3d.model.mesh)
selection.nsel_ext(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "_A", 0, 2)

# ---------------------------------------------------------------------------------------------------------------------
# 4. Terminal of the circuit-coupled conductor: couple its V DOFs
#    In 2D the V of a circuit-coupled conductor is the voltage across its out-of-plane length, a single value for the whole
#    cross-section. Select the cells of the ring (location="elem"), then their nodes (nsle), and give them the CP id 1
#    in the point array CP_0: dof_couplings["1"] = ["V"] ties all their V DOFs to the one of the master node,
#    the lowest node id of the group. Loads and Dirichlet values of the terminal are set on that master node.
# ---------------------------------------------------------------------------------------------------------------------
selection.select(pytec3d.model.mesh, "MATERIAL", 23, location="elem")
selection.nsle(pytec3d.model.mesh)
selection.set_value_on_selection(pytec3d.model.mesh, "CP_0", 1)
# Take the first selected node as master node (lowest ID)
master_node = vtknp.vtk_to_numpy(pytec3d.model.mesh.GetPointData().GetArray("selection_point")).nonzero()[0][0]

# ---------------------------------------------------------------------------------------------------------------------
# 5. Load the model into the solver
#    Validates the model (arrays, ids, model data against the element library) and creates the C++ solver session.
#    From now on the mesh and the model data are locked: only array VALUES (loads, Dirichlet values) may change.
# ---------------------------------------------------------------------------------------------------------------------
pytec3d.load_model_into_solver()

# ---------------------------------------------------------------------------------------------------------------------
# 6. Time integration: theta method with theta = 0.5 (Crank-Nicolson, second-order accurate)
#    time_step["time"] is the END time of each step; dt is the difference with the previous one (the first step starts
#    at t = 0 from zero fields). transient_effects stays True (default): the eddy-current terms are active.
#    select_all() only resets the selection; the loop selects the master node again before each change.
# ---------------------------------------------------------------------------------------------------------------------
selection.select_all(pytec3d.model.mesh)
pytec3d.time_step["theta"] = 0.5

# ---------------------------------------------------------------------------------------------------------------------
# 7. Time loop: step 1 current driven, step 2 voltage driven
#    * AMPS: nodal load of the V equation of a circuit-coupled conductor = total current through its cross-section [A].
#    * _V: Dirichlet value of the terminal voltage [V]. Fixing it changes the set of constrained DOFs, so the solver
#      rebuilds and refactorizes the system matrix at step 2 (it reuses it while only values change).
#    solve() writes the solution back into the mesh arrays; write_output() saves output/coil_2D_avolt_out_<n>.vtu.
# ---------------------------------------------------------------------------------------------------------------------
for i in range(1, 3):
    pytec3d.time_step["time"] = i * 0.00001
    if i == 1:
        # Current driven: 100 A through the lower ring, its voltage is computed.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "AMPS", 100)
    else:
        # Voltage driven: remove the current load and impose 10 kV; the current is computed.
        selection.select(pytec3d.model.mesh, "ID", master_node)
        selection.set_value_on_selection(pytec3d.model.mesh, "AMPS", 0)
        selection.set_value_on_selection(pytec3d.model.mesh, "_V", 10000)
    pytec3d.solve()
    pytec3d.write_output()

# ---------------------------------------------------------------------------------------------------------------------
# 8. Save the model (mesh, last solution and model data) to <mesh_name>.vtu in the working directory
#    Run from the case folder, this overwrites the input coil_2D_avolt.vtu with the solved model (the model data is
#    stored inside it, so the file can be read back with read_model() alone).
# ---------------------------------------------------------------------------------------------------------------------
pytec3d.save_model()
