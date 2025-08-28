import subprocess
import tempfile
import shutil
import os
import pytest
import numpy as np

import vtk


# ------------------------------------------------------------------------------
# Helper: run the solver once
# ------------------------------------------------------------------------------
def run_solver(input_vtu, nproc, output_dir):
    """Run the solver (serial or parallel) and put results in output_dir."""

    input_vtu = os.path.abspath(input_vtu)
    sim_dir = os.path.dirname(input_vtu)
    analysis_name = os.path.splitext(os.path.basename(input_vtu))[0]

    exe_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../bin/apec_m"))

    cmd = []
    if nproc > 1:
        cmd = ["mpirun", "-np", str(nproc)]

    cmd += [exe_path, "-b", input_vtu]

    print(f"\n[TEST] Running solver: {' '.join(cmd)} in {sim_dir}")
    result = subprocess.run(cmd, cwd=sim_dir, capture_output=True, text=True)
    print("[TEST] Solver stdout:\n", result.stdout)
    print("[TEST] Solver stderr:\n", result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"Solver failed with exit code {result.returncode}")

    # Move outputs from the simulation's output directory
    output_src_dir = os.path.join(sim_dir, "output")
    print(f"[TEST] Moving output files from {output_src_dir} to {output_dir}")
    for f in os.listdir(output_src_dir):
        if f.startswith(analysis_name + "_out_"):
            print(f"[TEST] Moving {f}")
            shutil.move(os.path.join(output_src_dir, f),
                        os.path.join(output_dir, f))


# ------------------------------------------------------------------------------
# Helper: read a specific timestep output and extract value
# ------------------------------------------------------------------------------
def read_value_by_libmesh_id(output_dir, analysis_name, timestep,
                             array_name, component, libmesh_id, is_point_data=True):
    """Read value from a given array/component at the node/cell with given libmesh id."""

    filename = os.path.join(output_dir, f"{analysis_name}_out_{timestep}.pvtu")

    reader = vtk.vtkXMLPUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    ugrid = reader.GetOutput()

    # Select correct id array
    if is_point_data:
        id_array = ugrid.GetPointData().GetArray("libmesh_node_id")
        data = ugrid.GetPointData().GetArray(array_name)
    else:
        id_array = ugrid.GetCellData().GetArray("libmesh_elem_id")
        data = ugrid.GetCellData().GetArray(array_name)

    if id_array is None:
        raise ValueError(f"ID array not found in {filename}")
    if data is None:
        raise ValueError(f"Array '{array_name}' not found in {filename}")

    # Convert VTK array to numpy for efficient masking
    id_np = np.array([id_array.GetTuple1(i) for i in range(id_array.GetNumberOfTuples())], dtype=int)
    mask = (id_np == libmesh_id)
    indices = np.where(mask)[0]
    if len(indices) == 0:
        raise ValueError(f"libmesh_id {libmesh_id} not found in {filename}")
    found_idx = indices[0]

    ncomp = data.GetNumberOfComponents()
    if component >= ncomp:
        raise IndexError(f"Array '{array_name}' has {ncomp} components, "
                         f"cannot access component {component}")

    return data.GetTuple(found_idx)[component]


# ------------------------------------------------------------------------------
# Define analysis cases with multiple checks
# ------------------------------------------------------------------------------
analysis_cases = [
    {
        "analysis_name": "3D_loop_cond_base",
        "nproc": 4,
        "checks": [
            {
                "timestep": 1,
                "array": "Jt",
                "component": 1,
                "libmesh_id": 6254,
                "expected": 0.875,
                "tol": 1e-5,
                "is_point": False,
            },
            {
                "timestep": 2,
                "array": "A",
                "component": 1,
                "libmesh_id": 1644,
                "expected": 1.24741e-11,
                "tol": 1e-5,
                "is_point": True,
            },
            {
                "timestep": 2,
                "array": "VOLT",
                "component": 0,
                "libmesh_id": 3837,
                "expected": -0.45695E-015,
                "tol": 1e-5,
                "is_point": True,
            },
            {
                "timestep": 3,
                "array": "Jt",
                "component": 1,
                "libmesh_id": 8754,
                "expected": 0.10152,
                "tol": 1e-5,
                "is_point": False,
            },
            {
                "timestep": 4,
                "array": "Be",
                "component": 0,
                "libmesh_id": 2530,
                "expected": -0.93074E-010,
                "tol": 1e-5,
                "is_point": False,
            },
        ],
    },
    {
        "analysis_name": "ITER_coil_inf_2D",
        "nproc": 4,
        "checks": [
            # {
            #     "timestep": 2,
            #     "array": "Jt_z",
            #     "component": 0,
            #     "libmesh_id": 331,
            #     "expected": 0.11956E-001,
            #     "tol": 1e-5,
            #     "is_point": False,
            # },
            {
                "timestep": 15,
                "array": "A_z",
                "component": 0,
                "libmesh_id": 152,
                "expected": 0.24491E-009,
                "tol": 1e-5,
                "is_point": True,
            },
            {
                "timestep": 8,
                "array": "A_z",
                "component": 0,
                "libmesh_id": 184,
                "expected": 0.22862E-009,
                "tol": 1e-5,
                "is_point": True,
            },
            # {
            #     "timestep": 3,
            #     "array": "Jt",
            #     "component": 1,
            #     "libmesh_id": 8754,
            #     "expected": 0.10152,
            #     "tol": 1e-5,
            #     "is_point": False,
            # },
            {
                "timestep": 8,
                "array": "Be_x",
                "component": 0,
                "libmesh_id": 234,
                "expected": 0.65660E-008,
                "tol": 1e-5,
                "is_point": False,
            },
            {
                "timestep": 6,
                "array": "Be_x",
                "component": 0,
                "libmesh_id": 265,
                "expected": -0.27643E-009,
                "tol": 1e-5,
                "is_point": False,
            },
        ],
    },
    # Add more cases as needed
]


# ------------------------------------------------------------------------------
# Parametrize analyses
# ------------------------------------------------------------------------------
@pytest.mark.parametrize("case", analysis_cases, ids=lambda c: c["analysis_name"])
def test_fem(case):
    analysis_name = case["analysis_name"]
    nproc = case["nproc"]
    input_vtu = os.path.abspath(
        os.path.join(os.path.dirname(__file__), f"../tests/{analysis_name}/{analysis_name}.vtu")
    )

    print(f"\n[TEST] Starting test for {analysis_name}")

    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"[TEST] Temporary output directory: {tmpdir}")
        run_solver(input_vtu, nproc, tmpdir)

        # perform all checks
        for i, chk in enumerate(case["checks"], 1):
            print(f"[TEST] Check {i}: {chk}")
            val = read_value_by_libmesh_id(
                tmpdir,
                analysis_name,
                chk["timestep"],
                chk["array"],
                chk["component"],
                chk["libmesh_id"],
                chk["is_point"],
            )
            print(f"[TEST] Value read: {val}, expected: {chk['expected']}")
            assert np.isclose(val, chk["expected"], atol=chk["tol"]), (
                f"Mismatch in {analysis_name}, "
                f"timestep={chk['timestep']}, array={chk['array']}, "
                f"comp={chk['component']}, libmesh_id={chk['libmesh_id']}: "
                f"got {val}, expected {chk['expected']}"
            )
        print(f"[TEST] All checks passed for {analysis_name}")