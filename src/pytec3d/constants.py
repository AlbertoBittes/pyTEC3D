import vtk
import numpy as np

# Default solver settings
DEFAULT_SOLVER_SETTINGS = {
    "ksp_type": "preonly",
    "pc_type": "lu",
    "initial_linear_tolerance": 1e-6,
    "max_linear_iterations": 100,
    "minimum_linear_tolerance": 1e-10,
}


# Default time parameters: theta for first-order formulations, Newmark beta and gamma (average acceleration)
# for second-order ones (the element types' "time_order").
DEFAULT_TIME_PARAMETERS = {
    "time": 1.0,
    "transient_effects": True,
    "theta": 1.0,
    "beta": 0.25,
    "gamma": 0.5,
}

# Default pyTEC3D algorithm parameters/status
DEFAULT_pyTEC3D_ALGORITHM_INPUT_PARAMETERS = {
    "target_central_axis_x": 0.0,
    "target_central_axis_y": 0.0,
    "target_central_axis_z": 0.0,
    "R0": 0.0,
    "max_picard_iterations": 1,
    "tolerance_value": 1e-6,
    "bradx_value": 0.0,
    "brady_value": 0.0,
    "bverx_value": 0.0,
    "bvery_value": 0.0,
    "alpha_mp": 1.0001,
    "alpha_np": 0.83197,
    "plasma_current": 0.0,
    "beta_p": 0.0,
    "AF0": 5.3 * 6.2,
    "equilibrium_mesh": "",
    "static_mesh": "",
    "evolution_mesh": "",
    "extrusion_angle": 360.0,
    "n_toroidal_divisions": -1,
    "halo_delta1": 0.0,
    "halo_delta2": 0.0,
    "t0_halo": 0.0,
    # Generic control-current mapping.
    # Per-coil correction current is:
    #   I_corr = radial_gain * radial_correction_current
    #          + vertical_gain * vertical_correction_current
    # and the applied current density increment is dJ = I_corr / CS.
    # Defaults preserve legacy behavior:
    #   PFC_2 += radial, PFC_3 += vertical,
    #   PFC_4 += vertical, PFC_5 -= radial.
    "control_coils": [
        {"component": "PFC_2", "radial_gain": 1.0, "vertical_gain": 0.0},
        {"component": "PFC_3", "radial_gain": 0.0, "vertical_gain": 1.0},
        {"component": "PFC_4", "radial_gain": 0.0, "vertical_gain": 1.0},
        {"component": "PFC_5", "radial_gain": -1.0, "vertical_gain": 0.0},
    ],
}

DEFAULT_pyTEC3D_ALGORITHM_STATUS = {
    "psi_magnetic_axis": 0.0,
    "psi_max": 0.0,
    "psi_min": 0.0,
    "psi_2_min": 0.0,
    "magnetic_axis_node_id": -1,
    "limiter_node_id": -1,
    "x_point_node_id": -1,
    "magnetic_axis": (0.0, 0.0, 0.0),
    "b_target_axis": (0.0, 0.0),
    "last_psi_limiter": 0.0,
    "last_psi_x_point": 0.0,
    "separatrix_configuration": "uninitialized",
    "limiter_point": (np.nan, np.nan, np.nan),
    "x_point": (np.nan, np.nan, np.nan),
    "flux_norm": np.nan,
    "flux_tolerance": np.nan,
    "radial_correction_current": 0.0,
    "vertical_correction_current": 0.0,
    "total_plasma_current_before_scaling": 0.0,
    "plasma_current_scaling_lambda": 0.0,
    "plasma_element_count": 0,
    "picard_iteration": 0,
    "converged": False,
}

# The grid arrays of a model follow from the element library (pytec3d.element_library.grid_arrays).
