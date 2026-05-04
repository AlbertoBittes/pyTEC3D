import vtk
import numpy as np

# Global variable list, index is the number
VARIABLES_GLOBAL_LIST = [
    "A",
    "Az",
    "VOLT",
    "EMF",
    "B",
    "H",
    "Be",
    "He",
    "Jt",
    "AMPS",
    "Js",
]

# Default solver settings
DEFAULT_SOLVER_SETTINGS = {
    "ksp_type": "kspreonly",
    "pc_type": "lu",
    "initial_linear_tolerance": 1e-6,
    "max_linear_iterations": 100,
    "minimum_linear_tolerance": 1e-10,
    "linear_tolerance_multiplier": 1.0,
    "write_full_matrix": 0,
    "timesolver_tolerance": 1e-6,
    "analytic_jacobians": 0,
    "verify_analytic_jacobians": 0,
    "numerical_jacobian_h": 1e-6,
    "quiet": 0,
    "max_nonlinear_it": 25,
    "min_steplength": 1e-6,
    "relative_step_tolerance": 1e-6,
    "require_residual_reduction": 0,
    "reduce_deltat_on_diffsolver_failure": 1,
    "continue_after_max_iterations": 0,
    "continue_after_backtrack_failure": 0,
    "timesolver_core": 1,
    "print_solution_norms": 0,
    "print_solutions": 0,
    "print_residual_norms": 0,
    "print_residuals": 0,
    "print_jacobian_norms": 0,
    "print_jacobians": 0,
    "print_element_solutions": 0,
    "print_element_residuals": 0,
    "print_element_jacobians": 0,
    "constrain_in_solver": 0,
    "use_petsc_snes": 0,
    "min_step_length": 1e-6,
}


# Add a dictionary with the name of the setting and its type
DEFAULT_SOLVER_SETTINGS_TYPES = {
    "ksp_type": str,
    "pc_type": str,
    "initial_linear_tolerance": float,
    "max_linear_iterations": int,
    "minimum_linear_tolerance": float,
    "linear_tolerance_multiplier": float,
    "write_full_matrix": int,
    "timesolver_tolerance": float,
    "analytic_jacobians": int,
    "verify_analytic_jacobians": int,
    "numerical_jacobian_h": float,
    "quiet": int,
    "max_nonlinear_it": int,
    "min_steplength": float,
    "relative_step_tolerance": float,
    "require_residual_reduction": int,
    "reduce_deltat_on_diffsolver_failure": int,
    "continue_after_max_iterations": int,
    "continue_after_backtrack_failure": int,
    "timesolver_core": int,
    "print_solution_norms": int,
    "print_solutions": int,
    "print_residual_norms": int,
    "print_residuals": int,
    "print_jacobian_norms": int,
    "print_jacobians": int,
    "print_element_solutions": int,
    "print_element_residuals": int,
    "print_element_jacobians": int,
    "constrain_in_solver": int,
    "use_petsc_snes": int,
    "min_step_length": float,
}

# Default time parameters
DEFAULT_TIME_PARAMETERS = {
    "time": 0.0,
    "transient_effects": True,
    "theta": 1,
}

# Add a dictionary with the name of the time parameter and its type
DEFAULT_TIME_PARAMETERS_TYPES = {
    "time": float,
    "transient_effects": bool,
    "theta": float,
}

# Default APEC algorithm parameters/status
DEFAULT_APEC_ALGORITHM_INPUT_PARAMETERS = {
    "target_central_axis_x": 0.0,
    "target_central_axis_y": 0.0,
    "target_central_axis_z": 0.0,
    "R0": 0.0,
    "max_picard_iterations": 1,
    "min_distance": 0.3,
    "cp_gain": 0.0,
    "tolerance_value": 1e-6,
    "tolerance": 1e-6,
    "bradx_value": 0.0,
    "brady_value": 0.0,
    "bverx_value": 0.0,
    "bvery_value": 0.0,
    "alpha_mp": 1.0001,
    "alpha_np": 0.83197,
    "plasma_current": 0.0,
    "beta_p": 0.0,
    "AF0": 5.3 * 6.2,
    "print_full_output": 0,
    "equilibrium_mesh": "",
    "static_mesh": "",
    "evolution_mesh": "",
    "n_slices": -1,
}

DEFAULT_APEC_ALGORITHM_STATUS = {
    "psi_magnetic_axis": 0.0,
    "psi_max": 0.0,
    "psi_min": 0.0,
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

# Assign types to materials columns
MATERIALS_COLUMNS_TYPES = {
    "MATERIAL_ID": int,
    "MUR_X": float,
    "MUR_Y": float,
    "MUR_Z": float,
    "RHO_X": float,
    "RHO_Y": float,
    "RHO_Z": float,
    "M_X": float,
    "M_Y": float,
    "M_Z": float,
}

# Assign types to DOF coupling columns
DOF_COUPLINGS_COLUMNS_TYPES = {
    "CP_ID": int,
    "VARIABLE": str,
}

# Assign types to real constants columns
CONSTANTS_COLUMNS_TYPES = {
    "CONSTANTS_ID": int,
    "THK": float,
    "CS": float,
    "TURNS": float,
    "VOL": float,
    "RAD": float,
    "T_X": float,
    "T_Y": float,
    "T_Z": float,
    "R": float,
    "S": float,
    "V_X": float,
    "V_Y": float,
    "V_Z": float,
    "OMEGA_X": float,
    "OMEGA_Y": float,
    "OMEGA_Z": float,
    "LOC_X": float,
    "LOC_Y": float,
    "LOC_Z": float,
}

# Assign types to element types columns
ELEMENT_TYPES_COLUMNS_TYPES = {
    "TYPE_ID": int,
    "ELEMENT_TYPE": int,
    "OPTION1": int,
    "OPTION2": int,
    "OPTION3": int,
    "OPTION4": int,
    "OPTION5": int,
}


# Assign types to loads and boundary conditions columns
CELL_POINT_DATA_COLUMNS_TYPES = {
    "VTK_ID": int,
    "LABEL": str,
    "COMPONENT": int,
    "VALUE": float,
}

DEFAULT_COLUMNS_ATTRIBUTES_TYPES = {
    "element_types": ELEMENT_TYPES_COLUMNS_TYPES,
    "materials": MATERIALS_COLUMNS_TYPES,
    "constants": CONSTANTS_COLUMNS_TYPES,
    "solver_settings": DEFAULT_SOLVER_SETTINGS_TYPES,
    "time_step": DEFAULT_TIME_PARAMETERS_TYPES,
    "node": CELL_POINT_DATA_COLUMNS_TYPES,
    "elem": CELL_POINT_DATA_COLUMNS_TYPES,
    "dof_couplings": DOF_COUPLINGS_COLUMNS_TYPES,
}

# Define the mapping of sheet names to processing methods
SHEET_PROCESSING_MAP = {
    "element_types": "read_excel_df",
    "materials": "read_excel_df",
    "constants": "read_excel_df",
    "solver_settings": "read_excel_dict",
    "time_step": "read_excel_dict",
    "apec_settings": "read_excel_dict",
    "node": ("update_vtk_grid_from_excel", "node"),
    "elem": ("update_vtk_grid_from_excel", "elem"),
    "dof_couplings": "read_excel_df",
}

VTK_ARRAYS_TO_INIT = {
    "angles": {
        "location": "node",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "A": {
        "location": "node",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "Az": {
        "location": "node",
        "components": 1,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "VOLT": {
        "location": "node",
        "components": 1,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "EMF": {
        "location": "node",
        "components": 1,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "_A": {
        "location": "node",
        "components": 3,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "_Az": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "_VOLT": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "_EMF": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "B": {
        "location": "node",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "H": {
        "location": "node",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "Be": {
        "location": "elem",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "He": {
        "location": "elem",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "Jt": {
        "location": "elem",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "AMPS": {
        "location": "node",
        "components": 1,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "Js": {
        "location": "elem",
        "components": 3,
        "fill": 0.0,
        "dtype": vtk.VTK_DOUBLE,
    },
    "CP_0": {
        "location": "node",
        "components": 1,
        "fill": 0,
        "dtype": vtk.VTK_INT,
    },
}
