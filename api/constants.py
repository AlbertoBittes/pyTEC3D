import vtk
import numpy as np

# Global variable list, index is the number
VARIABLES_GLOBAL_LIST = [
    "A",
    "Az",
    "VOLT",
    "V",
    "EMF",
    "E",
    "AMPS",
    "Js",
    "B",
    "H",
    "PSI",
    "Be",
    "He",
    "PSIe",
    "Jt",
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
    "time_step_number": 1,
    "time": 0.0,
    "number_of_substeps": 1,
    "transient_effects": 1,
    "ramped_loads": 1,
    "theta": 1,
    "write_substep_results": 0,
}

# Add a dictionary with the name of the time parameter and its type
DEFAULT_TIME_PARAMETERS_TYPES = {
    "time_step_number": int,
    "time": float,
    "number_of_substeps": int,
    "transient_effects": int,
    "ramped_loads": int,
    "theta": float,
    "write_substep_results": int,
}

DEFAULT_APEC_SETTINGS = {
    "target_central_axis_x": 0.0,
    "target_central_axis_y": 0.0,
    "target_central_axis_z": 0.0,
    "target_central_axis_element_id": 0,
    "R0": 0.0,
    "max_picard_iterations": 1,
    "min_distance": 0.3,
    "cp_gain": 0.0,
    "tolerance_value": 1e-6,
    "bradx_value": 0.0,
    "brady_value": 0.0,
    "bverx_value": 0.0,
    "bvery_value": 0.0,
    "alpha_value": 1.0,
    "plasma_current": 0.0,
    "beta_p": 0.0,
    "print_full_output": 0,
    "main_mesh": "",
    "static_mesh": "",
    "n_slices": 3,
    "toroidal_slices_nodes": np.zeros((0, 0), dtype=np.int64),
    "toroidal_slices_elems": np.zeros((0, 0), dtype=np.int64),
}

REQUIRED_APEC_SETTINGS = [
    "R0",
    "max_picard_iterations",
    "tolerance_value",
    "alpha_value",
    "plasma_current",
    "beta_p",
]

DEFAULT_APEC_SETTINGS_TYPES = {
    "target_central_axis_x": float,
    "target_central_axis_y": float,
    "target_central_axis_z": float,
    "target_central_axis_element_id": int,
    "R0": float,
    "max_picard_iterations": int,
    "min_distance": float,
    "cp_gain": float,
    "tolerance_value": float,
    "bradx_value": float,
    "brady_value": float,
    "bverx_value": float,
    "bvery_value": float,
    "alpha_value": float,
    "plasma_current": float,
    "beta_p": float,
    "print_full_output": int,
    "main_mesh": str,
    "toroidal_slices_nodes": np.ndarray,
    "toroidal_slices_elems": np.ndarray,
}
# Backwards compatibility with earlier naming
APEC_SETTINGS_TYPES = DEFAULT_APEC_SETTINGS_TYPES

# Assign types to materials columns
MATERIALS_COLUMNS_TYPES = {
    "ansys_material_type": int,
    "MURX": float,
    "MURY": float,
    "MURZ": float,
    "RSVX": float,
    "RSVY": float,
    "RSVZ": float,
    "PERX": float,
    "PERY": float,
    "PERZ": float,
    "MGXX": float,
    "MGYY": float,
    "MGZZ": float,
}

# Assign types to DOF coupling columns
DOF_COUPLINGS_COLUMNS_TYPES = {
    "cp_id": int,
    "variable": str,
}

# Assign types to real constants columns
REAL_CONSTANTS_COLUMNS_TYPES = {
    "ansys_real_constant": int,
    "THK": float,
    "SC": float,
    "NC": float,
    "VC": float,
    "RAD": float,
    "TX": float,
    "TY": float,
    "TZ": float,
    "R": float,
    "SYM": float,
    "VELOX": float,
    "VELOY": float,
    "VELOZ": float,
    "OMEGAX": float,
    "OMEGAY": float,
    "OMEGAZ": float,
    "XLOC": float,
    "YLOC": float,
    "ZLOC": float,
    "CPSIGN": float,
    "CURR": float,
}

# Assign types to element types columns
ELEMENT_TYPES_COLUMNS_TYPES = {
    "ansys_etype": int,
    "element_type": str,
    "KEYOPT1": int,
    "KEYOPT2": int,
    "KEYOPT3": int,
    "KEYOPT4": int,
    "KEYOPT5": int,
}


# Assign types to loads and boundary conditions columns
CELL_POINT_DATA_COLUMNS_TYPES = {
    "vtk_id": int,
    "label": str,
    "component": int,
    "value": float,
}

DEFAULT_COLUMNS_ATTRIBUTES_TYPES = {
    "element_types": ELEMENT_TYPES_COLUMNS_TYPES,
    "materials": MATERIALS_COLUMNS_TYPES,
    "real_constants": REAL_CONSTANTS_COLUMNS_TYPES,
    "solver_settings": DEFAULT_SOLVER_SETTINGS_TYPES,
    "time_parameters": DEFAULT_TIME_PARAMETERS_TYPES,
    "node": CELL_POINT_DATA_COLUMNS_TYPES,
    "elem": CELL_POINT_DATA_COLUMNS_TYPES,
    "apec_settings": DEFAULT_APEC_SETTINGS_TYPES,
    "dof_couplings": DOF_COUPLINGS_COLUMNS_TYPES,
}

# Define the mapping of sheet names to processing methods
SHEET_PROCESSING_MAP = {
    "element_types": "read_excel_df",
    "materials": "read_excel_df",
    "real_constants": "read_excel_df",
    "solver_settings": "read_excel_dict",
    "time_parameters": "read_excel_dict",
    "apec_settings": "read_excel_dict",
    "node": ("update_vtk_grid_from_excel", "node"),
    "elem": ("update_vtk_grid_from_excel", "elem"),
    "dof_couplings": "read_excel_df",
}

WRITING_FUNCTIONS_MAP = {
    "element_types": "_write_to_binary",
    "time_parameters": "_write_to_binary",
    "solver_settings": "_write_to_binary",
    "materials": "_write_to_binary",
    "real_constants": "_write_to_binary",
    "node": "_write_numpy_array_to_binary",
    "elem": "_write_numpy_array_to_binary",
    "dof_couplings": "_write_to_binary",
    "apec_settings": "_write_to_binary",
}

WRITING_LABELS = list(WRITING_FUNCTIONS_MAP.keys())
WRITING_LABEL_TO_ID = {label: idx for idx, label in enumerate(WRITING_LABELS)}

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
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "Az": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "VOLT": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "V": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "EMF": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "E": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "AMPS": {
        "location": "node",
        "components": 1,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "Js": {
        "location": "elem",
        "components": 3,
        "fill": np.nan,
        "dtype": vtk.VTK_DOUBLE,
    },
    "CP_0": {
        "location": "node",
        "components": 1,
        "fill": 0,
        "dtype": vtk.VTK_INT,
    },
}

ARRAYS_TO_REMOVE = [
    "selection_point",
    "selection_cell",
    "cell_centroids",
    "cell_volume",
]

# Labels to read for each load step from Excel
EXCEL_LOAD_STEP_LABELS = [
    "node",
    "elem",
]
