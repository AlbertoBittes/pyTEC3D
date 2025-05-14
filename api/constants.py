# Global variable mappings
VARIABLES_GLOBAL_NUMBERS = {"A": 0, "VOLT": 1, "Az": 2}
VARIABLES_LOADS_NUMBERS = {"Js": 0}

# Default solver settings
DEFAULT_SOLVER_SETTINGS = {
    "ksp_type": "kspreonly",
    "pc_type": "lu",
    "initial_linear_tolerance": 1e-6,
    "max_linear_iterations": 100,
    "minimum_linear_tolerance": 1e-10,
    "linear_tolerance_multiplier": 1.0,
    "theta": 0.5,
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
    "theta": float,
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
    "write_substep_results": 0,
}

# Add a dictionary with the name of the time parameter and its type
DEFAULT_TIME_PARAMETERS_TYPES = {
    "time_step_number": int,
    "time": float,
    "number_of_substeps": int,
    "transient_effects": int,
    "ramped_loads": int,
    "write_substep_results": int,
}

# Column definitions for materials, real constants, etc.
MATERIALS_COLUMNS = [
    "ansys_material_type",
    "MURX",
    "MURY",
    "MURZ",
    "RSVX",
    "RSVY",
    "RSVZ",
    "PERX",
    "PERY",
    "PERZ",
    "MGXX",
    "MGYY",
    "MGZZ",
    "MURX_B-H",
    "MURY_B-H",
    "MURZ_B-H",
]

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
    "MURX_B-H": int,
    "MURY_B-H": int,
    "MURZ_B-H": int,
}

REAL_CONSTANTS_COLUMNS = [
    "ansys_real_constant",
    "CARE",
    "TURN",
    "LENG",
    "DIRZ",
    "FILL",
    "SC",
    "NC",
    "RAD",
    "TX",
    "TY",
    "TZ",
    "R",
    "SYM",
    "VELOX",
    "VELOY",
    "VELOZ",
    "OMEGAX",
    "OMEGAY",
    "OMEGAZ",
    "XLOC",
    "YLOC",
    "ZLOC",
]

# Assign types to real constants columns
REAL_CONSTANTS_COLUMNS_TYPES = {
    "ansys_real_constant": int,
    "CARE": float,
    "TURN": float,
    "LENG": float,
    "DIRZ": float,
    "FILL": float,
    "SC": float,
    "NC": float,
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
}

ELEMENT_TYPES_COLUMNS = [
    "ansys_etype",
    "element_type",
    "KEYOPT1",
    "KEYOPT2",
    "KEYOPT3",
]

# Assign types to element types columns
ELEMENT_TYPES_COLUMNS_TYPES = {
    "ansys_etype": int,
    "element_type": str,
    "KEYOPT1": int,
    "KEYOPT2": int,
    "KEYOPT3": int,
}

LOADS_BCS_COLUMNS = [
    "vtk_id",
    "label",
    "component",
    "value",
]

# Assign types to loads and boundary conditions columns
LOADS_BCS_COLUMNS_TYPES = {
    "vtk_id": int,
    "label": str,
    "component": int,
    "value": float,
}

TRANSIENT_ANALYSIS_COLUMNS = [
    "time_step_number",
    "time",
    "number_of_substeps",
    "transient_effects",
    "ramped_loads",
    "write_substep_results",
    "load_step_input_file",
]

TRANSIENT_ANALYSIS_COLUMNS_TYPES = {
    "time_step_number": int,
    "time": float,
    "number_of_substeps": int,
    "transient_effects": int,
    "ramped_loads": int,
    "write_substep_results": int,
    "load_step_input_file": str,
}

DEFAULT_COLUMNS_ATTRIBUTES = {
    "element_types": ELEMENT_TYPES_COLUMNS,
    "materials": MATERIALS_COLUMNS,
    "real_constants": REAL_CONSTANTS_COLUMNS,
    "solver_settings": DEFAULT_SOLVER_SETTINGS,
    "time_parameters": DEFAULT_TIME_PARAMETERS,
    "loads": LOADS_BCS_COLUMNS,
    "boundary_conditions": LOADS_BCS_COLUMNS,
    "transient_analysis": TRANSIENT_ANALYSIS_COLUMNS,
}

DEFAULT_COLUMNS_ATTRIBUTES_TYPES = {
    "element_types": ELEMENT_TYPES_COLUMNS_TYPES,
    "materials": MATERIALS_COLUMNS_TYPES,
    "real_constants": REAL_CONSTANTS_COLUMNS_TYPES,
    "solver_settings": DEFAULT_SOLVER_SETTINGS_TYPES,
    "time_parameters": DEFAULT_TIME_PARAMETERS_TYPES,
    "loads": LOADS_BCS_COLUMNS_TYPES,
    "boundary_conditions": LOADS_BCS_COLUMNS_TYPES,
    "transient_analysis": TRANSIENT_ANALYSIS_COLUMNS_TYPES,
}

# Define the mapping of sheet names to processing methods
SHEET_PROCESSING_MAP = {
    "element_types": "read_excel_df",
    "materials": "read_excel_df",
    "real_constants": "read_excel_df",
    "solver_settings": "read_excel_dict",
    "time_parameters": "read_excel_dict",
    "loads": ("update_vtk_grid_from_excel", "cell"),
    "boundary_conditions": ("update_vtk_grid_from_excel", "point"),
    "transient_analysis": "read_excel_df",
}

WRITING_FUNCTIONS_MAP = {
    "time_parameters": "write_settings_to_binary",
    "solver_settings": "write_settings_to_binary",
    "materials": "write_dataframe_to_binary",
    "real_constants": "write_dataframe_to_binary",
    "loads": "write_numpy_array_to_binary",
    "boundary_conditions": "write_numpy_array_to_binary",
}
