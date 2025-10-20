#include "apec_m_periodic_boundary_factory.h"
#include "globals.h"
#include "libmesh/dof_map.h"

ApecMPeriodicBoundaryFactory::ApecMPeriodicBoundaryFactory(std::shared_ptr<vtkFieldData> field_data, EquationSystems *eq_sys)
    : FEDataObjectMap<vtkFieldData>(field_data),
      FEParentObjectBase<EquationSystems>(eq_sys)
{
    if (eq_sys)
        this->populate_periodic_boundaries_from_map();
}

void ApecMPeriodicBoundaryFactory::populate_periodic_boundaries_from_map()
{
    vtkFieldData *field_data = this->get_data_map();
    if (!field_data)
        libmesh_error_msg("No vtkFieldData stored in the factory.");

    periodic_boundaries.clear();

    // Arrays for boundaries
    auto *main_boundary_arr = field_data->GetAbstractArray("main_boundary");
    auto *paired_boundary_arr = field_data->GetAbstractArray("paired_boundary");

    // Arrays for periodic matrix (flattened)
    auto *periodic_matrix_arr = vtkDoubleArray::SafeDownCast(field_data->GetAbstractArray("periodic_matrix"));
    auto *periodic_matrix_size = vtkIntArray::SafeDownCast(field_data->GetAbstractArray("periodic_matrix_size"));
    auto *periodic_matrix_offset = vtkIntArray::SafeDownCast(field_data->GetAbstractArray("periodic_matrix_offset"));

    // Array for variables (single string per row)
    auto *variables_arr = vtkStringArray::SafeDownCast(field_data->GetAbstractArray("variables"));

    if (!main_boundary_arr || !paired_boundary_arr ||
        !periodic_matrix_arr || !periodic_matrix_size || !periodic_matrix_offset ||
        !variables_arr)
        libmesh_error_msg("Required arrays not found in vtkFieldData.");

    vtkIdType num_rows = main_boundary_arr->GetNumberOfTuples();

    for (vtkIdType i = 0; i < num_rows; ++i)
    {
        boundary_id_type main_boundary = static_cast<boundary_id_type>(main_boundary_arr->GetVariantValue(i).ToInt());
        boundary_id_type paired_boundary = static_cast<boundary_id_type>(paired_boundary_arr->GetVariantValue(i).ToInt());

        // --- PERIODIC MATRIX (flattened, with size/offset) ---
        int pm_size = periodic_matrix_size->GetValue(i); // should be 16 for 4x4
        int pm_offset = periodic_matrix_offset->GetValue(i);

        std::vector<std::vector<Real>> periodic_matrix(4, std::vector<Real>(4));
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                periodic_matrix[r][c] = periodic_matrix_arr->GetValue(pm_offset + r * 4 + c);

        // --- VARIABLES (single string per row) ---
        std::string variable = variables_arr->GetValue(i);

        // Construct boundary
        periodic_boundaries.emplace_back(
            main_boundary, paired_boundary, periodic_matrix,
            PeriodicBoundaryBase::FORWARD,
            variable,
            this->get_parent_object());
    }
}

void ApecMPeriodicBoundaryFactory::apply_periodic_boundaries()
{
    EquationSystems *eq_sys = this->get_parent_object();
    if (!eq_sys)
        libmesh_error_msg("No EquationSystems associated with ApecMPeriodicBoundaryFactory.");

    DofMap &dof_map = eq_sys->get_system(0).get_dof_map();

    for (const auto &pb : periodic_boundaries)
    {
        dof_map.add_periodic_boundary(pb);
    }
}