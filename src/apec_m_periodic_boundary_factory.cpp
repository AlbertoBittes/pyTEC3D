#include "apec_m_periodic_boundary_factory.h"
#include "globals.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh.h"

#include <cmath>
#include <vtkAbstractArray.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

ApecMPeriodicBoundaryFactory::ApecMPeriodicBoundaryFactory(std::shared_ptr<vtkFieldData> field_data, EquationSystems *eq_sys)
    : FEDataObjectMap<vtkFieldData>(field_data),
      FEParentObjectBase<EquationSystems>(eq_sys)
{
    if (eq_sys)
        this->initialize(eq_sys);
}

void ApecMPeriodicBoundaryFactory::populate_periodic_boundaries_from_map()
{
    vtkFieldData *field_data = this->get_data_map();
    if (!field_data)
    {
        libMesh::out << "No vtkFieldData stored in the factory. Skipping periodic boundary creation." << std::endl;
        return;
    }

    periodic_boundaries.clear();

    // Arrays for boundaries
    auto *main_boundary_arr = vtkDataArray::SafeDownCast(field_data->GetAbstractArray("main_boundary"));
    auto *paired_boundary_arr = vtkDataArray::SafeDownCast(field_data->GetAbstractArray("paired_boundary"));

    // Arrays for periodic matrix (flattened)
    auto *periodic_matrix_base = field_data->GetAbstractArray("periodic_matrix");
    auto *periodic_matrix_size = vtkIntArray::SafeDownCast(field_data->GetAbstractArray("periodic_matrix_size"));
    auto *periodic_matrix_offset = vtkIntArray::SafeDownCast(field_data->GetAbstractArray("periodic_matrix_offset"));

    // Array for variables (single string per row)
    auto *variables_arr = vtkStringArray::SafeDownCast(field_data->GetAbstractArray("variables"));

    if (!main_boundary_arr || !paired_boundary_arr ||
        !periodic_matrix_base || !periodic_matrix_size || !periodic_matrix_offset ||
        !variables_arr)
    {
        libMesh::out << "Periodic boundary arrays not found in field data. Skipping periodic boundary creation." << std::endl;
        return;
    }

    vtkIdType num_rows = main_boundary_arr->GetNumberOfTuples();
    if (num_rows == 0)
        return;

    vtkIdType paired_rows = paired_boundary_arr->GetNumberOfTuples();
    vtkIdType variables_rows = variables_arr->GetNumberOfTuples();
    vtkIdType size_rows = periodic_matrix_size->GetNumberOfTuples();
    vtkIdType offset_rows = periodic_matrix_offset->GetNumberOfTuples();

    if (paired_rows != num_rows || variables_rows != num_rows ||
        size_rows != num_rows || offset_rows != num_rows)
    {
        libMesh::out << "Mismatch in periodic boundary array sizes. Skipping periodic boundary creation." << std::endl;
        return;
    }

    auto *periodic_matrix_arr_double = vtkDoubleArray::SafeDownCast(periodic_matrix_base);
    auto *periodic_matrix_arr_int = vtkIntArray::SafeDownCast(periodic_matrix_base);
    if (!periodic_matrix_arr_double && !periodic_matrix_arr_int)
    {
        libMesh::out << "Unsupported periodic matrix array type. Skipping periodic boundary creation." << std::endl;
        return;
    }

    for (vtkIdType i = 0; i < num_rows; ++i)
    {
        boundary_id_type main_boundary = static_cast<boundary_id_type>(main_boundary_arr->GetComponent(i, 0));
        boundary_id_type paired_boundary = static_cast<boundary_id_type>(paired_boundary_arr->GetComponent(i, 0));

        // --- PERIODIC MATRIX (flattened, with size/offset) ---
        int pm_size = periodic_matrix_size->GetValue(i);
        if (pm_size <= 0)
            continue;

        int pm_offset = periodic_matrix_offset->GetValue(i);

        std::vector<Real> flattened(pm_size);
        for (int k = 0; k < pm_size; ++k)
        {
            if (periodic_matrix_arr_double)
                flattened[k] = periodic_matrix_arr_double->GetValue(pm_offset + k);
            else
                flattened[k] = static_cast<Real>(periodic_matrix_arr_int->GetValue(pm_offset + k));
        }

        int rows = static_cast<int>(std::sqrt(static_cast<double>(pm_size)));
        int cols = rows > 0 ? pm_size / rows : pm_size;
        if (rows * cols != pm_size || rows == 0)
        {
            rows = 1;
            cols = pm_size;
        }

        std::vector<std::vector<Real>> periodic_matrix(rows, std::vector<Real>(cols));
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c)
                periodic_matrix[r][c] = flattened[r * cols + c];

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

void ApecMPeriodicBoundaryFactory::initialize(EquationSystems *eq_sys)
{
    libMesh::out << "Applying periodic boundaries..." << std::endl;
    this->set_parent_object(eq_sys);
    this->populate_periodic_boundaries_from_map();
    this->apply_periodic_boundaries();
}
