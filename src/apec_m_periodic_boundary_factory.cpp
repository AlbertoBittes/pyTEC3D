#include "apecm_periodic_boundary_factory.h"
#include "libmesh/libmesh_logging.h"
#include <vtkAbstractArray.h>
#include <vtkStringArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>

void ApecMPeriodicBoundaryFactory::populate_data_map_from_vtk(vtkFieldData *field_data)
{
    if (!field_data)
        libmesh_error_msg("No vtkFieldData provided to populate ApecMPeriodicBoundaryFactory.");

    vtkIntArray *main_boundary_arr = vtkIntArray::SafeDownCast(field_data->GetArray("main_boundary"));
    vtkIntArray *paired_boundary_arr = vtkIntArray::SafeDownCast(field_data->GetArray("paired_boundary"));
    vtkDoubleArray *periodic_matrix_arr = vtkDoubleArray::SafeDownCast(field_data->GetArray("periodic_matrix"));
    vtkDoubleArray *variables_transformation_arr = vtkDoubleArray::SafeDownCast(field_data->GetArray("variables_transformation_matrix"));
    vtkStringArray *variables_arr = vtkStringArray::SafeDownCast(field_data->GetAbstractArray("variables"));

    vtkDoubleArray *periodic_matrix_shape = vtkDoubleArray::SafeDownCast(field_data->GetArray("periodic_matrix_shape"));
    vtkDoubleArray *variables_transformation_shape = vtkDoubleArray::SafeDownCast(field_data->GetArray("variables_transformation_matrix_shape"));

    if (!main_boundary_arr || !paired_boundary_arr || !periodic_matrix_arr || !variables_transformation_arr || !variables_arr)
        libmesh_error_msg("Missing one or more required vtkFieldData arrays for periodic boundaries.");

    vtkIdType n_tuples = main_boundary_arr->GetNumberOfTuples();
    if (paired_boundary_arr->GetNumberOfTuples() != n_tuples)
        libmesh_error_msg("Mismatch in tuple count between main_boundary and paired_boundary.");

    // --- Parse and populate data_map ---
    data_map.clear();
    data_map.reserve(n_tuples);

    vtkIdType matrix_index = 0;
    vtkIdType transform_index = 0;

    for (vtkIdType i = 0; i < n_tuples; ++i)
    {
        ApecMPeriodicBoundaryData entry;
        entry.main_boundary = static_cast<unsigned int>(main_boundary_arr->GetValue(i));
        entry.paired_boundary = static_cast<unsigned int>(paired_boundary_arr->GetValue(i));

        // Parse periodic_matrix (assumes 2D matrices concatenated row-wise)
        if (periodic_matrix_shape)
        {
            // Each entry: ndim, dim0, dim1, ...
            int ndim = static_cast<int>(periodic_matrix_shape->GetValue(3 * i));
            int rows = static_cast<int>(periodic_matrix_shape->GetValue(3 * i + 1));
            int cols = static_cast<int>(periodic_matrix_shape->GetValue(3 * i + 2));
            std::vector<std::vector<Real>> mat(rows, std::vector<Real>(cols));
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c)
                    mat[r][c] = static_cast<Real>(periodic_matrix_arr->GetValue(matrix_index++));
            entry.periodic_matrix = std::move(mat);
        }

        // Parse variable list (semicolon-separated strings)
        std::string var_str = variables_arr->GetValue(i);
        std::stringstream ss(var_str);
        std::string token;
        while (std::getline(ss, token, ';'))
            entry.variables.push_back(token);

        // Parse variables_transformation_matrix (2D)
        if (variables_transformation_shape)
        {
            int ndim = static_cast<int>(variables_transformation_shape->GetValue(3 * i));
            int rows = static_cast<int>(variables_transformation_shape->GetValue(3 * i + 1));
            int cols = static_cast<int>(variables_transformation_shape->GetValue(3 * i + 2));
            std::vector<std::vector<Real>> mat(rows, std::vector<Real>(cols));
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c)
                    mat[r][c] = static_cast<Real>(variables_transformation_arr->GetValue(transform_index++));
            entry.variables_transformation_matrix = std::move(mat);
        }

        data_map.push_back(std::move(entry));
    }

    libMesh::out << "ApecMPeriodicBoundaryFactory: Parsed " << n_tuples << " periodic boundaries from vtkFieldData." << std::endl;
}

void ApecMPeriodicBoundaryFactory::create_boundaries_from_data()
{
    boundaries.clear();
    for (const auto &entry : data_map)
        boundaries.push_back(std::make_unique<ApecMPeriodicBoundary>(entry));

    libMesh::out << "ApecMPeriodicBoundaryFactory: Created " << boundaries.size() << " ApecMPeriodicBoundary objects." << std::endl;
}