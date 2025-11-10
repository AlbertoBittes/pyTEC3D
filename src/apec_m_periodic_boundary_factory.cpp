#include "apec_m_periodic_boundary_factory.h"
#include "globals.h"
#include "libmesh/dof_map.h"
#include "libmesh/libmesh.h"

#include <map>
#include <set>

namespace
{
const std::vector<int> &extract_int_array(const Globals::NDArray &array)
{
    if (!std::holds_alternative<std::vector<int>>(array.values))
        libmesh_error_msg("Periodic boundary factory expected integer array.");
    return std::get<std::vector<int>>(array.values);
}

const std::vector<Real> &extract_real_array(const Globals::NDArray &array)
{
    if (!std::holds_alternative<std::vector<Real>>(array.values))
        libmesh_error_msg("Periodic boundary factory expected real array.");
    return std::get<std::vector<Real>>(array.values);
}

const std::vector<std::string> &extract_string_array(const Globals::NDArray &array)
{
    if (!std::holds_alternative<std::vector<std::string>>(array.values))
        libmesh_error_msg("Periodic boundary factory expected string array.");
    return std::get<std::vector<std::string>>(array.values);
}
}

ApecMPeriodicBoundaryFactory::ApecMPeriodicBoundaryFactory(std::ifstream &file, EquationSystems *eq_sys)
    : BinaryDataFEDataObject(file),
      FEParentObjectBase<EquationSystems>(eq_sys)
{
    if (eq_sys)
        this->initialize(eq_sys);
}

const std::vector<ApecMPeriodicBoundary> &ApecMPeriodicBoundaryFactory::get_periodic_boundaries() const
{
    return periodic_boundaries;
}

void ApecMPeriodicBoundaryFactory::populate_periodic_boundaries_from_map()
{
    periodic_boundaries.clear();

    if (empty() || this->get_data_map()->columns.empty())
        return;

    if (!has_column("main_boundary") || !has_column("paired_boundary") ||
        !has_column("periodic_matrix") || !has_column("variables"))
        libmesh_error_msg("Periodic boundary factory: required keys missing from settings data.");

    const auto &main_array = get_column("main_boundary");
    const auto &paired_array = get_column("paired_boundary");
    const auto &matrix_array = get_column("periodic_matrix");
    const auto &variables_array = get_column("variables");

    const auto &main_values = extract_int_array(main_array);
    const auto &paired_values = extract_int_array(paired_array);
    const auto &matrix_values = extract_real_array(matrix_array);
    const auto &variable_values = extract_string_array(variables_array);

    if (main_values.size() != paired_values.size())
        libmesh_error_msg("Periodic boundary factory: boundary array size mismatch.");

    std::size_t count = main_values.size();
    if (count == 0)
        return;

    std::size_t matrix_stride = 0;
    if (matrix_array.shape.size() < 2)
        libmesh_error_msg("Periodic boundary factory: periodic matrix must have at least 2 dimensions.");

    if (matrix_array.shape[0] != count)
        libmesh_error_msg("Periodic boundary factory: periodic matrix first dimension mismatch.");

    std::size_t rows = matrix_array.shape[matrix_array.shape.size() - 2];
    std::size_t cols = matrix_array.shape.back();
    matrix_stride = rows * cols;

    if (matrix_values.size() != count * matrix_stride)
        libmesh_error_msg("Periodic boundary factory: periodic matrix flattened size mismatch.");

    std::size_t variable_stride = 1;
    if (!variables_array.shape.empty())
    {
        if (variables_array.shape[0] != count)
            libmesh_error_msg("Periodic boundary factory: variables first dimension mismatch.");
        variable_stride = 1;
        for (std::size_t dim = 1; dim < variables_array.shape.size(); ++dim)
            variable_stride *= variables_array.shape[dim];
    }
    if (variable_values.size() != count * variable_stride)
        libmesh_error_msg("Periodic boundary factory: variables flattened size mismatch.");

    periodic_boundaries.reserve(count);

    for (std::size_t i = 0; i < count; ++i)
    {
        boundary_id_type main_boundary = static_cast<boundary_id_type>(main_values[i]);
        boundary_id_type paired_boundary = static_cast<boundary_id_type>(paired_values[i]);

        std::vector<std::vector<Real>> periodic_matrix(rows, std::vector<Real>(cols));
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t c = 0; c < cols; ++c)
            {
                std::size_t flat = i * matrix_stride + r * cols + c;
                periodic_matrix[r][c] = matrix_values[flat];
            }

        std::vector<std::string> variable_labels;
        variable_labels.reserve(variable_stride);
        for (std::size_t v = 0; v < variable_stride; ++v)
        {
            std::size_t flat = i * variable_stride + v;
            if (flat < variable_values.size())
                variable_labels.push_back(variable_values[flat]);
        }
        if (variable_labels.empty())
            variable_labels.emplace_back("A");

        periodic_boundaries.emplace_back(
            main_boundary, paired_boundary, periodic_matrix,
            PeriodicBoundaryBase::FORWARD,
            variable_labels,
            this->get_parent_object());
    }

    this->compress_periodic_boundaries();
}

void ApecMPeriodicBoundaryFactory::compress_periodic_boundaries()
{
    if (periodic_boundaries.empty())
        return;

    EquationSystems *eq_sys = this->get_parent_object();
    if (!eq_sys)
        return;

    std::map<std::pair<boundary_id_type, boundary_id_type>, std::vector<ApecMPeriodicBoundary>> grouped;
    for (const auto &pb : periodic_boundaries)
        grouped[{pb.myboundary, pb.pairedboundary}].push_back(pb);

    std::vector<ApecMPeriodicBoundary> compressed;
    compressed.reserve(grouped.size());

    for (auto &entry : grouped)
    {
        auto &group = entry.second;
        const ApecMPeriodicBoundary &base = group.front();

        const TensorValue<Real> &base_rotation = base.get_rotation_matrix();
        const VectorValue<Real> &base_translation = base.get_translation_vector();

        std::set<unsigned int> combined_vars = base.get_variables();
        std::set<std::string> combined_labels(base.get_variable_labels().begin(), base.get_variable_labels().end());

        for (std::size_t idx = 1; idx < group.size(); ++idx)
        {
            const auto &candidate = group[idx];

            const TensorValue<Real> &rot_candidate = candidate.get_rotation_matrix();
            const VectorValue<Real> &trans_candidate = candidate.get_translation_vector();

            for (unsigned int r = 0; r < 3; ++r)
                for (unsigned int c = 0; c < 3; ++c)
                    if (base_rotation(r, c) != rot_candidate(r, c))
                        libmesh_error_msg("Periodic boundary rotation mismatch while compressing.");

            for (unsigned int comp = 0; comp < 3; ++comp)
                if (base_translation(comp) != trans_candidate(comp))
                    libmesh_error_msg("Periodic boundary translation mismatch while compressing.");

            const auto &candidate_vars = candidate.get_variables();
            for (unsigned int var : candidate_vars)
            {
                if (!combined_vars.insert(var).second)
                    libmesh_error_msg("Periodic boundary variables overlap detected while compressing.");
            }

            const auto &candidate_labels = candidate.get_variable_labels();
            combined_labels.insert(candidate_labels.begin(), candidate_labels.end());
        }

        std::vector<std::vector<Real>> periodic_matrix = {
            {base_rotation(0, 0), base_rotation(0, 1), base_rotation(0, 2), base_translation(0)},
            {base_rotation(1, 0), base_rotation(1, 1), base_rotation(1, 2), base_translation(1)},
            {base_rotation(2, 0), base_rotation(2, 1), base_rotation(2, 2), base_translation(2)},
            {0.0, 0.0, 0.0, 1.0}};

        std::vector<std::string> merged_labels(combined_labels.begin(), combined_labels.end());

        ApecMPeriodicBoundary merged_boundary(entry.first.first,
                                              entry.first.second,
                                              periodic_matrix,
                                              PeriodicBoundaryBase::FORWARD,
                                              merged_labels,
                                              eq_sys);

        compressed.push_back(std::move(merged_boundary));
    }

    periodic_boundaries.swap(compressed);

    // for (auto &pb : periodic_boundaries)
    // {
    //     if (pb.has_transformation_matrix())
    //     {
    //         const DenseMatrix<Real> &matrix = pb.get_transformation_matrix();
    //         if (matrix.m() == 1 && matrix.n() == 1)
    //             pb.clear_transformation_matrix();
    //     }
    // }
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

void ApecMPeriodicBoundaryFactory::print_data_info() const
{
    if (empty() || this->get_data_map()->columns.empty())
    {
        libMesh::out << "  Periodic boundaries: none configured." << std::endl;
        return;
    }

    if (!has_column("main_boundary") || !has_column("paired_boundary"))
    {
        libMesh::out << "  Periodic boundary data present but required columns are missing." << std::endl;
        BinaryDataFEDataObject::print_data_info();
        return;
    }

    const auto &main_array = get_column("main_boundary");
    const auto &paired_array = get_column("paired_boundary");
    const auto &matrix_array = get_column("periodic_matrix");

    const auto &main_values = extract_int_array(main_array);
    const auto &paired_values = extract_int_array(paired_array);
    const auto &matrix_values = extract_real_array(matrix_array);

    const std::size_t count = std::min(main_values.size(), paired_values.size());
    libMesh::out << "  Periodic boundary pairs: " << count << std::endl;
    if (!count)
        return;

    std::size_t rows = 0;
    std::size_t cols = 0;
    if (matrix_array.shape.size() >= 2)
    {
        rows = matrix_array.shape[matrix_array.shape.size() - 2];
        cols = matrix_array.shape.back();
    }
    const std::size_t matrix_stride = rows && cols ? rows * cols : 0;

    const std::vector<std::string> *variables_ptr = nullptr;
    std::size_t variable_stride = 0;
    if (has_column("variables"))
    {
        const auto &variables_array = get_column("variables");
        variables_ptr = &extract_string_array(variables_array);
        if (!variables_array.shape.empty())
        {
            variable_stride = 1;
            for (std::size_t dim = 1; dim < variables_array.shape.size(); ++dim)
                variable_stride *= variables_array.shape[dim];
        }
        else
            variable_stride = variables_ptr->empty() ? 0 : 1;
    }

    for (std::size_t i = 0; i < count; ++i)
    {
        libMesh::out << "    " << main_values[i] << " <-> " << paired_values[i];
        if (matrix_stride && matrix_values.size() >= (i + 1) * matrix_stride)
        {
            libMesh::out << " translation=(";
            for (std::size_t r = 0; r < rows; ++r)
            {
                if (r)
                    libMesh::out << ", ";
                const std::size_t idx = i * matrix_stride + r * cols + (cols ? cols - 1 : 0);
                if (idx < matrix_values.size())
                    libMesh::out << matrix_values[idx];
                else
                    libMesh::out << "?";
            }
            libMesh::out << ")";
        }

        if (variables_ptr && variable_stride)
        {
            const std::size_t start = i * variable_stride;
            const std::size_t end = std::min(start + variable_stride, variables_ptr->size());
            if (start < end)
            {
                libMesh::out << " variables=[";
                for (std::size_t idx = start; idx < end; ++idx)
                {
                    if (idx > start)
                        libMesh::out << ", ";
                    libMesh::out << (*variables_ptr)[idx];
                }
                libMesh::out << "]";
            }
        }

        libMesh::out << std::endl;
    }
}
