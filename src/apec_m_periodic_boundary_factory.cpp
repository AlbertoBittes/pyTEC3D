#include "apec_m_periodic_boundary_factory.h"
#include "globals.h"

ApecMPeriodicBoundaryFactory::ApecMPeriodicBoundaryFactory(EquationSystems *eq_sys)
    : FEParentObjectBase<EquationSystems>(eq_sys)
{
}

// --------------------------------------------------------
// Stage 1: Populate intermediate map
// --------------------------------------------------------
void ApecMPeriodicBoundaryFactory::populate_intermediate_map_from_vtk(vtkFieldData *field_data)
{
    if (!field_data)
        libmesh_error_msg("No vtkFieldData provided.");

    std::vector<std::string> array_names = {
        "main_boundary",
        "paired_boundary",
        "periodic_matrix",
        "variables",
        "variables_transformation_matrix"};

    this->read_from_vtk_field_data(field_data, array_names);
}

// --------------------------------------------------------
// Stage 2: Populate final periodic boundary objects
// --------------------------------------------------------
void ApecMPeriodicBoundaryFactory::populate_final_map()
{
    auto *intermediate_map = this->get_data_map();
    if (!intermediate_map)
        libmesh_error_msg("Intermediate map is empty. Call populate_intermediate_map_from_vtk first.");

    periodic_boundaries.clear();

    for (const auto &[row_id, fields] : *intermediate_map)
    {
        if (fields.size() < 4)
        {
            libmesh_warning("Row " << row_id << " incomplete, skipping.");
            continue;
        }

        boundary_id_type main_boundary = static_cast<boundary_id_type>(std::get<int>(fields[0]));
        boundary_id_type paired_boundary = static_cast<boundary_id_type>(std::get<int>(fields[1]));
        const auto &periodic_matrix = std::get<std::vector<std::vector<Real>>>(fields[2]);
        const auto &variables = std::get<std::vector<std::vector<std::string>>>(fields[3])[0];
        std::vector<std::vector<Real>> variables_transformation_matrix;

        if (fields.size() > 4 && std::holds_alternative<std::vector<std::vector<Real>>>(fields[4]))
            variables_transformation_matrix = std::get<std::vector<std::vector<Real>>>(fields[4]);

        ApecMPeriodicBoundary boundary(main_boundary, paired_boundary, periodic_matrix,
                                       PeriodicBoundaryBase::FORWARD,
                                       variables, this->get_parent_object());

        periodic_boundaries.push_back(std::move(boundary));
    }
}

// --------------------------------------------------------
// Broadcast data
// --------------------------------------------------------
// void ApecMPeriodicBoundaryFactory::broadcast(const libMesh::Parallel::Communicator &comm)
// {
//     comm.broadcast<std::unordered_map<unsigned int,
//                                       std::vector<std::variant<
//                                           int, Real, std::string,
//                                           std::vector<std::vector<int>>,
//                                           std::vector<std::vector<Real>>,
//                                           std::vector<std::vector<std::string>>>>>>(*get_data_map());
// }

// --------------------------------------------------------
// Getter
// --------------------------------------------------------
const std::vector<ApecMPeriodicBoundary> &ApecMPeriodicBoundaryFactory::get_periodic_boundaries() const
{
    return periodic_boundaries;
}

template class vtkFieldDataFEDataObject<
    std::unordered_map<unsigned int,
                       std::vector<std::variant<
                           int,
                           Real,
                           std::string,
                           std::vector<std::vector<int>>,
                           std::vector<std::vector<Real>>,
                           std::vector<std::vector<std::string>>>>>>;