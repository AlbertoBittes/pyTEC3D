#ifndef APEC_M_PERIODIC_BOUNDARY_FACTORY_H
#define APEC_M_PERIODIC_BOUNDARY_FACTORY_H

#include "apec_m_periodic_boundary.h" // your ApecMPeriodicBoundary class
#include "fe_data_object.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/equation_systems.h"
#include "fe_parent_object.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <tuple>

using namespace libMesh;

class ApecMPeriodicBoundaryFactory
    : public vtkFieldDataFEDataObject<std::unordered_map<unsigned int,
                                                         std::tuple<
                                                             unsigned int,                   // main_boundary
                                                             unsigned int,                   // paired_boundary
                                                             std::vector<std::vector<Real>>, // periodic_matrix
                                                             std::vector<std::string>,       // variables
                                                             std::vector<std::vector<Real>>  // variables_transformation_matrix
                                                             >>>,
      public FEParentObjectBase<EquationSystems>
{
public:
    ApecMPeriodicBoundaryFactory(EquationSystems *eq_sys = nullptr)
        : FEParentObjectBase<EquationSystems>(eq_sys),
          vtkFieldDataFEDataObject<std::unordered_map<unsigned int, std::tuple<unsigned int, unsigned int, std::vector<std::vector<Real>>, std::vector<std::string>, std::vector<std::vector<Real>>>>>() {
          };

    // -------------------------------
    // Stage 1: Populate intermediate map from VTK field_data
    // -------------------------------
    void populate_intermediate_map_from_vtk(vtkFieldData *field_data) override
    {
        if (!field_data)
            libmesh_error_msg("No vtkFieldData provided.");

        // Specify array names in the field_data
        std::vector<std::string> array_names = {
            "main_boundary",
            "paired_boundary",
            "periodic_matrix",
            "variables",
        };
        // "variables_transformation_matrix"};

        // Base class handles the reading, including _shape and _offset for iterables
        this->read_from_vtk_field_data(field_data, array_names);
    }

    // -------------------------------
    // Stage 2: Populate actual ApecMPeriodicBoundary objects
    // -------------------------------
    void populate_final_map() override
    {
        auto *intermediate_map = this->get_data_map();
        if (!intermediate_map)
            libmesh_error_msg("Intermediate map is empty. Call populate_intermediate_map_from_vtk first.");

        periodic_boundaries.clear();

        for (const auto &[row_id, row_tuple] : *intermediate_map)
        {
            boundary_id_type main_boundary = static_cast<boundary_id_type>(std::get<0>(row_tuple));
            boundary_id_type paired_boundary = static_cast<boundary_id_type>(std::get<1>(row_tuple));
            const std::vector<std::vector<Real>> &periodic_matrix = std::get<2>(row_tuple);
            const std::vector<std::string> &variables = std::get<3>(row_tuple);
            // const std::vector<std::vector<Real>> &variables_transformation_matrix = std::get<4>(row_tuple);

            // Construct the ApecMPeriodicBoundary object
            ApecMPeriodicBoundary boundary(main_boundary, paired_boundary, periodic_matrix,
                                           PeriodicBoundaryBase::FORWARD, variables, this->get_parent_object());

            periodic_boundaries.push_back(std::move(boundary));
        }
    }

    // Getter for created periodic boundaries
    inline const std::vector<ApecMPeriodicBoundary> &get_periodic_boundaries() const
    {
        return periodic_boundaries;
    }

private:
    std::vector<ApecMPeriodicBoundary> periodic_boundaries;
};

#endif // APEC_M_PERIODIC_BOUNDARY_FACTORY_H