#ifndef APEC_M_PERIODIC_BOUNDARY_FACTORY_H
#define APEC_M_PERIODIC_BOUNDARY_FACTORY_H

#include "apec_m_periodic_boundary.h"
#include "fe_data_object.h"
#include "fe_parent_object.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/equation_systems.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <variant>

using namespace libMesh;

// ------------------------------------------------------
// Factory to build ApecMPeriodicBoundary from vtkFieldData
// ------------------------------------------------------
class ApecMPeriodicBoundaryFactory
    : public vtkFieldDataFEDataObject<
          std::unordered_map<unsigned int,
                             std::vector<std::variant<
                                 int, Real, std::string,
                                 std::vector<std::vector<int>>,
                                 std::vector<std::vector<Real>>,
                                 std::vector<std::vector<std::string>>>>>>,
      public FEParentObjectBase<EquationSystems>
{
public:
    explicit ApecMPeriodicBoundaryFactory(EquationSystems *eq_sys = nullptr);

    void populate_intermediate_map_from_vtk(vtkFieldData *field_data) override;
    void populate_final_map() override;
    // void broadcast(const libMesh::Parallel::Communicator &comm) override;

    const std::vector<ApecMPeriodicBoundary> &get_periodic_boundaries() const;

private:
    std::vector<ApecMPeriodicBoundary> periodic_boundaries;
};

#endif // APEC_M_PERIODIC_BOUNDARY_FACTORY_H