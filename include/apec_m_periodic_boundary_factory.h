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
    : public FEDataObjectMap<vtkFieldData>,
      public FEParentObjectBase<EquationSystems>
{
public:
    explicit ApecMPeriodicBoundaryFactory(std::shared_ptr<vtkFieldData> field_data, EquationSystems *eq_sys = nullptr);
    // NEW: Populate periodic boundaries directly from stored vtkFieldData
    void populate_periodic_boundaries_from_map();

    const std::vector<ApecMPeriodicBoundary> &get_periodic_boundaries() const;

    void apply_periodic_boundaries();
    void initialize(EquationSystems *eq_sys);

private:
    std::vector<ApecMPeriodicBoundary> periodic_boundaries;
};

#endif // APEC_M_PERIODIC_BOUNDARY_FACTORY_H
