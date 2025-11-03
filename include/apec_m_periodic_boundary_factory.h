#ifndef APEC_M_PERIODIC_BOUNDARY_FACTORY_H
#define APEC_M_PERIODIC_BOUNDARY_FACTORY_H

#include "apec_m_periodic_boundary.h"
#include "fe_data_object.h"
#include "fe_parent_object.h"
#include "globals.h"
#include "libmesh/equation_systems.h"
#include "libmesh/periodic_boundary_base.h"
#include <memory>
#include <string>
#include <vector>

using namespace libMesh;

// ------------------------------------------------------
// Factory to build ApecMPeriodicBoundary instances from serialized settings data
// ------------------------------------------------------
class ApecMPeriodicBoundaryFactory
    : public BinaryDataFEDataObject,
      public FEParentObjectBase<EquationSystems>
{
public:
    explicit ApecMPeriodicBoundaryFactory(std::shared_ptr<BinaryColumnarData> settings = nullptr, EquationSystems *eq_sys = nullptr);
    // NEW: Populate periodic boundaries directly from stored vtkFieldData
    void populate_periodic_boundaries_from_map();

    const std::vector<ApecMPeriodicBoundary> &get_periodic_boundaries() const;

    void apply_periodic_boundaries();
    void initialize(EquationSystems *eq_sys);
    void compress_periodic_boundaries();

private:
    std::vector<ApecMPeriodicBoundary> periodic_boundaries;
};

#endif // APEC_M_PERIODIC_BOUNDARY_FACTORY_H
