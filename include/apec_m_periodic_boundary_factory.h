#ifndef APEC_M_PERIODIC_BOUNDARY_FACTORY_H
#define APEC_M_PERIODIC_BOUNDARY_FACTORY_H

#include "apec_m_periodic_boundary.h"
#include "fe_data_object.h"
#include "fe_parent_object.h"

#include "libmesh/system.h"
#include "libmesh/periodic_boundary_base.h"
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace libMesh;

// ------------------------------------------------------
// Factory to build ApecMPeriodicBoundary instances from serialized settings data
// ------------------------------------------------------
class ApecMPeriodicBoundaryFactory
    : public BinaryDataFEDataObject,
      public FEParentObjectBase<System>
{
public:
    explicit ApecMPeriodicBoundaryFactory(std::ifstream &file);
    // NEW: Populate periodic boundaries directly from stored vtkFieldData

    const std::vector<ApecMPeriodicBoundary> &get_periodic_boundaries() const;

    void initialize(System *sys);

    void print_data_info() const override;

private:
    void apply_periodic_boundaries();
    void compress_periodic_boundaries();
    void populate_periodic_boundaries_from_map();

    std::vector<ApecMPeriodicBoundary> periodic_boundaries;
};

#endif // APEC_M_PERIODIC_BOUNDARY_FACTORY_H
