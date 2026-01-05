#ifndef TOROIDAL_PROJECTION_H
#define TOROIDAL_PROJECTION_H

#include <string>
#include <vector>

#include <libmesh/mesh_base.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/point.h>

#include "coordinate_system.h"

class ApecMEquationSystems;
class ToroidalProjection
{
public:
    using Map = std::vector<std::vector<libMesh::dof_id_type>>;

    void set_maps(const Map &nodal_map, const Map &elem_map);

    void project(const ApecMEquationSystems &input_eq,
                 ApecMEquationSystems &target_eq,
                 const std::string &input_variable_label,
                 const std::string &target_variable_label,
                 const CylindricalCoordinateSystem &coordinate_system = CylindricalCoordinateSystem());

    void toroidal_average(const ApecMEquationSystems &input_eq,
                          ApecMEquationSystems &target_eq,
                          const std::string &input_variable_label,
                          const std::string &target_variable_label,
                          const CylindricalCoordinateSystem &coordinate_system = CylindricalCoordinateSystem());

private:
    Map nodal_map_fwd_;
    Map nodal_map_rev_;
    Map elem_map_fwd_;
    Map elem_map_rev_;
};

#endif // TOROIDAL_PROJECTION_H
