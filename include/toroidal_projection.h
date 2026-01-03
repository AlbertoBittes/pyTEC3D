#ifndef TOROIDAL_PROJECTION_H
#define TOROIDAL_PROJECTION_H

#include <string>
#include <unordered_map>
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

    void project(const ApecMEquationSystems &input_eq,
                 const ApecMEquationSystems &target_eq,
                 const std::string &input_variable_label,
                 const libMesh::NumericVector<libMesh::Number> &input_vector,
                 const std::string &target_variable_label,
                 libMesh::NumericVector<libMesh::Number> &target_vector,
                 const CylindricalCoordinateSystem &coordinate_system = CylindricalCoordinateSystem(),
                 libMesh::Real tolerance = 1e-6);

    void project_inverse(const ApecMEquationSystems &input_eq,
                         const ApecMEquationSystems &target_eq,
                         const std::string &input_variable_label,
                         const libMesh::NumericVector<libMesh::Number> &input_vector,
                         const std::string &target_variable_label,
                         libMesh::NumericVector<libMesh::Number> &target_vector,
                         const CylindricalCoordinateSystem &coordinate_system = CylindricalCoordinateSystem(),
                         libMesh::Real tolerance = 1e-6);

private:
    struct CachedEntry
    {
        Map map;
    };

    struct QuantKey
    {
        long long r = 0;
        long long a = 0;

        bool operator==(const QuantKey &other) const
        {
            return r == other.r && a == other.a;
        }
    };

    struct QuantKeyHash
    {
        std::size_t operator()(const QuantKey &key) const
        {
            std::size_t h1 = std::hash<long long>{}(key.r);
            std::size_t h2 = std::hash<long long>{}(key.a);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };

    const CachedEntry &get_or_build_cache(const libMesh::MeshBase &input_mesh,
                                          const libMesh::MeshBase &target_mesh,
                                          bool input_is_nodal,
                                          bool target_is_nodal,
                                          const CylindricalCoordinateSystem &coordinate_system,
                                          libMesh::Real tolerance);

    CachedEntry build_cache(const libMesh::MeshBase &input_mesh,
                            const libMesh::MeshBase &target_mesh,
                            bool input_is_nodal,
                            bool target_is_nodal,
                            const CylindricalCoordinateSystem &coordinate_system,
                            libMesh::Real tolerance) const;

    std::string make_cache_key(const libMesh::MeshBase &input_mesh,
                               const libMesh::MeshBase &target_mesh,
                               bool input_is_nodal,
                               bool target_is_nodal,
                               const CylindricalCoordinateSystem &coordinate_system,
                               libMesh::Real tolerance,
                               std::size_t input_count,
                               std::size_t target_count) const;

    std::unordered_map<std::string, CachedEntry> cache_;
};

#endif // TOROIDAL_PROJECTION_H
