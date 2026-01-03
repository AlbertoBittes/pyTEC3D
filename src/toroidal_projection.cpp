#include "toroidal_projection.h"
#include "coordinate_system.h"
#include "apec_m_mesh.h"
#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "variables_map.h"

#include <array>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>

#include <libmesh/elem.h>
#include <libmesh/libmesh.h>
#include <libmesh/mesh_base.h>
#include <libmesh/numeric_vector.h>

namespace
{
    const VariableInfo &get_info_by_global_label(const std::string &global_label)
    {
        const auto key = VariablesMap::get_key_by_predicate_one(
            [&global_label](const VariableInfo &info)
            {
                return info.global_label == global_label;
            });
        return VariablesMap::get_info_by_key(key);
    }

    void fill_nodal_values_from_local_solution(const MeshBase &mesh,
                                               const System &system,
                                               unsigned int var_number,
                                               const NumericVector<Number> &solution_vector,
                                               std::vector<Real> &values)
    {
        const auto &dof_map = system.get_dof_map();
        std::vector<dof_id_type> dof_indices;
        for (const auto *node : mesh.local_node_ptr_range())
        {
            const auto id = node->id();
            dof_map.dof_indices(node, dof_indices, var_number);
            if (!dof_indices.empty())
                values[id] = Real(solution_vector.el(dof_indices.front()));
        }
    }

    void fill_element_values_from_local_solution(const MeshBase &mesh,
                                                 const System &system,
                                                 unsigned int var_number,
                                                 const NumericVector<Number> &solution_vector,
                                                 std::vector<Real> &values)
    {
        const auto &dof_map = system.get_dof_map();
        std::vector<dof_id_type> dof_indices;
        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            if (elem->infinite() || elem->type() == NODEELEM)
                continue;
            const auto id = elem->id();
            dof_map.dof_indices(elem, dof_indices, var_number);
            if (!dof_indices.empty())
                values[id] = Real(solution_vector.el(dof_indices.front()));
        }
    }

    std::map<dof_id_type, Point> map_mesh_local_coordinates(const ApecMMesh &mesh,
                                                            bool is_nodal,
                                                            const CylindricalCoordinateSystem &coordinate_system)
    {
        if (is_nodal)
        {
            const auto node_range = mesh.local_node_ptr_range();
            return mesh.map_local_coordinates(node_range.begin(), node_range.end(), coordinate_system);
        }

        std::vector<const Elem *> elements;
        elements.reserve(mesh.n_active_elem());
        for (const auto &elem : mesh.active_local_element_ptr_range())
        {
            if (elem->infinite() || elem->type() == NODEELEM)
                continue;
            elements.push_back(elem);
        }
        return mesh.map_local_coordinates(elements.begin(), elements.end(), coordinate_system);
    }

    struct ProjectionContext
    {
        const MeshBase &input_mesh;
        const MeshBase &target_mesh;
        const System &input_system;
        const System &target_system;
        bool input_is_nodal;
        bool target_is_nodal;
        unsigned int input_var_number;
        unsigned int target_var_number;
        std::size_t input_count;
        std::size_t target_count;
    };

    std::vector<Real> build_input_values(const ProjectionContext &context)
    {
        std::vector<Real> input_values(context.input_count, 0.0);
        if (!context.input_system.current_local_solution)
            libmesh_error_msg("ToroidalProjection: input system current_local_solution is null.");
        const NumericVector<Number> &solution_vector = *context.input_system.current_local_solution;
        if (context.input_is_nodal)
        {
            fill_nodal_values_from_local_solution(context.input_mesh,
                                                  context.input_system,
                                                  context.input_var_number,
                                                  solution_vector,
                                                  input_values);
        }
        else
        {
            fill_element_values_from_local_solution(context.input_mesh,
                                                    context.input_system,
                                                    context.input_var_number,
                                                    solution_vector,
                                                    input_values);
        }
        context.input_mesh.comm().sum(input_values);
        return input_values;
    }

    template <typename Func>
    void for_row_range(const MeshBase &mesh, std::size_t row_count, Func func)
    {
        const auto rank = mesh.comm().rank();
        const auto size = mesh.comm().size();
        const std::size_t row_begin = (row_count * static_cast<std::size_t>(rank)) / static_cast<std::size_t>(size);
        const std::size_t row_end = (row_count * static_cast<std::size_t>(rank + 1)) / static_cast<std::size_t>(size);
        for (std::size_t row_index = row_begin; row_index < row_end; ++row_index)
            func(row_index);
    }

    ProjectionContext build_projection_context(const ApecMEquationSystems &input_eq,
                                               const ApecMEquationSystems &target_eq,
                                               const std::string &input_variable_label,
                                               const std::string &target_variable_label)
    {
        const auto &input_mesh = input_eq.get_mesh();
        const auto &target_mesh = target_eq.get_mesh();

        const auto &input_info = get_info_by_global_label(input_variable_label);
        const auto input_system_type = input_info.system_type;
        const auto &input_system = input_eq.get_system(ApecMSystemBase::system_name_from_type(input_system_type));
        const bool input_is_nodal = ApecMSystemBase::is_nodal_system(input_system_type);
        if (!input_system.has_variable(input_variable_label))
            libmesh_error_msg("ToroidalProjection: system does not contain variable: " + input_variable_label);
        const unsigned int input_var_number = input_system.variable_number(input_variable_label);

        const auto &target_info = get_info_by_global_label(target_variable_label);
        const auto target_system_type = target_info.system_type;
        const auto &target_system = target_eq.get_system(ApecMSystemBase::system_name_from_type(target_system_type));
        const bool target_is_nodal = ApecMSystemBase::is_nodal_system(target_system_type);
        if (!target_system.has_variable(target_variable_label))
            libmesh_error_msg("ToroidalProjection: system does not contain variable: " + target_variable_label);
        const unsigned int target_var_number = target_system.variable_number(target_variable_label);

        const std::size_t input_count = input_is_nodal ? input_mesh.n_nodes() : input_mesh.n_elem();
        const std::size_t target_count = target_is_nodal ? target_mesh.n_nodes() : target_mesh.n_elem();

        return ProjectionContext{
            input_mesh,
            target_mesh,
            input_system,
            target_system,
            input_is_nodal,
            target_is_nodal,
            input_var_number,
            target_var_number,
            input_count,
            target_count};
    }

    void apply_target_values(const ProjectionContext &context,
                             const std::vector<Real> &target_values,
                             const std::vector<Real> &target_mask,
                             NumericVector<Number> &target_vector)
    {
        const auto &target_dof_map = context.target_system.get_dof_map();
        std::vector<dof_id_type> dof_indices;

        if (context.target_is_nodal)
        {
            for (const auto *node : context.target_mesh.local_node_ptr_range())
            {
                const auto id = node->id();
                if (id >= context.target_count || target_mask[id] == 0.0)
                    continue;
                target_dof_map.dof_indices(node, dof_indices, context.target_var_number);
                if (!dof_indices.empty())
                    target_vector.set(dof_indices.front(), target_values[id]);
            }
        }
        else
        {
            for (const auto &elem : context.target_mesh.active_local_element_ptr_range())
            {
                if (elem->infinite() || elem->type() == NODEELEM)
                    continue;
                const auto id = elem->id();
                if (id >= context.target_count || target_mask[id] == 0.0)
                    continue;
                target_dof_map.dof_indices(elem, dof_indices, context.target_var_number);
                if (!dof_indices.empty())
                    target_vector.set(dof_indices.front(), target_values[id]);
            }
        }
        target_vector.close();
    }

}

const ToroidalProjection::CachedEntry &ToroidalProjection::get_or_build_cache(
    const MeshBase &input_mesh,
    const MeshBase &target_mesh,
    bool input_is_nodal,
    bool target_is_nodal,
    const CylindricalCoordinateSystem &coordinate_system,
    Real tolerance)
{
    const std::size_t input_count = input_is_nodal ? input_mesh.n_nodes() : input_mesh.n_elem();
    const std::size_t target_count = target_is_nodal ? target_mesh.n_nodes() : target_mesh.n_elem();

    const auto key = make_cache_key(
        input_mesh,
        target_mesh,
        input_is_nodal,
        target_is_nodal,
        coordinate_system,
        tolerance,
        input_count,
        target_count);

    auto it = cache_.find(key);
    if (it != cache_.end())
        return it->second;

    CachedEntry entry = build_cache(
        input_mesh,
        target_mesh,
        input_is_nodal,
        target_is_nodal,
        coordinate_system,
        tolerance);

    auto inserted = cache_.emplace(key, std::move(entry));
    return inserted.first->second;
}

ToroidalProjection::CachedEntry ToroidalProjection::build_cache(
    const MeshBase &input_mesh,
    const MeshBase &target_mesh,
    bool input_is_nodal,
    bool target_is_nodal,
    const CylindricalCoordinateSystem &coordinate_system,
    Real tolerance) const
{
    if (tolerance <= 0.0)
        libmesh_error_msg("tolerance must be positive.");

    const std::size_t input_count = input_is_nodal ? input_mesh.n_nodes() : input_mesh.n_elem();
    const std::size_t target_count = target_is_nodal ? target_mesh.n_nodes() : target_mesh.n_elem();

    const auto *input_apec_mesh = dynamic_cast<const ApecMMesh *>(&input_mesh);
    if (!input_apec_mesh)
        libmesh_error_msg("ToroidalProjection: input mesh is not an ApecMMesh.");
    const auto *target_apec_mesh = dynamic_cast<const ApecMMesh *>(&target_mesh);
    if (!target_apec_mesh)
        libmesh_error_msg("ToroidalProjection: target mesh is not an ApecMMesh.");

    auto input_map = map_mesh_local_coordinates(*input_apec_mesh, input_is_nodal, coordinate_system);
    auto target_map = map_mesh_local_coordinates(*target_apec_mesh, target_is_nodal, coordinate_system);

    std::map<dof_id_type, std::array<Real, 3>> input_local_map;
    std::map<dof_id_type, std::array<Real, 3>> target_local_map;
    for (const auto &entry : input_map)
        input_local_map.emplace(entry.first, std::array<Real, 3>{{entry.second(0), entry.second(1), entry.second(2)}});
    for (const auto &entry : target_map)
        target_local_map.emplace(entry.first, std::array<Real, 3>{{entry.second(0), entry.second(1), entry.second(2)}});

    input_mesh.comm().set_union(input_local_map);
    target_mesh.comm().set_union(target_local_map);

    const Real scale = 1.0 / tolerance;
    std::unordered_map<QuantKey, std::vector<dof_id_type>, QuantKeyHash> target_index;
    target_index.reserve(target_count);

    for (const auto &entry : target_local_map)
    {
        if (entry.first >= target_count)
            continue;
        const auto &local = entry.second;
        const Real radius = local[0];
        const Real axial = local[2];
        const auto key = QuantKey{
            static_cast<long long>(std::llround(radius * scale)),
            static_cast<long long>(std::llround(axial * scale))};
        target_index[key].push_back(entry.first);
    }

    CachedEntry entry;
    entry.map.resize(input_count);

    for (std::size_t i = 0; i < input_count; ++i)
    {
        auto &row = entry.map[i];
        row.push_back(static_cast<dof_id_type>(i));
        const auto input_id = static_cast<dof_id_type>(i);
        const auto map_it = input_local_map.find(input_id);
        if (map_it == input_local_map.end())
            continue;
        const auto &local = map_it->second;
        const Real radius = local[0];
        const Real axial = local[2];
        const auto key = QuantKey{
            static_cast<long long>(std::llround(radius * scale)),
            static_cast<long long>(std::llround(axial * scale))};
        const auto it = target_index.find(key);
        if (it == target_index.end())
            continue;
        const auto &matches = it->second;
        row.insert(row.end(), matches.begin(), matches.end());
    }

    return entry;
}

std::string ToroidalProjection::make_cache_key(
    const MeshBase &input_mesh,
    const MeshBase &target_mesh,
    bool input_is_nodal,
    bool target_is_nodal,
    const CylindricalCoordinateSystem &coordinate_system,
    Real tolerance,
    std::size_t input_count,
    std::size_t target_count) const
{
    const Point origin = coordinate_system.origin();
    const auto rotation = coordinate_system.local_to_global_rotation(origin);
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << reinterpret_cast<std::uintptr_t>(&input_mesh) << '|'
        << reinterpret_cast<std::uintptr_t>(&target_mesh) << '|'
        << input_is_nodal << '|'
        << target_is_nodal << '|'
        << input_count << '|'
        << target_count << '|'
        << std::setprecision(12)
        << origin(0) << ',' << origin(1) << ',' << origin(2) << '|'
        << rotation(0, 0) << ',' << rotation(0, 1) << ',' << rotation(0, 2) << ','
        << rotation(1, 0) << ',' << rotation(1, 1) << ',' << rotation(1, 2) << ','
        << rotation(2, 0) << ',' << rotation(2, 1) << ',' << rotation(2, 2) << '|'
        << tolerance;
    return oss.str();
}

void ToroidalProjection::project(
    const ApecMEquationSystems &input_eq,
    const ApecMEquationSystems &target_eq,
    const std::string &input_variable_label,
    const std::string &target_variable_label,
    const CylindricalCoordinateSystem &coordinate_system,
    Real tolerance)
{
    const auto context = build_projection_context(input_eq,
                                                  target_eq,
                                                  input_variable_label,
                                                  target_variable_label);

    const auto &cached = get_or_build_cache(
        context.input_mesh,
        context.target_mesh,
        context.input_is_nodal,
        context.target_is_nodal,
        coordinate_system,
        tolerance);
    const auto &active_map = cached.map;
    // target_vector.zero();

    std::vector<Real> input_values = build_input_values(context);

    std::vector<Real> target_values(context.target_count, 0.0);
    std::vector<Real> target_mask(context.target_count, 0.0);

    const std::size_t row_count = active_map.size();
    for_row_range(context.target_mesh, row_count, [&](std::size_t row_index)
                  {
        const auto &row = active_map[row_index];
        if (row.empty())
            return;
        const auto input_id = row.front();
        if (input_id >= input_values.size())
            return;
        const Real value = input_values[input_id];
        for (std::size_t k = 1; k < row.size(); ++k)
        {
            const auto target_id = row[k];
            if (target_id >= context.target_count)
                continue;
            target_values[target_id] = value;
            target_mask[target_id] = 1.0;
        } });

    context.target_mesh.comm().sum(target_values);
    context.target_mesh.comm().sum(target_mask);

    NumericVector<Number> &target_vector = *context.target_system.current_local_solution;
    apply_target_values(context, target_values, target_mask, target_vector);
}

void ToroidalProjection::toroidal_average(
    const ApecMEquationSystems &input_eq,
    const ApecMEquationSystems &target_eq,
    const std::string &input_variable_label,
    const std::string &target_variable_label,
    const CylindricalCoordinateSystem &coordinate_system,
    Real tolerance)
{
    const auto context = build_projection_context(input_eq,
                                                  target_eq,
                                                  input_variable_label,
                                                  target_variable_label);

    const auto &cached = get_or_build_cache(
        context.target_mesh,
        context.input_mesh,
        context.target_is_nodal,
        context.input_is_nodal,
        coordinate_system,
        tolerance);
    const auto &active_map = cached.map;
    // target_vector.close();
    // target_vector.zero();

    std::vector<Real> input_values = build_input_values(context);

    std::vector<Real> target_values(context.target_count, 0.0);
    std::vector<Real> target_mask(context.target_count, 0.0);

    const std::size_t row_count = active_map.size();
    for_row_range(context.target_mesh, row_count, [&](std::size_t row_index)
                  {
        const auto &row = active_map[row_index];
        if (row.empty())
            return;
        const auto target_id = row.front();
        if (target_id >= context.target_count)
            return;
        if (row.size() == 1)
        {
            target_values[target_id] = std::numeric_limits<Real>::quiet_NaN();
            target_mask[target_id] = 1.0;
            return;
        }
        Real accum = 0.0;
        Real count = 0.0;
        for (std::size_t k = 1; k < row.size(); ++k)
        {
            const auto input_id = row[k];
            if (input_id >= context.input_count)
                continue;
            accum += input_values[input_id];
            count += 1.0;
        }
        if (count > 0.0)
            target_values[target_id] = accum / count;
        else
            target_values[target_id] = std::numeric_limits<Real>::quiet_NaN();
        target_mask[target_id] = 1.0; });

    context.target_mesh.comm().sum(target_values);
    context.target_mesh.comm().sum(target_mask);

    NumericVector<Number> &target_vector = *context.target_system.current_local_solution;
    apply_target_values(context, target_values, target_mask, target_vector);
}
