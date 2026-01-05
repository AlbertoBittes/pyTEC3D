#include "toroidal_projection.h"
#include "coordinate_system.h"
#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "variables_map.h"

#include <limits>

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

    void fill_nodal_values_from_solution(const MeshBase &mesh,
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

    void fill_element_values_from_solution(const MeshBase &mesh,
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


    struct ProjectionContext
    {
        const MeshBase &input_mesh;
        const MeshBase &target_mesh;
        const System &input_system;
        System &target_system;
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

        const NumericVector<Number> &solution_vector = *context.input_system.solution;
        if (context.input_is_nodal)
        {
            fill_nodal_values_from_solution(context.input_mesh,
                                            context.input_system,
                                            context.input_var_number,
                                            solution_vector,
                                            input_values);
        }
        else
        {
            fill_element_values_from_solution(context.input_mesh,
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
                                               ApecMEquationSystems &target_eq,
                                               const std::string &input_variable_label,
                                               const std::string &target_variable_label)
    {
        const auto &input_mesh = input_eq.get_mesh();
        const auto &target_mesh = target_eq.get_mesh();

        const auto &input_info = get_info_by_global_label(input_variable_label);
        const auto input_system_type = VariablesMap::get_system_type_for_variable(input_eq,
                                                                                  input_info.global_name);
        const auto &input_system = input_eq.get_system(ApecMSystemBase::system_name_from_type(input_system_type));
        const bool input_is_nodal = ApecMSystemBase::is_nodal_system(input_system_type);
        if (!input_system.has_variable(input_variable_label))
            libmesh_error_msg("ToroidalProjection: system does not contain variable: " + input_variable_label);
        const unsigned int input_var_number = input_system.variable_number(input_variable_label);

        const auto &target_info = get_info_by_global_label(target_variable_label);
        const auto target_system_type = VariablesMap::get_system_type_for_variable(target_eq,
                                                                                   target_info.global_name);
        auto &target_system = target_eq.get_system(ApecMSystemBase::system_name_from_type(target_system_type));
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

    ToroidalProjection::Map build_reverse_map(const ToroidalProjection::Map &forward)
    {
        std::size_t max_target = 0;
        for (const auto &row : forward)
        {
            if (row.size() <= 1)
                continue;
            for (std::size_t k = 1; k < row.size(); ++k)
                max_target = std::max<std::size_t>(max_target, row[k]);
        }

        ToroidalProjection::Map reverse;
        reverse.resize(max_target + 1);
        for (std::size_t id = 0; id <= max_target; ++id)
            reverse[id].push_back(static_cast<dof_id_type>(id));

        for (const auto &row : forward)
        {
            if (row.empty())
                continue;
            const auto source_id = row.front();
            for (std::size_t k = 1; k < row.size(); ++k)
            {
                const auto target_id = row[k];
                if (target_id >= reverse.size())
                    continue;
                reverse[target_id].push_back(source_id);
            }
        }
        return reverse;
    }

}

void ToroidalProjection::set_maps(const Map &nodal_map, const Map &elem_map)
{
    if (!nodal_map.empty())
    {
        nodal_map_fwd_ = nodal_map;
        nodal_map_rev_ = build_reverse_map(nodal_map);
    }
    if (!elem_map.empty())
    {
        elem_map_fwd_ = elem_map;
        elem_map_rev_ = build_reverse_map(elem_map);
    }
}

void ToroidalProjection::project(
    const ApecMEquationSystems &input_eq,
    ApecMEquationSystems &target_eq,
    const std::string &input_variable_label,
    const std::string &target_variable_label,
    const CylindricalCoordinateSystem &coordinate_system)
{
    const auto context = build_projection_context(input_eq,
                                                  target_eq,
                                                  input_variable_label,
                                                  target_variable_label);
    (void)coordinate_system;

    const auto &active_map = context.input_is_nodal ? nodal_map_fwd_ : elem_map_fwd_;
    if (active_map.empty())
        libmesh_error_msg("ToroidalProjection: projection map is empty.");
    if (active_map.size() != context.input_count)
        libmesh_error_msg("ToroidalProjection: projection map row count does not match input mesh.");

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

    NumericVector<Number> &target_vector = *context.target_system.solution;
    apply_target_values(context, target_values, target_mask, target_vector);
    context.target_system.update();
}

void ToroidalProjection::toroidal_average(
    const ApecMEquationSystems &input_eq,
    ApecMEquationSystems &target_eq,
    const std::string &input_variable_label,
    const std::string &target_variable_label,
    const CylindricalCoordinateSystem &coordinate_system)
{
    const auto context = build_projection_context(input_eq,
                                                  target_eq,
                                                  input_variable_label,
                                                  target_variable_label);
    (void)coordinate_system;

    const auto &active_map = context.target_is_nodal ? nodal_map_rev_ : elem_map_rev_;
    if (active_map.empty())
        libmesh_error_msg("ToroidalProjection: projection map is empty.");
    if (active_map.size() != context.target_count)
        libmesh_error_msg("ToroidalProjection: projection map row count does not match target mesh.");

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

    NumericVector<Number> &target_vector = *context.target_system.solution;
    apply_target_values(context, target_values, target_mask, target_vector);
    context.target_system.update();
}
