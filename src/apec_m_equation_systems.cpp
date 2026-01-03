#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "apec_m_mesh.h"
#include "coordinate_system.h"
#include "time_step.h"
#include "variables_map.h"
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <utility>

namespace
{
    VectorValue<Real> read_vector(const NumericVector<Number> &vec,
                                  const std::array<dof_id_type, 3> &dofs)
    {
        return VectorValue<Real>(vec(dofs[0]), vec(dofs[1]), vec(dofs[2]));
    }

    void write_vector(NumericVector<Number> &vec,
                      const std::array<dof_id_type, 3> &dofs,
                      const VectorValue<Real> &values)
    {
        vec.set(dofs[0], values(0));
        vec.set(dofs[1], values(1));
        vec.set(dofs[2], values(2));
    }

    template <typename TransformFn>
    void apply_vector_transform(System &system,
                                bool is_nodal,
                                const std::array<unsigned int, 3> &var_numbers,
                                TransformFn &&transform)
    {
        const auto &dof_map = system.get_dof_map();
        std::vector<dof_id_type> dof_indices;
        const auto &mesh = system.get_mesh();

        if (is_nodal)
        {
            for (const auto *node : mesh.local_node_ptr_range())
            {
                std::array<dof_id_type, 3> dofs{};
                bool complete = true;
                for (unsigned int i = 0; i < 3; ++i)
                {
                    dof_map.dof_indices(node, dof_indices, var_numbers[i]);
                    if (dof_indices.empty())
                    {
                        complete = false;
                        break;
                    }
                    dofs[i] = dof_indices.front();
                }
                if (!complete)
                    continue;

                transform(*node, dofs);
            }
        }
        else
        {
            for (const auto &elem : mesh.active_local_element_ptr_range())
            {
                if (elem->infinite() || elem->type() == NODEELEM)
                    continue;
                std::array<dof_id_type, 3> dofs{};
                bool complete = true;
                for (unsigned int i = 0; i < 3; ++i)
                {
                    dof_map.dof_indices(elem, dof_indices, var_numbers[i]);
                    if (dof_indices.empty())
                    {
                        complete = false;
                        break;
                    }
                    dofs[i] = dof_indices.front();
                }
                if (!complete)
                    continue;

                transform(elem->centroid(), dofs);
            }
        }
    }
} // namespace

void ApecMEquationSystems::apply_ramped_effects(const Substep &time_substep)
{
    auto ramp_sets = collect_ramp_vectors();

    for (auto &[old_vec, current_vec, ramped_vec] : ramp_sets)
        ApecMSystem<System>::ramp_vectors(time_substep, *old_vec, *current_vec, *ramped_vec);
}

void ApecMEquationSystems::configure_systems()
{
    libMesh::out << "Building system of equations..." << std::endl;

    this->clear_preserving_solver_settings();

    libMesh::out << "Checking system options..." << std::endl;

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        ApecMSystemBase::SystemType type = static_cast<ApecMSystemBase::SystemType>(i);
        const std::string sys_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            auto &linear_system = this->add_system<ApecMLinearSystem>(sys_name);
            if (auto assembler = linear_system.get_assembler())
                assembler->set_parent_object(this);
        }
        else
            this->add_system<ApecMExplicitSystem>(sys_name);
    }
}

void ApecMEquationSystems::initialize_system_contexts()
{
    auto system_builder_view = this->get_context_view();
    if (!system_builder_view)
        return;

    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            if (auto *system = dynamic_cast<ApecMLinearSystem *>(&base_system))
            {
                system->set_context_view(system_builder_view);
                system->init_variables();
            }
        }
        else if (auto *system = dynamic_cast<ApecMExplicitSystem *>(&base_system))
        {
            system->set_context_view(system_builder_view);
            system->init_variables();
        }
    }
}

void ApecMEquationSystems::post_initialize_systems()
{
    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (type == ApecMSystemBase::DOF_SOLUTION_SYSTEM)
        {
            if (auto *system = dynamic_cast<ApecMLinearSystem *>(&base_system))
                system->post_init();
        }
        else if (auto *system = dynamic_cast<ApecMExplicitSystem *>(&base_system))
        {
            system->post_init();
        }
    }
}

void ApecMEquationSystems::initialize_systems(const std::shared_ptr<FEMDataRepository> &model_data)
{
    if (!model_data)
    {
        libMesh::out << "ApecMEquationSystems: model repository not provided for physics initialization." << std::endl;
        return;
    }

    auto system_view = std::make_shared<ContextView>(ContextView::from(model_data));
    this->set_context_view(system_view);

    libMesh::out << "Distributing degrees of freedom..." << std::endl;
    initialize_system_contexts();
}

void ApecMEquationSystems::initialize_physics()
{
    this->init();
    post_initialize_systems();
    this->print_info();
}

void ApecMEquationSystems::clear_preserving_solver_settings()
{
    Real tol = parameters.get<Real>("linear solver tolerance");
    unsigned int max_iter = parameters.get<unsigned int>("linear solver maximum iterations");

    this->clear();

    parameters.set<Real>("linear solver tolerance") = tol;
    parameters.set<unsigned int>("linear solver maximum iterations") = max_iter;
}

void ApecMEquationSystems::zero_all_systems()
{
    for (int i = 0; i < ApecMSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<ApecMSystemBase::SystemType>(i);
        const auto &info = ApecMSystemBase::SYSTEM_INFO_MAP.at(type);
        auto &base_system = this->get_system<System>(info);

        if (auto *zeroable = dynamic_cast<ApecMSystemBase *>(&base_system))
            zeroable->zero_all();
    }
}

std::shared_ptr<ApecMEquationSystems> ApecMEquationSystems::build(ApecMMesh &mesh, APEC_M *program)
{
    if (program == nullptr)
        libmesh_error_msg("APEC_M program pointer not provided to ApecMEquationSystems::build.");
    auto equation_systems = std::make_shared<ApecMEquationSystems>(mesh);
    equation_systems->set_parent_object(program);
    equation_systems->configure_systems();
    return equation_systems;
}

bool ApecMEquationSystems::set_vector_rhs_in_coordinate_system(ApecMSystemBase::SystemType system_type,
                                                               const CoordinateSystem &coordinate_system,
                                                               const std::string &variable_label)
{
    const auto &sys_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(system_type);
    auto &system = this->get_system<ExplicitSystem>(sys_name);

    std::array<unsigned int, 3> var_numbers{};
    if (!VariablesMap::resolve_vector_var_numbers(system, variable_label, var_numbers))
        return false;

    if (!system.rhs)
        libmesh_error_msg("System RHS vector is null.");

    NumericVector<Number> &rhs = *system.rhs;

    const bool is_nodal = ApecMSystemBase::is_nodal_system(system_type);

    auto transform_values = [&](const Point &point,
                                const std::array<dof_id_type, 3> &dofs)
    {
        VectorValue<Real> global = read_vector(rhs, dofs);
        VectorValue<Real> local = coordinate_system.global_to_local_rotation(point) * global;
        write_vector(rhs, dofs, local);
    };

    apply_vector_transform(system, is_nodal, var_numbers, transform_values);

    rhs.close();
    return true;
}

bool ApecMEquationSystems::set_vector_solution_from_rhs_global(ApecMSystemBase::SystemType system_type,
                                                               const CoordinateSystem &coordinate_system,
                                                               const std::string &variable_label)
{
    const auto &sys_name = ApecMSystemBase::SYSTEM_INFO_MAP.at(system_type);
    auto &system = this->get_system<ExplicitSystem>(sys_name);

    std::array<unsigned int, 3> var_numbers{};
    if (!VariablesMap::resolve_vector_var_numbers(system, variable_label, var_numbers))
        return false;

    if (!system.rhs)
        libmesh_error_msg("System RHS vector is null.");
    if (!system.current_local_solution)
        libmesh_error_msg("System current_local_solution is null.");

    NumericVector<Number> &rhs = *system.rhs;
    NumericVector<Number> &solution = *system.current_local_solution;

    const bool is_nodal = ApecMSystemBase::is_nodal_system(system_type);

    auto transform_values = [&](const Point &point,
                                const std::array<dof_id_type, 3> &dofs)
    {
        VectorValue<Real> local = read_vector(rhs, dofs);
        VectorValue<Real> global = coordinate_system.local_to_global_rotation(point) * local;
        write_vector(solution, dofs, global);
    };

    apply_vector_transform(system, is_nodal, var_numbers, transform_values);

    solution.close();
    return true;
}

std::vector<ApecMEquationSystems::RampVectorTuple> ApecMEquationSystems::collect_ramp_vectors()
{
    std::vector<RampVectorTuple> ramp_sets;

    System &main_system = this->get_system<System>("ApecMLinearSystem");
    ramp_sets.emplace_back(
        main_system.request_vector("old_node_dof_constraints"),
        main_system.request_vector("node_dof_constraints"),
        main_system.request_vector("ramped_node_dof_constraints"));

    System &element_loads_system = this->get_system<System>("ElementLoadsSystem");
    ramp_sets.emplace_back(
        element_loads_system.request_vector("old_loads"),
        element_loads_system.request_vector("new_loads"),
        element_loads_system.current_local_solution.get());

    System &nodal_loads_system = this->get_system<System>("NodalLoadsSystem");
    ramp_sets.emplace_back(
        nodal_loads_system.request_vector("old_loads"),
        nodal_loads_system.request_vector("new_loads"),
        nodal_loads_system.current_local_solution.get());

    return ramp_sets;
}
