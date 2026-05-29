#include "apec_m_system.h"

#include <array>
#include <iomanip>
#include <sstream>
#include <set>
#include <unordered_set>
#include <vector>

#include "system_matrix_reduction.h"
#include "element_types_factory.h"
#include "global_constraints_object.h"
#include "apec_m_logging.h"

#include "libmesh/libmesh.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/numeric_vector.h"
#include "variables_map.h"

namespace
{
    std::string to_scientific_string(Real value)
    {
        std::ostringstream out;
        out << std::scientific << std::setprecision(6) << value;
        return out.str();
    }

    std::string format_time_step_table(const ApecMLinearSystem::TimeStep &ts, const std::string &title)
    {
        std::ostringstream out;
        const auto append_row = [&out](const std::string &label, const std::string &value)
        {
            out << "| " << std::left << std::setw(20) << label
                << " | " << std::left << std::setw(20) << value << " |\n";
        };

        out << "+---------------------------------------------+\n";
        out << "| " << std::left << std::setw(44) << title << "|\n";
        out << "+----------------------+----------------------+\n";
        append_row("index", std::to_string(ts.time_step_idx));
        append_row("time", to_scientific_string(ts.time));
        append_row("dt", to_scientific_string(ts.dt));
        append_row("theta", to_scientific_string(ts.theta));
        append_row("transient_effects", ts.transient_effects ? "true" : "false");
        out << "+----------------------+----------------------+";
        return out.str();
    }
} // namespace

const std::map<ApecMSystemBase::SystemType, std::string> ApecMSystemBase::SYSTEM_INFO_MAP = {
    {ApecMSystemBase::DOF_SOLUTION_SYSTEM, "DOF Solution"},
    {ApecMSystemBase::NODAL_SOLUTION, "Nodal Solution"},
    {ApecMSystemBase::ELEMENT_SOLUTION, "Element Solution"},
};

const std::string &ApecMSystemBase::system_name_from_type(ApecMSystemBase::SystemType system_type)
{
    auto it = ApecMSystemBase::SYSTEM_INFO_MAP.find(system_type);
    if (it == ApecMSystemBase::SYSTEM_INFO_MAP.end())
        libmesh_error_msg("System type not found.");
    return it->second;
}

void ApecMLinearSystem::release_subset_solve()
{
    if (this->_subset != nullptr)
    {
        delete this->_subset;
        this->_subset = nullptr;
    }
}

void ApecMLinearSystem::prepare_time_step_solve(const pybind11::dict &time_step_dict)
{
    const Real requested_time = pybind11::cast<Real>(time_step_dict["time"]);
    const Real requested_theta = pybind11::cast<Real>(time_step_dict["theta"]);
    const bool requested_transient_effect = pybind11::cast<bool>(time_step_dict["transient_effects"]);

    if (!time_step)
    {
        time_step = std::make_unique<TimeStep>();
        time_step->time = requested_time;
        time_step->theta = requested_theta;
        time_step->transient_effects = requested_transient_effect;
        time_step->time_step_idx = 0;
        time_step->dt = requested_time;
        time_step->updated = true;
        *this->old_local_solution = *this->current_local_solution;
        this->old_local_solution->close();
        apec::logging::info("Advanced to new time step");
        apec::logging::table_root(format_time_step_table(*time_step, "Time Step Information"));
        return;
    }

    const Real current_time = time_step->time;
    const Real time_scale = std::max(std::abs(requested_time), std::abs(current_time));
    const Real time_tol = std::numeric_limits<Real>::epsilon() * time_scale;
    const Real time_delta = requested_time - current_time;

    if (time_delta > time_tol)
    {
        old_time_step = std::move(time_step);
        time_step = std::make_unique<TimeStep>();
        time_step->time = requested_time;
        time_step->theta = requested_theta;
        time_step->transient_effects = requested_transient_effect;
        time_step->time_step_idx = old_time_step->time_step_idx + 1;
        time_step->dt = requested_time - current_time;
        time_step->updated = true;
        *this->old_local_solution = *this->current_local_solution;
        this->old_local_solution->close();
        apec::logging::info("Advanced to new time step");
        apec::logging::table_root(format_time_step_table(*time_step, "Time Step Information"));
        return;
    }

    if (std::abs(time_delta) <= time_tol)
    {
        const Real current_theta = time_step->theta;
        const bool current_transient_effect = time_step->transient_effects;
        const Real theta_tol = std::numeric_limits<Real>::epsilon() * std::abs(current_theta);

        if (std::abs(requested_theta - current_theta) > theta_tol ||
            requested_transient_effect != current_transient_effect)
        {
            libmesh_error_msg("ApecMLinearSystem::prepare_time_step_solve: theta/transient_effects changed while time is unchanged.");
        }

        time_step->time = requested_time;
        time_step->theta = requested_theta;
        time_step->transient_effects = requested_transient_effect;
        time_step->updated = false;
        return;
    }

    libmesh_error_msg("ApecMLinearSystem::prepare_time_step_solve: received time smaller than current time step.");
}

void ApecMLinearSystem::compute_solution_time_derivative()
{
    if (!time_step)
        libmesh_error_msg("ApecMLinearSystem::compute_solution_time_derivative: time_step is null.");

    const bool transient_effect = time_step->transient_effects;
    const Real dt = time_step->dt;
    const Real theta = time_step->theta;
    NumericVector<Number> *time_derivative_vector = this->request_vector("solution_time_derivative");
    if (!transient_effect)
    {
        // If transient effects are not considered, the solution time derivative is zero
        this->older_local_solution->zero();
        time_derivative_vector->zero();
    }
    else
    {
        Real theta_x_dt = theta * dt;
        Real theta_minus_1_dt = (theta - 1) * dt;

        this->older_local_solution->scale(theta_minus_1_dt);
        *this->older_local_solution += *this->current_local_solution;
        *this->older_local_solution -= *this->old_local_solution;
        this->older_local_solution->scale(1 / theta_x_dt);

        time_derivative_vector->zero();
        *time_derivative_vector += *this->current_local_solution;
        *time_derivative_vector -= *this->old_local_solution;
        time_derivative_vector->scale(1 / dt);
    }
    this->older_local_solution->close();
    time_derivative_vector->close();
}

void ApecMLinearSystem::delete_rows_cols_constrained_dofs()
{
    apec::logging::info("Removing constrained DOFs from system matrix and RHS");
    SystemMatrixReduction *subset_ptr = new SystemMatrixReduction(*this);
    this->release_subset_solve();
    this->restrict_solve_to(subset_ptr, SUBSET_ZERO);
}

void ApecMLinearSystem::zero_all()
{
    ApecMSystem<TransientSystem<LinearImplicitSystem>>::zero_all();
    this->get_system_matrix().zero();
}

template <typename SysType>
void ApecMSystem<SysType>::init_variables()
{
    std::string fe_fam;
    std::string order;

    std::string vars_set;

    auto *equation_systems = dynamic_cast<ApecMEquationSystems *>(&this->get_equation_systems());
    auto *solver_parent = equation_systems ? equation_systems->get_parent_object() : nullptr;
    if (!solver_parent)
        libmesh_error_msg("ApecMSystem::init_variables: solver parent not available.");

    auto element_types_collector = solver_parent->element_types_factory;
    // auto vars_to_subdomain_ids = element_types_collector->associateVariablesWithElementIds();

    // Do not do the loop here, attacht the initialization object to the other systems too
    // for (unsigned int k = 0; k < equation_systems.n_systems(); k++)
    // {
    // Use new helper functions to map system name to SystemType
    ApecMSystemBase::SystemType type = ApecMSystemBase::NUM_SYSTEMS;
    const std::string sys_name = this->name();
    for (const auto &pair : ApecMSystemBase::SYSTEM_INFO_MAP)
    {
        if (pair.second == sys_name)
        {
            type = pair.first;
            break;
        }
    }
    if (type == ApecMSystemBase::NUM_SYSTEMS)
        libmesh_error_msg("System type not found for system name: " + sys_name);

    if (!element_types_collector || !element_types_collector->get_element_types_map())
        libmesh_error_msg("Element types map not available for system variables.");
    auto vars_to_subdomain_ids = VariablesMap::associate_variables_with_element_IDs(
        *element_types_collector->get_element_types_map(),
        type);

    std::set<subdomain_id_type> mesh_subdomain_ids;
    this->get_mesh().subdomain_ids(mesh_subdomain_ids);

    for (auto &pair : vars_to_subdomain_ids)
    {
        const auto &variable_string = VariablesMap::get_info_by_key(pair.first).global_label;
        auto &subdomain_ids = pair.second;
        for (auto it = subdomain_ids.begin(); it != subdomain_ids.end();)
        {
            if (!mesh_subdomain_ids.count(*it))
                it = subdomain_ids.erase(it);
            else
                ++it;
        }
        if (subdomain_ids.empty())
            continue;
        auto elem_type = element_types_collector->get_element_type(*subdomain_ids.begin());
        auto fe_family = VariablesMap::get_variable_fe_family(*elem_type, pair.first);
        auto order = VariablesMap::get_variable_order(*elem_type, pair.first);

        // Add the variable to the system, specifying subdomains
        this->add_variable(
            variable_string,
            order,
            fe_family,
            &subdomain_ids);
    }
    // Add more conditions as needed for other combinations
}

template class ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>;
template class ApecMSystem<libMesh::ExplicitSystem>;
