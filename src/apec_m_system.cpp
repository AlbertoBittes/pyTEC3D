#include "apec_m_system.h"

#include <array>
#include <set>
#include <unordered_set>
#include <vector>

#include "system_matrix_reduction.h"
#include "element_types_factory.h"
#include "global_constraints_object.h"
#include "apec_m_mesh.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh_serializer.h"
#include "libmesh/numeric_vector.h"
#include "variables_map.h"

namespace
{
    template <typename SysLike>
    bool transform_coordinate_system_solution(ApecMSystem<SysLike> &system,
                                              NumericVector<Number> &solution,
                                              const CoordinateSystem &coordinate_system,
                                              const std::string &variable_label,
                                              const std::vector<dof_id_type> &selected_entity_ids,
                                              bool use_global_to_local_rotation)
    {
        const bool is_nodal = ApecMSystemBase::is_nodal_system(system.system_type);
        const ApecMMesh *mesh = dynamic_cast<const ApecMMesh *>(&system.get_mesh());

        const std::vector<std::string> vector_components = VariablesMap::get_global_labels_for_label(variable_label);
        if (vector_components.size() < 3 || selected_entity_ids.empty())
            return true;

        std::vector<Point> positions(selected_entity_ids.size());
        for (std::size_t i = 0; i < selected_entity_ids.size(); ++i)
        {
            const dof_id_type object_id = selected_entity_ids[i];
            if (is_nodal)
            {
                const Node *node_ptr = mesh->query_node_ptr(object_id);
                positions[i] = ApecMMesh::entity_point(*node_ptr);
            }
            else
            {
                const Elem *elem_ptr = mesh->query_elem_ptr(object_id);
                positions[i] = ApecMMesh::entity_point(*elem_ptr);
            }
        }

        std::vector<std::vector<dof_id_type>> dofs(3);
        system.dof_numbers_vars(vector_components[0], !is_nodal, selected_entity_ids, dofs[0], false);
        system.dof_numbers_vars(vector_components[1], !is_nodal, selected_entity_ids, dofs[1], false);
        system.dof_numbers_vars(vector_components[2], !is_nodal, selected_entity_ids, dofs[2], false);

        std::vector<ApecMTensorValue> rotations(selected_entity_ids.size());
        for (std::size_t i = 0; i < selected_entity_ids.size(); ++i)
        {
            if (use_global_to_local_rotation)
                rotations[i] = coordinate_system.global_to_local_rotation(positions[i]);
            else
                rotations[i] = coordinate_system.local_to_global_rotation(positions[i]);
        }

        std::vector<std::vector<Real>> values(3);
        solution.get(dofs[0], values[0]);
        solution.get(dofs[1], values[1]);
        solution.get(dofs[2], values[2]);

        for (std::size_t i = 0; i < selected_entity_ids.size(); ++i)
        {
            const auto &rotation = rotations[i];
            const Point input_vec(values[0][i], values[1][i], values[2][i]);
            const Point output_vec = rotation * input_vec;
            values[0][i] = output_vec(0);
            values[1][i] = output_vec(1);
            values[2][i] = output_vec(2);
        }

        solution.insert(values[0], dofs[0]);
        solution.insert(values[1], dofs[1]);
        solution.insert(values[2], dofs[2]);

        solution.close();
        system.update();
        return true;
    }
} // namespace

const std::map<ApecMSystemBase::SystemType, std::string> ApecMSystemBase::SYSTEM_INFO_MAP = {
    {ApecMSystemBase::DOF_SOLUTION_SYSTEM, "ApecMLinearSystem"},
    {ApecMSystemBase::NODAL_LOADS, "NodalLoadsSystem"},
    {ApecMSystemBase::ELEMENT_LOADS, "ElementLoadsSystem"},
    {ApecMSystemBase::NODAL_SOLUTION, "NodalSolutionSystem"},
    {ApecMSystemBase::ELEMENT_SOLUTION, "ElementSolutionSystem"},
    {ApecMSystemBase::SUPPORT_SOLUTION, "SolutionSupportSystem"},
    {ApecMSystemBase::NODAL_B, "NodalBSolutionSystem"},
};

bool ApecMSystemBase::is_nodal_system(ApecMSystemBase::SystemType system_type)
{
    return system_type == ApecMSystemBase::DOF_SOLUTION_SYSTEM ||
           system_type == ApecMSystemBase::NODAL_LOADS ||
           system_type == ApecMSystemBase::NODAL_SOLUTION ||
           system_type == ApecMSystemBase::SUPPORT_SOLUTION ||
           system_type == ApecMSystemBase::NODAL_B;
}

const std::string &ApecMSystemBase::system_name_from_type(ApecMSystemBase::SystemType system_type)
{
    auto it = ApecMSystemBase::SYSTEM_INFO_MAP.find(system_type);
    if (it == ApecMSystemBase::SYSTEM_INFO_MAP.end())
        libmesh_error_msg("System type not found.");
    return it->second;
}

VectorFieldCoordinateTransformation::VectorFieldCoordinateTransformation(ApecMSystemBase &system_in,
                                                                         const CoordinateSystem &coordinate_system_in,
                                                                         std::vector<std::string> variable_labels_in)
    : FEParentObjectBase<ApecMSystemBase>(&system_in),
      coordinate_system(&coordinate_system_in),
      variable_labels(std::move(variable_labels_in))
{
    for (const auto &label : variable_labels)
        get_parent_object()->global_to_local_coordinate_system_solution(*coordinate_system, label);
}

VectorFieldCoordinateTransformation::~VectorFieldCoordinateTransformation()
{
    if (!get_parent_object() || !coordinate_system)
        return;

    for (const auto &label : variable_labels)
        get_parent_object()->local_to_global_coordinate_system_solution(*coordinate_system, label);
}

template <typename SysLike>
bool ApecMSystem<SysLike>::global_to_local_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                                      const std::string &variable_label,
                                                                      const std::vector<dof_id_type> &selected_dof_ids)
{
    return transform_coordinate_system_solution(
        *this,
        *this->solution,
        coordinate_system,
        variable_label,
        selected_dof_ids,
        true);
}

template <typename SysLike>
bool ApecMSystem<SysLike>::local_to_global_coordinate_system_solution(const CoordinateSystem &coordinate_system,
                                                                      const std::string &variable_label,
                                                                      const std::vector<dof_id_type> &selected_dof_ids)
{
    return transform_coordinate_system_solution(
        *this,
        *this->solution,
        coordinate_system,
        variable_label,
        selected_dof_ids,
        false);
}

void ApecMLinearSystem::release_subset_solve()
{
    if (this->_subset != nullptr)
    {
        delete this->_subset;
        this->_subset = nullptr;
    }
}

void ApecMLinearSystem::compute_solution_time_derivative(const Substep &time_substep)
{
    bool transient_effect = time_substep.get_transient_effect();
    NumericVector<Number> *time_derivative_vector = this->request_vector("solution_time_derivative");
    if (!transient_effect)
    {
        // If transient effects are not considered, the solution time derivative is zero
        this->older_local_solution->zero();
        time_derivative_vector->zero();
    }
    else
    {
        Real theta = time_substep.get_theta();
        Real dt = time_substep.get_dt();
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
    libMesh::out << "Removing constrained DOFs from the system..." << std::endl;
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

    auto element_types_collector = context<ElementTypesFactory>();
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

template <typename SysLike>
void ApecMSystem<SysLike>::dof_numbers_vars(const std::string &variable_name,
                                            bool select_elements,
                                            const std::vector<dof_id_type> &object_ids,
                                            std::vector<dof_id_type> &dof_ids,
                                            bool global)
{
    dof_ids.assign(object_ids.size(), libMesh::invalid_uint);

    EquationSystems &equation_systems = this->get_equation_systems();
    MeshBase &mesh = equation_systems.get_mesh();
    const unsigned int var_number = this->variable_number(variable_name);

    auto append_dofs = [&](const MeshBase &mesh_view)
    {
        for (std::size_t i = 0; i < object_ids.size(); ++i)
        {
            const dof_id_type object_id = object_ids[i];
            if (select_elements)
            {
                const Elem *elem_ptr = mesh_view.query_elem_ptr(object_id);
                dof_ids[i] = elem_ptr->dof_number(this->number(), var_number, 0);
            }
            else
            {
                const Node *node_ptr = mesh_view.query_node_ptr(object_id);
                dof_ids[i] = node_ptr->dof_number(this->number(), var_number, 0);
            }
        }
    };

    MeshBase &serial_mesh = mesh;
    std::unique_ptr<MeshSerializer> serializer;
    if (global)
        serializer = std::make_unique<MeshSerializer>(serial_mesh, true);

    append_dofs(serial_mesh);
}

template class ApecMSystem<libMesh::TransientSystem<libMesh::LinearImplicitSystem>>;
template class ApecMSystem<libMesh::ExplicitSystem>;
