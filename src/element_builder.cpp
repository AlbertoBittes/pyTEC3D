#include "element_builder.h"
#include "binary_reader.h"
#include "apec_m.h"
// #include <libmesh/equation_systems.h>
#include <libmesh/dof_map.h>
#include <libmesh/numeric_vector.h>
// #include <libmesh/string_to_enum.h>
// Define the Finite Element object.
#include <libmesh/fe.h>
#include <libmesh/petsc_linear_solver.h>

// Define Gauss quadrature rules.
#include "session.h"       // Ensure the Session class is included
#include "element_types.h" // Include the header defining ElementType
#include "apec_m_system.h"
#include "element_types_factory.h"

#include "variables_map.h"

//
//
// constructor and init
//
//

ElementBuilder::ElementBuilder(EquationSystems *equation_systems,
                               std::shared_ptr<ContextView> ctx_view)
    : FEParentObjectBase<EquationSystems>(equation_systems),
      FEMService<ElementTypesFactory, MaterialsObject, RealConstantsObject>(ctx_view),
      matrices(MatricesSet::NUM_MATRICES_SETS),
      vectors(VectorsSet::NUM_VECTORS_SETS),
      cylindrical_coordinate_system_()
{
    equation_systems_ = dynamic_cast<ApecMEquationSystems *>(equation_systems);
    if (!equation_systems_)
        libmesh_error_msg("ElementBuilder: equation systems is not ApecMEquationSystems.");

    dof_solution_system_ = dynamic_cast<ApecMLinearSystem *>(&equation_systems_->get_system(ApecMSystemBase::DOF_SOLUTION_SYSTEM));
    if (!dof_solution_system_)
        libmesh_error_msg("ElementBuilder: DOF solution system is not ApecMLinearSystem.");

    nodal_solution_system_ = &equation_systems_->get_system(ApecMSystemBase::NODAL_SOLUTION);
    element_solution_system_ = &equation_systems_->get_system(ApecMSystemBase::ELEMENT_SOLUTION);
    support_solution_system_ = &equation_systems_->get_system(ApecMSystemBase::SUPPORT_SOLUTION);
    nodal_b_system_ = &equation_systems_->get_system(ApecMSystemBase::NODAL_B);
    nodal_loads_system_ = &equation_systems_->get_system(ApecMSystemBase::NODAL_LOADS);
    element_loads_system_ = &equation_systems_->get_system(ApecMSystemBase::ELEMENT_LOADS);
    dof_solution_dof_map_ = &dof_solution_system_->get_dof_map();

    materials_object_ = this->context<MaterialsObject>();
    real_constants_object_ = this->context<RealConstantsObject>();
}

std::unique_ptr<ElementBuilder> ElementBuilder::make_element_builder(ApecMEquationSystems *eq_sys)
{
    auto program = eq_sys->get_parent_object();

    auto analysis = program->get_analysis();
    auto analysis_ctx = analysis->get_context_view();
    auto builder_view = analysis_ctx->template project_view<ElementTypesFactory, MaterialsObject, RealConstantsObject>();

    return std::make_unique<ElementBuilder>(eq_sys, builder_view);
}

//
//
// Methods to size element matrices and vectors
//
//

void ElementBuilder::build_vectors(std::vector<std::vector<unsigned int>> *sizes_map, unsigned int vectors_set)
{
    this->vectors[vectors_set].resize(sizes_map->size());
    for (unsigned int i = 0; i < sizes_map->size(); ++i)
    {
        const auto &size_pair = (*sizes_map)[i];
        unsigned int rows = size_pair[0];
        this->vectors[vectors_set][i].resize(rows);
    }
}

void ElementBuilder::build_matrices(std::vector<std::vector<unsigned int>> *sizes_map, unsigned int matrices_set)
{

    unsigned int num_matrices = static_cast<unsigned int>(sizes_map->size());
    this->matrices[matrices_set].resize(num_matrices);
    for (unsigned int i = 0; i < num_matrices; ++i)
    {
        const auto &size_pair = (*sizes_map)[i];
        unsigned int rows = size_pair[0];
        unsigned int cols = size_pair[1];
        this->matrices[matrices_set][i].resize(rows, cols);
    }
}

void ElementBuilder::prepare_element_matrix_building()
{
    // reinit to be done inside each element
    // as well as the reception of shape functions
    // this->elem_data.element_type->reinit_fe(this->elem);
    std::vector<std::vector<unsigned int>> matrix_sizes = {
        {n_dofs, n_dofs},
        {n_dofs, n_dofs}};

    // One vector of size n_dofs
    std::vector<std::vector<unsigned int>> vector_sizes = {
        {n_dofs}};

    // Usage:
    this->build_matrices(&matrix_sizes, MatricesSet::DOF_SOLUTION_MATRICES);
    this->build_vectors(&vector_sizes, VectorsSet::RHS_VECTORS);
    this->build_matrices(elem_data.element_type->get_solution_matrix_idx_size_map(true).get(), MatricesSet::NODAL_SOLUTION_MATRICES);
    this->build_matrices(elem_data.element_type->get_solution_matrix_idx_size_map(false).get(), MatricesSet::ELEMENT_SOLUTION_MATRICES);
    // this->build_vectors(vector_idx_size_map, 1);
}

//
//
// Setters
//
//

void ElementBuilder::set_element_data(Elem *elem)
{
    auto element_types = this->context<ElementTypesFactory>();

    elem_data.element_type = element_types->get_element_type(elem->subdomain_id());
    if (!elem->infinite())
    {
        elem_data.material_num = elem->get_extra_integer(1);
        elem_data.real_constant_num = elem->get_extra_integer(2);
    }
    elem_data.elem = elem;

    std::vector<dof_id_type> dof_indices;
    dof_solution_dof_map_->dof_indices(elem, dof_indices);
    this->n_dofs = static_cast<unsigned int>(dof_indices.size());
}

//
//
// Methods to build element matrices and rhs
//
//

void ElementBuilder::build_element()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_element_rhs()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_rhs();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_rhs();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_rhs();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

//
//
// solution calculation methods
//
//

void ElementBuilder::build_nodal_solution()
{
    if (!this->compute_nodal_solution_vectors)
        return;

    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_nodal_solution();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_nodal_solution();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_nodal_solution();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_bn_solution()
{
    if (!element_solution_system_ || !nodal_b_system_)
        return;

    const auto api_vars = VariablesMap::get_api_variables_for_element(*elem_data.element_type, ApecMSystemBase::NODAL_B);
    if (api_vars.empty())
        return;

    std::vector<DenseVector<Number>> nodal_b_vectors;
    nodal_b_vectors.resize(api_vars.size());
    DenseVector<Number> bn_values;
    bool have_bn = false;

    auto build_bn = [&]()
    {
        bn_values.resize(0);

        const auto append_component_average = [&](const std::string &be_label, const std::string &bn_label)
        {
            if (!element_solution_system_->has_variable(be_label) || !nodal_b_system_->has_variable(bn_label))
                return;

            std::vector<dof_id_type> dof_indices;
            unsigned int be_var_num = element_solution_system_->variable_number(be_label);
            unsigned int bn_var_num = nodal_b_system_->variable_number(bn_label);
            element_solution_system_->get_dof_map().dof_indices(elem_data.elem, dof_indices, be_var_num);
            if (dof_indices.empty())
                return;

            Real be_val = element_solution_system_->current_local_solution->el(dof_indices[0]);

            for (unsigned int i = 0; i < elem_data.elem->n_nodes(); i++)
            {
                dof_id_type dof_index = libMesh::invalid_uint;
                try
                {
                    dof_index = elem_data.elem->node_ptr(i)->dof_number(ApecMSystemBase::NODAL_B, bn_var_num, 0);
                }
                catch (...)
                {
                    continue;
                }
                if (dof_index < nodal_b_system_->solution->first_local_index() || dof_index >= nodal_b_system_->solution->last_local_index())
                    continue;

                const unsigned int previous_count = bn_node_counts[dof_index];
                const unsigned int next_count = previous_count + 1;
                const Number previous_avg = nodal_b_system_->solution->el(dof_index);
                const Number next_avg = (previous_avg * previous_count + be_val) / static_cast<Real>(next_count);

                bn_node_counts[dof_index] = next_count;
                DenseVector<Number> bn_entry(1);
                bn_entry(0) = next_avg;
                bn_values.append(bn_entry);
            }
        };

        append_component_average("Be_x", "Bn_x");
        append_component_average("Be_y", "Bn_y");
        append_component_average("Be_z", "Bn_z");
        have_bn = true;
    };

    this->fill_vectors_by_api(api_vars, nodal_b_vectors, [&](APIVariableName api_var) -> const DenseVector<Number> *
                              {
                                  switch (api_var)
                                  {
                                  case APIVariableName::Bn:
                                      build_bn();
                                      return have_bn ? &bn_values : nullptr;
                                  default:
                                      return nullptr;
                                  } });

    this->append_and_write_solution(*nodal_b_system_, nodal_b_vectors);
}

void ElementBuilder::build_element_solution()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_element_solution();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_element_solution();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_element_solution();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

void ElementBuilder::build_support_solution()
{
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_support_solution();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_support_solution();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_support_solution();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

std::vector<DenseVector<Number>> ElementBuilder::get_variables_solution(const NumericVector<Number> *solution,
                                                                        int dof_map_index)
{
    System *sys = &equation_systems_->get_system(dof_map_index);
    DofMap *dof_map = (dof_map_index == ApecMSystemBase::DOF_SOLUTION_SYSTEM) ? dof_solution_dof_map_ : &sys->get_dof_map();
    // Cast dof_map_index to SystemType for SYSTEM_INFO_MAP lookup
    ApecMSystemBase::SystemType sys_type = static_cast<ApecMSystemBase::SystemType>(dof_map_index);

    std::vector<DenseVector<Number>> sol_vectors;
    if (!solution || !dof_map)
        return sol_vectors;

    // This is now maybe too slow, but it works for now
    // TODO: Check if there is room for optimization later
    const auto vars = VariablesMap::get_active_variables_for_element(*elem_data.element_type, sys_type);
    sol_vectors.reserve(vars.size());

    for (const auto &var_enum : vars)
    {
        unsigned int system_var_num = sys->variable_number(VariablesMap::get_info_by_key(var_enum).global_label);
        std::vector<dof_id_type> dof_indices_variable;
        dof_map->dof_indices(elem_data.elem, dof_indices_variable, system_var_num);
        DenseVector<Number> local_sol(dof_indices_variable.size());
        for (unsigned int i = 0; i < dof_indices_variable.size(); ++i)
            local_sol(i) = solution->el(dof_indices_variable[i]);
        sol_vectors.emplace_back(std::move(local_sol));
    }

    return sol_vectors;
}

void ElementBuilder::append_and_write_solution(System &system,
                                               const std::vector<DenseVector<Number>> &ordered_vectors)
{
    DenseVector<Number> values;
    for (const auto &vec : ordered_vectors)
    {
        values.append(vec);
    }

    if (values.size() == 0)
        return;

    std::vector<dof_id_type> dof_indices;
    DofMap *dof_map = &system.get_dof_map();
    dof_map->dof_indices(elem_data.elem, dof_indices);

    for (unsigned int i = 0; i < dof_indices.size() && i < values.size(); i++)
    {
        dof_id_type dof_index = dof_indices[i];
        if ((system.solution->first_local_index() <= dof_index) &&
            (dof_index < system.solution->last_local_index()))
            system.solution->set(dof_index, values(i));
    }
}
