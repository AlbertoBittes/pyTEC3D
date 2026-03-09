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
    nodal_loads_system_ = &equation_systems_->get_system(ApecMSystemBase::NODAL_LOADS);
    element_loads_system_ = &equation_systems_->get_system(ApecMSystemBase::ELEMENT_LOADS);
    dof_solution_dof_map_ = &dof_solution_system_->get_dof_map();

    materials_object_ = this->context<MaterialsObject>();
    real_constants_object_ = this->context<RealConstantsObject>();

    const auto &mesh = equation_systems_->get_mesh();
    mat_index = mesh.get_elem_integer_index("ansys_material_type");
    real_constant_index = mesh.get_elem_integer_index("ansys_real_constant");
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
        elem_data.material_num = elem->get_extra_integer(mat_index);
        elem_data.real_constant_num = elem->get_extra_integer(real_constant_index);
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
