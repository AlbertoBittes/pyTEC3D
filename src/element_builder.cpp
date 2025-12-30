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
      FEMService<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject>(ctx_view),
      matrices(MatricesSet::NUM_MATRICES_SETS),
      vectors(VectorsSet::NUM_VECTORS_SETS)
{
}

std::unique_ptr<ElementBuilder> ElementBuilder::make_element_builder(ApecMEquationSystems *eq_sys)
{
    auto program = eq_sys->get_parent_object();

    auto analysis = program->get_analysis();
    auto analysis_ctx = analysis->get_context_view();
    auto builder_view = analysis_ctx->template project_view<ElementTypesFactory, TimeStep, MaterialsObject, RealConstantsObject>();

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
    case ElementTypeEnum::INFIN110:
        this->build_infin110();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112();
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
    case ElementTypeEnum::INFIN110:
        this->build_infin110_rhs();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111_rhs();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112_rhs();
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

void ElementBuilder::build_solution_vectors()
{
    // Build the nodal solution vector
    // This is a placeholder for the actual implementation
    auto etype = static_cast<ElementTypeEnum>(elem_data.element_type->get_keyopts()[0]);

    switch (etype)
    {
    case ElementTypeEnum::SOLID97:
        this->build_solid97_solution_vectors();
        break;
    case ElementTypeEnum::SOLID236:
        this->build_solid236_solution_vectors();
        break;
    case ElementTypeEnum::PLANE233:
        this->build_plane233_solution_vectors();
        break;
    case ElementTypeEnum::INFIN110:
        this->build_infin110_solution_vectors();
        break;
    case ElementTypeEnum::INFIN111:
        this->build_infin111_solution_vectors();
        break;
    case ElementTypeEnum::INFIN112:
        this->build_infin112_solution_vectors();
        break;
    default:
        libmesh_error_msg("Element type not supported: " + std::to_string(static_cast<unsigned int>(etype)));
    }
}

bool ElementBuilder::build_solution_vectors_from_emat_files(dof_id_type element_id,
                                                            const std::filesystem::path &emat_folderpath)
{
    auto eq_sys = this->get_parent_object();

    const auto elem_id_uint = static_cast<unsigned int>(element_id);
    auto &mesh = eq_sys->get_mesh();
    Elem *elem_ptr = mesh.query_elem_ptr(elem_id_uint);
    if (!elem_ptr)
        return false;

    this->set_element_data(elem_ptr);

    BinaryReader reader;
    if (!reader.read_solution_matrices_for_element(*this, element_id, emat_folderpath))
        return false;

    this->build_solution_vectors();
    return true;
}

void ElementBuilder::fill_solution_vector_support(const NumericVector<Number> *solution, int dof_map_index)
{
    solution_vector_support.clear();

    auto eq_sys = this->get_parent_object();

    System *sys = &eq_sys->get_system(dof_map_index);
    DofMap *dof_map = &sys->get_dof_map();
    // Cast dof_map_index to SystemType for SYSTEM_INFO_MAP lookup
    ApecMSystemBase::SystemType sys_type = static_cast<ApecMSystemBase::SystemType>(dof_map_index);

    // This is now maybe too slow, but it works for now
    // TODO: Check if there is room for optimization later
    const auto vars = VariablesMap::get_active_variables_for_element(*elem_data.element_type, sys_type);

    for (const auto &var_enum : vars)
    {
        unsigned int system_var_num = sys->variable_number(VariablesMap::get_info_by_key(var_enum).global_label);
        std::vector<dof_id_type> dof_indices_variable;
        if (dof_map)
        {
            dof_map->dof_indices(elem_data.elem, dof_indices_variable, system_var_num);
        }
        solution_vector_support.dof_indices_var.emplace_back(std::move(dof_indices_variable));

        if (solution && dof_map)
        {
            const auto &dof_indices_ref = solution_vector_support.dof_indices_var.back();
            DenseVector<Number> local_sol(dof_indices_ref.size());
            for (unsigned int i = 0; i < dof_indices_ref.size(); ++i)
            {
                local_sol(i) = solution->el(dof_indices_ref[i]);
            }
            solution_vector_support.sol_vectors.emplace_back(std::move(local_sol));
        }
    }
}
