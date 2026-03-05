#include "apec_m_evolution_analysis.h"
#include "apec_m_mesh.h"
#include "apec_m.h"

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

#include <libmesh/libmesh.h>
#include <libmesh/mesh_serializer.h>

namespace
{
    constexpr const char *KEY_TAVG_EQ_PSI_DOFS = "tavg.eq.psi.dofs";
    constexpr const char *KEY_TAVG_EQ_B_RAD_DOFS = "tavg.eq.b_rad.dofs";
    constexpr const char *KEY_TAVG_EQ_B_VER_DOFS = "tavg.eq.b_ver.dofs";
    constexpr const char *KEY_TAVG_EQ_PSIE_DOFS = "tavg.eq.psie.dofs";
    constexpr const char *KEY_TAVG_EVOL_PSI_DOF_GROUPS = "tavg.evol.psi.dof_groups";
    constexpr const char *KEY_TAVG_EVOL_B_RAD_DOF_GROUPS = "tavg.evol.b_rad.dof_groups";
    constexpr const char *KEY_TAVG_EVOL_B_VER_DOF_GROUPS = "tavg.evol.b_ver.dof_groups";
    constexpr const char *KEY_TAVG_EVOL_PSIE_DOF_GROUPS = "tavg.evol.psie.dof_groups";

    constexpr const char *KEY_TPROJ_EQ_JS_MAIN_DOFS = "tproj.eq.js.main.dofs";
    constexpr const char *KEY_TPROJ_STATIC_JS_MAIN_DOF_GROUPS = "tproj.static.js.main.dof_groups";
    constexpr const char *KEY_TPROJ_STATIC_JS_C1_DOF_GROUPS = "tproj.static.js.c1.dof_groups";
    constexpr const char *KEY_TPROJ_STATIC_JS_C2_DOF_GROUPS = "tproj.static.js.c2.dof_groups";
    constexpr const char *KEY_TPROJ_STATIC_LOCAL_ELEM_ID_GROUPS = "tproj.static.local.elem_id_groups";

    constexpr const char *KEY_IMPRESSED_EVOL_JT_Z = "impressed.evolution.jt.z";
    constexpr const char *KEY_IMPRESSED_STATIC_JS_Z = "impressed.static.js.z";
    constexpr const char *KEY_IMPRESSED_EVOL_JT_X = "impressed.evolution.jt.x";
    constexpr const char *KEY_IMPRESSED_STATIC_JS_X = "impressed.static.js.x";
    constexpr const char *KEY_IMPRESSED_EVOL_JT_Y = "impressed.evolution.jt.y";
    constexpr const char *KEY_IMPRESSED_STATIC_JS_Y = "impressed.static.js.y";

    constexpr const char *KEY_SET_PLASMA_STATIC_JS_Z = "set_plasma.static.js.z";
    constexpr const char *KEY_SET_PLASMA_EVOL_JS_Z = "set_plasma.evolution.js.z";
    constexpr const char *KEY_SET_PLASMA_STATIC_JS_X = "set_plasma.static.js.x";
    constexpr const char *KEY_SET_PLASMA_EVOL_JS_X = "set_plasma.evolution.js.x";
    constexpr const char *KEY_SET_PLASMA_STATIC_JS_Y = "set_plasma.static.js.y";
    constexpr const char *KEY_SET_PLASMA_EVOL_JS_Y = "set_plasma.evolution.js.y";

    std::vector<std::string> get_j_b_labels(unsigned int mesh_dim)
    {
        if (mesh_dim == 2)
        {
            return {"_z", "_x", "_y"};
        }
        else
        {
            return {"_y", "_x", "_z"};
        }
    }
} // namespace

APECMEvolutionAnalysis::APECMEvolutionAnalysis(APEC_M *program,
                                               std::shared_ptr<AnalysisBase::ContextView> context_view)
    : APECMEquilibriumAnalysis(program, context_view),
      cylindrical_coordinate_system(CylindricalCoordinateSystem())
{
}

void APECMEvolutionAnalysis::set_poloidal_currents()
{
    // libMesh::out << "Projection from main mesh to static mesh..." << std::endl;
    // // Rotate the J and B to cylindrical coordinates
    // auto settings = this->context<ApecMSettingsObject>();

    // APEC_M *program = this->get_parent_object();

    // const std::string &evolution_mesh = settings->evolution_mesh();
    // const std::string &equilibrium_mesh = settings->equilibrium_mesh();

    // auto main_eq_sys = program->get_equation_system(evolution_mesh);
    // auto static_eq_sys = program->get_equation_system(equilibrium_mesh);

    // ApecMExplicitSystem &main_sys_element_sol = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    // ApecMExplicitSystem &static_sys_element_sol = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    // ApecMExplicitSystem &main_sys_element_loads = main_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    // ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    // auto toroidal_proj = toroidal_projection.get();

    // auto j_b_labels_evol = get_j_b_labels(main_eq_sys->get_mesh().mesh_dimension());
    // auto j_b_labels_stat = get_j_b_labels(static_eq_sys->get_mesh().mesh_dimension());
}

void APECMEvolutionAnalysis::set_toroidal_slices()
{
    auto settings = this->context<ApecMSettingsObject>();
    if (!settings)
        libmesh_error_msg("APECMEvolutionAnalysis: APEC settings are not available.");

    toroidal_slices_nodes_ = settings->toroidal_slices_nodes();
    toroidal_slices_elems_ = settings->toroidal_slices_elems();
    dofs_cache_.is_initialized = false;
    dofs_cache_.vectors.clear();
    dofs_cache_.vector_groups.clear();
}

const std::vector<dof_id_type> &APECMEvolutionAnalysis::cached_dofs(const std::string &key) const
{
    const auto it = dofs_cache_.vectors.find(key);
    if (it == dofs_cache_.vectors.end())
        libmesh_error_msg("APECMEvolutionAnalysis: missing cached DOFs for key '" + key + "'.");
    return it->second;
}

const std::vector<std::vector<dof_id_type>> &APECMEvolutionAnalysis::cached_dof_groups(const std::string &key) const
{
    const auto it = dofs_cache_.vector_groups.find(key);
    if (it == dofs_cache_.vector_groups.end())
        libmesh_error_msg("APECMEvolutionAnalysis: missing cached DOF group for key '" + key + "'.");
    return it->second;
}

void APECMEvolutionAnalysis::cache_slice_dofs()
{
    if (dofs_cache_.is_initialized)
        return;

    auto settings = this->context<ApecMSettingsObject>();
    if (!settings)
        libmesh_error_msg("APECMEvolutionAnalysis: APEC settings are not available.");

    APEC_M *program = this->get_parent_object();

    const std::string &equilibrium_mesh_name = settings->equilibrium_mesh();
    const std::string &evolution_mesh_name = settings->evolution_mesh();
    const std::string &static_mesh_name = settings->static_mesh();

    auto equilibrium_eq_sys = program->get_equation_system(equilibrium_mesh_name);
    auto evolution_eq_sys = program->get_equation_system(evolution_mesh_name);
    auto static_eq_sys = program->get_equation_system(static_mesh_name);

    auto equilibrium_mesh = program->get_mesh(equilibrium_mesh_name);
    auto evolution_mesh = program->get_mesh(evolution_mesh_name);
    auto static_mesh = program->get_mesh(static_mesh_name);
    MeshBase &serial_equilibrium_mesh = *equilibrium_mesh;
    MeshBase &serial_evolution_mesh = *evolution_mesh;
    MeshBase &serial_static_mesh = *static_mesh;

    std::unique_ptr<MeshSerializer> equilibrium_serializer = std::make_unique<MeshSerializer>(serial_equilibrium_mesh, true);
    std::unique_ptr<MeshSerializer> evolution_serializer = std::make_unique<MeshSerializer>(serial_evolution_mesh, true);
    std::unique_ptr<MeshSerializer> static_serializer = std::make_unique<MeshSerializer>(serial_static_mesh, true);

    ApecMExplicitSystem &equilibrium_sys_support_sol = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::SUPPORT_SOLUTION);
    ApecMExplicitSystem &equilibrium_sys_element_sol = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &equilibrium_sys_element_loads = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    ApecMExplicitSystem &evolution_sys_support_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::SUPPORT_SOLUTION);
    ApecMExplicitSystem &evolution_sys_element_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &evolution_sys_element_loads = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    const processor_id_type proc_id = evolution_mesh->comm().rank();
    const dof_id_type invalid_dof = libMesh::invalid_uint;

    dofs_cache_.vectors.clear();
    dofs_cache_.vector_groups.clear();

    std::vector<dof_id_type> passive_structure_elems_evolution;
    std::vector<dof_id_type> passive_structure_elems_static;
    std::vector<dof_id_type> coils_elems_evolution;
    std::vector<dof_id_type> coils_elems_static;
    std::vector<dof_id_type> elems_evolution_impressed;
    std::vector<dof_id_type> elems_static_impressed;
    evolution_mesh->select("PASSIVE_STRUCTURES", 1, true, passive_structure_elems_evolution);
    static_mesh->select("PASSIVE_STRUCTURES", 1, true, passive_structure_elems_static);
    evolution_mesh->select("COILS", 1, true, coils_elems_evolution);
    static_mesh->select("COILS", 1, true, coils_elems_static);
    elems_evolution_impressed.insert(elems_evolution_impressed.end(), passive_structure_elems_evolution.begin(), passive_structure_elems_evolution.end());
    elems_evolution_impressed.insert(elems_evolution_impressed.end(), coils_elems_evolution.begin(), coils_elems_evolution.end());
    elems_static_impressed.insert(elems_static_impressed.end(), passive_structure_elems_static.begin(), passive_structure_elems_static.end());
    elems_static_impressed.insert(elems_static_impressed.end(), coils_elems_static.begin(), coils_elems_static.end());

    evolution_sys_element_sol.dof_numbers_vars("Jt_z", elems_evolution_impressed, dofs_cache_.vectors[KEY_IMPRESSED_EVOL_JT_Z], false);
    static_sys_element_loads.dof_numbers_vars("Js_z", elems_static_impressed, dofs_cache_.vectors[KEY_IMPRESSED_STATIC_JS_Z], false);
    if (evolution_mesh->mesh_dimension() == 3)
    {
        evolution_sys_element_sol.dof_numbers_vars("Jt_x", elems_evolution_impressed, dofs_cache_.vectors[KEY_IMPRESSED_EVOL_JT_X], false);
        static_sys_element_loads.dof_numbers_vars("Js_x", elems_static_impressed, dofs_cache_.vectors[KEY_IMPRESSED_STATIC_JS_X], false);
        evolution_sys_element_sol.dof_numbers_vars("Jt_y", elems_evolution_impressed, dofs_cache_.vectors[KEY_IMPRESSED_EVOL_JT_Y], false);
        static_sys_element_loads.dof_numbers_vars("Js_y", elems_static_impressed, dofs_cache_.vectors[KEY_IMPRESSED_STATIC_JS_Y], false);
    }

    std::vector<dof_id_type> elems_evolution_vacuum;
    std::vector<dof_id_type> elems_static_vacuum;
    evolution_mesh->select("VACUUM_CHAMBER", 1, true, elems_evolution_vacuum);
    static_mesh->select("VACUUM_CHAMBER", 1, true, elems_static_vacuum);

    evolution_sys_element_loads.dof_numbers_vars("Js_z", elems_evolution_vacuum, dofs_cache_.vectors[KEY_SET_PLASMA_EVOL_JS_Z], false);
    static_sys_element_loads.dof_numbers_vars("Js_z", elems_static_vacuum, dofs_cache_.vectors[KEY_SET_PLASMA_STATIC_JS_Z], false);
    if (evolution_mesh->mesh_dimension() == 3)
    {
        evolution_sys_element_loads.dof_numbers_vars("Js_x", elems_evolution_vacuum, dofs_cache_.vectors[KEY_SET_PLASMA_EVOL_JS_X], false);
        static_sys_element_loads.dof_numbers_vars("Js_x", elems_static_vacuum, dofs_cache_.vectors[KEY_SET_PLASMA_STATIC_JS_X], false);
        evolution_sys_element_loads.dof_numbers_vars("Js_y", elems_evolution_vacuum, dofs_cache_.vectors[KEY_SET_PLASMA_EVOL_JS_Y], false);
        static_sys_element_loads.dof_numbers_vars("Js_y", elems_static_vacuum, dofs_cache_.vectors[KEY_SET_PLASMA_STATIC_JS_Y], false);
    }

    const auto j_b_labels_equilibrium = get_j_b_labels(equilibrium_mesh->mesh_dimension());
    const auto j_b_labels_evolution = get_j_b_labels(evolution_mesh->mesh_dimension());
    const auto j_b_labels_static = get_j_b_labels(static_mesh->mesh_dimension());

    const unsigned int eq_psi_var = equilibrium_sys_support_sol.variable_number("PSI");
    const unsigned int evol_psi_var = evolution_sys_support_sol.variable_number("PSI");
    const unsigned int eq_b_rad_var = equilibrium_sys_element_sol.variable_number("Be" + j_b_labels_equilibrium[1]);
    const unsigned int eq_b_ver_var = equilibrium_sys_element_sol.variable_number("Be" + j_b_labels_equilibrium[2]);
    const unsigned int eq_psie_var = equilibrium_sys_element_sol.variable_number("PSIe");
    const unsigned int evol_b_rad_var = evolution_sys_element_sol.variable_number("Be" + j_b_labels_evolution[1]);
    const unsigned int evol_b_ver_var = evolution_sys_element_sol.variable_number("Be" + j_b_labels_evolution[2]);
    const unsigned int evol_psie_var = evolution_sys_element_sol.variable_number("PSIe");

    const std::string js_eq_main_name = "Js" + j_b_labels_equilibrium[0];
    const std::string js_static_main_name = "Js" + j_b_labels_static[0];
    const std::string js_static_c1_name = "Js" + j_b_labels_static[1];
    const std::string js_static_c2_name = "Js" + j_b_labels_static[2];
    const unsigned int eq_js_main_var = equilibrium_sys_element_loads.variable_number(js_eq_main_name);
    const unsigned int static_js_main_var = static_sys_element_loads.variable_number(js_static_main_name);
    const bool has_js_static_c1 = static_sys_element_loads.has_variable(js_static_c1_name);
    const bool has_js_static_c2 = static_sys_element_loads.has_variable(js_static_c2_name);
    const unsigned int static_js_c1_var = has_js_static_c1 ? static_sys_element_loads.variable_number(js_static_c1_name) : libMesh::invalid_uint;
    const unsigned int static_js_c2_var = has_js_static_c2 ? static_sys_element_loads.variable_number(js_static_c2_name) : libMesh::invalid_uint;

    const numeric_index_type eq_support_first = equilibrium_sys_support_sol.solution->first_local_index();
    const numeric_index_type eq_support_last = equilibrium_sys_support_sol.solution->last_local_index();
    const numeric_index_type eq_element_first = equilibrium_sys_element_sol.solution->first_local_index();
    const numeric_index_type eq_element_last = equilibrium_sys_element_sol.solution->last_local_index();
    const numeric_index_type evol_support_first = evolution_sys_support_sol.current_local_solution->first_local_index();
    const numeric_index_type evol_support_last = evolution_sys_support_sol.current_local_solution->last_local_index();
    const numeric_index_type evol_element_first = evolution_sys_element_sol.current_local_solution->first_local_index();
    const numeric_index_type evol_element_last = evolution_sys_element_sol.current_local_solution->last_local_index();
    const numeric_index_type eq_loads_first = equilibrium_sys_element_loads.current_local_solution->first_local_index();
    const numeric_index_type eq_loads_last = equilibrium_sys_element_loads.current_local_solution->last_local_index();
    const numeric_index_type static_loads_first = static_sys_element_loads.solution->first_local_index();
    const numeric_index_type static_loads_last = static_sys_element_loads.solution->last_local_index();

    const std::size_t n_node_rows = toroidal_slices_nodes_.empty() ? 0 : toroidal_slices_nodes_[0].size();
    const std::size_t n_elem_rows = toroidal_slices_elems_.empty() ? 0 : toroidal_slices_elems_[0].size();

    auto &eq_psi_dofs = dofs_cache_.vectors[KEY_TAVG_EQ_PSI_DOFS];
    auto &eq_b_rad_dofs = dofs_cache_.vectors[KEY_TAVG_EQ_B_RAD_DOFS];
    auto &eq_b_ver_dofs = dofs_cache_.vectors[KEY_TAVG_EQ_B_VER_DOFS];
    auto &eq_psie_dofs = dofs_cache_.vectors[KEY_TAVG_EQ_PSIE_DOFS];
    eq_psi_dofs.assign(n_node_rows, invalid_dof);
    eq_b_rad_dofs.assign(n_elem_rows, invalid_dof);
    eq_b_ver_dofs.assign(n_elem_rows, invalid_dof);
    eq_psie_dofs.assign(n_elem_rows, invalid_dof);

    auto &evol_psi_groups = dofs_cache_.vector_groups[KEY_TAVG_EVOL_PSI_DOF_GROUPS];
    auto &evol_b_rad_groups = dofs_cache_.vector_groups[KEY_TAVG_EVOL_B_RAD_DOF_GROUPS];
    auto &evol_b_ver_groups = dofs_cache_.vector_groups[KEY_TAVG_EVOL_B_VER_DOF_GROUPS];
    auto &evol_psie_groups = dofs_cache_.vector_groups[KEY_TAVG_EVOL_PSIE_DOF_GROUPS];
    evol_psi_groups.assign(toroidal_slices_nodes_.size() > 0 ? toroidal_slices_nodes_.size() - 1 : 0,
                           std::vector<dof_id_type>(n_node_rows, invalid_dof));
    evol_b_rad_groups.assign(toroidal_slices_elems_.size() > 0 ? toroidal_slices_elems_.size() - 1 : 0,
                             std::vector<dof_id_type>(n_elem_rows, invalid_dof));
    evol_b_ver_groups.assign(toroidal_slices_elems_.size() > 0 ? toroidal_slices_elems_.size() - 1 : 0,
                             std::vector<dof_id_type>(n_elem_rows, invalid_dof));
    evol_psie_groups.assign(toroidal_slices_elems_.size() > 0 ? toroidal_slices_elems_.size() - 1 : 0,
                            std::vector<dof_id_type>(n_elem_rows, invalid_dof));

    if (!toroidal_slices_nodes_.empty())
    {
        const auto &base_node_ids = toroidal_slices_nodes_[0];
        const std::size_t n_base_rows = std::min(n_node_rows, base_node_ids.size());
        for (std::size_t row = 0; row < n_base_rows; ++row)
        {
            const Node *node_ptr = equilibrium_mesh->query_node_ptr(base_node_ids[row]);
            if (!node_ptr || node_ptr->processor_id() != proc_id)
                continue;

            const dof_id_type psi_dof = node_ptr->dof_number(ApecMSystemBase::SUPPORT_SOLUTION, eq_psi_var, 0);
            if (psi_dof >= eq_support_first && psi_dof < eq_support_last)
                eq_psi_dofs[row] = psi_dof;
        }
    }

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_nodes_.size(); ++slice_idx)
    {
        const auto &slice_node_ids = toroidal_slices_nodes_[slice_idx];
        const std::size_t n_rows = std::min(n_node_rows, slice_node_ids.size());
        auto &slice_dofs = evol_psi_groups[slice_idx - 1];
        for (std::size_t row = 0; row < n_rows; ++row)
        {
            const Node *node_ptr = evolution_mesh->query_node_ptr(slice_node_ids[row]);
            if (!node_ptr || node_ptr->processor_id() != proc_id)
                continue;

            const dof_id_type psi_dof = node_ptr->dof_number(ApecMSystemBase::SUPPORT_SOLUTION, evol_psi_var, 0);
            if (psi_dof >= evol_support_first && psi_dof < evol_support_last)
                slice_dofs[row] = psi_dof;
        }
    }

    if (!toroidal_slices_elems_.empty())
    {
        const auto &base_elem_ids = toroidal_slices_elems_[0];
        const std::size_t n_base_rows = std::min(n_elem_rows, base_elem_ids.size());
        for (std::size_t row = 0; row < n_base_rows; ++row)
        {
            const Elem *elem_ptr = equilibrium_mesh->query_elem_ptr(base_elem_ids[row]);
            if (!elem_ptr || elem_ptr->processor_id() != proc_id)
                continue;

            const dof_id_type b_rad_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, eq_b_rad_var, 0);
            const dof_id_type b_ver_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, eq_b_ver_var, 0);
            const dof_id_type psie_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, eq_psie_var, 0);
            if (b_rad_dof >= eq_element_first && b_rad_dof < eq_element_last)
                eq_b_rad_dofs[row] = b_rad_dof;
            if (b_ver_dof >= eq_element_first && b_ver_dof < eq_element_last)
                eq_b_ver_dofs[row] = b_ver_dof;
            if (psie_dof >= eq_element_first && psie_dof < eq_element_last)
                eq_psie_dofs[row] = psie_dof;
        }
    }

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_elems_.size(); ++slice_idx)
    {
        const auto &slice_elem_ids = toroidal_slices_elems_[slice_idx];
        const std::size_t n_rows = std::min(n_elem_rows, slice_elem_ids.size());
        auto &slice_b_rad_dofs = evol_b_rad_groups[slice_idx - 1];
        auto &slice_b_ver_dofs = evol_b_ver_groups[slice_idx - 1];
        auto &slice_psie_dofs = evol_psie_groups[slice_idx - 1];
        for (std::size_t row = 0; row < n_rows; ++row)
        {
            const Elem *elem_ptr = evolution_mesh->query_elem_ptr(slice_elem_ids[row]);
            if (!elem_ptr || elem_ptr->processor_id() != proc_id)
                continue;

            const dof_id_type b_rad_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, evol_b_rad_var, 0);
            const dof_id_type b_ver_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, evol_b_ver_var, 0);
            const dof_id_type psie_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_SOLUTION, evol_psie_var, 0);
            if (b_rad_dof >= evol_element_first && b_rad_dof < evol_element_last)
                slice_b_rad_dofs[row] = b_rad_dof;
            if (b_ver_dof >= evol_element_first && b_ver_dof < evol_element_last)
                slice_b_ver_dofs[row] = b_ver_dof;
            if (psie_dof >= evol_element_first && psie_dof < evol_element_last)
                slice_psie_dofs[row] = psie_dof;
        }
    }

    auto &eq_js_main_dofs = dofs_cache_.vectors[KEY_TPROJ_EQ_JS_MAIN_DOFS];
    eq_js_main_dofs.assign(n_elem_rows, invalid_dof);
    auto &static_js_main_groups = dofs_cache_.vector_groups[KEY_TPROJ_STATIC_JS_MAIN_DOF_GROUPS];
    auto &static_local_elem_id_groups = dofs_cache_.vector_groups[KEY_TPROJ_STATIC_LOCAL_ELEM_ID_GROUPS];
    static_js_main_groups.assign(toroidal_slices_elems_.size() > 0 ? toroidal_slices_elems_.size() - 1 : 0,
                                 std::vector<dof_id_type>(n_elem_rows, invalid_dof));
    static_local_elem_id_groups.assign(toroidal_slices_elems_.size() > 0 ? toroidal_slices_elems_.size() - 1 : 0,
                                       std::vector<dof_id_type>());

    if (has_js_static_c1)
        dofs_cache_.vector_groups[KEY_TPROJ_STATIC_JS_C1_DOF_GROUPS].assign(toroidal_slices_elems_.size() > 0 ? toroidal_slices_elems_.size() - 1 : 0,
                                                                            std::vector<dof_id_type>(n_elem_rows, invalid_dof));
    if (has_js_static_c2)
        dofs_cache_.vector_groups[KEY_TPROJ_STATIC_JS_C2_DOF_GROUPS].assign(toroidal_slices_elems_.size() > 0 ? toroidal_slices_elems_.size() - 1 : 0,
                                                                            std::vector<dof_id_type>(n_elem_rows, invalid_dof));

    if (!toroidal_slices_elems_.empty())
    {
        const auto &base_elem_ids = toroidal_slices_elems_[0];
        const std::size_t n_base_rows = std::min(n_elem_rows, base_elem_ids.size());
        for (std::size_t row = 0; row < n_base_rows; ++row)
        {
            const Elem *elem_ptr = equilibrium_mesh->query_elem_ptr(base_elem_ids[row]);
            if (!elem_ptr || elem_ptr->processor_id() != proc_id)
                continue;

            const dof_id_type js_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_LOADS, eq_js_main_var, 0);
            if (js_dof >= eq_loads_first && js_dof < eq_loads_last)
                eq_js_main_dofs[row] = js_dof;
        }
    }

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_elems_.size(); ++slice_idx)
    {
        const auto &slice_elem_ids = toroidal_slices_elems_[slice_idx];
        const std::size_t n_rows = std::min(n_elem_rows, slice_elem_ids.size());
        auto &slice_js_main_dofs = static_js_main_groups[slice_idx - 1];
        auto &slice_local_elem_ids = static_local_elem_id_groups[slice_idx - 1];
        auto *slice_js_c1_dofs = has_js_static_c1 ? &dofs_cache_.vector_groups[KEY_TPROJ_STATIC_JS_C1_DOF_GROUPS] : nullptr;
        auto *slice_js_c2_dofs = has_js_static_c2 ? &dofs_cache_.vector_groups[KEY_TPROJ_STATIC_JS_C2_DOF_GROUPS] : nullptr;
        for (std::size_t row = 0; row < n_rows; ++row)
        {
            const Elem *elem_ptr = static_mesh->query_elem_ptr(slice_elem_ids[row]);
            if (!elem_ptr || elem_ptr->processor_id() != proc_id)
                continue;

            slice_local_elem_ids.push_back(slice_elem_ids[row]);

            const dof_id_type js_main_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_LOADS, static_js_main_var, 0);
            if (js_main_dof >= static_loads_first && js_main_dof < static_loads_last)
                slice_js_main_dofs[row] = js_main_dof;

            if (has_js_static_c1)
            {
                const dof_id_type js_c1_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_LOADS, static_js_c1_var, 0);
                if (js_c1_dof >= static_loads_first && js_c1_dof < static_loads_last)
                    (*slice_js_c1_dofs)[slice_idx - 1][row] = js_c1_dof;
            }

            if (has_js_static_c2)
            {
                const dof_id_type js_c2_dof = elem_ptr->dof_number(ApecMSystemBase::ELEMENT_LOADS, static_js_c2_var, 0);
                if (js_c2_dof >= static_loads_first && js_c2_dof < static_loads_last)
                    (*slice_js_c2_dofs)[slice_idx - 1][row] = js_c2_dof;
            }
        }
    }

    dofs_cache_.is_initialized = true;
}

void APECMEvolutionAnalysis::toroidal_average()
{
    libMesh::out << "Toroidal average from evolution mesh to equilibrium mesh..." << std::endl;
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &equilibrium_mesh_name = settings->equilibrium_mesh();
    const std::string &evolution_mesh_name = settings->evolution_mesh();

    auto equilibrium_eq_sys = program->get_equation_system(equilibrium_mesh_name);
    auto evolution_eq_sys = program->get_equation_system(evolution_mesh_name);

    auto equilibrium_mesh = program->get_mesh(equilibrium_mesh_name);
    auto evolution_mesh = program->get_mesh(evolution_mesh_name);

    ApecMExplicitSystem &equilibrium_sys_support_sol = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::SUPPORT_SOLUTION);
    ApecMExplicitSystem &equilibrium_sys_element_sol = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &evolution_sys_support_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::SUPPORT_SOLUTION);
    ApecMExplicitSystem &evolution_sys_element_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);

    this->cache_slice_dofs();

    auto &comm = evolution_mesh->comm();

    const dof_id_type invalid_dof = libMesh::invalid_uint;
    const auto &eq_psi_dofs = this->cached_dofs(KEY_TAVG_EQ_PSI_DOFS);
    const auto &eq_b_rad_dofs = this->cached_dofs(KEY_TAVG_EQ_B_RAD_DOFS);
    const auto &eq_b_ver_dofs = this->cached_dofs(KEY_TAVG_EQ_B_VER_DOFS);
    const auto &eq_psie_dofs = this->cached_dofs(KEY_TAVG_EQ_PSIE_DOFS);
    const auto &evol_psi_groups = this->cached_dof_groups(KEY_TAVG_EVOL_PSI_DOF_GROUPS);
    const auto &evol_b_rad_groups = this->cached_dof_groups(KEY_TAVG_EVOL_B_RAD_DOF_GROUPS);
    const auto &evol_b_ver_groups = this->cached_dof_groups(KEY_TAVG_EVOL_B_VER_DOF_GROUPS);
    const auto &evol_psie_groups = this->cached_dof_groups(KEY_TAVG_EVOL_PSIE_DOF_GROUPS);

    const std::size_t n_node_rows = eq_psi_dofs.size();
    const std::size_t n_elem_rows = eq_b_rad_dofs.size();

    std::vector<Real> psi_sums(n_node_rows, 0.0);
    std::vector<Real> psi_counts(n_node_rows, 0.0);
    std::vector<Real> b_rad_sums(n_elem_rows, 0.0);
    std::vector<Real> b_ver_sums(n_elem_rows, 0.0);
    std::vector<Real> psie_sums(n_elem_rows, 0.0);
    std::vector<Real> elem_counts(n_elem_rows, 0.0);

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_nodes_.size(); ++slice_idx)
    {
        const auto &slice_psi_dofs = evol_psi_groups[slice_idx - 1];
        const std::size_t n_rows = std::min(n_node_rows, slice_psi_dofs.size());
        for (std::size_t row = 0; row < n_rows; ++row)
        {
            const dof_id_type psi_dof = slice_psi_dofs[row];
            if (psi_dof == invalid_dof)
                continue;

            psi_sums[row] += evolution_sys_support_sol.current_local_solution->el(psi_dof);
            psi_counts[row] += 1.0;
        }
    }

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_elems_.size(); ++slice_idx)
    {
        const auto &slice_b_rad_dofs = evol_b_rad_groups[slice_idx - 1];
        const auto &slice_b_ver_dofs = evol_b_ver_groups[slice_idx - 1];
        const auto &slice_psie_dofs = evol_psie_groups[slice_idx - 1];
        const std::size_t n_rows = std::min({n_elem_rows, slice_b_rad_dofs.size(), slice_b_ver_dofs.size(), slice_psie_dofs.size()});
        for (std::size_t row = 0; row < n_rows; ++row)
        {
            const dof_id_type b_rad_dof = slice_b_rad_dofs[row];
            const dof_id_type b_ver_dof = slice_b_ver_dofs[row];
            const dof_id_type psie_dof = slice_psie_dofs[row];
            if (b_rad_dof == invalid_dof || b_ver_dof == invalid_dof || psie_dof == invalid_dof)
                continue;

            b_rad_sums[row] += evolution_sys_element_sol.current_local_solution->el(b_rad_dof);
            b_ver_sums[row] += evolution_sys_element_sol.current_local_solution->el(b_ver_dof);
            psie_sums[row] += evolution_sys_element_sol.current_local_solution->el(psie_dof);
            elem_counts[row] += 1.0;
        }
    }

    comm.sum(psi_sums);
    comm.sum(psi_counts);
    comm.sum(b_rad_sums);
    comm.sum(b_ver_sums);
    comm.sum(psie_sums);
    comm.sum(elem_counts);

    equilibrium_sys_element_sol.global_to_local_coordinate_system_solution(cylindrical_coordinate_system, "Be");

    for (std::size_t row = 0; row < n_node_rows; ++row)
    {
        if (psi_counts[row] <= 0.0)
            continue;

        const dof_id_type psi_dof = eq_psi_dofs[row];
        if (psi_dof == invalid_dof)
            continue;

        equilibrium_sys_support_sol.solution->set(psi_dof, psi_sums[row] / psi_counts[row]);
    }

    for (std::size_t row = 0; row < n_elem_rows; ++row)
    {
        if (elem_counts[row] <= 0.0)
            continue;

        const dof_id_type b_rad_dof = eq_b_rad_dofs[row];
        const dof_id_type b_ver_dof = eq_b_ver_dofs[row];
        const dof_id_type psie_dof = eq_psie_dofs[row];
        if (b_rad_dof == invalid_dof || b_ver_dof == invalid_dof || psie_dof == invalid_dof)
            continue;

        const Real inv_count = 1.0 / elem_counts[row];
        equilibrium_sys_element_sol.solution->set(b_rad_dof, b_rad_sums[row] * inv_count);
        equilibrium_sys_element_sol.solution->set(b_ver_dof, b_ver_sums[row] * inv_count);
        equilibrium_sys_element_sol.solution->set(psie_dof, psie_sums[row] * inv_count);
    }

    equilibrium_sys_support_sol.solution->close();
    equilibrium_sys_support_sol.update();
    equilibrium_sys_element_sol.solution->close();
    equilibrium_sys_element_sol.update();
}

void APECMEvolutionAnalysis::impressed_eddy_currents_coils_currents()
{
    libMesh::out << "Projection from main mesh to static mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &evolution_mesh_name = settings->evolution_mesh();
    const std::string &static_mesh_name = settings->static_mesh();

    auto evolution_eq_sys = program->get_equation_system(evolution_mesh_name);
    auto static_eq_sys = program->get_equation_system(static_mesh_name);

    auto evolution_mesh = program->get_mesh(evolution_mesh_name);

    ApecMExplicitSystem &evolution_sys_element_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    this->cache_slice_dofs();

    unsigned int dim = evolution_mesh->mesh_dimension();

    std::vector<Real> jt_z;
    std::vector<Real> jt_x;
    std::vector<Real> jt_y;

    const auto &dofs_jt_z = this->cached_dofs(KEY_IMPRESSED_EVOL_JT_Z);
    const auto &dofs_js_z = this->cached_dofs(KEY_IMPRESSED_STATIC_JS_Z);

    evolution_sys_element_sol.solution->get(dofs_jt_z, jt_z);

    static_sys_element_loads.solution->insert(jt_z, dofs_js_z);

    if (dim == 3)
    {
        const auto &dofs_jt_x = this->cached_dofs(KEY_IMPRESSED_EVOL_JT_X);
        const auto &dofs_js_x = this->cached_dofs(KEY_IMPRESSED_STATIC_JS_X);
        const auto &dofs_jt_y = this->cached_dofs(KEY_IMPRESSED_EVOL_JT_Y);
        const auto &dofs_js_y = this->cached_dofs(KEY_IMPRESSED_STATIC_JS_Y);

        evolution_sys_element_sol.solution->get(dofs_jt_x, jt_x);
        evolution_sys_element_sol.solution->get(dofs_jt_y, jt_y);

        static_sys_element_loads.solution->insert(jt_x, dofs_js_x);
        static_sys_element_loads.solution->insert(jt_y, dofs_js_y);
    }

    static_sys_element_loads.solution->close();
    static_sys_element_loads.update();
}

void APECMEvolutionAnalysis::set_plasma_current()
{
    libMesh::out << "Projection from static mesh to main mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &static_mesh_name = settings->static_mesh();
    const std::string &evolution_mesh_name = settings->evolution_mesh();

    auto static_eq_sys = program->get_equation_system(static_mesh_name);
    auto evolution_eq_sys = program->get_equation_system(evolution_mesh_name);

    auto static_mesh = program->get_mesh(static_mesh_name);
    auto evolution_mesh = program->get_mesh(evolution_mesh_name);

    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    ApecMExplicitSystem &evolution_sys_element_loads = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    this->cache_slice_dofs();

    unsigned int dim = evolution_mesh->mesh_dimension();

    std::vector<Real> js_z;
    std::vector<Real> js_x;
    std::vector<Real> js_y;

    const auto &dofs_js_z_evol = this->cached_dofs(KEY_SET_PLASMA_EVOL_JS_Z);
    const auto &dofs_js_z_static = this->cached_dofs(KEY_SET_PLASMA_STATIC_JS_Z);

    static_sys_element_loads.solution->get(dofs_js_z_static, js_z);

    evolution_sys_element_loads.solution->insert(js_z, dofs_js_z_evol);

    if (dim == 3)
    {
        const auto &dofs_js_x_evol = this->cached_dofs(KEY_SET_PLASMA_EVOL_JS_X);
        const auto &dofs_js_x_static = this->cached_dofs(KEY_SET_PLASMA_STATIC_JS_X);
        const auto &dofs_js_y_evol = this->cached_dofs(KEY_SET_PLASMA_EVOL_JS_Y);
        const auto &dofs_js_y_static = this->cached_dofs(KEY_SET_PLASMA_STATIC_JS_Y);

        static_sys_element_loads.solution->get(dofs_js_x_static, js_x);
        static_sys_element_loads.solution->get(dofs_js_y_static, js_y);

        evolution_sys_element_loads.solution->insert(js_x, dofs_js_x_evol);
        evolution_sys_element_loads.solution->insert(js_y, dofs_js_y_evol);
    }

    evolution_sys_element_loads.solution->close();
    evolution_sys_element_loads.update();
}

void APECMEvolutionAnalysis::toroidal_projection_plasma_current()
{
    libMesh::out << "Projection from static mesh to main mesh..." << std::endl;
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &static_mesh_name = settings->static_mesh();
    const std::string &equilibrium_mesh_name = settings->equilibrium_mesh();

    auto static_eq_sys = program->get_equation_system(static_mesh_name);
    auto equilibrium_eq_sys = program->get_equation_system(equilibrium_mesh_name);

    auto static_mesh = program->get_mesh(static_mesh_name);
    auto equilibrium_mesh = program->get_mesh(equilibrium_mesh_name);

    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    ApecMExplicitSystem &equilibrium_sys_element_loads = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    this->cache_slice_dofs();

    auto &comm = static_mesh->comm();
    auto j_b_labels_static = get_j_b_labels(static_mesh->mesh_dimension());
    const std::string js_component_1 = "Js" + j_b_labels_static[1];
    const std::string js_component_2 = "Js" + j_b_labels_static[2];
    const bool has_js_component_1 = static_sys_element_loads.has_variable(js_component_1);
    const bool has_js_component_2 = static_sys_element_loads.has_variable(js_component_2);
    const dof_id_type invalid_dof = libMesh::invalid_uint;
    const auto &eq_js_main_dofs = this->cached_dofs(KEY_TPROJ_EQ_JS_MAIN_DOFS);
    const auto &static_js_main_groups = this->cached_dof_groups(KEY_TPROJ_STATIC_JS_MAIN_DOF_GROUPS);
    const auto &static_local_elem_id_groups = this->cached_dof_groups(KEY_TPROJ_STATIC_LOCAL_ELEM_ID_GROUPS);
    const std::vector<std::vector<dof_id_type>> *static_js_c1_groups = has_js_component_1 ? &this->cached_dof_groups(KEY_TPROJ_STATIC_JS_C1_DOF_GROUPS) : nullptr;
    const std::vector<std::vector<dof_id_type>> *static_js_c2_groups = has_js_component_2 ? &this->cached_dof_groups(KEY_TPROJ_STATIC_JS_C2_DOF_GROUPS) : nullptr;

    const std::size_t n_rows = eq_js_main_dofs.size();
    std::vector<Real> js_main_sums(n_rows, 0.0);
    std::vector<Real> js_main_counts(n_rows, 0.0);

    for (std::size_t row = 0; row < n_rows; ++row)
    {
        const dof_id_type eq_dof = eq_js_main_dofs[row];
        if (eq_dof == invalid_dof)
            continue;
        js_main_sums[row] += equilibrium_sys_element_loads.current_local_solution->el(eq_dof);
        js_main_counts[row] += 1.0;
    }

    comm.sum(js_main_sums);
    comm.sum(js_main_counts);

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_elems_.size(); ++slice_idx)
    {
        const auto &slice_js_main_dofs = static_js_main_groups[slice_idx - 1];
        const std::size_t n_slice_rows = std::min(n_rows, slice_js_main_dofs.size());

        for (std::size_t row = 0; row < n_slice_rows; ++row)
        {
            const dof_id_type js_main_dof = slice_js_main_dofs[row];
            if (js_main_dof == invalid_dof)
                continue;

            const Real js_value = (js_main_counts[row] > 0.0) ? (js_main_sums[row] / js_main_counts[row]) : 0.0;
            static_sys_element_loads.solution->set(js_main_dof, js_value);

            if (has_js_component_1)
            {
                const dof_id_type js1_dof = (*static_js_c1_groups)[slice_idx - 1][row];
                if (js1_dof != invalid_dof)
                    static_sys_element_loads.solution->set(js1_dof, 0.0);
            }

            if (has_js_component_2)
            {
                const dof_id_type js2_dof = (*static_js_c2_groups)[slice_idx - 1][row];
                if (js2_dof != invalid_dof)
                    static_sys_element_loads.solution->set(js2_dof, 0.0);
            }
        }

        const auto &local_slice_elem_ids = static_local_elem_id_groups[slice_idx - 1];
        if (!local_slice_elem_ids.empty())
            static_sys_element_loads.local_to_global_coordinate_system_solution(cylindrical_coordinate_system, "Js", local_slice_elem_ids);
    }

    static_sys_element_loads.solution->close();
    static_sys_element_loads.update();
}

bool APECMEvolutionAnalysis::picard_loop(Substep &substep)
{
    APEC_M *program = this->get_parent_object();

    auto settings = this->context<ApecMSettingsObject>();
    bool print_full_output = settings->print_full_output();

    const std::string &evolution_mesh = settings->evolution_mesh();
    const std::string &static_mesh = settings->static_mesh();
    const std::string &equilibrium_mesh = settings->equilibrium_mesh();

    if (first_absolute_iteration)
        set_toroidal_slices();

    if (substep.get_parent_object()->get_previous_time_step() == nullptr)
    {
        program->set_active_mesh(equilibrium_mesh);
        bool conv = APECMEquilibriumAnalysis::picard_loop(substep);
        program->set_active_mesh(evolution_mesh);
        // set_plasma_current();
        // set_poloidal_currents();
        solve_equation_systems(substep, false);
        return conv;
    }

    program->set_active_mesh(evolution_mesh);
    set_plasma_current();
    set_poloidal_currents();

    libMesh::out << "Solving evolution step..." << std::endl;
    solve_equation_systems(substep, false);

    if (print_full_output)
        print_output(substep, true);

    program->set_active_mesh(equilibrium_mesh);

    toroidal_average(); // With this you set the averaged psi, eddy currents, and B fields on the static mesh

    find_separatrix();

    apply_jplasma_gradshafranov();

    program->set_active_mesh(static_mesh);

    toroidal_projection_plasma_current();

    impressed_eddy_currents_coils_currents();

    libMesh::out << "Solving equilibrium step..." << std::endl;
    solve_equation_systems(substep, false);

    if (print_full_output)
        print_output(substep, true);

    bool converged = check_poloidal_flux_convergence();

    append_iteration_diagnostics(substep, converged);

    if (converged)
        return true;

    return false;
}
