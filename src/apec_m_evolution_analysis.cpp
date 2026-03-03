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
}

void APECMEvolutionAnalysis::toroidal_average()
{
    libMesh::out << "Projection from main mesh to static mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &equilibrium_mesh_name = settings->equilibrium_mesh();
    const std::string &evolution_mesh_name = settings->evolution_mesh();

    auto equilibrium_eq_sys = program->get_equation_system(equilibrium_mesh_name);
    auto evolution_eq_sys = program->get_equation_system(evolution_mesh_name);

    auto equilibrium_mesh = program->get_mesh(equilibrium_mesh_name);
    auto evolution_mesh = program->get_mesh(evolution_mesh_name);
    MeshBase &serial_equilibrium_mesh = *equilibrium_mesh;
    MeshBase &serial_evolution_mesh = *evolution_mesh;
    std::unique_ptr<MeshSerializer> equilibrium_serializer = std::make_unique<MeshSerializer>(serial_equilibrium_mesh, true);
    std::unique_ptr<MeshSerializer> evolution_serializer = std::make_unique<MeshSerializer>(serial_evolution_mesh, true);

    ApecMExplicitSystem &equilibrium_sys_nodal_sol = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::NODAL_B);
    ApecMExplicitSystem &equilibrium_sys_support_sol = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::SUPPORT_SOLUTION);
    ApecMExplicitSystem &equilibrium_sys_element_sol = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &evolution_sys_nodal_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::NODAL_B);
    ApecMExplicitSystem &evolution_sys_support_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::SUPPORT_SOLUTION);
    ApecMExplicitSystem &evolution_sys_element_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);

    equilibrium_sys_nodal_sol.solution->localize(*equilibrium_sys_nodal_sol.current_local_solution);
    equilibrium_sys_support_sol.solution->localize(*equilibrium_sys_support_sol.current_local_solution);
    equilibrium_sys_element_sol.solution->localize(*equilibrium_sys_element_sol.current_local_solution);
    evolution_sys_nodal_sol.solution->localize(*evolution_sys_nodal_sol.current_local_solution);
    evolution_sys_support_sol.solution->localize(*evolution_sys_support_sol.current_local_solution);
    evolution_sys_element_sol.solution->localize(*evolution_sys_element_sol.current_local_solution);

    auto j_b_labels_equilibrium = get_j_b_labels(equilibrium_mesh->mesh_dimension());
    auto j_b_labels_evolution = get_j_b_labels(evolution_mesh->mesh_dimension());

    std::vector<std::vector<dof_id_type>> dofs_b_rad_evolution(toroidal_slices_nodes_.size() - 1);
    std::vector<std::vector<dof_id_type>> dofs_b_ver_evolution(toroidal_slices_nodes_.size() - 1);
    std::vector<std::vector<dof_id_type>> dofs_psi_evolution(toroidal_slices_nodes_.size() - 1);
    std::vector<std::vector<dof_id_type>> dofs_psie_evolution(toroidal_slices_elems_.size() - 1);

    std::vector<dof_id_type> dofs_b_rad_equilibrium;
    std::vector<dof_id_type> dofs_b_ver_equilibrium;
    std::vector<dof_id_type> dofs_psi_equilibrium;
    std::vector<dof_id_type> dofs_psie_equilibrium;

    equilibrium_sys_nodal_sol.global_to_local_coordinate_system_solution(cylindrical_coordinate_system, "Bn");

    equilibrium_sys_nodal_sol.dof_numbers_vars("Bn" + j_b_labels_equilibrium[1], false, toroidal_slices_nodes_[0], dofs_b_rad_equilibrium, true);
    equilibrium_sys_nodal_sol.dof_numbers_vars("Bn" + j_b_labels_equilibrium[2], false, toroidal_slices_nodes_[0], dofs_b_ver_equilibrium, true);
    equilibrium_sys_support_sol.dof_numbers_vars("PSI", false, toroidal_slices_nodes_[0], dofs_psi_equilibrium, true);
    equilibrium_sys_element_sol.dof_numbers_vars("PSIe", true, toroidal_slices_elems_[0], dofs_psie_equilibrium, true);

    std::vector<Real> b_rad_avg(dofs_b_rad_equilibrium.size(), 0.0);
    std::vector<Real> b_ver_avg(dofs_b_ver_equilibrium.size(), 0.0);
    std::vector<Real> psi_avg(dofs_psi_equilibrium.size(), 0.0);
    std::vector<Real> psie_avg(dofs_psie_equilibrium.size(), 0.0);

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_nodes_.size(); ++slice_idx)
    {
        const auto &slice_node_ids = toroidal_slices_nodes_[slice_idx];
        evolution_sys_nodal_sol.dof_numbers_vars("Bn" + j_b_labels_evolution[1], false, slice_node_ids, dofs_b_rad_evolution[slice_idx - 1], true);
        evolution_sys_nodal_sol.dof_numbers_vars("Bn" + j_b_labels_evolution[2], false, slice_node_ids, dofs_b_ver_evolution[slice_idx - 1], true);
        evolution_sys_support_sol.dof_numbers_vars("PSI", false, slice_node_ids, dofs_psi_evolution[slice_idx - 1], true);

        std::vector<Real> b_rad_slice;
        std::vector<Real> b_ver_slice;
        std::vector<Real> psi_slice;

        evolution_sys_nodal_sol.current_local_solution->get(dofs_b_rad_evolution[slice_idx - 1], b_rad_slice);
        evolution_sys_nodal_sol.current_local_solution->get(dofs_b_ver_evolution[slice_idx - 1], b_ver_slice);
        evolution_sys_support_sol.current_local_solution->get(dofs_psi_evolution[slice_idx - 1], psi_slice);

        for (unsigned int i = 0; i < b_rad_avg.size(); ++i)
            b_rad_avg[i] += b_rad_slice[i];
        for (unsigned int i = 0; i < b_ver_avg.size(); ++i)
            b_ver_avg[i] += b_ver_slice[i];
        for (unsigned int i = 0; i < psi_avg.size(); ++i)
            psi_avg[i] += psi_slice[i];
    }

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_elems_.size(); ++slice_idx)
    {
        const auto &slice_elem_ids = toroidal_slices_elems_[slice_idx];
        evolution_sys_element_sol.dof_numbers_vars("PSIe", true, slice_elem_ids, dofs_psie_evolution[slice_idx - 1], true);

        std::vector<Real> psie_slice;
        evolution_sys_element_sol.current_local_solution->get(dofs_psie_evolution[slice_idx - 1], psie_slice);

        for (unsigned int i = 0; i < psie_avg.size(); ++i)
            psie_avg[i] += psie_slice[i];
    }

    for (unsigned int i = 0; i < b_rad_avg.size(); ++i)
        b_rad_avg[i] /= static_cast<Real>(toroidal_slices_nodes_.size() - 1);
    for (unsigned int i = 0; i < b_ver_avg.size(); ++i)
        b_ver_avg[i] /= static_cast<Real>(toroidal_slices_nodes_.size() - 1);
    for (unsigned int i = 0; i < psi_avg.size(); ++i)
        psi_avg[i] /= static_cast<Real>(toroidal_slices_nodes_.size() - 1);
    for (unsigned int i = 0; i < psie_avg.size(); ++i)
        psie_avg[i] /= static_cast<Real>(toroidal_slices_elems_.size() - 1);

    equilibrium_sys_nodal_sol.solution->insert(b_rad_avg, dofs_b_rad_equilibrium);
    equilibrium_sys_nodal_sol.solution->insert(b_ver_avg, dofs_b_ver_equilibrium);
    equilibrium_sys_support_sol.solution->insert(psi_avg, dofs_psi_equilibrium);
    equilibrium_sys_element_sol.solution->insert(psie_avg, dofs_psie_equilibrium);

    equilibrium_sys_nodal_sol.solution->close();
    equilibrium_sys_nodal_sol.update();
    equilibrium_sys_support_sol.solution->close();
    equilibrium_sys_support_sol.update();
    equilibrium_sys_element_sol.solution->close();
    equilibrium_sys_element_sol.update();
}

void APECMEvolutionAnalysis::impressed_eddy_currents_coils_currents()
{
    libMesh::out << "Projection from static mesh to main mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &evolution_mesh_name = settings->evolution_mesh();
    const std::string &static_mesh_name = settings->static_mesh();

    auto evolution_eq_sys = program->get_equation_system(evolution_mesh_name);
    auto static_eq_sys = program->get_equation_system(static_mesh_name);

    auto evolution_mesh = program->get_mesh(evolution_mesh_name);
    auto static_mesh = program->get_mesh(static_mesh_name);

    ApecMExplicitSystem &evolution_sys_element_sol = evolution_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_SOLUTION);
    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    unsigned int dim = evolution_mesh->mesh_dimension();

    std::vector<dof_id_type> passive_structure_elems_evolution;
    std::vector<dof_id_type> passive_structure_elems_static;

    std::vector<dof_id_type> coils_elems_evolution;
    std::vector<dof_id_type> coils_elems_static;

    std::vector<dof_id_type> elems_evolution;
    std::vector<dof_id_type> elems_static;

    std::vector<dof_id_type> dofs_jt_z;
    std::vector<dof_id_type> dofs_js_z;
    std::vector<dof_id_type> dofs_jt_x;
    std::vector<dof_id_type> dofs_js_x;
    std::vector<dof_id_type> dofs_jt_y;
    std::vector<dof_id_type> dofs_js_y;

    std::vector<Real> jt_z;
    std::vector<Real> jt_x;
    std::vector<Real> jt_y;

    evolution_mesh->select("PASSIVE_STRUCTURES", 1, true, passive_structure_elems_evolution);
    static_mesh->select("PASSIVE_STRUCTURES", 1, true, passive_structure_elems_static);

    evolution_mesh->select("COILS", 1, true, coils_elems_evolution);
    static_mesh->select("COILS", 1, true, coils_elems_static);

    elems_evolution.insert(elems_evolution.end(), passive_structure_elems_evolution.begin(), passive_structure_elems_evolution.end());
    elems_evolution.insert(elems_evolution.end(), coils_elems_evolution.begin(), coils_elems_evolution.end());
    elems_static.insert(elems_static.end(), passive_structure_elems_static.begin(), passive_structure_elems_static.end());
    elems_static.insert(elems_static.end(), coils_elems_static.begin(), coils_elems_static.end());

    evolution_sys_element_sol.dof_numbers_vars("Jt_z", true, elems_evolution, dofs_jt_z, true);
    static_sys_element_loads.dof_numbers_vars("Js_z", true, elems_static, dofs_js_z, true);

    evolution_sys_element_sol.solution->get(dofs_jt_z, jt_z);

    static_sys_element_loads.solution->insert(jt_z, dofs_js_z);

    if (dim == 3)
    {
        evolution_sys_element_sol.dof_numbers_vars("Jt_x", true, elems_evolution, dofs_jt_x, true);
        static_sys_element_loads.dof_numbers_vars("Js_x", true, elems_static, dofs_js_x, true);
        evolution_sys_element_sol.dof_numbers_vars("Jt_y", true, elems_evolution, dofs_jt_y, true);
        static_sys_element_loads.dof_numbers_vars("Js_y", true, elems_static, dofs_js_y, true);

        evolution_sys_element_sol.solution->get(dofs_jt_x, jt_x);
        evolution_sys_element_sol.solution->get(dofs_jt_y, jt_y);

        static_sys_element_loads.solution->insert(jt_x, dofs_js_x);
        static_sys_element_loads.solution->insert(jt_y, dofs_js_y);
    }

    static_sys_element_loads.solution->close();
    static_sys_element_loads.update();
}

void APECMEvolutionAnalysis::toroidal_projection_plasma_current()
{
    libMesh::out << "Projection from static mesh to main mesh..." << std::endl;
    // Rotate the J and B to cylindrical coordinates
    auto settings = this->context<ApecMSettingsObject>();

    APEC_M *program = this->get_parent_object();

    const std::string &static_mesh_name = settings->static_mesh();
    const std::string &equilibrium_mesh_name = settings->equilibrium_mesh();

    auto static_eq_sys = program->get_equation_system(static_mesh_name);
    auto equilibrium_eq_sys = program->get_equation_system(equilibrium_mesh_name);

    auto static_mesh = program->get_mesh(static_mesh_name);
    auto equilibrium_mesh = program->get_mesh(equilibrium_mesh_name);
    MeshBase &serial_static_mesh = *static_mesh;
    MeshBase &serial_equilibrium_mesh = *equilibrium_mesh;
    std::unique_ptr<MeshSerializer> static_serializer = std::make_unique<MeshSerializer>(serial_static_mesh, true);
    std::unique_ptr<MeshSerializer> equilibrium_serializer = std::make_unique<MeshSerializer>(serial_equilibrium_mesh, true);

    ApecMExplicitSystem &static_sys_element_loads = static_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);
    ApecMExplicitSystem &equilibrium_sys_element_loads = equilibrium_eq_sys->get_system<ApecMExplicitSystem>(ApecMSystemBase::ELEMENT_LOADS);

    static_sys_element_loads.solution->localize(*static_sys_element_loads.current_local_solution);
    equilibrium_sys_element_loads.solution->localize(*equilibrium_sys_element_loads.current_local_solution);

    auto j_b_labels_static = get_j_b_labels(static_mesh->mesh_dimension());
    auto j_b_labels_equilibrium = get_j_b_labels(equilibrium_mesh->mesh_dimension());

    std::vector<dof_id_type> dofs_js_z_equilibrium;

    std::vector<std::vector<dof_id_type>> dofs_js_tor_static(toroidal_slices_elems_.size() - 1);

    std::vector<Real> js_z;

    equilibrium_sys_element_loads.dof_numbers_vars("Js" + j_b_labels_equilibrium[0], true, toroidal_slices_elems_[0], dofs_js_z_equilibrium, true);

    equilibrium_sys_element_loads.current_local_solution->get(dofs_js_z_equilibrium, js_z);

    for (unsigned int slice_idx = 1; slice_idx < toroidal_slices_elems_.size(); ++slice_idx)
    {
        const auto &slice_elem_ids = toroidal_slices_elems_[slice_idx];
        static_sys_element_loads.dof_numbers_vars("Js" + j_b_labels_static[0], true, slice_elem_ids, dofs_js_tor_static[slice_idx - 1], true);
        static_sys_element_loads.solution->insert(js_z, dofs_js_tor_static[slice_idx - 1]);

        const std::string js_component_1 = "Js" + j_b_labels_static[1];
        if (static_sys_element_loads.has_variable(js_component_1))
        {
            std::vector<dof_id_type> dofs_js_component_1;
            static_sys_element_loads.dof_numbers_vars(js_component_1, true, slice_elem_ids, dofs_js_component_1, true);
            std::vector<Real> zeros_component_1(dofs_js_component_1.size(), 0.0);
            static_sys_element_loads.solution->insert(zeros_component_1, dofs_js_component_1);
        }

        const std::string js_component_2 = "Js" + j_b_labels_static[2];
        if (static_sys_element_loads.has_variable(js_component_2))
        {
            std::vector<dof_id_type> dofs_js_component_2;
            static_sys_element_loads.dof_numbers_vars(js_component_2, true, slice_elem_ids, dofs_js_component_2, true);
            std::vector<Real> zeros_component_2(dofs_js_component_2.size(), 0.0);
            static_sys_element_loads.solution->insert(zeros_component_2, dofs_js_component_2);
        }

        static_sys_element_loads.local_to_global_coordinate_system_solution(cylindrical_coordinate_system, "Js", slice_elem_ids);
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
        if (conv)
        {
            program->set_active_mesh(evolution_mesh);
            toroidal_projection_plasma_current();
            set_poloidal_currents();
            solve_equation_systems(substep, false);
        }
        return conv;
    }

    program->set_active_mesh(evolution_mesh);
    toroidal_projection_plasma_current();
    set_poloidal_currents();

    libMesh::out << "Solving evolution step..." << std::endl;
    solve_equation_systems(substep, false);

    if (print_full_output)
        print_output(substep, true);

    program->set_active_mesh(equilibrium_mesh);

    toroidal_average(); // With this you set the averaged psi, eddy currents, and B fields on the static mesh

    find_separatrix();

    apply_jplasma_gradshafranov();

    bool converged = check_poloidal_flux_convergence();

    program->set_active_mesh(static_mesh);

    toroidal_projection_plasma_current();

    impressed_eddy_currents_coils_currents();

    libMesh::out << "Solving equilibrium step..." << std::endl;
    solve_equation_systems(substep, false);

    append_iteration_diagnostics(substep, converged);

    if (converged)
        return true;

    if (print_full_output)
        print_output(substep, true);

    return false;
}
