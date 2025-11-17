#include "transient_analysis.h"

#include <algorithm>
#include <filesystem>
#include <set>
#include <utility>

#include <libmesh/libmesh.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/dof_map.h>

#include "apec_m_eq_system_object.h"
#include "apec_m_system_objects.h"
#include "binary_reader.h"
#include "element_types_factory.h"
#include "fem_data_repository.h"
#include "globals.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "time_step.h"
#include "vtkio_updated.h"

namespace
{
    bool has_bin_extension(const std::filesystem::path &path)
    {
        return path.has_extension() && path.extension() == ".bin";
    }
}

TransientAnalysis::TransientAnalysis(ApecMEquationSystems *eq_systems)
    : TransientAnalysis(AnalysisBase::Type::Transient, eq_systems, "TransientAnalysis")
{
}

TransientAnalysis::TransientAnalysis(AnalysisBase::Type type,
                                     ApecMEquationSystems *eq_systems,
                                     std::string analysis_name)
    : AnalysisBase(type, eq_systems, std::move(analysis_name))
{
}

std::filesystem::path TransientAnalysis::load_steps_folder() const
{
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto &folder = session->get_path("load_steps_folderpath");
    return folder;
}

std::vector<std::filesystem::path> TransientAnalysis::list_load_step_files() const
{
    std::vector<std::filesystem::path> files;
    const auto folder = load_steps_folder();
    if (folder.empty() || !std::filesystem::exists(folder))
        return files;

    for (const auto &entry : std::filesystem::directory_iterator(folder))
    {
        if (!entry.is_regular_file())
            continue;
        const auto &path = entry.path();
        if (!has_bin_extension(path))
            continue;
        files.push_back(path);
    }

    std::sort(files.begin(), files.end());
    return files;
}

std::size_t TransientAnalysis::load_step_file_count() const
{
    return list_load_step_files().size();
}

std::shared_ptr<TimeStep> TransientAnalysis::read_time_step(const std::filesystem::path &filepath)
{
    if (filepath.empty())
    {
        libMesh::out << "Load step path is empty." << std::endl;
        return nullptr;
    }

    if (!std::filesystem::exists(filepath))
    {
        libMesh::out << "Load step file '" << filepath << "' not found." << std::endl;
        return nullptr;
    }

    auto *equation_systems = this->get_parent_object();
    if (!equation_systems)
    {
        libMesh::out << "Equation systems not available for load step analysis." << std::endl;
        return nullptr;
    }

    try
    {
        auto load_step_repo = std::make_shared<FEMDataRepository>();
        BinaryReader reader;
        reader.read_binary_fem_data(filepath, load_step_repo, std::shared_ptr<ApecMEquationSystems>(equation_systems, [](ApecMEquationSystems *) {}));
        auto time_step = load_step_repo->get_object<TimeStep>();
        if (!time_step)
        {
            libMesh::out << "Time step data missing in load step binary '" << filepath << "'." << std::endl;
            return nullptr;
        }

        auto previous_time_step = this->context<TimeStep>();
        time_step->populate_time_step(previous_time_step);
        update_time_step_context(time_step);
        return time_step;
    }
    catch (const std::exception &ex)
    {
        libMesh::out << "Failed to read load step '" << filepath << "': " << ex.what() << std::endl;
        return nullptr;
    }
    catch (...)
    {
        libMesh::out << "Failed to read load step '" << filepath << "' due to an unknown error." << std::endl;
        return nullptr;
    }

    return nullptr;
}

void TransientAnalysis::configure_services() const
{
    auto *equation_systems = this->get_parent_object();

    auto analysis_ctx = this->get_context_view();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto bc_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep>() : nullptr;
    auto global_constraints_object = system.get_global_constraints_object();
    global_constraints_object->update_context_views(bc_view);

    auto assembler_view = analysis_ctx ? analysis_ctx->template project_view<TimeStep, Session>() : nullptr;

    auto assembler = system.get_assembler();
    assembler->update_context_views(assembler_view);
}

bool TransientAnalysis::execute_solve()
{
    auto *equation_systems = this->get_parent_object();
    if (!equation_systems)
        return false;

    auto time_step = this->context<TimeStep>();
    if (!time_step)
    {
        libMesh::out << "Load step does not contain time stepping information." << std::endl;
        return false;
    }

    configure_services();

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");

    for (auto &substep : *time_step)
    {
        prepare_time_substep_solve(substep);

        solve_equation_systems(substep);

        print_output(substep);

        *system.old_local_solution = *system.current_local_solution;
    }

    return true;
}

void TransientAnalysis::prepare_time_substep_solve(const Substep &substep)
{
    auto *equation_systems = this->get_parent_object();
    if (!equation_systems)
        return;

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();
    if (!assembler)
    {
        libMesh::out << "Assembler not available for linear system." << std::endl;
        return;
    }

    equation_systems->apply_ramped_effects(substep);

    system.get_dof_map().swap_dof_constraints();

    system.reinit_constraints();

    assembler->write_element_matrices_and_rhs();
}

void TransientAnalysis::solve_equation_systems(const Substep &substep)
{
    auto *equation_systems = this->get_parent_object();
    if (!equation_systems)
        return;

    auto &system = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem");
    auto assembler = system.get_assembler();
    if (!assembler)
    {
        libMesh::out << "Assembler not available for linear system." << std::endl;
        return;
    }

    system.solve();
    system.get_dof_map().enforce_constraints_exactly(system);
    system.compute_solution_time_derivative(substep);
    assembler->write_solution_vectors(true);
}

void TransientAnalysis::print_output(const Substep &substep)
{
    auto *equation_systems = this->get_parent_object();

    auto time_step = this->context<TimeStep>();

    auto session = this->context<Session>();

    auto assembler = equation_systems->get_system<ApecMLinearSystem>("ApecMLinearSystem").get_assembler();

    const bool should_write = substep.get_write_results() ||
                              substep.get_index() == substep.get_parent_object()->get_n_substeps() - 1;

    if (!should_write || !session || !equation_systems || !assembler)
        return;

    std::set<std::string> system_names = {"ApecMLinearSystem", "NodalSolutionSystem"};

    session->set_output_filepath(cumulative_substep_index);
    VTKIO_Updated vtk_output(equation_systems->get_mesh(),
                             &equation_systems->get_system<ApecMExplicitSystem>("ElementSolutionSystem"));
    vtk_output.write_equation_systems(session->get_path("final_output_path"),
                                      *equation_systems,
                                      &system_names);
    ++cumulative_substep_index;
}

std::filesystem::path TransientAnalysis::resolve_load_step_path(const std::string &filename) const
{
    auto folder = load_steps_folder();
    if (folder.empty())
        return {};

    if (filename.empty())
        return {};

    auto candidate = std::filesystem::path(filename);
    if (!candidate.is_absolute())
        candidate = folder / candidate;

    if (!candidate.has_extension())
        candidate.replace_extension(".bin");

    return candidate;
}

std::filesystem::path TransientAnalysis::resolve_default_base_name() const
{
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto folder = load_steps_folder();
    if (folder.empty())
        return {};

    auto token = session->get_mesh_filename_without_extension();
    if (token.empty())
        token = "load_step";

    std::filesystem::path base = folder / (token + "_loadstep");
    return base;
}

bool TransientAnalysis::solve(const std::string &filename)
{
    const auto path = resolve_load_step_path(filename);
    if (path.empty())
    {
        libMesh::out << "Provide a load step filename to solve." << std::endl;
        return false;
    }

    if (!read_time_step(path))
        return false;

    return execute_solve();
}

bool TransientAnalysis::solve_transient(const std::string &base_name)
{
    std::filesystem::path base;
    if (base_name.empty())
        base = resolve_default_base_name();
    else
    {
        base = resolve_load_step_path(base_name);
        base = base.parent_path() / base.stem();
    }

    if (base.empty())
    {
        libMesh::out << "Unable to determine base name for transient solve." << std::endl;
        return false;
    }

    unsigned int step_index = 1;
    bool solved_any = false;
    while (true)
    {
        auto next_path = base;
        next_path.replace_filename(base.stem().string() + "_" + std::to_string(step_index));
        next_path.replace_extension(".bin");

        if (!std::filesystem::exists(next_path))
            break;

        libMesh::out << "Solving transient load step '" << next_path << "'." << std::endl;
        if (!read_time_step(next_path))
            return false;

        if (!execute_solve())
            return false;

        solved_any = true;
        ++step_index;
    }

    if (!solved_any)
    {
        libMesh::out << "No transient load steps found for base '" << base << "'." << std::endl;
        return false;
    }

    return true;
}
