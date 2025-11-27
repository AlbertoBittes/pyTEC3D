#include "transient_analysis.h"

#include <algorithm>
#include <filesystem>
#include <set>
#include <utility>

#include <libmesh/libmesh.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>
#include <libmesh/dof_map.h>

#include "apec_m.h"
#include "apec_m_eq_system_object.h"
#include "apec_m_system_objects.h"
#include "binary_reader.h"
#include "element_types_factory.h"
#include "fem_data_repository.h"
#include "globals.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "apec_m_settings.h"
#include "time_step.h"
#include "vtkio_updated.h"

namespace
{
    bool has_bin_extension(const std::filesystem::path &path)
    {
        return path.has_extension() && path.extension() == ".bin";
    }
}

TransientAnalysis::TransientAnalysis(APEC_M *program,
                                     std::shared_ptr<AnalysisBase::ContextView> context_view)
    : TransientAnalysis(AnalysisBase::Type::TRANSIENT, program, "TransientAnalysis", std::move(context_view))
{
}

TransientAnalysis::TransientAnalysis(AnalysisBase::Type type,
                                     APEC_M *program,
                                     std::string analysis_name,
                                     std::shared_ptr<AnalysisBase::ContextView> context_view)
    : StaticAnalysis(type, program, std::move(context_view), std::move(analysis_name))
{
}

std::vector<std::filesystem::path> TransientAnalysis::list_load_step_files() const
{
    std::vector<std::filesystem::path> files;
    auto session = this->context<Session>();
    if (!session)
        return {};

    const auto &folder = session->get_path("load_steps_folderpath");
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

bool TransientAnalysis::solve_transient()
{
    auto program = this->get_parent_object();
    if (!program || !this->active_equation_system_ptr())
    {
        libMesh::out << "No load step analysis available for the active mesh." << std::endl;
        return false;
    }

    const auto grid_base = resolve_default_grid_base_name();
    const auto time_base = resolve_default_time_base_name();
    if (grid_base.empty() || time_base.empty())
    {
        libMesh::out << "Unable to determine default base names for transient solve." << std::endl;
        return false;
    }

    unsigned int step_index = static_cast<unsigned int>(load_steps_completed_count + 1);

    bool solved_any = false;
    for (; true; ++step_index)
    {
        auto grid_path = grid_base;
        grid_path.replace_filename(grid_base.stem().string() + "_" + std::to_string(step_index));
        grid_path.replace_extension(".bin");

        auto time_path = time_base;
        time_path.replace_filename(time_base.stem().string() + "_" + std::to_string(step_index));
        time_path.replace_extension(".bin");

        if (!std::filesystem::exists(grid_path) || !std::filesystem::exists(time_path))
            break;

        libMesh::out << "Solving transient step " << step_index << " using grid '" << grid_path << "' and timestep '" << time_path << "'." << std::endl;

        if (!program->load_data({grid_path.string(), time_path.string()}))
            return false;

        if (!execute_solve())
            return false;

        load_steps_completed_count++;
        solved_any = true;
    }

    if (!solved_any)
    {
        libMesh::out << "No transient grid/timestep pairs found." << std::endl;
        return false;
    }

    return true;
}
