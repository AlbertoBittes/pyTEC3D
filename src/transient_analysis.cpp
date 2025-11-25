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

bool TransientAnalysis::solve_transient(const std::string &base_name)
{
    std::filesystem::path base;
    unsigned int step_index = 1;
    if (base_name.empty())
    {
        base = resolve_default_base_name();
        if (base.empty())
        {
            libMesh::out << "Unable to determine base name for transient solve." << std::endl;
            return false;
        }
        step_index = static_cast<unsigned int>(load_steps_completed_count + 1);
    }
    else
    {
        auto resolved_path = resolve_load_step_path(base_name);
        if (resolved_path.empty())
        {
            libMesh::out << "Unable to resolve load step path for '" << base_name << "'." << std::endl;
            return false;
        }
        base = resolved_path.parent_path() / resolved_path.stem();
    }

    if (base.empty())
    {
        libMesh::out << "Unable to determine base name for transient solve." << std::endl;
        return false;
    }

    bool solved_any = false;
    for (; true; ++step_index)
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

        load_steps_completed_count++;
        solved_any = true;
    }

    if (!solved_any)
    {
        libMesh::out << "No transient load steps found for base '" << base << "'." << std::endl;
        return false;
    }

    return true;
}
