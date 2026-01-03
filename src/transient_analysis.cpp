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
#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "binary_reader.h"
#include "element_types_factory.h"
#include "fem_data_repository.h"

#include "materials_object.h"
#include "real_constants_object.h"
#include "apec_m_settings.h"
#include "time_step.h"
#include "vtkio_updated.h"

bool TransientAnalysis::solve_transient()
{
    auto program = this->get_parent_object();

    auto session = this->context<Session>();

    const auto grid_base = session->resolve_default_step_base("load_steps_folderpath", "grid", "_grid");
    const auto time_base = session->resolve_default_step_base("time_steps_folderpath", "timestep", "_timestep");
    if (grid_base.empty() || time_base.empty())
    {
        libMesh::out << "Unable to determine default base names for transient solve." << std::endl;
        return false;
    }

    auto time_step = this->context<TimeStep>();
    unsigned int step_index = static_cast<unsigned int>((time_step ? time_step->count_steps_including_previous() : 0) + 1);

    bool solved_any = false;
    for (; true; ++step_index)
    {
        auto grid_path = grid_base;
        grid_path.replace_filename(grid_base.stem().string() + "_" + std::to_string(step_index));
        grid_path.replace_extension(".bin");

        auto time_path = time_base;
        time_path.replace_filename(time_base.stem().string() + "_" + std::to_string(step_index));
        time_path.replace_extension(".bin");

        const bool has_grid = std::filesystem::exists(grid_path);
        const bool has_time = std::filesystem::exists(time_path);
        if (!has_time)
            break;

        if (!has_grid)
        {
            libMesh::out << "Transient step " << step_index
                         << ": load/bcs binary not found, using previous loads and bcs." << std::endl;
        }

        libMesh::out << "Solving transient step " << step_index << " using "
                     << (has_grid ? "grid '" + grid_path.string() + "' and " : "")
                     << "timestep '" << time_path << "'." << std::endl;

        std::vector<std::string> paths;
        paths.reserve(2);
        if (has_grid)
            paths.push_back(grid_path.string());
        paths.push_back(time_path.string());
        if (!program->load_data(paths))
            return false;

        if (!execute_solve())
            return false;

        solved_any = true;
    }

    if (!solved_any)
    {
        libMesh::out << "No transient grid/timestep pairs found." << std::endl;
        return false;
    }

    return true;
}
