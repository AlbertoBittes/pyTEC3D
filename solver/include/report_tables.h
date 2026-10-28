/**
 * @file report_tables.h
 * @brief Text tables of the solver log: banner, mesh, model, systems, solver configuration, time step, constraints and
 * the report of each solved step.
 *
 * Printed on the root rank with pytec3d::logging::table_root(). The functions that gather data over the ranks (mesh,
 * model, systems) are collective: every rank calls them.
 */
#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

class pyTEC3DEquationSystems;
struct SolverSettings;
struct StepReport;
struct TimeStep;

namespace libMesh
{
    class MeshBase;
}

namespace pytec3d
{
    struct ShapeQuality;
}

namespace pytec3d::report
{
    /**
     * @brief A text table: a title, optional column headers, rows and section titles; the columns fit their content.
     *
     *     +--------------------------------+
     *     | Title                          |
     *     +-------------+------------------+
     *     | key         | value            |
     *     +-------------+------------------+
     *     | Section                        |
     *     +-------------+------------------+
     *     | key         | value            |
     *     +-------------+------------------+
     */
    class Table
    {
    public:
        explicit Table(std::string title);

        /** @brief Column headers, printed under the title. */
        Table &header(std::vector<std::string> cells);
        /** @brief A row; numbers are right-aligned. */
        Table &row(std::vector<std::string> cells);
        /** @brief A title spanning the whole width, between rules. */
        Table &section(std::string title);

        std::string str() const;

    private:
        struct Line
        {
            enum class Kind { Row, Section } kind;
            std::vector<std::string> cells;
        };
        std::string title_;
        std::vector<std::string> header_;
        std::vector<Line> lines_;
    };

    /** @brief What the banner shows besides the versions of the libraries. */
    struct BannerInfo
    {
        std::string_view version;
        std::string_view build_mode;
        std::string_view build_stamp;
        std::string_view git_commit;
        std::string_view host;
        int ranks = 1;
        int threads = 1;
    };

    struct ConstraintVariableStat
    {
        std::string variable_name;
        std::size_t total_count = 0;
        std::size_t applied_count = 0;
    };

    std::string format_solver_banner(const BannerInfo &info);
    /** @brief Collective. */
    std::string format_mesh_summary_table(const libMesh::MeshBase &mesh);
    /** @brief @p flipped: the elements reoriented when the mesh was imported. */
    std::string format_mesh_quality_table(const std::vector<ShapeQuality> &quality, std::size_t flipped);
    /** @brief The element types by TYPE id (with their cell counts), materials, constants, laws, couplings. Collective. */
    std::string format_model_table(const pyTEC3DEquationSystems &equation_systems);
    /** @brief Collective. */
    std::string format_equation_systems_summary_table(const pyTEC3DEquationSystems &equation_systems);
    std::string format_solver_configuration_table(const SolverSettings &settings, bool nonlinear, std::string_view ksp_name,
                                                  std::string_view pc_name, std::string_view factor_package);
    std::string format_time_step_table(const TimeStep &time_step, unsigned int time_order,
                                       std::string_view title = "Time Step");
    std::string format_step_report_table(const StepReport &report);
    /** @brief The residual and step of each Newton iteration of the step. */
    std::string format_newton_history_table(const StepReport &report);
    std::string format_constraints_summary_table(std::size_t time_integration_count,
                                                 std::size_t dof_constraint_count,
                                                 std::size_t coupling_count,
                                                 std::size_t total_count,
                                                 std::size_t applied_count,
                                                 const std::vector<ConstraintVariableStat> &variable_stats,
                                                 std::size_t unknown_total_count,
                                                 std::size_t unknown_applied_count);
} // namespace pytec3d::report
