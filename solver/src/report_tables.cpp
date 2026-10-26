/**
 * @file report_tables.cpp
 * @brief Text tables of the solver log (report_tables.h).
 */
#include "report_tables.h"

#include "mesh_quality.h"
#include "model_element_types.h"
#include "pytec3d_equation_systems.h"
#include "pytec3d_system.h"
#include "solver_configuration.h"
#include "time_step.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/libmesh_config.h>
#include <libmesh/mesh_base.h>
#include <libmesh/mesh_tools.h>
#include <petscversion.h>

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
    using pytec3d::report::Table;

    std::string number(double value) { return fmt::format("{:.6g}", value); }
    std::string scientific(double value) { return fmt::format("{:.3e}", value); }
    std::string yes_no(bool value) { return value ? "yes" : "no"; }

    template <typename T>
    std::string join(const T &values, std::string_view separator = ", ")
    {
        std::string out;
        for (const auto &value : values)
        {
            if (!out.empty())
                out += separator;
            out += fmt::format("{}", value);
        }
        return out;
    }

    /** @brief "p0=12, p1=10" from per-rank counts. */
    std::string by_rank(const std::vector<libMesh::dof_id_type> &counts)
    {
        std::vector<std::string> parts;
        for (std::size_t rank = 0; rank < counts.size(); ++rank)
            parts.push_back("p" + std::to_string(rank) + "=" + std::to_string(counts[rank]));
        return join(parts);
    }

    /** @brief Ids as compact ranges: "1-4, 7, 9-10". */
    template <typename T>
    std::string ranges(const std::set<T> &ids)
    {
        std::vector<std::string> parts;
        for (auto it = ids.begin(); it != ids.end();)
        {
            const T first = *it;
            T last = first;
            for (++it; it != ids.end() && *it == last + 1; ++it)
                last = *it;
            parts.push_back(first == last ? std::to_string(first) : std::to_string(first) + "-" + std::to_string(last));
        }
        return join(parts);
    }

    bool looks_numeric(const std::string &text)
    {
        if (text.empty() || !(std::isdigit(static_cast<unsigned char>(text[0])) || text[0] == '-' || text[0] == '+'))
            return false;
        if (std::none_of(text.begin(), text.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); }))
            return false;
        return std::all_of(text.begin(), text.end(), [](char c)
                           { return std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == 'e' || c == 'E' ||
                                    c == '+' || c == '-'; });
    }

    /** @brief A parameter value: "1.5", "(1, 0, 0)", "2 I", "diag(1, 2, 3)" or the full tensor. */
    std::string value_text(const pytec3d::ParameterValue &value)
    {
        switch (value.type())
        {
        case pytec3d::ParameterType::Scalar:
            return number(value.scalar());
        case pytec3d::ParameterType::Vector:
        {
            const auto &v = value.vector();
            return "(" + number(v(0)) + ", " + number(v(1)) + ", " + number(v(2)) + ")";
        }
        case pytec3d::ParameterType::Tensor:
        {
            const auto &t = value.tensor();
            const bool diagonal = t(0, 1) == 0 && t(0, 2) == 0 && t(1, 0) == 0 && t(1, 2) == 0 && t(2, 0) == 0 && t(2, 1) == 0;
            if (diagonal && t(0, 0) == t(1, 1) && t(1, 1) == t(2, 2))
                return number(t(0, 0)) + " I";
            if (diagonal)
                return "diag(" + number(t(0, 0)) + ", " + number(t(1, 1)) + ", " + number(t(2, 2)) + ")";
            std::vector<std::string> rows;
            for (unsigned int i = 0; i < 3; ++i)
                rows.push_back("[" + number(t(i, 0)) + ", " + number(t(i, 1)) + ", " + number(t(i, 2)) + "]");
            return "[" + join(rows) + "]";
        }
        }
        return "?";
    }

    std::string law_names(pytec3d::registry::LawMask mask)
    {
        std::vector<std::string> names;
        for (std::size_t l = 0; l < pytec3d::registry::kNumLawId; ++l)
            if (mask & pytec3d::registry::bit(static_cast<pytec3d::registry::LawId>(l)))
                names.emplace_back(pytec3d::registry::name(static_cast<pytec3d::registry::LawId>(l)));
        return names.empty() ? "-" : join(names);
    }

    /** @brief The defined values of a materials or constants row, "name = value; ...". */
    template <typename Id, std::size_t N>
    std::string row_values(const pytec3d::PropertyRow<Id, N> &row)
    {
        std::vector<std::string> parts;
        for (std::size_t i = 0; i < N; ++i)
            if (const pytec3d::ParameterValue *value = row.find(static_cast<Id>(i)))
                parts.push_back(std::string(pytec3d::registry::name(static_cast<Id>(i))) + " = " + value_text(*value));
        return parts.empty() ? "-" : join(parts, "; ");
    }
} // namespace

namespace pytec3d::report
{
    // ------------------------------------------------------------------------------------------------- Table

    Table::Table(std::string title) : title_(std::move(title)) {}

    Table &Table::header(std::vector<std::string> cells)
    {
        header_ = std::move(cells);
        return *this;
    }

    Table &Table::row(std::vector<std::string> cells)
    {
        lines_.push_back({Line::Kind::Row, std::move(cells)});
        return *this;
    }

    Table &Table::section(std::string title)
    {
        lines_.push_back({Line::Kind::Section, {std::move(title)}});
        return *this;
    }

    std::string Table::str() const
    {
        std::size_t columns = header_.size();
        for (const Line &line : lines_)
            if (line.kind == Line::Kind::Row)
                columns = std::max(columns, line.cells.size());
        columns = std::max<std::size_t>(columns, 1);

        std::vector<std::size_t> widths(columns, 0);
        const auto measure = [&widths](const std::vector<std::string> &cells)
        {
            for (std::size_t c = 0; c < cells.size(); ++c)
                widths[c] = std::max(widths[c], cells[c].size());
        };
        measure(header_);
        for (const Line &line : lines_)
            if (line.kind == Line::Kind::Row)
                measure(line.cells);

        // The spanning lines (title, sections) widen the last column if they do not fit.
        std::size_t inner = std::accumulate(widths.begin(), widths.end(), std::size_t(0)) + 3 * (columns - 1);
        std::size_t spanning = title_.size();
        for (const Line &line : lines_)
            if (line.kind == Line::Kind::Section)
                spanning = std::max(spanning, line.cells[0].size());
        if (spanning > inner)
        {
            widths.back() += spanning - inner;
            inner = spanning;
        }

        std::string rule = "+";
        for (const std::size_t width : widths)
            rule += std::string(width + 2, '-') + "+";
        const std::string full_rule = "+" + std::string(inner + 2, '-') + "+";
        const auto spanning_line = [inner](const std::string &text)
        { return "| " + text + std::string(inner - text.size(), ' ') + " |"; };
        const auto row_line = [&widths, columns](const std::vector<std::string> &cells)
        {
            std::string out = "|";
            for (std::size_t c = 0; c < columns; ++c)
            {
                const std::string cell = c < cells.size() ? cells[c] : "";
                const std::string pad(widths[c] - cell.size(), ' ');
                out += " " + (looks_numeric(cell) ? pad + cell : cell + pad) + " |";
            }
            return out;
        };

        std::ostringstream out;
        out << full_rule << "\n" << spanning_line(title_) << "\n" << rule << "\n";
        if (!header_.empty())
            out << row_line(header_) << "\n" << rule << "\n";
        for (std::size_t i = 0; i < lines_.size(); ++i)
        {
            const Line &line = lines_[i];
            if (line.kind == Line::Kind::Row)
                out << row_line(line.cells) << "\n";
            else
            {
                if (i > 0)
                    out << rule << "\n";
                out << spanning_line(line.cells[0]) << "\n" << rule << "\n";
            }
        }
        if (lines_.empty() || lines_.back().kind == Line::Kind::Row)
            out << rule;
        std::string text = out.str();
        if (!text.empty() && text.back() == '\n')
            text.pop_back();
        return text;
    }

    // ------------------------------------------------------------------------------------------------- banner

    std::string format_solver_banner(const BannerInfo &info)
    {
        return Table("A P E C - M   S O L V E R  -  low-frequency electromagnetic finite elements (pyTEC3D)")
            .row({"Version", std::string(info.version)})
            .row({"Build mode", std::string(info.build_mode)})
            .row({"Build stamp", std::string(info.build_stamp)})
            .row({"Git commit", std::string(info.git_commit)})
#if defined(__clang__)
            .row({"Compiler", std::string("Clang ") + __clang_version__})
#elif defined(__GNUC__)
            .row({"Compiler", std::string("GCC ") + __VERSION__})
#endif
            .row({"libMesh", LIBMESH_LIB_VERSION})
            .row({"PETSc", fmt::format("{}.{}.{}", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR, PETSC_VERSION_SUBMINOR)})
            .row({"Host", std::string(info.host)})
            .row({"MPI ranks", std::to_string(info.ranks)})
            .row({"Threads per rank", std::to_string(info.threads)})
            .str();
    }

    // ------------------------------------------------------------------------------------------------- mesh

    std::string format_mesh_summary_table(const libMesh::MeshBase &mesh)
    {
        const auto &comm = mesh.comm();
        std::vector<libMesh::dof_id_type> nodes_by_rank, elements_by_rank;
        comm.allgather(mesh.n_local_nodes(), nodes_by_rank);
        comm.allgather(mesh.n_local_elem(), elements_by_rank);

        // Cell types and sizes over the active local elements, summed / reduced over the ranks.
        std::vector<libMesh::dof_id_type> type_counts(static_cast<std::size_t>(libMesh::INVALID_ELEM) + 1, 0);
        double h_min = std::numeric_limits<double>::max(), h_max = 0.0;
        for (const auto *elem : mesh.active_local_element_ptr_range())
        {
            ++type_counts[static_cast<std::size_t>(elem->type())];
            h_min = std::min(h_min, static_cast<double>(elem->hmin()));
            h_max = std::max(h_max, static_cast<double>(elem->hmax()));
        }
        comm.sum(type_counts);
        comm.min(h_min);
        comm.max(h_max);
        std::vector<std::string> types;
        for (std::size_t t = 0; t < type_counts.size(); ++t)
            if (type_counts[t])
                types.push_back(libMesh::Utility::enum_to_string(static_cast<libMesh::ElemType>(t)) + " " +
                                std::to_string(type_counts[t]));

        const libMesh::BoundingBox box = libMesh::MeshTools::create_bounding_box(mesh);
        std::set<libMesh::subdomain_id_type> subdomains;
        mesh.subdomain_ids(subdomains);
        std::vector<unsigned int> dimensions(mesh.elem_dimensions().begin(), mesh.elem_dimensions().end());
        std::vector<std::string> orders;
        for (const auto order : mesh.elem_default_orders())
            orders.push_back(libMesh::Utility::enum_to_string(order));

        Table table("Mesh");
        table.row({"Spatial dimension", std::to_string(mesh.spatial_dimension())})
            .row({"Element dimensions", join(dimensions)})
            .row({"Element orders", join(orders)})
            .row({"Nodes", std::to_string(mesh.n_nodes())})
            .row({"Elements (active)", std::to_string(mesh.n_elem()) + " (" + std::to_string(mesh.n_active_elem()) + ")"})
            .row({"Element types", join(types)})
            .row({"Subdomains (TYPE ids)", ranges(subdomains)})
            .row({"Element size h", "min " + scientific(h_min) + ", max " + scientific(h_max)})
            .section("Bounding box");
        for (unsigned int axis = 0; axis < 3; ++axis)
            table.row({std::string(1, "xyz"[axis]),
                       "[" + number(box.min()(axis)) + ", " + number(box.max()(axis)) + "]"});
        table.section("Distribution")
            .row({"Nodes by rank", by_rank(nodes_by_rank)})
            .row({"Elements by rank", by_rank(elements_by_rank)});
        return table.str();
    }

    std::string format_mesh_quality_table(const std::vector<ShapeQuality> &quality, std::size_t flipped)
    {
        Table table("Mesh quality");
        table.header({"Shape", "Elements", "Max edge ratio", "Min angle", "Max angle", "Min scaled Jacobian",
                      "Worst cell", "Distorted", "Degenerate"});
        for (const ShapeQuality &shape : quality)
        {
            const bool angles = shape.max_angle > 0.0;
            table.row({shape.shape, std::to_string(shape.elements), fmt::format("{:.3g}", shape.max_edge_ratio),
                       angles ? fmt::format("{:.1f}", shape.min_angle) : "-",
                       angles ? fmt::format("{:.1f}", shape.max_angle) : "-",
                       fmt::format("{:.3f}", shape.min_scaled_jacobian), std::to_string(shape.worst_element),
                       std::to_string(shape.distorted), std::to_string(shape.degenerate)});
        }
        table.section(fmt::format("Angles in degrees; distorted: scaled Jacobian < {}", kDistortedScaledJacobian))
            .row({"Reoriented on import", std::to_string(flipped)});
        return table.str();
    }

    // ------------------------------------------------------------------------------------------------- model

    std::string format_model_table(const pyTEC3DEquationSystems &equation_systems)
    {
        const pytec3d::ModelData &data = equation_systems.get_model_data();
        const ModelElementTypes &element_types = equation_systems.element_types();
        const libMesh::MeshBase &mesh = equation_systems.get_mesh();

        // Cells of each TYPE id, over all the ranks.
        std::vector<libMesh::subdomain_id_type> ids;
        for (const auto &[id, entry] : element_types.entries())
            ids.push_back(id);
        std::sort(ids.begin(), ids.end());
        std::vector<libMesh::dof_id_type> cells(ids.size(), 0);
        for (const auto *elem : mesh.active_local_element_ptr_range())
        {
            const auto it = std::lower_bound(ids.begin(), ids.end(), elem->subdomain_id());
            if (it != ids.end() && *it == elem->subdomain_id())
                ++cells[static_cast<std::size_t>(it - ids.begin())];
        }
        mesh.comm().sum(cells);

        const auto &laws = data.law_selection();
        const auto selected = [&laws](registry::LawSource source, unsigned int id)
        {
            const auto &by_id = laws[static_cast<std::size_t>(source)];
            const auto it = by_id.find(id);
            return law_names(it == by_id.end() ? 0 : it->second);
        };

        // The analysis: the dimensions of the element types and the 2D kinematics.
        std::set<unsigned int> dimensions;
        bool planar = false, axisymmetric = false;
        for (const auto &[id, entry] : element_types.entries())
        {
            dimensions.insert(entry.type->dimension());
            if (entry.type->dimension() == 2)
                (entry.type->options().axisymmetric ? axisymmetric : planar) = true;
        }
        std::vector<std::string> kinds;
        if (planar)
            kinds.emplace_back("2D planar");
        if (axisymmetric)
            kinds.emplace_back("2D axisymmetric (about the y axis)");
        if (dimensions.count(3))
            kinds.emplace_back("3D");
        Table analysis("Model: analysis");
        analysis.row({"Analysis", kinds.empty() ? "-" : join(kinds, " + ")})
            .row({"Units", "SI: m, kg, s, A, V, T; the model data and the mesh must use them"})
            .row({"Time order", std::to_string(element_types.time_order()) + " (0 static, 1 first order, 2 second order)"});
        if (planar || axisymmetric)
            analysis.row({"2D out-of-plane length", "per TYPE below (x 2 pi r if axisymmetric)"});

        // The out-of-plane length of a 2D element type: its parameter, or 1.
        const auto out_of_plane = [](const ElementType &type) -> std::string
        {
            if (type.dimension() != 2)
                return "-";
            const auto parameter = type.spec().out_of_plane_length;
            std::string length = parameter < 0 ? "1" : std::string(registry::name(static_cast<registry::ParameterId>(parameter)));
            return type.options().axisymmetric ? length + " x 2 pi r" : length;
        };

        Table types("Model: element types (TYPE ids)");
        types.header({"TYPE", "Element type", "Domain", "Dim.", "Order", "Axisym.", "Out-of-plane length", "Time order",
                      "Cells", "Laws selected", "Laws it can take"});
        for (std::size_t i = 0; i < ids.size(); ++i)
        {
            const auto &entry = element_types.get(ids[i]);
            const ElementType &type = *entry.type;
            types.row({std::to_string(ids[i]), std::string(type.name()), std::string(registry::name(entry.domain)),
                       std::to_string(type.dimension()), std::to_string(type.options().order),
                       yes_no(type.options().axisymmetric), out_of_plane(type), std::to_string(type.spec().time_order),
                       std::to_string(cells[i]), selected(registry::LawSource::Cell, ids[i]), law_names(type.laws())});
        }

        Table materials("Model: materials (MATERIAL ids)");
        materials.header({"MATERIAL", "Properties", "Nonlinear laws"});
        for (const auto &[id, row] : data.material_rows())
            materials.row({std::to_string(id), row_values(row), selected(registry::LawSource::Material, id)});

        Table constants("Model: real constants (CONSTANT ids)");
        constants.header({"CONSTANT", "Real constants", "Nonlinear laws"});
        for (const auto &[id, row] : data.constant_rows())
            constants.row({std::to_string(id), row_values(row), selected(registry::LawSource::Constant, id)});

        // DOF couplings, grouped by the fields they couple.
        std::map<int, std::vector<std::string>> fields_by_id;
        for (const auto &coupling : data.dof_couplings())
            fields_by_id[coupling.cp_id].push_back(coupling.variable);
        std::map<std::string, std::set<int>> ids_by_fields;
        for (const auto &[cp_id, fields] : fields_by_id)
            ids_by_fields[join(fields)].insert(cp_id);
        Table couplings("Model: DOF couplings (CP ids)");
        couplings.header({"Fields", "CP ids", "Count"});
        for (const auto &[fields, cp_ids] : ids_by_fields)
            couplings.row({fields, ranges(cp_ids), std::to_string(cp_ids.size())});
        if (ids_by_fields.empty())
            couplings.row({"-", "-", "0"});

        std::string out = analysis.str() + "\n" + types.str() + "\n" + materials.str();
        if (!data.constant_rows().empty())
            out += "\n" + constants.str();
        return out + "\n" + couplings.str();
    }

    // ------------------------------------------------------------------------------------------------- systems

    std::string format_equation_systems_summary_table(const pyTEC3DEquationSystems &equation_systems)
    {
        Table table("Systems");
        table.header({"System", "Variable", "FE family", "Order", "Subdomains", "Variable DOFs", "System DOFs",
                      "Local DOFs by rank"});
        libMesh::dof_id_type all_dofs = 0;
        for (int i = 0; i < pyTEC3DSystemBase::NUM_SYSTEMS; ++i)
        {
            const auto type = static_cast<pyTEC3DSystemBase::SystemType>(i);
            const auto &name = pyTEC3DSystemBase::SYSTEM_INFO_MAP.at(type);
            const auto &system = equation_systems.get_system<libMesh::System>(name);
            std::vector<libMesh::dof_id_type> local_by_rank;
            system.comm().allgather(system.n_local_dofs(), local_by_rank);
            all_dofs += system.n_dofs();
            if (system.n_vars() == 0)
            {
                table.row({name, "-", "-", "-", "-", "-", std::to_string(system.n_dofs()), by_rank(local_by_rank)});
                continue;
            }
            for (unsigned int var = 0; var < system.n_vars(); ++var)
            {
                const libMesh::FEType fe = system.variable_type(var);
                std::vector<libMesh::dof_id_type> local_indices;
                system.get_dof_map().local_variable_indices(local_indices, var);
                libMesh::dof_id_type variable_dofs = local_indices.size();
                system.comm().sum(variable_dofs);
                const auto &subdomains = system.variable(var).active_subdomains();
                table.row({var == 0 ? name : "", system.variable_name(var),
                           libMesh::Utility::enum_to_string(fe.family),
                           libMesh::Utility::enum_to_string(static_cast<libMesh::Order>(fe.order)),
                           subdomains.empty() ? "all" : ranges(subdomains), std::to_string(variable_dofs),
                           var == 0 ? std::to_string(system.n_dofs()) : "", var == 0 ? by_rank(local_by_rank) : ""});
            }
        }
        table.section("All systems: " + std::to_string(all_dofs) + " DOFs");
        return table.str();
    }

    // ------------------------------------------------------------------------------------------------- solver

    std::string format_solver_configuration_table(const SolverSettings &settings, bool nonlinear, std::string_view ksp_name,
                                                  std::string_view pc_name, std::string_view factor_package)
    {
        const LinearSolverSettings &linear = settings.linear;
        Table table("Solver configuration");
        table.section("Linear solver (PETSc KSP" + std::string(nonlinear ? ", in each Newton iteration)" : ")"))
            .row({"Type", std::string(ksp_name) == "preonly" ? "direct" : "iterative"})
            .row({"KSP", std::string(ksp_name)})
            .row({"Preconditioner", std::string(pc_name)})
            .row({"Factorization package", std::string(factor_package)})
            .row({"Relative tolerance", scientific(linear.relative_tolerance)})
            .row({"Absolute tolerance", scientific(linear.absolute_tolerance)})
            .row({"Divergence tolerance", linear.divergence_tolerance > 0 ? scientific(linear.divergence_tolerance)
                                                                          : "PETSc default"})
            .row({"Maximum iterations", std::to_string(linear.max_iterations)})
            .row({"Initial guess", linear.initial_guess_nonzero ? "previous solution" : "zero"});
        if (nonlinear)
        {
            const NewtonSettings &newton = settings.nonlinear;
            table.section("Nonlinear solver (libMesh Newton)")
                .row({"Maximum iterations", std::to_string(newton.max_iterations)})
                .row({"Residual tolerance", "relative " + scientific(newton.relative_residual_tolerance) + ", absolute " +
                                                scientific(newton.absolute_residual_tolerance)})
                .row({"Step tolerance", "relative " + scientific(newton.relative_step_tolerance) + ", absolute " +
                                            scientific(newton.absolute_step_tolerance)})
                .row({"Line search", newton.brent_line_search ? "Brent" : "backtracking"})
                .row({"Minimum step length", scientific(newton.minimum_step_length)})
                .row({"Linear tolerance multiplier", scientific(newton.linear_tolerance_multiplier)})
                .row({"Require residual reduction", yes_no(newton.require_residual_reduction)})
                .row({"Require finite residual", yes_no(newton.require_finite_residual)})
                .row({"Continue after max iterations", yes_no(newton.continue_after_max_iterations)})
                .row({"Continue after line search failure", yes_no(newton.continue_after_backtrack_failure)});
        }
        table.section("PETSc options");
        for (const auto &[name, value] : settings.petsc_options)
            table.row({name, value.empty() ? "(flag)" : value});
        if (settings.petsc_options.empty())
            table.row({"-", "-"});
        return table.str();
    }

    // ------------------------------------------------------------------------------------------------- time step

    std::string format_time_step_table(const TimeStep &time_step, unsigned int time_order, std::string_view title)
    {
        const TimeScheme scheme = TimeScheme::make(time_step, time_order);
        const std::string kind = !scheme.transient() ? "static"
                                 : time_order == 1   ? "theta method (theta = " + number(time_step.theta) + ")"
                                                     : "Newmark (beta = " + number(time_step.beta) + ", gamma = " +
                                                         number(time_step.gamma) + ")";
        Table table{std::string(title)};
        table.row({"Index", std::to_string(time_step.time_step_idx)})
            .row({"Time", scientific(time_step.time)})
            .row({"dt", scientific(time_step.dt)})
            .row({"Scheme", kind})
            .row({"Formulation time order", std::to_string(time_order)});
        if (scheme.transient())
            table.row({"d/dt coefficient c_v", scientific(scheme.c_v)});
        if (scheme.c_a != 0.0)
            table.row({"d2/dt2 coefficient c_a", scientific(scheme.c_a)});
        return table.str();
    }

    std::string format_step_report_table(const StepReport &report)
    {
        const double total = report.seconds_element_matrices + report.seconds_solve + report.seconds_nodal_fields +
                             report.seconds_element_fields + report.seconds_output;
        const auto seconds = [](double s) { return fmt::format("{:.3f}", s); };
        Table table("Step " + std::to_string(report.index) + " solved: t = " + scientific(report.time));
        table.row({"Kind", std::string(report.transient ? "transient" : "static") +
                               (report.advanced ? "" : ", same time solved again")})
            .row({"dt", scientific(report.dt)})
            .row({"Elements", std::to_string(report.elements)})
            .row({report.nonlinear ? "DOFs" : "DOFs (unconstrained)",
                  std::to_string(report.total_dofs) +
                      (report.nonlinear ? "" : " (" + std::to_string(report.unconstrained_dofs) + ")")})
            .row({"Element matrices", report.element_matrices_computed ? "computed" : "reused"})
            .row({"System matrix", report.system_matrix_updated ? "assembled and factorized" : "reused"});
        if (report.nonlinear)
            table.row({"Newton iterations", std::to_string(report.newton_iterations)})
                .row({"Newton result", report.newton_result})
                .row({"Linear iterations (all Newton steps)", std::to_string(report.linear_iterations)});
        else
            table.row({"Linear iterations", std::to_string(report.linear_iterations)})
                .row({"Linear residual", scientific(report.linear_residual)})
                .row({"True residual ||b - Ax|| / ||b||", report.true_residual < 0 ? "-" : scientific(report.true_residual)})
                .row({"Linear convergence", report.linear_reason});
        const FactorizationReport &f = report.factorization;
        if (f.available)
        {
            table.section(std::string("Factorization (MUMPS)") + (report.system_matrix_updated ? "" : ", reused"))
                .row({"Status INFOG(1), INFOG(2)", std::to_string(f.status) + ", " + std::to_string(f.status_detail) +
                                                       (f.status == 0 ? " (ok)" : f.status < 0 ? " (error)" : " (warning)")});
            if (f.status >= 0)
            {
                table.row({"Memory (MB)", fmt::format("{} total, {} on the largest rank", f.memory_total_mb, f.memory_max_mb)})
                    .row({"Entries of the matrix", fmt::format("{:.4g}", f.matrix_entries)})
                    .row({"Entries of the factors", fmt::format("{:.4g}", f.factor_entries)})
                    .row({"Fill-in", f.matrix_entries > 0 ? fmt::format("{:.2f}", f.factor_entries / f.matrix_entries) : "-"})
                    .row({"Floating-point operations", fmt::format("{:.4g}", f.flops)})
                    .row({"Delayed pivots", std::to_string(f.delayed_pivots)})
                    .row({"Null pivots", f.null_pivots < 0 ? "not checked (-mat_mumps_icntl_24 1)" : std::to_string(f.null_pivots)});
            }
        }
        if (!report.law_ranges.empty())
        {
            table.section("Nonlinear laws: operating range (largest element mean of the input)");
            for (const LawRange &range : report.law_ranges)
                table.row({range.law + " of " + range.source + " " + std::to_string(range.id),
                           fmt::format("max |{}| = {:.4g}", range.input, range.max_input) +
                               (range.curve_max > 0 ? fmt::format(", curve up to {:.4g}{}", range.curve_max,
                                                                  range.extrapolated ? ": EXTRAPOLATED" : "")
                                                    : "")});
        }
        table.section("Wall time (s)")
            .row({"Element matrices", seconds(report.seconds_element_matrices)})
            .row({"Assembly and solve", seconds(report.seconds_solve)})
            .row({"Nodal fields", seconds(report.seconds_nodal_fields)})
            .row({"Element fields", seconds(report.seconds_element_fields)})
            .row({"Output to the grid", seconds(report.seconds_output)})
            .row({"Total", seconds(total)});
        return table.str();
    }

    std::string format_newton_history_table(const StepReport &report)
    {
        Table table("Newton iterations: t = " + scientific(report.time) + " (" + report.newton_result + ")");
        table.header({"Iterate", "Residual ||R||", "Reduction", "Step ||du||", "||du|| / ||u||", "Line search step"});
        const double initial = report.newton_history.empty() ? 0.0 : report.newton_history.front().residual;
        for (std::size_t i = 0; i < report.newton_history.size(); ++i)
        {
            const NewtonIteration &it = report.newton_history[i];
            table.row({std::to_string(i), scientific(it.residual), initial > 0 ? scientific(it.residual / initial) : "-",
                       it.step < 0 ? "-" : scientific(it.step), it.relative_step < 0 ? "-" : scientific(it.relative_step),
                       it.step_length < 0 ? "-" : fmt::format("{:.3g}", it.step_length)});
        }
        return table.str();
    }

    // ------------------------------------------------------------------------------------------------- constraints

    std::string format_constraints_summary_table(std::size_t time_integration_count,
                                                 std::size_t dof_constraint_count,
                                                 std::size_t coupling_count,
                                                 std::size_t total_count,
                                                 std::size_t applied_count,
                                                 const std::vector<ConstraintVariableStat> &variable_stats,
                                                 std::size_t unknown_total_count,
                                                 std::size_t unknown_applied_count)
    {
        Table table("Constraints");
        table.row({"Static-step constraints (no stiffness)", std::to_string(time_integration_count)})
            .row({"Dirichlet constraints", std::to_string(dof_constraint_count)})
            .row({"Coupling constraints", std::to_string(coupling_count)})
            .row({"Total", std::to_string(total_count)})
            .row({"Applied", std::to_string(applied_count)})
            .section("By variable: total / applied");
        for (const auto &stat : variable_stats)
            table.row({stat.variable_name, std::to_string(stat.total_count) + " / " + std::to_string(stat.applied_count)});
        if (unknown_total_count > 0 || unknown_applied_count > 0)
            table.row({"(unknown)", std::to_string(unknown_total_count) + " / " + std::to_string(unknown_applied_count)});
        return table.str();
    }
} // namespace pytec3d::report
