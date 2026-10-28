/**
 * @file array_io.cpp
 * @brief Mesh import from the caller's mesh view and output of the solution systems into its arrays.
 */
#include "array_io.h"

#include "element_type.h"
#include "logging.h"
#include "model_element_types.h"
#include "pytec3d_system_base.h"
#include "registries.h"
#include "variables_map.h"

#include <libmesh/elem.h>
#include <libmesh/enum_io_package.h>
#include <libmesh/equation_systems.h>
#include <libmesh/fe_interface.h>
#include <libmesh/mesh_base.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/system.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <numeric>
#include <set>
#include <tuple>

namespace pytec3d
{
    namespace
    {
        using libMesh::dof_id_type;
        using libMesh::ElemType;
        using libMesh::FEInterface;
        using libMesh::Number;
        using libMesh::System;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort(message);
        }

        /** @brief The libMesh element of a VTK cell type (Lagrange mapping). */
        ElemType elem_type_of(std::int64_t vtk_type, std::size_t cell)
        {
            switch (vtk_type)
            {
            case cell_type::LINE: return libMesh::EDGE2;
            case cell_type::QUADRATIC_EDGE: return libMesh::EDGE3;
            case cell_type::TRIANGLE: return libMesh::TRI3;
            case cell_type::QUADRATIC_TRIANGLE: return libMesh::TRI6;
            case cell_type::QUAD: return libMesh::QUAD4;
            case cell_type::QUADRATIC_QUAD: return libMesh::QUAD8;
            case cell_type::BIQUADRATIC_QUAD: return libMesh::QUAD9;
            case cell_type::TETRA: return libMesh::TET4;
            case cell_type::QUADRATIC_TETRA: return libMesh::TET10;
            case cell_type::HEXAHEDRON: return libMesh::HEX8;
            case cell_type::QUADRATIC_HEXAHEDRON: return libMesh::HEX20;
            case cell_type::TRIQUADRATIC_HEXAHEDRON: return libMesh::HEX27;
            case cell_type::WEDGE: return libMesh::PRISM6;
            case cell_type::QUADRATIC_WEDGE: return libMesh::PRISM15;
            case cell_type::BIQUADRATIC_QUADRATIC_WEDGE: return libMesh::PRISM18;
            case cell_type::PYRAMID: return libMesh::PYRAMID5;
            case cell_type::QUADRATIC_PYRAMID: return libMesh::PYRAMID13;
            default:
                fail("Mesh import: cell " + std::to_string(cell) + " has the VTK cell type " + std::to_string(vtk_type) +
                     ", which the solver does not support.");
            }
        }

        struct WriteEntry
        {
            std::string array;           ///< Target array.
            unsigned int component = 0;  ///< Component of the target array.
            std::size_t source_offset = 0;
        };

        struct WritePlan
        {
            std::vector<WriteEntry> entries;
            std::size_t tuple_stride = 0;
            bool nodal = true; ///< Point arrays, else cell arrays.
        };

        std::string strip_trailing_lowercase(std::string name)
        {
            while (!name.empty() && std::islower(static_cast<unsigned char>(name.back())))
                name.pop_back();
            return name;
        }

        /**
         * @brief Where each array component is in the solution vector built by libMesh: per node, the variables in
         * system order (nodal layout); per variable block, the variables sorted by name (element layout, as
         * EquationSystems::build_parallel_elemental_solution_vector() orders them).
         */
        WritePlan write_plan(const System &system, bool nodal)
        {
            WritePlan plan;
            plan.nodal = nodal;

            const libMesh::MeshBase &mesh = system.get_mesh();
            std::vector<unsigned int> variables(system.n_vars());
            std::iota(variables.begin(), variables.end(), 0u);
            if (!nodal)
                std::sort(variables.begin(), variables.end(), [&system](unsigned int a, unsigned int b)
                          { return system.variable_name(a) < system.variable_name(b); });

            std::size_t running_offset = 0;
            for (const unsigned int var : variables)
            {
                const std::string &var_name = system.variable_name(var);
                // A vector variable (nodal, or MONOMIAL_VEC element data) gives one value per component.
                unsigned int n_components = 1u;
                const auto n_vec_dim = FEInterface::n_vec_dim(mesh, system.variable_type(var));
                if (n_vec_dim > 1)
                    n_components = libMesh::cast_int<unsigned int>(n_vec_dim);

                const auto keys = VariablesMap::get_keys_by_predicate(
                    [&var_name](const VariableInfo &info)
                    { return info.global_label == var_name; });

                if (keys.size() == 1 && n_components == 1)
                {
                    const auto &var_info = VariablesMap::get_info_by_key(keys.front());
                    plan.entries.push_back({std::string(var_info.api_label), var_info.component_index, running_offset});
                }
                else if (n_components > 1)
                {
                    // The array of a vector variable is named by its api label (the field name).
                    const std::string base_name = keys.size() == 1 ? std::string(VariablesMap::get_info_by_key(keys.front()).api_label)
                                                                   : strip_trailing_lowercase(var_name);
                    if (base_name.empty())
                        fail("Solution output: cannot derive the array name of the vector variable '" + var_name + "'.");
                    for (unsigned int comp = 0; comp < n_components; ++comp)
                        plan.entries.push_back({base_name, comp, running_offset + comp});
                }

                running_offset += n_components;
            }

            plan.tuple_stride = running_offset;
            return plan;
        }

        WritePlan write_plan(const libMesh::EquationSystems &systems, const std::string &system_name)
        {
            if (!systems.has_system(system_name))
                fail("Solution output: system '" + system_name + "' not found.");
            const bool nodal = pyTEC3DSystemBase::system_type_from_name(system_name) != pyTEC3DSystemBase::ELEMENT_SOLUTION;
            return write_plan(systems.get_system(system_name), nodal);
        }

        /** @brief The field that receives @p entry, or nullptr when the caller does not provide it. */
        const RealArray *target(const ModelArrays &arrays, const WriteEntry &entry, bool nodal)
        {
            const RealArray *view = nodal ? arrays.point_field(entry.array) : arrays.cell_field(entry.array);
            return view && entry.component < view->n_components() ? view : nullptr;
        }

        std::size_t components_of(registry::FieldId field)
        {
            return registry::is_vector(field) ? 3 : 1;
        }
    } // namespace

    std::vector<FieldArrayRequirement> required_field_arrays(const ModelElementTypes &element_types)
    {
        std::vector<FieldArrayRequirement> out;
        std::set<std::pair<bool, std::string>> seen;
        const auto add = [&](registry::FieldId field, bool point, const std::string &prefix = "")
        {
            const std::string name = prefix + std::string(registry::name(field));
            if (seen.insert({point, name}).second)
                out.push_back({name, point, components_of(field)});
        };

        std::set<const ElementType *> types;
        for (const auto &[subdomain, entry] : element_types.entries())
            types.insert(entry.type);
        for (const ElementType *type : types)
        {
            for (std::size_t index = 0; index < type->spec().fields.size(); ++index)
            {
                const formulation::FieldSpec &field = type->field(index);
                switch (field.kind)
                {
                case FieldKind::Dof:
                    add(field.id, true);
                    add(field.id, true, "_");
                    if (field.nodal_load != registry::FieldId::COUNT)
                        add(field.nodal_load, true);
                    break;
                case FieldKind::Nodal:
                    add(field.id, true);
                    break;
                case FieldKind::Element:
                    add(field.id, false);
                    break;
                }
            }
            for (std::size_t p = 0; p < registry::kNumParameterId; ++p)
            {
                const auto id = static_cast<registry::ParameterId>(p);
                if (!registry::defined(type->id(), id))
                    continue;
                const registry::ParameterInfo &info = registry::info(id);
                if (info.source == registry::ParameterSource::Cell)
                    add(info.array, false);
                else if (info.source == registry::ParameterSource::Point)
                    add(info.array, true);
            }
        }
        return out;
    }

    void check_field_arrays(const ModelArrays &arrays, const std::vector<FieldArrayRequirement> &required,
                            std::size_t n_points, std::size_t n_cells, bool warn_unused)
    {
        std::string problems;
        std::set<std::pair<bool, std::string>> used;
        for (const FieldArrayRequirement &field : required)
        {
            used.insert({field.point, field.name});
            const char *where = field.point ? "point" : "cell";
            const RealArray *view = field.point ? arrays.point_field(field.name) : arrays.cell_field(field.name);
            const std::size_t n_tuples = field.point ? n_points : n_cells;
            if (!view)
                problems += "\n  - the " + std::string(where) + " field '" + field.name + "' is missing";
            else if (view->n_tuples() != n_tuples || view->n_components() != field.n_components)
                problems += "\n  - the " + std::string(where) + " field '" + field.name + "' has shape (" +
                            std::to_string(view->n_tuples()) + ", " + std::to_string(view->n_components()) +
                            "), expected (" + std::to_string(n_tuples) + ", " + std::to_string(field.n_components) + ")";
        }
        if (!problems.empty())
            logging::error_and_abort("Model input: the field arrays of the element types of the model:" + problems);

        if (!warn_unused)
            return;
        std::string unused;
        for (const auto &[name, view] : arrays.point_fields)
            if (!used.count({true, name}))
                unused += (unused.empty() ? "" : ", ") + name + " (point)";
        for (const auto &[name, view] : arrays.cell_fields)
            if (!used.count({false, name}))
                unused += (unused.empty() ? "" : ", ") + name + " (cell)";
        if (!unused.empty())
            logging::warn("Field arrays not used by the element types of the model, ignored: {}.", unused);
    }

    void import_mesh(libMesh::MeshBase &mesh, const MeshView &view, const IdArray &type_ids)
    {
        if (mesh.default_mapping_type() != libMesh::LAGRANGE_MAP)
            fail("Mesh import: only the Lagrange mapping is supported.");
        mesh.clear();

        const std::size_t n_points = view.n_points();
        for (std::size_t i = 0; i < n_points; ++i)
            mesh.add_point(libMesh::Point(view.points(i, 0), view.points(i, 1), view.points(i, 2)),
                           static_cast<dof_id_type>(i));

        std::array<bool, 4> elems_of_dimension{};
        std::vector<dof_id_type> conn;
        const std::size_t n_cells = view.n_cells();
        for (std::size_t i = 0; i < n_cells; ++i)
        {
            auto elem = libMesh::Elem::build(elem_type_of(view.cell_types(i), i));
            const std::int64_t begin = view.offsets(i);
            const std::int64_t n_cell_points = view.offsets(i + 1) - begin;
            if (n_cell_points != static_cast<std::int64_t>(elem->n_nodes()))
                fail("Mesh import: cell " + std::to_string(i) + " has " + std::to_string(n_cell_points) +
                     " points, its cell type " + std::to_string(elem->n_nodes()) + ".");

            for (auto j : elem->node_index_range())
            {
                const std::int64_t point = view.connectivity(static_cast<std::size_t>(begin) + j);
                if (point < 0 || static_cast<std::size_t>(point) >= n_points)
                    fail("Mesh import: cell " + std::to_string(i) + " refers to the point " + std::to_string(point) +
                         ", out of range.");
                elem->set_node(j) = mesh.node_ptr(static_cast<dof_id_type>(point));
            }

            // The node order of the VTK cell, in libMesh order.
            elem->connectivity(0, libMesh::VTK, conn);
            for (unsigned int j = 0, n_conn = libMesh::cast_int<unsigned int>(conn.size()); j != n_conn; ++j)
                elem->set_node(j) = mesh.node_ptr(conn[j]);

            elem->set_id(static_cast<dof_id_type>(i));
            elem->subdomain_id() = static_cast<libMesh::subdomain_id_type>(type_ids(i));
            elems_of_dimension[elem->dim()] = true;
            mesh.add_elem(std::move(elem));
        }

        for (unsigned char i = 0; i != 4; ++i)
            if (elems_of_dimension[i])
                mesh.set_mesh_dimension(i);

#if LIBMESH_DIM < 3
        if (mesh.mesh_dimension() > LIBMESH_DIM)
            fail("Cannot read a " + std::to_string(mesh.mesh_dimension()) + "D mesh: libMesh is configured without " +
                 std::to_string(mesh.mesh_dimension()) + "D support.");
#endif
    }

    std::vector<std::string> missing_output_arrays(const libMesh::EquationSystems &systems, const std::string &system_name,
                                                   const ModelArrays &arrays)
    {
        const WritePlan plan = write_plan(systems, system_name);
        std::vector<std::string> missing;
        std::set<std::string> seen;
        for (const WriteEntry &entry : plan.entries)
            if (!target(arrays, entry, plan.nodal) && seen.insert(entry.array).second)
                missing.push_back(entry.array + (plan.nodal ? " (point)" : " (cell)"));
        return missing;
    }

    void write_system_to_arrays(libMesh::EquationSystems &systems, const std::string &system_name, const ModelArrays &arrays)
    {
        const WritePlan plan = write_plan(systems, system_name);
        const System &system = systems.get_system(system_name);

        std::unique_ptr<libMesh::NumericVector<Number>> parallel_solution;
        if (plan.nodal)
        {
            const std::set<std::string> system_names = {system_name};
            parallel_solution = systems.build_parallel_solution_vector(&system_names);
        }
        else
        {
            if (plan.tuple_stride == 0)
                fail("Solution output: the element system has no variable.");
            // libMesh selects a vector variable by the names of its components, n_vec_dim of them.
            static const std::array<const char *, 3> suffix = {"_x", "_y", "_z"};
            std::vector<std::string> var_names;
            for (unsigned int i = 0; i < system.n_vars(); ++i)
            {
                const libMesh::FEType &type = system.variable_type(i);
                if (FEInterface::field_type(type) == libMesh::TYPE_VECTOR)
                    for (unsigned int c = 0; c < FEInterface::n_vec_dim(system.get_mesh(), type); ++c)
                        var_names.push_back(system.variable_name(i) + suffix[c]);
                else
                    var_names.push_back(system.variable_name(i));
            }
            parallel_solution = systems.build_parallel_elemental_solution_vector(var_names);
        }
        if (!parallel_solution)
            fail("Solution output: libMesh did not build the solution vector of '" + system_name + "'.");

        std::vector<Number> solution;
        parallel_solution->localize(solution);

        for (const WriteEntry &entry : plan.entries)
        {
            const RealArray *array = target(arrays, entry, plan.nodal);
            if (!array)
                continue;
            const RealArray &out = *array;
            const std::size_t n_tuples = out.n_tuples();
            for (std::size_t tuple = 0; tuple < n_tuples; ++tuple)
            {
                const std::size_t idx = plan.nodal ? tuple * plan.tuple_stride + entry.source_offset
                                                   : entry.source_offset * n_tuples + tuple;
                if (idx >= solution.size())
                    fail("Solution output: solution index out of range for the array '" + entry.array + "'.");
                out(tuple, entry.component) = libMesh::libmesh_real(solution[idx]);
            }
        }
    }
} // namespace pytec3d
