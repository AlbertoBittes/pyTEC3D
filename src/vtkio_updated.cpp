// The libMesh Finite Element Library. Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson,
// Roy H. Stogner. Licensed under the GNU LGPL v2.1 or later.
/**
 * @file vtkio_updated.cpp
 * @brief VTKIOUpdated implementation: VTK mesh read/write with APEC-M solution-vector extensions.
 *
 * Derived from the libMesh VTKIO class; adds ability to read/write APEC-M solution vectors
 * stored as VTK point- and cell-data arrays alongside the standard mesh topology.
 */

// Local includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_base.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/system.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/dof_map.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/mesh_communication.h"

#include "vtkio_updated.h"

// I get a lot of "warning: extra ';' inside a class [-Wextra-semi]" from clang
// on VTK header files.
#include "libmesh/ignore_warnings.h"
#include "libmesh/fe_interface.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataArray.h"
#include "vtkIntArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkStringArray.h"
#include "vtkVersion.h"
#include "vtkVersionMacros.h"

#include "libmesh/restore_warnings.h"
#include "variables_map.h"
#include "apec_m_system.h"
#include "apec_m_logging.h"
// C++ includes
#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>

using libMesh::cast_int;
using libMesh::dof_id_type;
using libMesh::EDGE2;
using libMesh::EDGE3;
using libMesh::Elem;
using libMesh::ElemMappingType;
using libMesh::ElemType;
using libMesh::EquationSystems;
using libMesh::FEFamily;
using libMesh::FEInterface;
using libMesh::HEX20;
using libMesh::HEX27;
using libMesh::HEX8;
using libMesh::MeshBase;
using libMesh::MeshInput;
using libMesh::MeshOutput;
using libMesh::Node;
using libMesh::Number;
using libMesh::Order;
using libMesh::Point;
using libMesh::PRISM15;
using libMesh::PRISM18;
using libMesh::PRISM6;
using libMesh::PYRAMID5;
using libMesh::QUAD4;
using libMesh::QUAD8;
using libMesh::QUAD9;
using libMesh::Real;
using libMesh::System;
using libMesh::TET10;
using libMesh::TET4;
using libMesh::TRI3;
using libMesh::TRI3SUBDIVISION;
using libMesh::TRI6;
using libMesh::VTK;

namespace Utility = libMesh::Utility;

// A convenient macro for comparing VTK versions.  Returns 1 if the
// // current VTK version is < major.minor.subminor and zero otherwise.
// #define VTK_VERSION_LESS_THAN(major, minor, subminor)                                                            \
//     (VTK_VERSION_NUMBER < VTK_VERSION_CHECK((major), (minor), (subminor)))

namespace
{
    struct GridWritePlan
    {
        std::vector<VariableInfo> vars_to_write;
        std::vector<std::size_t> source_offsets;
        std::size_t tuple_stride = 0;
    };

    GridWritePlan build_grid_write_plan(const System &system, const bool nodal_layout)
    {
        GridWritePlan plan;

        const MeshBase &mesh = system.get_mesh();
        std::size_t running_offset = 0;
        for (unsigned int var = 0; var < system.n_vars(); ++var)
        {
            const std::string &var_name = system.variable_name(var);
            unsigned int n_components = 1u;
            if (nodal_layout)
            {
                const auto n_vec_dim = FEInterface::n_vec_dim(mesh, system.variable_type(var));
                if (n_vec_dim > 1)
                    n_components = cast_int<unsigned int>(n_vec_dim);
            }

            const auto keys = VariablesMap::get_keys_by_predicate(
                [&var_name](const VariableInfo &info)
                { return info.global_label == var_name; });

            // Vector-valued libMesh variables are present in storage, but not written to VTK arrays.
            if (keys.size() == 1 && n_components == 1)
            {
                plan.vars_to_write.push_back(VariablesMap::get_info_by_key(keys.front()));
                plan.source_offsets.push_back(running_offset);
            }

            running_offset += n_components;
        }

        plan.tuple_stride = running_offset;
        return plan;
    }

    void write_solution_to_grid(vtkDataSetAttributes *data_attributes,
                                const vtkIdType n_tuples,
                                const GridWritePlan &plan,
                                const std::vector<Number> &localized_solution,
                                const bool nodal_layout)
    {
        if (!nodal_layout && plan.tuple_stride == 0)
            apec::logging::error_and_abort("VTKIO_Updated::write_solution_to_grid: invalid elemental tuple stride (zero).");

        const std::size_t n_vars = plan.vars_to_write.size();
        for (std::size_t var = 0; var < n_vars; ++var)
        {
            const auto &var_info = plan.vars_to_write[var];
            vtkDataArray *array =
                vtkDataArray::SafeDownCast(data_attributes->GetArray(var_info.api_label.c_str()));
            if (!array)
                apec::logging::error_and_abort("VTKIO_Updated::write_solution_to_grid: VTK array '" +
                                               var_info.api_label + "' not found.");

            const std::size_t source_offset = plan.source_offsets[var];

            for (vtkIdType tuple_id = 0; tuple_id < n_tuples; ++tuple_id)
            {
                const std::size_t idx = nodal_layout
                                            ? static_cast<std::size_t>(tuple_id) * plan.tuple_stride + source_offset
                                            : source_offset * static_cast<std::size_t>(n_tuples) + static_cast<std::size_t>(tuple_id);
                if (idx >= localized_solution.size())
                    apec::logging::error_and_abort("VTKIO_Updated::write_solution_to_grid: computed solution index out of range.");
                array->SetComponent(tuple_id, var_info.component_index, localized_solution[idx]);
            }
        }
    }
} // namespace

VTKIO_Updated::VTKIO_Updated(MeshBase &mesh, vtkUnstructuredGrid *vtk_grid) : MeshInput<MeshBase>(mesh, /*is_parallel_format=*/true),
                                                                              MeshOutput<MeshBase>(mesh, /*is_parallel_format=*/true),
                                                                              _vtk_grid(vtk_grid), _compress(false), _eq_sys(nullptr)
{
}

VTKIO_Updated::VTKIO_Updated(const MeshBase &mesh, EquationSystems *eq_systems, vtkUnstructuredGrid *vtk_grid) : MeshOutput<MeshBase>(mesh, /*is_parallel_format=*/true),
                                                                                                                 _vtk_grid(vtk_grid), _compress(false), _eq_sys(eq_systems)
{
}

// Static initialisation of the element type maps.
std::map<ElemMappingType, VTKIO_Updated::ElementMaps> VTKIO_Updated::_element_maps = VTKIO_Updated::build_element_maps();

/// Constructs the map from libMesh element types to VTK cell types for each mapping strategy.
std::map<ElemMappingType, VTKIO_Updated::ElementMaps> VTKIO_Updated::build_element_maps()
{
    std::map<ElemMappingType, VTKIO_Updated::ElementMaps> all_maps;
    ElementMaps em; // Lagrange element maps

    em.associate(EDGE2, VTK_LINE);
    em.associate(EDGE3, VTK_QUADRATIC_EDGE);
    em.associate(TRI3, VTK_TRIANGLE);
    em.associate(TRI6, VTK_QUADRATIC_TRIANGLE);
    em.associate(QUAD4, VTK_QUAD);
    em.associate(QUAD8, VTK_QUADRATIC_QUAD);
    em.associate(TET4, VTK_TETRA);
    em.associate(TET10, VTK_QUADRATIC_TETRA);
    em.associate(HEX8, VTK_HEXAHEDRON);
    em.associate(HEX20, VTK_QUADRATIC_HEXAHEDRON);
    em.associate(HEX27, VTK_TRIQUADRATIC_HEXAHEDRON);
    em.associate(PRISM6, VTK_WEDGE);
    em.associate(PRISM15, VTK_QUADRATIC_WEDGE);
    em.associate(PRISM18, VTK_BIQUADRATIC_QUADRATIC_WEDGE);
    em.associate(PYRAMID5, VTK_PYRAMID);
    em.associate(QUAD9, VTK_BIQUADRATIC_QUAD);

    // TRI3SUBDIVISION is for writing only
    em.writing_map[TRI3SUBDIVISION] = VTK_TRIANGLE;
    all_maps[ElemMappingType::LAGRANGE_MAP] = em;

    // VTK_BEZIER_* element types require VTK >= 9.0; fall back to Lagrange maps on older builds.
    all_maps[ElemMappingType::RATIONAL_BERNSTEIN_MAP] = em;

    ElementMaps bem; // Rational Bernstein element maps (VTK 9.0+)
    bem.associate(EDGE2, VTK_LINE);
    bem.associate(EDGE3, VTK_BEZIER_CURVE);
    bem.associate(TRI3, VTK_TRIANGLE);
    bem.associate(TRI6, VTK_BEZIER_TRIANGLE);
    bem.associate(QUAD4, VTK_QUAD);
    bem.associate(QUAD8, VTK_QUADRATIC_QUAD);
    bem.associate(QUAD9, VTK_BEZIER_QUADRILATERAL);
    bem.associate(TET4, VTK_TETRA);
    bem.associate(TET10, VTK_QUADRATIC_TETRA);
    bem.associate(HEX8, VTK_HEXAHEDRON);
    bem.associate(HEX20, VTK_QUADRATIC_HEXAHEDRON);
    bem.associate(HEX27, VTK_BEZIER_HEXAHEDRON);
    bem.associate(PRISM6, VTK_WEDGE);
    bem.associate(PRISM15, VTK_QUADRATIC_WEDGE);
    bem.associate(PRISM18, VTK_BEZIER_WEDGE);
    bem.associate(PYRAMID5, VTK_PYRAMID);
    bem.writing_map[TRI3SUBDIVISION] = VTK_TRIANGLE;
    all_maps[ElemMappingType::RATIONAL_BERNSTEIN_MAP] = bem;

    return all_maps;
}

void VTKIO_Updated::read(vtkUnstructuredGrid *grid)
{
    if (!grid)
        apec::logging::error_and_abort("VTKIO_Updated::read(vtkUnstructuredGrid*): input grid is null.");

    _vtk_grid = grid;
    build_mesh_from_active_vtk_grid();
}

void VTKIO_Updated::set_vtk_grid(vtkUnstructuredGrid *grid)
{
    if (!grid)
        apec::logging::error_and_abort("VTKIO_Updated::set_vtk_grid: input grid is null.");
    _vtk_grid = grid;
}

void VTKIO_Updated::write_system_grid(const std::string &system_name)
{
    if (!_eq_sys)
        apec::logging::error_and_abort("VTKIO_Updated::write_system_grid: equation systems pointer is null.");
    if (!_vtk_grid)
        apec::logging::error_and_abort("VTKIO_Updated::write_system_grid: active vtk grid is null.");
    if (!_eq_sys->has_system(system_name))
        apec::logging::error_and_abort("VTKIO_Updated::write_system_grid: system '" + system_name + "' not found.");

    const auto system_type = ApecMSystemBase::system_type_from_name(system_name);
    const System &system = _eq_sys->get_system(system_name);
    if (system_type != ApecMSystemBase::ELEMENT_SOLUTION)
    {
        const std::set<std::string> system_names = {system_name};
        const auto write_plan = build_grid_write_plan(system, true);

        if ("Nodal Solution" == system_name)
            apec::logging::info("Averaging nodal solution at nodes");
        auto parallel_solution = _eq_sys->build_parallel_solution_vector(&system_names);
        if (!parallel_solution)
            apec::logging::error_and_abort("VTKIO_Updated::write_system_grid: failed to build parallel nodal solution vector.");

        std::vector<Number> localized_solution;
        parallel_solution->localize(localized_solution);
        write_solution_to_grid(_vtk_grid->GetPointData(),
                               _vtk_grid->GetNumberOfPoints(),
                               write_plan,
                               localized_solution,
                               true);
    }
    else
    {
        std::vector<std::string> var_names;
        var_names.reserve(system.n_vars());
        for (unsigned int i = 0; i < system.n_vars(); ++i)
            var_names.push_back(system.variable_name(i));
        const auto write_plan = build_grid_write_plan(system, false);

        auto parallel_solution = _eq_sys->build_parallel_elemental_solution_vector(var_names);
        if (!parallel_solution)
            apec::logging::error_and_abort("VTKIO_Updated::write_system_grid: failed to build parallel elemental solution vector.");

        std::vector<Number> localized_solution;
        parallel_solution->localize(localized_solution);
        write_solution_to_grid(_vtk_grid->GetCellData(),
                               _vtk_grid->GetNumberOfCells(),
                               write_plan,
                               localized_solution,
                               false);
    }

    _vtk_grid->Modified();
}

void VTKIO_Updated::build_mesh_from_active_vtk_grid()
{
    if (!_vtk_grid)
        apec::logging::error_and_abort("VTKIO_Updated: active vtk grid is null.");

    elems_of_dimension.clear();
    elems_of_dimension.resize(4, false);

    MeshBase &mesh = MeshInput<MeshBase>::mesh();
    mesh.clear();

    vtkAbstractArray *abstract_elemtype_id =
        _vtk_grid->GetCellData()->GetAbstractArray("TYPE");

    vtkIntArray *elemtype_id = vtkIntArray::SafeDownCast(abstract_elemtype_id);

    if (abstract_elemtype_id && !elemtype_id)
        libmesh_warning("Found non-integral ansys_elem_type array; forced to ignore it.\n"
                        "This is technically valid but probably broken.");

    const unsigned int vtk_num_points = static_cast<unsigned int>(_vtk_grid->GetNumberOfPoints());

    for (unsigned int i = 0; i < vtk_num_points; ++i)
    {
        double pnt[3];
        _vtk_grid->GetPoint(static_cast<vtkIdType>(i), pnt);
        Point xyz(pnt[0], pnt[1], pnt[2]);
        mesh.add_point(xyz, i);
    }

    const unsigned int vtk_num_cells = static_cast<unsigned int>(_vtk_grid->GetNumberOfCells());
    auto &element_map = libmesh_map_find(_element_maps, mesh.default_mapping_type());

    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    for (unsigned int i = 0; i < vtk_num_cells; ++i)
    {
        _vtk_grid->GetCell(i, cell);

        ElemType libmesh_elem_type = element_map.find(cell->GetCellType());
        auto elem = Elem::build(libmesh_elem_type);

        for (auto j : elem->node_index_range())
        {
            const auto vtk_point_id = cell->GetPointId(j);
            const dof_id_type libmesh_node_id = vtk_point_id;
            elem->set_node(j) = mesh.node_ptr(libmesh_node_id);
        }

        std::vector<dof_id_type> conn;
        elem->connectivity(0, VTK, conn);

        for (unsigned int j = 0, n_conn = cast_int<unsigned int>(conn.size());
             j != n_conn; ++j)
            elem->set_node(j) = mesh.node_ptr(conn[j]);

        elem->set_id(i);

        if (elemtype_id)
        {
            auto elemtype = elemtype_id->GetValue(i);
            elem->subdomain_id() = elemtype;
        }
        else
            apec::logging::error_and_abort("No TYPE array found in cell data; all subdomain ids set to zero.");

        elems_of_dimension[elem->dim()] = true;
        mesh.add_elem(std::move(elem));
    }

    for (unsigned char i = 0; i != 4; ++i)
        if (elems_of_dimension[i])
            mesh.set_mesh_dimension(i);

#if LIBMESH_DIM < 3
    if (mesh.mesh_dimension() > LIBMESH_DIM)
        apec::logging::error_and_abort("Cannot open dimension " + std::to_string(mesh.mesh_dimension()) +
                                       " mesh file when configured without " +
                                       std::to_string(mesh.mesh_dimension()) + "D support.");
#endif // LIBMESH_DIM < 3
}

vtkUnstructuredGrid *VTKIO_Updated::get_vtk_grid()
{
    return _vtk_grid;
}

void VTKIO_Updated::set_compression(bool b)
{
    this->_compress = b;
}

void VTKIO_Updated::write_nodal_data(const std::string &,
                                     const std::vector<Number> &,
                                     const std::vector<std::string> &)
{
}
