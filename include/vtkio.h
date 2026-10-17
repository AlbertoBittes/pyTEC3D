// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
// Distributed under the GNU LGPL v2.1 or later.  See COPYING for the full licence.

/**
 * @file vtkio.h
 * @brief Extended libMesh VTK I/O that handles pyTEC3D-specific data arrays.
 *
 * Derives from libMesh::MeshInput and MeshOutput.  The read() path imports a
 * vtkUnstructuredGrid into a libMesh DistributedMesh, honouring infinite-element
 * flags and element-type subdomain IDs.  The write_system_grid() path exports a
 * named solution system's fields back onto the VTK grid so Python can access them.
 */

#pragma once

// Local includes
#include <libmesh/libmesh_common.h>
#include <libmesh/mesh_input.h>
#include <libmesh/mesh_output.h>
#include <libmesh/utility.h>
#include <libmesh/equation_systems.h>
#include <libmesh/ignore_warnings.h>
#include <libmesh/restore_warnings.h>

// Ignore "deprecated...header" warning from strstream
#include <vtkType.h>
#include <vtkSmartPointer.h>

// C++ includes
#include <cstddef>
#include <map>
#include <string>
#include <vector>
// Forward declarations
class vtkUnstructuredGrid;

/**
 * @brief Extended libMesh VTK mesh reader/writer for the pyTEC3D solver.
 *
 * In read mode, the constructor accepts an existing vtkUnstructuredGrid and
 * translates it into a libMesh::MeshBase.  In write mode (write_system_grid()),
 * a named equation system's solution is scattered back onto the VTK grid so that
 * Python can inspect it without extra data copies.
 */
class VTKIO : public libMesh::MeshInput<libMesh::MeshBase>,
              public libMesh::MeshOutput<libMesh::MeshBase>
{
public:
    explicit VTKIO(libMesh::MeshBase &mesh, vtkUnstructuredGrid *vtk_grid);

    explicit VTKIO(const libMesh::MeshBase &mesh, libMesh::EquationSystems *eq_systems, vtkUnstructuredGrid *vtk_grid);

    using libMesh::MeshOutput<libMesh::MeshBase>::write_nodal_data;

    virtual void write_nodal_data(const std::string &,
                                  const std::vector<libMesh::Number> &,
                                  const std::vector<std::string> &) override;

    /** @brief Reads the given @p grid into the libMesh mesh. */
    virtual void read(vtkUnstructuredGrid *grid);
    /** @brief No-op override satisfying the libMesh::MeshInput<> pure-virtual requirement. */
    virtual void read(const std::string &) override {}
    virtual void write(const std::string &) override {}

    /** @brief Exports the named solution system's fields onto the active VTK grid. */
    void write_system_grid(const std::string &system_name);

    /** @brief Replaces the active VTK grid pointer. */
    void set_vtk_grid(vtkUnstructuredGrid *grid);

    /** @brief Toggles zlib compression for any VTK output. */
    void set_compression(bool b);

    /** @brief Returns the active VTK grid pointer. */
    vtkUnstructuredGrid *get_vtk_grid();

private:
    /** @brief Translates the active VTK grid into the libMesh mesh. */
    void build_mesh_from_active_vtk_grid();

    vtkSmartPointer<vtkUnstructuredGrid> _vtk_grid; ///< Shared VTK grid (read: allocated internally; write: supplied externally).

    bool _compress; ///< When true, VTK output is zlib-compressed.
    libMesh::EquationSystems *_eq_sys;
    std::map<libMesh::dof_id_type, libMesh::dof_id_type> _local_node_map;

    struct ElementMaps
    {
        /** @brief Associates @p libmesh_type with @p vtk_type (searchable in both directions). */
        void associate(libMesh::ElemType libmesh_type, vtkIdType vtk_type);

        /** @brief The VTK type of @p libmesh_type; aborts via libMesh diagnostics if unknown. */
        vtkIdType find(libMesh::ElemType libmesh_type);

        /** @brief The libMesh type of @p vtk_type; aborts via libMesh diagnostics if unknown. */
        libMesh::ElemType find(vtkIdType vtk_type);

        std::map<libMesh::ElemType, vtkIdType> writing_map;
        std::map<vtkIdType, libMesh::ElemType> reading_map;
    };

    static std::map<libMesh::ElemMappingType, ElementMaps> _element_maps;

    static std::map<libMesh::ElemMappingType, ElementMaps> build_element_maps();

    // void broadcast_all_fe_data_objects();
};
