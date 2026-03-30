// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef VTK_IO_UPDATED_H
#define VTK_IO_UPDATED_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/mesh_input.h"
#include "libmesh/mesh_output.h"
#include "libmesh/utility.h"
#include "libmesh/equation_systems.h"

#ifdef LIBMESH_HAVE_VTK
// Ignore "deprecated...header" warning from strstream
#include "libmesh/ignore_warnings.h"
#include "vtkType.h"
#include "vtkSmartPointer.h"
#include "libmesh/restore_warnings.h"
#endif

// C++ includes
#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include "element_types.h"

using namespace libMesh;
// Forward declarations
class vtkUnstructuredGrid;

// Extended VTK I/O class from libMesh (see libMesh documentation), it allows to read also
// other data of the model (ElementTypes, Couplings, infinite boundaries, ...)
class VTKIO_Updated : public MeshInput<MeshBase>,
                      public MeshOutput<MeshBase>
{
public:
    explicit VTKIO_Updated(MeshBase &mesh, vtkUnstructuredGrid *vtk_grid);

    explicit VTKIO_Updated(const MeshBase &mesh, EquationSystems *eq_systems, vtkUnstructuredGrid *vtk_grid);

    using MeshOutput<MeshBase>::write_nodal_data;

    virtual void write_nodal_data(const std::string &,
                                  const std::vector<Number> &,
                                  const std::vector<std::string> &) override;

    virtual void read(const std::string &) override {}
    virtual void read(vtkUnstructuredGrid *grid);
    virtual void write(const std::string &) override {}

    void write_system_grid(const std::string &system_name);
    void set_vtk_grid(vtkUnstructuredGrid *grid);

#ifdef LIBMESH_HAVE_VTK

    void set_compression(bool b);

    vtkUnstructuredGrid *get_vtk_grid();

    // void read_apec_m_model(const std::string &mesh_filepath);

private:
    void build_mesh_from_active_vtk_grid();

    // void system_vectors_to_vtk(const EquationSystems & es,
    //                            vtkUnstructuredGrid * & grid);

    vtkSmartPointer<vtkUnstructuredGrid> _vtk_grid;

    bool _compress;
    EquationSystems *_eq_sys;
    std::map<dof_id_type, dof_id_type> _local_node_map;

    struct ElementMaps
    {
        // Associate libmesh_type with vtk_type (searchable in both directions).
        void associate(ElemType libmesh_type, vtkIdType vtk_type)
        {
            writing_map[libmesh_type] = vtk_type;
            reading_map[vtk_type] = libmesh_type;
        }

        // Find an entry in the writing map, or throw an error.
        vtkIdType find(ElemType libmesh_type)
        {
            return libmesh_map_find(writing_map, libmesh_type);
        }

        // Find an entry in the reading map, or throw an error.
        ElemType find(vtkIdType vtk_type)
        {
            return libmesh_map_find(reading_map, vtk_type);
        }

        std::map<ElemType, vtkIdType> writing_map;
        std::map<vtkIdType, ElemType> reading_map;
    };

    static std::map<ElemMappingType, ElementMaps> _element_maps;

    static std::map<ElemMappingType, ElementMaps> build_element_maps();

    // void broadcast_all_fe_data_objects();

#endif
};

#endif // VTKIO_UPDATED_H
