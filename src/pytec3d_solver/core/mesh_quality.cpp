/**
 * @file mesh_quality.cpp
 * @brief Edge ratio, angles and scaled Jacobian of the mesh elements.
 */
#include "mesh_quality.h"

#include "logging.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_quality.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_to_string.h>
#include <libmesh/fe_base.h>
#include <libmesh/libmesh_exceptions.h>
#include <libmesh/mesh_base.h>
#include <libmesh/parallel_implementation.h>
#include <libmesh/quadrature_gauss.h>

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <utility>

namespace pytec3d
{
    namespace
    {
        /** @brief An FE map and quadrature rule per (dimension, order). */
        class JacobianSampler
        {
        public:
            /** @brief min / max of det(J) over the quadrature points of @p elem; 0 or less if degenerate. */
            libMesh::Real scaled_jacobian(const libMesh::Elem &elem)
            {
                Entry &entry = entry_for(elem);
                try
                {
                    entry.fe->reinit(&elem);
                }
                catch (const libMesh::DegenerateMap &)
                {
                    return 0.0;
                }
                const std::vector<libMesh::Real> &JxW = entry.fe->get_JxW();
                const std::vector<libMesh::Real> &weights = entry.qrule->get_weights();
                libMesh::Real lowest = std::numeric_limits<libMesh::Real>::max();
                libMesh::Real highest = -lowest;
                for (std::size_t q = 0; q < JxW.size(); ++q)
                {
                    const libMesh::Real det = JxW[q] / weights[q];
                    lowest = std::min(lowest, det);
                    highest = std::max(highest, det);
                }
                return highest > 0.0 ? lowest / highest : 0.0;
            }

        private:
            struct Entry
            {
                std::unique_ptr<libMesh::FEBase> fe;
                std::unique_ptr<libMesh::QGauss> qrule;
            };

            Entry &entry_for(const libMesh::Elem &elem)
            {
                const auto key = std::make_pair(elem.dim(), static_cast<int>(elem.default_order()));
                Entry &entry = entries_[key];
                if (!entry.fe)
                {
                    entry.fe = libMesh::FEBase::build(elem.dim(), libMesh::FEType(elem.default_order(), libMesh::LAGRANGE));
                    entry.qrule = std::make_unique<libMesh::QGauss>(elem.dim(),
                                                                    static_cast<libMesh::Order>(2 * elem.default_order()));
                    entry.fe->attach_quadrature_rule(entry.qrule.get());
                    entry.fe->get_JxW();
                }
                return entry;
            }

            std::map<std::pair<unsigned short, int>, Entry> entries_;
        };
    } // namespace

    std::vector<ShapeQuality> compute_mesh_quality(const libMesh::MeshBase &mesh)
    {
        // Indexed by ElemType, so that the reductions line up on every rank.
        constexpr std::size_t n_types = libMesh::INVALID_ELEM;
        const libMesh::Real huge = std::numeric_limits<libMesh::Real>::max();
        std::vector<libMesh::Real> max_edge_ratio(n_types, 0.0), min_angle(n_types, huge), max_angle(n_types, 0.0),
            min_scaled(n_types, huge);
        std::vector<std::size_t> count(n_types, 0), distorted(n_types, 0), degenerate(n_types, 0);
        std::vector<libMesh::dof_id_type> worst(n_types, 0);

        JacobianSampler sampler;
        for (const libMesh::Elem *elem : mesh.active_local_element_ptr_range())
        {
            const auto t = static_cast<std::size_t>(elem->type());
            if (elem->dim() == 0 || t >= n_types)
                continue;
            ++count[t];
            // The generic measures of libMesh::Elem, the same for every shape.
            max_edge_ratio[t] = std::max(max_edge_ratio[t], elem->libMesh::Elem::quality(libMesh::EDGE_LENGTH_RATIO));
            if (elem->dim() >= 2)
            {
                min_angle[t] = std::min(min_angle[t], elem->libMesh::Elem::quality(libMesh::MIN_ANGLE));
                max_angle[t] = std::max(max_angle[t], elem->libMesh::Elem::quality(libMesh::MAX_ANGLE));
            }
            const libMesh::Real scaled = sampler.scaled_jacobian(*elem);
            if (scaled <= 0.0)
                ++degenerate[t];
            else if (scaled < kDistortedScaledJacobian)
                ++distorted[t];
            if (scaled < min_scaled[t])
            {
                min_scaled[t] = scaled;
                worst[t] = elem->id();
            }
        }

        const auto &comm = mesh.comm();
        comm.sum(count);
        comm.sum(distorted);
        comm.sum(degenerate);
        comm.max(max_edge_ratio);
        comm.min(min_angle);
        comm.max(max_angle);

        std::vector<ShapeQuality> out;
        for (std::size_t t = 0; t < n_types; ++t)
        {
            if (count[t] == 0)
                continue;
            libMesh::Real scaled = min_scaled[t];
            unsigned int rank = 0;
            comm.minloc(scaled, rank);
            libMesh::dof_id_type worst_id = worst[t];
            comm.broadcast(worst_id, rank);

            ShapeQuality quality;
            quality.shape = libMesh::Utility::enum_to_string(static_cast<libMesh::ElemType>(t));
            quality.elements = count[t];
            quality.max_edge_ratio = max_edge_ratio[t];
            quality.min_angle = min_angle[t] == huge ? 0.0 : min_angle[t];
            quality.max_angle = max_angle[t];
            quality.min_scaled_jacobian = scaled;
            quality.worst_element = worst_id;
            quality.distorted = distorted[t];
            quality.degenerate = degenerate[t];
            out.push_back(std::move(quality));
        }
        return out;
    }

    void check_mesh_quality(const std::vector<ShapeQuality> &quality)
    {
        for (const ShapeQuality &shape : quality)
        {
            if (shape.degenerate > 0)
                logging::warn("{} {} element(s) are degenerate or inverted (det(J) <= 0 at a quadrature point; e.g. cell "
                              "{}): the solution there is meaningless.",
                              shape.degenerate, shape.shape, shape.worst_element);
            else if (shape.distorted > 0)
                logging::warn("{} {} element(s) are distorted (scaled Jacobian < {}; the worst, {:.3f}, at cell {}).",
                              shape.distorted, shape.shape, kDistortedScaledJacobian, shape.min_scaled_jacobian,
                              shape.worst_element);
        }
    }
} // namespace pytec3d
