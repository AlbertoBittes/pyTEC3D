/**
 * @file quadrature.cpp
 * @brief Batched quadrature data of the mesh and reference shape tables (quadrature.h).
 */
#include "quadrature.h"

#include <cmath>
#include <string>

#include "logging.h"

namespace pytec3d::fem
{
    QuadratureData::QuadratureData(const FEMesh &mesh, const std::vector<int> &orders)
    {
        if (static_cast<int>(orders.size()) != mesh.n_local_elements())
            logging::error_and_abort("QuadratureData: one order per local element.");
        for (int e = 0; e < mesh.n_local_elements(); ++e)
            rules_.push_back(&mfem::IntRules.Get(mesh.mesh().GetElementGeometry(e), orders[static_cast<std::size_t>(e)]));
        compute(mesh);
    }

    QuadratureData::QuadratureData(const FEMesh &mesh, const std::vector<const mfem::IntegrationRule *> &rules)
        : rules_(rules)
    {
        if (static_cast<int>(rules.size()) != mesh.n_local_elements())
            logging::error_and_abort("QuadratureData: one rule per local element.");
        compute(mesh);
    }

    void QuadratureData::compute(const FEMesh &mesh)
    {
        const mfem::ParMesh &pmesh = mesh.mesh();
        const int n_elements = mesh.n_local_elements();
        const int dim = pmesh.Dimension();
        const int sdim = pmesh.SpaceDimension();
        offsets_.assign(static_cast<std::size_t>(n_elements) + 1, 0);
        for (int e = 0; e < n_elements; ++e)
            offsets_[static_cast<std::size_t>(e) + 1] = offsets_[static_cast<std::size_t>(e)] +
                                                       static_cast<std::size_t>(rules_[static_cast<std::size_t>(e)]->GetNPoints());
        const std::size_t n_points = offsets_.back();
        X_.assign(3 * n_points, 0.0);
        J_.assign(9 * n_points, 0.0);
        det_.assign(n_points, 0.0);

        // The geometry: the nodes of a curved mesh, else the vertices with the linear Lagrange elements.
        const mfem::GridFunction *nodes = pmesh.GetNodes();
        const mfem::FiniteElementSpace *nodes_fes = nodes ? nodes->FESpace() : nullptr;
        mfem::H1_FECollection linear(1, dim);

        // Batches of the elements of one (geometry, order): one table of the reference shapes of the geometry.
        std::map<std::pair<int, const mfem::IntegrationRule *>, std::vector<int>> batches;
        for (int e = 0; e < n_elements; ++e)
            batches[{pmesh.GetElementGeometry(e), rules_[static_cast<std::size_t>(e)]}].push_back(e);

        mfem::Array<int> dofs;
        std::vector<double> x;
        for (const auto &[key, elements] : batches)
        {
            const auto geometry = static_cast<mfem::Geometry::Type>(key.first);
            const mfem::IntegrationRule &rule = *rules_[static_cast<std::size_t>(elements.front())];
            const mfem::FiniteElement &fe = nodes_fes ? *nodes_fes->GetFE(elements.front())
                                                      : *linear.FiniteElementForGeometry(geometry);
            const int nq = rule.GetNPoints();
            const int nd = fe.GetDof();
            std::vector<double> B(static_cast<std::size_t>(nq * nd)), G(static_cast<std::size_t>(nq * nd * dim));
            mfem::Vector shape(nd);
            mfem::DenseMatrix dshape(nd, dim);
            for (int q = 0; q < nq; ++q)
            {
                fe.CalcShape(rule.IntPoint(q), shape);
                fe.CalcDShape(rule.IntPoint(q), dshape);
                for (int i = 0; i < nd; ++i)
                {
                    B[static_cast<std::size_t>(q * nd + i)] = shape(i);
                    for (int k = 0; k < dim; ++k)
                        G[static_cast<std::size_t>((q * nd + i) * dim + k)] = dshape(i, k);
                }
            }

            // The node coordinates of the batch, [element][node][3], then J, X, det for all of them.
            x.assign(elements.size() * static_cast<std::size_t>(nd) * 3, 0.0);
            for (std::size_t b = 0; b < elements.size(); ++b)
            {
                const int e = elements[b];
                double *xe = x.data() + b * static_cast<std::size_t>(nd) * 3;
                if (nodes_fes)
                {
                    nodes_fes->GetElementDofs(e, dofs);
                    for (int i = 0; i < nd; ++i)
                        for (int c = 0; c < sdim; ++c)
                            xe[3 * i + c] = (*nodes)(nodes_fes->DofToVDof(dofs[i], c));
                }
                else
                {
                    pmesh.GetElementVertices(e, dofs);
                    for (int i = 0; i < nd; ++i)
                        for (int c = 0; c < sdim; ++c)
                            xe[3 * i + c] = pmesh.GetVertex(dofs[i])[c];
                }
            }
            for (std::size_t b = 0; b < elements.size(); ++b)
            {
                const int e = elements[b];
                const double *xe = x.data() + b * static_cast<std::size_t>(nd) * 3;
                const std::size_t first = offsets_[static_cast<std::size_t>(e)];
                for (int q = 0; q < nq; ++q)
                {
                    double *X = X_.data() + 3 * (first + static_cast<std::size_t>(q));
                    double *J = J_.data() + 9 * (first + static_cast<std::size_t>(q));
                    const double *Bq = B.data() + static_cast<std::size_t>(q * nd);
                    const double *Gq = G.data() + static_cast<std::size_t>(q * nd * dim);
                    for (int i = 0; i < nd; ++i)
                        for (int c = 0; c < sdim; ++c)
                        {
                            X[c] += Bq[i] * xe[3 * i + c];
                            for (int k = 0; k < dim; ++k)
                                J[3 * c + k] += xe[3 * i + c] * Gq[i * dim + k];
                        }
                    double d = 0.0;
                    if (dim == 3)
                        d = J[0] * (J[4] * J[8] - J[5] * J[7]) - J[1] * (J[3] * J[8] - J[5] * J[6]) +
                            J[2] * (J[3] * J[7] - J[4] * J[6]);
                    else if (dim == 2)
                        d = J[0] * J[4] - J[1] * J[3];
                    det_[first + static_cast<std::size_t>(q)] = d;
                }
            }
        }
    }

    const ReferenceTable &ReferenceTables::get(const mfem::FiniteElement &fe, const mfem::IntegrationRule &rule) const
    {
        const auto key = std::make_pair(&fe, &rule);
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = tables_.find(key);
        if (it != tables_.end())
            return *it->second;

        auto table = std::make_unique<ReferenceTable>();
        const int nq = rule.GetNPoints();
        const int nd = fe.GetDof();
        const int dim = fe.GetDim();
        table->n_points = nq;
        table->n_dofs = nd;
        table->dim = dim;
        if (fe.GetRangeType() == mfem::FiniteElement::SCALAR)
        {
            table->B.resize(static_cast<std::size_t>(nq * nd));
            table->G.resize(static_cast<std::size_t>(nq * nd * dim));
            mfem::Vector shape(nd);
            mfem::DenseMatrix dshape(nd, dim);
            for (int q = 0; q < nq; ++q)
            {
                fe.CalcShape(rule.IntPoint(q), shape);
                fe.CalcDShape(rule.IntPoint(q), dshape);
                for (int i = 0; i < nd; ++i)
                {
                    table->B[static_cast<std::size_t>(q * nd + i)] = shape(i);
                    for (int k = 0; k < dim; ++k)
                        table->G[static_cast<std::size_t>((q * nd + i) * dim + k)] = dshape(i, k);
                }
            }
        }
        else
        {
            const int cdim = dim == 3 ? 3 : 1;
            table->curl_dim = cdim;
            table->B.resize(static_cast<std::size_t>(nq * nd * dim));
            table->G.resize(static_cast<std::size_t>(nq * nd * cdim));
            mfem::DenseMatrix vshape(nd, dim), curl(nd, cdim);
            for (int q = 0; q < nq; ++q)
            {
                fe.CalcVShape(rule.IntPoint(q), vshape);
                fe.CalcCurlShape(rule.IntPoint(q), curl);
                for (int i = 0; i < nd; ++i)
                {
                    for (int k = 0; k < dim; ++k)
                        table->B[static_cast<std::size_t>((q * nd + i) * dim + k)] = vshape(i, k);
                    for (int k = 0; k < cdim; ++k)
                        table->G[static_cast<std::size_t>((q * nd + i) * cdim + k)] = curl(i, k);
                }
            }
        }
        return *tables_.emplace(key, std::move(table)).first->second;
    }
} // namespace pytec3d::fem
