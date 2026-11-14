/**
 * @file dof_space.cpp
 * @brief Variables, element blocks and the DOFs of the DOF system (dof_space.h).
 */
#include "dof_space.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>

#include "fe_collections.h"
#include "logging.h"

namespace pytec3d::fem
{
    namespace
    {
        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("DOFs: " + message);
        }

        /** @brief The edges of a reference geometry: vertex pairs. */
        const int *reference_edges(mfem::Geometry::Type geometry, int &n)
        {
            using G = mfem::Geometry;
            switch (geometry)
            {
            case G::TRIANGLE: n = 3; return &G::Constants<G::TRIANGLE>::Edges[0][0];
            case G::SQUARE: n = 4; return &G::Constants<G::SQUARE>::Edges[0][0];
            case G::TETRAHEDRON: n = 6; return &G::Constants<G::TETRAHEDRON>::Edges[0][0];
            case G::CUBE: n = 12; return &G::Constants<G::CUBE>::Edges[0][0];
            case G::PRISM: n = 9; return &G::Constants<G::PRISM>::Edges[0][0];
            case G::PYRAMID: n = 8; return &G::Constants<G::PYRAMID>::Edges[0][0];
            default: n = 0; return nullptr;
            }
        }

        /** @brief Lexicographic order of coordinates (x, then y, then z), exact, as libMesh's Point comparison. */
        bool lexicographically_less(const double *a, const double *b, int dim)
        {
            for (int i = 0; i < dim; ++i)
            {
                if (a[i] < b[i])
                    return true;
                if (a[i] > b[i])
                    return false;
            }
            return false;
        }
    } // namespace

    // ---- FECatalog -------------------------------------------------------------------------------------------------

    FECatalog::FECatalog(int dimension, bool serendipity) : dimension_(dimension)
    {
        h1_[1] = std::make_unique<mfem::H1_FECollection>(1, dimension);
        h1_[2] = std::make_unique<QuadraticH1Collection>(dimension, serendipity);
        nd_ = std::make_unique<mfem::ND_FECollection>(1, dimension);
    }

    FECatalog::~FECatalog() = default;

    const mfem::FiniteElementCollection &FECatalog::collection(FEFamily family, int order) const
    {
        switch (family)
        {
        case FEFamily::LAGRANGE:
        case FEFamily::LAGRANGE_VEC:
        case FEFamily::L2_LAGRANGE:
        case FEFamily::L2_LAGRANGE_VEC:
            if (order < 1 || order > 2)
                fail("Lagrange fields of order " + std::to_string(order) + ": the orders are 1 and 2.");
            return *h1_[order];
        case FEFamily::NEDELEC_ONE:
            if (order != 1)
                fail("Nedelec fields of order " + std::to_string(order) + ": the order is 1.");
            return *nd_;
        default:
            fail(std::string("no finite element collection for the family ") + formulation::family_name(family) + ".");
        }
    }

    const mfem::FiniteElement *FECatalog::element(FEFamily family, int order, mfem::Geometry::Type geometry) const
    {
        if (family == FEFamily::MONOMIAL || family == FEFamily::MONOMIAL_VEC)
            return nullptr;
        const mfem::FiniteElement *fe = collection(family, order).FiniteElementForGeometry(geometry);
        if (!fe)
            fail(std::string("the family ") + formulation::family_name(family) + " of order " + std::to_string(order) +
                 " has no element on the geometry " + mfem::Geometry::Name[geometry] + ".");
        return fe;
    }

    // ---- SystemLayout ----------------------------------------------------------------------------------------------

    SystemLayout::SystemLayout(const FEMesh &mesh, const FECatalog &catalog, std::vector<Variable> variables)
        : mesh_(&mesh), catalog_(&catalog), variables_(std::move(variables))
    {
    }

    int SystemLayout::variable_number(const std::string &label) const
    {
        for (unsigned int v = 0; v < n_vars(); ++v)
            if (variables_[v].label == label)
                return static_cast<int>(v);
        return -1;
    }

    unsigned int SystemLayout::n_element_dofs(unsigned int v, int e) const
    {
        const Variable &variable = variables_[v];
        if (!variable.active_on(mesh_->subdomain(e)))
            return 0;
        if (variable.monomial())
            return variable.n_components;
        const mfem::FiniteElement *fe =
            catalog_->element(variable.family, variable.order, mesh_->mesh().GetElementGeometry(e));
        return static_cast<unsigned int>(fe->GetDof() * variable.vdim);
    }

    unsigned int SystemLayout::n_element_dofs(int e) const
    {
        unsigned int n = 0;
        for (unsigned int v = 0; v < n_vars(); ++v)
            n += n_element_dofs(v, e);
        return n;
    }

    // ---- DofSpace --------------------------------------------------------------------------------------------------

    DofSpace::DofSpace(const FEMesh &mesh, const FECatalog &catalog, std::vector<Variable> variables)
        : SystemLayout(mesh, catalog, std::move(variables))
    {
        auto &pmesh = const_cast<mfem::ParMesh &>(mesh.mesh());
        std::size_t offset = 0;
        for (const Variable &variable : variables_)
        {
            if (variable.monomial())
                fail("the DOF variable '" + variable.label + "' is of a monomial family.");
            // The DOFs of an edge field are addressed on the midside points (Dirichlet values, loads, couplings, the
            // tree gauge): a straight mesh has none.
            if (variable.edge() && mesh.geometric_order() < 2)
                fail("the edge (NEDELEC_ONE) field '" + variable.label + "' needs a quadratic mesh: its DOFs are on the "
                     "midside points of the cells.");
            const mfem::FiniteElementCollection &fec = catalog.collection(variable.family, variable.order);
            spaces_.push_back(std::make_unique<mfem::ParFiniteElementSpace>(
                &pmesh, const_cast<mfem::FiniteElementCollection *>(&fec), variable.vdim, mfem::Ordering::byNODES));
            offsets_.push_back(offset);
            offset += static_cast<std::size_t>(spaces_.back()->GetVSize());
        }
        offsets_.push_back(offset);
        global_.assign(offset, -1);
        owned_.assign(offset, 0);
        variable_of_.resize(offset);
        edge_sign_.assign(offset, 1.0);
        for (unsigned int v = 0; v < n_vars(); ++v)
            std::fill(variable_of_.begin() + static_cast<std::ptrdiff_t>(offsets_[v]),
                      variable_of_.begin() + static_cast<std::ptrdiff_t>(offsets_[v + 1]), v);

        build_numbering();
        build_edge_signs();
        build_element_dofs();

        // The points of the DOFs.
        point_dofs_.resize(n_vars());
        for (unsigned int v = 0; v < n_vars(); ++v)
        {
            const mfem::ParFiniteElementSpace &fes = *spaces_[v];
            const std::vector<std::int64_t> points = mesh.dof_points(fes);
            for (int d = 0; d < fes.GetNDofs(); ++d)
            {
                const std::int64_t p = points[static_cast<std::size_t>(d)];
                if (p < 0)
                    continue;
                bool unknown = false;
                for (int c = 0; c < fes.GetVDim(); ++c)
                    unknown = unknown || global_[local_index(v, fes.DofToVDof(d, c))] >= 0;
                if (unknown)
                    point_dofs_[v].emplace(p, d);
            }
        }

        // P_L: one 1 per local index that is an unknown, in the column of its global index.
        const MPI_Comm comm = mesh.comm();
        long long n_local = static_cast<long long>(global_.size()), first_local = 0, n_local_global = 0;
        MPI_Exscan(&n_local, &first_local, 1, MPI_LONG_LONG, MPI_SUM, comm);
        MPI_Allreduce(&n_local, &n_local_global, 1, MPI_LONG_LONG, MPI_SUM, comm);
        n_local_global_ = static_cast<HYPRE_BigInt>(n_local_global);
        local_partition_[0] = static_cast<HYPRE_BigInt>(first_local);
        local_partition_[1] = static_cast<HYPRE_BigInt>(first_local + n_local);
        true_partition_[0] = static_cast<HYPRE_BigInt>(first_owned_);
        true_partition_[1] = static_cast<HYPRE_BigInt>(first_owned_ + n_owned_);
        std::vector<int> I(global_.size() + 1, 0);
        std::vector<HYPRE_BigInt> J;
        std::vector<mfem::real_t> data;
        for (std::size_t l = 0; l < global_.size(); ++l)
        {
            if (global_[l] >= 0)
            {
                J.push_back(static_cast<HYPRE_BigInt>(global_[l]));
                data.push_back(1.0);
            }
            I[l + 1] = static_cast<int>(J.size());
        }
        prolongation_ = std::make_unique<mfem::HypreParMatrix>(
            comm, static_cast<int>(global_.size()), static_cast<HYPRE_BigInt>(n_local_global),
            static_cast<HYPRE_BigInt>(n_global_), I.data(), J.data(), data.data(), local_partition_, true_partition_);
    }

    DofSpace::~DofSpace() = default;

    void DofSpace::build_numbering()
    {
        const MPI_Comm comm = mesh_->comm();
        const int n_ranks = mesh_->n_ranks();
        mfem::Array<int> vdofs;

        // The unknowns: the DOFs of the elements where the variable is active, on any rank.
        std::vector<mfem::Array<int>> active(n_vars());
        for (unsigned int v = 0; v < n_vars(); ++v)
        {
            mfem::ParFiniteElementSpace &fes = *spaces_[v];
            active[v].SetSize(fes.GetVSize());
            active[v] = 0;
            for (int e = 0; e < mesh_->n_local_elements(); ++e)
            {
                if (!variables_[v].active_on(mesh_->subdomain(e)))
                    continue;
                fes.GetElementVDofs(e, vdofs);
                for (const int d : vdofs)
                    active[v][d >= 0 ? d : -1 - d] = 1;
            }
            fes.GroupComm().Reduce<int>(active[v].GetData(), mfem::GroupCommunicator::BitOR<int>);
            fes.GroupComm().Bcast(active[v].GetData());
        }

        // The owned unknowns of this rank, variable by variable, numbered after those of the lower ranks.
        std::vector<PetscInt> owned_count(n_vars(), 0);
        for (unsigned int v = 0; v < n_vars(); ++v)
            for (int d = 0; d < spaces_[v]->GetVSize(); ++d)
                if (active[v][d] && spaces_[v]->GetLocalTDofNumber(d) >= 0)
                    ++owned_count[v];
        n_owned_ = 0;
        for (const PetscInt c : owned_count)
            n_owned_ += c;
        owned_by_rank_.assign(static_cast<std::size_t>(n_ranks), 0);
        MPI_Allgather(&n_owned_, 1, MPIU_INT, owned_by_rank_.data(), 1, MPIU_INT, comm);
        first_owned_ = 0;
        for (int r = 0; r < mesh_->rank(); ++r)
            first_owned_ += owned_by_rank_[static_cast<std::size_t>(r)];
        n_global_ = 0;
        for (const PetscInt c : owned_by_rank_)
            n_global_ += c;
        n_global_by_variable_.assign(n_vars(), 0);
        MPI_Allreduce(owned_count.data(), n_global_by_variable_.data(), static_cast<int>(n_vars()), MPIU_INT, MPI_SUM, comm);

        PetscInt next = first_owned_;
        for (unsigned int v = 0; v < n_vars(); ++v)
        {
            mfem::ParFiniteElementSpace &fes = *spaces_[v];
            // Global indices as doubles (exact below 2^53) for the broadcast from the owner to the sharing ranks.
            mfem::Array<double> numbers(fes.GetVSize());
            numbers = -1.0;
            for (int d = 0; d < fes.GetVSize(); ++d)
                if (active[v][d] && fes.GetLocalTDofNumber(d) >= 0)
                {
                    numbers[d] = static_cast<double>(next++);
                    owned_[local_index(v, d)] = 1;
                }
            fes.GroupComm().Bcast(numbers.GetData());
            for (int d = 0; d < fes.GetVSize(); ++d)
                global_[local_index(v, d)] = active[v][d] ? static_cast<PetscInt>(numbers[d]) : -1;
        }
        if (next != first_owned_ + n_owned_)
            fail("inconsistent numbering of the owned DOFs.");
        for (std::size_t l = 0; l < global_.size(); ++l)
            if (owned_[l] == 0 && global_[l] >= 0 &&
                (global_[l] >= first_owned_ && global_[l] < first_owned_ + n_owned_))
                fail("a shared DOF has an index of this rank.");
    }

    void DofSpace::build_edge_signs()
    {
        // The sign of each Nedelec DOF: +1 when the global shape function of MFEM has a positive circulation along its
        // edge from the lexicographically smaller end point, else -1. On the reference element, the circulation of a
        // covariant shape function is that of the physical one along the image of the reference edge.
        const mfem::ParMesh &pmesh = mesh_->mesh();
        const int dim = mesh_->dimension();
        mfem::Array<int> vdofs;
        mfem::DenseMatrix shape;
        const mfem::IntegrationRule &line = mfem::IntRules.Get(mfem::Geometry::SEGMENT, 7);
        for (unsigned int v = 0; v < n_vars(); ++v)
        {
            if (!variables_[v].edge())
                continue;
            const mfem::ParFiniteElementSpace &fes = *spaces_[v];
            std::vector<char> done(static_cast<std::size_t>(fes.GetVSize()), 0);
            for (int e = 0; e < mesh_->n_local_elements(); ++e)
            {
                fes.GetElementVDofs(e, vdofs);
                const mfem::FiniteElement &fe = *fes.GetFE(e);
                const mfem::Geometry::Type geometry = fe.GetGeomType();
                const mfem::IntegrationRule &corners = *mfem::Geometries.GetVertices(geometry);
                const int *vertices = pmesh.GetElement(e)->GetVertices();
                int n_edges = 0;
                const int *edges = reference_edges(geometry, n_edges);
                // Circulation of every shape function along every reference edge, a -> b.
                mfem::DenseMatrix circulation(fe.GetDof(), n_edges);
                circulation = 0.0;
                shape.SetSize(fe.GetDof(), dim);
                for (int k = 0; k < n_edges; ++k)
                {
                    const mfem::IntegrationPoint &a = corners.IntPoint(edges[2 * k]);
                    const mfem::IntegrationPoint &b = corners.IntPoint(edges[2 * k + 1]);
                    const double t[3] = {b.x - a.x, b.y - a.y, b.z - a.z};
                    for (int q = 0; q < line.GetNPoints(); ++q)
                    {
                        const mfem::IntegrationPoint &s = line.IntPoint(q);
                        mfem::IntegrationPoint ip;
                        ip.Set3(a.x + s.x * t[0], a.y + s.x * t[1], a.z + s.x * t[2]);
                        fe.CalcVShape(ip, shape);
                        for (int i = 0; i < fe.GetDof(); ++i)
                        {
                            double dot = 0.0;
                            for (int c = 0; c < dim; ++c)
                                dot += shape(i, c) * t[c];
                            circulation(i, k) += s.weight * dot;
                        }
                    }
                }
                for (int i = 0; i < fe.GetDof(); ++i)
                {
                    const int d = vdofs[i] >= 0 ? vdofs[i] : -1 - vdofs[i];
                    if (done[static_cast<std::size_t>(d)])
                        continue;
                    int k_best = 0;
                    for (int k = 1; k < n_edges; ++k)
                        if (std::abs(circulation(i, k)) > std::abs(circulation(i, k_best)))
                            k_best = k;
                    const double c = circulation(i, k_best);
                    if (std::abs(c) < 0.5)
                        fail("a Nedelec shape function has no circulation along an edge of its element.");
                    const double *pa = pmesh.GetVertex(vertices[edges[2 * k_best]]);
                    const double *pb = pmesh.GetVertex(vertices[edges[2 * k_best + 1]]);
                    const double direction = lexicographically_less(pa, pb, dim) ? 1.0 : -1.0;
                    const double mfem_sign = vdofs[i] >= 0 ? 1.0 : -1.0;
                    edge_sign_[local_index(v, d)] = mfem_sign * c * direction > 0.0 ? 1.0 : -1.0;
                    done[static_cast<std::size_t>(d)] = 1;
                }
            }
        }
    }

    void DofSpace::build_element_dofs()
    {
        mfem::Array<int> vdofs;
        element_offsets_.assign(1, 0);
        for (int e = 0; e < mesh_->n_local_elements(); ++e)
        {
            for (unsigned int v = 0; v < n_vars(); ++v)
            {
                if (!variables_[v].active_on(mesh_->subdomain(e)))
                    continue;
                spaces_[v]->GetElementVDofs(e, vdofs);
                for (const int d : vdofs)
                {
                    const std::size_t l = local_index(v, d >= 0 ? d : -1 - d);
                    element_local_.push_back(static_cast<int>(l));
                    element_sign_.push_back((d >= 0 ? 1.0 : -1.0) * edge_sign_[l]);
                }
            }
            element_offsets_.push_back(element_local_.size());
        }
    }

    long long DofSpace::point_dof(unsigned int v, std::int64_t point, unsigned int c) const
    {
        const auto it = point_dofs_[v].find(point);
        if (it == point_dofs_[v].end())
            return -1;
        const mfem::ParFiniteElementSpace &fes = *spaces_[v];
        if (static_cast<int>(c) >= fes.GetVDim())
            return -1;
        const std::size_t l = local_index(v, fes.DofToVDof(it->second, static_cast<int>(c)));
        return global_[l] >= 0 ? static_cast<long long>(l) : -1;
    }

    void DofSpace::element_dofs(int e, const int *&local, const double *&sign, unsigned int &n) const
    {
        const std::size_t begin = element_offsets_[static_cast<std::size_t>(e)];
        n = static_cast<unsigned int>(element_offsets_[static_cast<std::size_t>(e) + 1] - begin);
        local = element_local_.data() + begin;
        sign = element_sign_.data() + begin;
    }

    void DofSpace::reduce_or(std::vector<int> &marker) const
    {
        for (unsigned int v = 0; v < n_vars(); ++v)
        {
            int *data = marker.data() + offsets_[v];
            spaces_[v]->GroupComm().Reduce<int>(data, mfem::GroupCommunicator::BitOR<int>);
            spaces_[v]->GroupComm().Bcast(data);
        }
    }

    std::unique_ptr<mfem::HypreParMatrix> DofSpace::assemble(mfem::SparseMatrix &A) const
    {
        // The L matrices of the ranks as one block-diagonal parallel matrix, then P_L^T A P_L (MFEM's ParallelAssemble).
        mfem::HypreParMatrix local(mesh_->comm(), n_local_global_, local_partition(), &A);
        return std::unique_ptr<mfem::HypreParMatrix>(mfem::RAP(&local, prolongation_.get()));
    }

    void DofSpace::gather(const mfem::Vector &x, std::vector<double> &local) const
    {
        mfem::Vector l(static_cast<int>(global_.size()));
        prolongation_->Mult(x, l);
        local.assign(l.GetData(), l.GetData() + l.Size());
    }

    void DofSpace::element_values(int e, const std::vector<double> &local, std::vector<double> &values) const
    {
        const int *dofs = nullptr;
        const double *sign = nullptr;
        unsigned int n = 0;
        element_dofs(e, dofs, sign, n);
        values.resize(n);
        for (unsigned int i = 0; i < n; ++i)
            values[i] = sign[i] * local[static_cast<std::size_t>(dofs[i])];
    }
} // namespace pytec3d::fem
