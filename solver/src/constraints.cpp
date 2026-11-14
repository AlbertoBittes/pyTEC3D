/**
 * @file constraints.cpp
 * @brief The constraints of the DOF system (constraints.h).
 */
#include "constraints.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "element_type.h"
#include "logging.h"
#include "model_element_types.h"
#include "registries.h"
#include "report_tables.h"

namespace pytec3d
{
    namespace
    {
        using la::Real;
        using la::RealTensorValue;
        using la::RealVectorValue;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort(message);
        }

        /** @brief True when the rotation @p R maps the xy plane onto itself (a rotation about z): the only local frame in
         * which a planar field (x, y only) can be given. */
        bool is_rotation_about_z(const RealTensorValue &R)
        {
            const Real tol = 1e-12;
            return std::abs(R(0, 2)) < tol && std::abs(R(1, 2)) < tol && std::abs(R(2, 0)) < tol && std::abs(R(2, 1)) < tol;
        }

        RealTensorValue rotation_of(const ConstraintInputs &in, std::int64_t point)
        {
            if (!in.node_rotations || in.node_rotations->empty())
                return RealTensorValue(1, 0, 0, 0, 1, 0, 0, 0, 1);
            return (*in.node_rotations)[static_cast<std::size_t>(point)];
        }

        /** @brief The constraint array of a variable: "_" + the field label of its grid arrays. */
        std::string constraint_array_name(const fem::Variable &variable)
        {
            return "_" + std::string(registry::info(variable.id).api_label);
        }

        /** @brief Component of the constraint array of an edge field holding its value (its midside point). */
        constexpr unsigned int kEdgeConstraintComponent = 2;

        /** @brief Rows being built on this rank, for the owned DOFs. */
        using RowMap = std::map<PetscInt, ConstraintRow>;

        void sort_row(ConstraintRow &row)
        {
            std::sort(row.masters.begin(), row.masters.end());
        }

        /** @brief The global DOF of (variable, point, component) on this rank, or -1. */
        PetscInt global_dof(const fem::DofSpace &dofs, unsigned int v, std::int64_t point, unsigned int c)
        {
            const long long l = dofs.point_dof(v, point, c);
            return l < 0 ? -1 : dofs.global(static_cast<std::size_t>(l));
        }

        bool owned_dof(const fem::DofSpace &dofs, unsigned int v, std::int64_t point, unsigned int c, PetscInt &g)
        {
            const long long l = dofs.point_dof(v, point, c);
            if (l < 0 || !dofs.owned(static_cast<std::size_t>(l)))
                return false;
            g = dofs.global(static_cast<std::size_t>(l));
            return true;
        }

        // ---- Dirichlet values -------------------------------------------------------------------------------------

        RowMap point_constraints(const ConstraintInputs &in)
        {
            const fem::DofSpace &dofs = *in.dofs;
            RowMap rows;
            const Real tol = std::numeric_limits<Real>::epsilon();
            for (unsigned int v = 0; v < dofs.n_vars(); ++v)
            {
                const fem::Variable &variable = dofs.variable(v);
                const RealArray *array = in.arrays->point_field(constraint_array_name(variable));
                if (!array)
                    continue;

                if (!variable.vector_lagrange())
                {
                    // A scalar, one component of a vector field, or an edge field (component 2 of its midside point).
                    const unsigned int component = variable.edge() ? kEdgeConstraintComponent
                                                                   : (variable.n_components == 1 ? variable.components[0] : 0u);
                    if (component >= array->n_components())
                        continue;
                    for (const auto &[point, dof] : dofs.point_dofs(v))
                    {
                        (void)dof;
                        PetscInt g = -1;
                        if (!owned_dof(dofs, v, point, 0, g))
                            continue;
                        const Real value = (*array)(static_cast<std::size_t>(point), component);
                        if (std::isnan(value))
                            continue;
                        ConstraintRow row;
                        row.value = value;
                        rows.emplace(g, row);
                    }
                    continue;
                }

                // A vector field: its components in the node frame. A planar field (2D elements) has x, y only.
                const bool planar = variable.vdim < 3;
                for (const auto &[point, dof] : dofs.point_dofs(v))
                {
                    (void)dof;
                    std::array<PetscInt, 3> g = {-1, -1, -1};
                    std::array<bool, 3> owned = {false, false, false};
                    std::array<Real, 3> values = {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN(),
                                                  std::numeric_limits<Real>::quiet_NaN()};
                    for (unsigned int c = 0; c < static_cast<unsigned int>(variable.vdim); ++c)
                    {
                        const long long l = dofs.point_dof(v, point, c);
                        if (l < 0 || c >= array->n_components())
                            continue;
                        g[c] = dofs.global(static_cast<std::size_t>(l));
                        owned[c] = dofs.owned(static_cast<std::size_t>(l));
                        values[c] = (*array)(static_cast<std::size_t>(point), c);
                    }
                    if (!owned[0] && !owned[1] && !owned[2])
                        continue;
                    if (std::isnan(values[0]) && std::isnan(values[1]) && std::isnan(values[2]))
                        continue;

                    const RealTensorValue local_to_global = rotation_of(in, point);
                    const RealTensorValue global_to_local = local_to_global.transpose();
                    if (planar && !is_rotation_about_z(local_to_global))
                    {
                        std::ostringstream msg;
                        msg << "Point constraint: node " << point << " has a local frame that is not a rotation about z, "
                            << "but the vector field '" << variable.label << "' has only the in-plane components x, y there "
                            << "(2D elements).";
                        fail(msg.str());
                    }

                    std::unordered_set<PetscInt> constrained;
                    const auto add = [&](unsigned int c, ConstraintRow row)
                    {
                        if (!owned[c])
                        {
                            constrained.insert(g[c]);
                            return;
                        }
                        sort_row(row);
                        if (rows.emplace(g[c], row).second)
                            constrained.insert(g[c]);
                    };

                    // All the components given: the global values directly.
                    const unsigned int n_field = planar ? 2 : 3;
                    bool full = true;
                    for (unsigned int c = 0; c < n_field && full; ++c)
                        full = g[c] >= 0 && !std::isnan(values[c]);
                    if (full)
                    {
                        const RealVectorValue global = local_to_global * RealVectorValue(values[0], values[1], planar ? 0.0 : values[2]);
                        for (unsigned int c = 0; c < n_field; ++c)
                        {
                            ConstraintRow row;
                            row.value = global(c);
                            add(c, row);
                        }
                        continue;
                    }

                    // Some components: rows on the global DOFs, from independent equations of the local ones.
                    std::array<unsigned int, 3> given = {0, 0, 0};
                    unsigned int n_given = 0;
                    for (unsigned int c = 0; c < 3; ++c)
                        if (!std::isnan(values[c]))
                            given[n_given++] = c;

                    bool handled = false;
                    if (n_given == 2)
                    {
                        const unsigned int eq0 = given[0], eq1 = given[1];
                        std::array<unsigned int, 3> candidates = {0, 0, 0};
                        unsigned int n_candidates = 0;
                        for (unsigned int j = 0; j < 3; ++j)
                            if (g[j] >= 0 && !constrained.count(g[j]))
                                candidates[n_candidates++] = j;
                        if (n_candidates >= 2)
                        {
                            Real best = 0.0;
                            unsigned int s0 = 0, s1 = 0;
                            for (unsigned int a = 0; a < n_candidates; ++a)
                                for (unsigned int b = a + 1; b < n_candidates; ++b)
                                {
                                    const unsigned int j0 = candidates[a], j1 = candidates[b];
                                    const Real det = global_to_local(eq0, j0) * global_to_local(eq1, j1) -
                                                     global_to_local(eq0, j1) * global_to_local(eq1, j0);
                                    if (std::abs(det) > best)
                                    {
                                        best = std::abs(det);
                                        s0 = j0;
                                        s1 = j1;
                                    }
                                }
                            if (best >= tol)
                            {
                                const Real a00 = global_to_local(eq0, s0), a01 = global_to_local(eq0, s1);
                                const Real a10 = global_to_local(eq1, s0), a11 = global_to_local(eq1, s1);
                                const Real det = a00 * a11 - a01 * a10;
                                const Real inv[2][2] = {{a11 / det, -a01 / det}, {-a10 / det, a00 / det}};
                                const unsigned int slaves[2] = {s0, s1};
                                for (unsigned int k = 0; k < 2; ++k)
                                {
                                    ConstraintRow row;
                                    for (unsigned int j = 0; j < 3; ++j)
                                    {
                                        if (j == s0 || j == s1 || g[j] < 0)
                                            continue;
                                        row.masters.emplace_back(
                                            g[j], -(inv[k][0] * global_to_local(eq0, j) + inv[k][1] * global_to_local(eq1, j)));
                                    }
                                    row.value = inv[k][0] * values[eq0] + inv[k][1] * values[eq1];
                                    add(slaves[k], row);
                                }
                                handled = true;
                            }
                        }
                    }
                    if (handled)
                        continue;

                    for (unsigned int i = 0; i < 3; ++i)
                    {
                        if (std::isnan(values[i]))
                            continue;
                        unsigned int pivot = i;
                        Real pivot_coeff = global_to_local(i, pivot);
                        if (g[pivot] < 0 || std::abs(pivot_coeff) < tol)
                        {
                            Real best = 0.0;
                            for (unsigned int j = 0; j < 3; ++j)
                                if (g[j] >= 0 && std::abs(global_to_local(i, j)) > best)
                                {
                                    best = std::abs(global_to_local(i, j));
                                    pivot_coeff = global_to_local(i, j);
                                    pivot = j;
                                }
                        }
                        if (g[pivot] < 0 || std::abs(pivot_coeff) < tol)
                            continue;
                        ConstraintRow row;
                        for (unsigned int j = 0; j < 3; ++j)
                            if (j != pivot && g[j] >= 0)
                                row.masters.emplace_back(g[j], -global_to_local(i, j) / pivot_coeff);
                        row.value = values[i] / pivot_coeff;
                        add(pivot, row);
                    }
                }
            }
            return rows;
        }

        // ---- DOF couplings ----------------------------------------------------------------------------------------

        /** @brief The variables a coupling of @p label couples: the variable of that label, or those of the field. */
        std::vector<unsigned int> coupled_variables(const fem::DofSpace &dofs, const std::string &label)
        {
            const int exact = dofs.variable_number(label);
            if (exact >= 0)
                return {static_cast<unsigned int>(exact)};
            std::vector<std::pair<unsigned int, unsigned int>> matches; // (component, variable)
            for (unsigned int v = 0; v < dofs.n_vars(); ++v)
                if (registry::info(dofs.variable(v).id).api_label == label)
                    matches.emplace_back(dofs.variable(v).components[0], v);
            if (matches.empty())
                fail("Dof couplings: '" + label + "' is not a variable of the DOF system.");
            std::sort(matches.begin(), matches.end());
            std::vector<unsigned int> out;
            for (const auto &m : matches)
                out.push_back(m.second);
            return out;
        }

        /**
         * @brief Whether the edge DOFs of two coupled midside points have the same orientation: the CP labels of the end
         * points (in the orientation of the DOFs, lexicographic) must match directly or reversed.
         */
        bool same_edge_orientation(const fem::FEMesh &mesh, const std::unordered_map<std::int64_t, int> &cp_of_point,
                                   int cp_id, std::int64_t master, std::int64_t slave)
        {
            std::int64_t m0, m1, s0, s1;
            if (!mesh.edge_corners(master, m0, m1) || !mesh.edge_corners(slave, s0, s1))
                return true;
            const auto greater = [&mesh](std::int64_t a, std::int64_t b)
            {
                const double *p = mesh.point(a), *q = mesh.point(b);
                for (int i = 0; i < 3; ++i)
                {
                    if (p[i] > q[i])
                        return true;
                    if (p[i] < q[i])
                        return false;
                }
                return false;
            };
            // The corners ordered larger first: the order in which the CP labels of the two edges are compared.
            const auto order = [&](std::int64_t a, std::int64_t b) { return greater(a, b) ? std::array<std::int64_t, 2>{a, b}
                                                                                           : std::array<std::int64_t, 2>{b, a}; };
            const auto mc = order(m0, m1), sc = order(s0, s1);
            const auto cp_of = [&](std::int64_t p)
            {
                const auto it = cp_of_point.find(p);
                return it == cp_of_point.end() ? 0 : it->second;
            };
            const std::array<int, 3> mt = {cp_of(mc[0]), cp_of(master), cp_of(mc[1])};
            const std::array<int, 3> st = {cp_of(sc[0]), cp_of(slave), cp_of(sc[1])};
            std::ostringstream where;
            where << "CP_ID=" << cp_id << " (master midside node " << master << ", slave midside node " << slave << ")";
            if ((mt[0] != 0) + (mt[2] != 0) == 0 && (st[0] != 0) + (st[2] != 0) == 0)
                fail("Coupling constraint: unable to infer the edge orientation because no corner CP label was found for "
                     "either edge in " + where.str() + ".");
            std::array<int, 3> ms = mt, ss = st;
            std::sort(ms.begin(), ms.end());
            std::sort(ss.begin(), ss.end());
            if (ms != ss)
                fail("Coupling constraint: incompatible edge CP token triplets for " + where.str() + ".");
            const bool direct = mt[0] == st[0] && mt[1] == st[1] && mt[2] == st[2];
            const bool reversed = mt[0] == st[2] && mt[1] == st[1] && mt[2] == st[0];
            const bool direct_ok = !((mt[0] == 0 && st[0] == 0 && mc[0] != sc[0]) || (mt[2] == 0 && st[2] == 0 && mc[1] != sc[1]));
            const bool reversed_ok = !((mt[0] == 0 && st[2] == 0 && mc[0] != sc[1]) || (mt[2] == 0 && st[0] == 0 && mc[1] != sc[0]));
            if (!(direct && direct_ok) && !(reversed && reversed_ok))
                fail("Coupling constraint: unable to infer the edge orientation from CP tokens/coordinates for " + where.str() +
                     ". If CP is missing on matched corners, the missing-CP corners must be the same node.");
            return direct && direct_ok;
        }

        RowMap coupling_constraints(const ConstraintInputs &in)
        {
            RowMap rows;
            const fem::DofSpace &dofs = *in.dofs;
            if (!in.couplings || in.couplings->empty() || !in.pairs || in.pairs->empty())
                return rows;
            const MPI_Comm comm = in.mesh->comm();

            // The global DOFs of the coupled points, from the ranks that have them: (point, variable, component, dof).
            std::set<std::int64_t> points;
            for (const DofCoupling &coupling : *in.couplings)
            {
                const auto it = in.pairs->find(coupling.cp_id);
                if (it == in.pairs->end())
                    continue;
                for (const auto &[a, b] : it->second)
                {
                    points.insert(a);
                    points.insert(b);
                }
            }
            std::vector<long long> known;
            for (const std::int64_t p : points)
                for (unsigned int v = 0; v < dofs.n_vars(); ++v)
                    for (unsigned int c = 0; c < static_cast<unsigned int>(dofs.variable(v).vdim); ++c)
                    {
                        const PetscInt g = global_dof(dofs, v, p, c);
                        if (g >= 0)
                            known.insert(known.end(), {static_cast<long long>(p), v, c, static_cast<long long>(g)});
                    }
            int n_ranks = 1;
            MPI_Comm_size(comm, &n_ranks);
            int count = static_cast<int>(known.size());
            std::vector<int> counts(static_cast<std::size_t>(n_ranks)), displs(static_cast<std::size_t>(n_ranks), 0);
            MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
            for (int r = 1; r < n_ranks; ++r)
                displs[static_cast<std::size_t>(r)] = displs[static_cast<std::size_t>(r) - 1] + counts[static_cast<std::size_t>(r) - 1];
            std::vector<long long> all(static_cast<std::size_t>(displs.back() + counts.back()));
            MPI_Allgatherv(known.data(), count, MPI_LONG_LONG, all.data(), counts.data(), displs.data(), MPI_LONG_LONG, comm);
            std::map<std::tuple<long long, long long, long long>, PetscInt> dof_of;
            for (std::size_t i = 0; i + 3 < all.size(); i += 4)
                dof_of.emplace(std::make_tuple(all[i], all[i + 1], all[i + 2]), static_cast<PetscInt>(all[i + 3]));
            const auto lookup = [&](std::int64_t p, unsigned int v, unsigned int c) -> PetscInt
            {
                const auto it = dof_of.find(std::make_tuple(static_cast<long long>(p), static_cast<long long>(v),
                                                            static_cast<long long>(c)));
                return it == dof_of.end() ? -1 : it->second;
            };
            const auto owned = [&](PetscInt g) { return g >= dofs.first_owned() && g < dofs.first_owned() + dofs.n_owned(); };

            // The CP label of each point of a pair (the first set, by id), for the orientation of edge DOFs.
            std::unordered_map<std::int64_t, int> cp_of_point;
            for (const auto &[cp, pairs] : *in.pairs)
                for (const auto &[a, b] : pairs)
                {
                    cp_of_point.emplace(a, cp);
                    cp_of_point.emplace(b, cp);
                }

            for (const DofCoupling &coupling : *in.couplings)
            {
                if (coupling.variable.empty())
                    continue;
                const auto it = in.pairs->find(coupling.cp_id);
                if (it == in.pairs->end() || it->second.empty())
                    continue;
                const std::vector<unsigned int> variables = coupled_variables(dofs, coupling.variable);
                for (const auto &[master, slave] : it->second)
                    for (const unsigned int v : variables)
                    {
                        const fem::Variable &variable = dofs.variable(v);
                        if (variable.vector_lagrange())
                        {
                            const unsigned int n = static_cast<unsigned int>(variable.vdim);
                            std::array<PetscInt, 3> first = {-1, -1, -1}, second = {-1, -1, -1};
                            bool complete = true;
                            for (unsigned int c = 0; c < n; ++c)
                            {
                                first[c] = lookup(master, v, c);
                                second[c] = lookup(slave, v, c);
                                complete = complete && first[c] >= 0 && second[c] >= 0;
                            }
                            if (!complete || !owned(second[0]))
                                continue;
                            const RealTensorValue transform = rotation_of(in, slave) * rotation_of(in, master).transpose();
                            if (n == 2 && !is_rotation_about_z(transform))
                            {
                                std::ostringstream msg;
                                msg << "Coupling constraint: the local frames of master node " << master << " and slave node "
                                    << slave << " differ by more than a rotation about z, but the vector field '"
                                    << variable.label << "' has only the in-plane components x, y there (2D elements).";
                                fail(msg.str());
                            }
                            for (unsigned int i = 0; i < n; ++i)
                            {
                                ConstraintRow row;
                                for (unsigned int j = 0; j < n; ++j)
                                    row.masters.emplace_back(first[j], transform(i, j));
                                sort_row(row);
                                if (owned(second[i]))
                                    rows.emplace(second[i], row);
                            }
                            continue;
                        }
                        const PetscInt first = lookup(master, v, 0);
                        const PetscInt second = lookup(slave, v, 0);
                        if (first < 0 || second < 0 || !owned(second))
                            continue;
                        const double sign =
                            variable.edge() && !same_edge_orientation(*in.mesh, cp_of_point, coupling.cp_id, master, slave) ? -1.0 : 1.0;
                        ConstraintRow row;
                        row.masters.emplace_back(first, sign);
                        rows.emplace(second, row);
                    }
            }
            return rows;
        }

        // ---- the pins of static steps -----------------------------------------------------------------------------

        RowMap static_constraints(const ConstraintInputs &in)
        {
            RowMap rows;
            if (!in.static_step)
                return rows;
            const fem::DofSpace &dofs = *in.dofs;
            const fem::FEMesh &mesh = *in.mesh;
            fem::FECatalog const &catalog = dofs.catalog();
            std::vector<int> with(dofs.n_local(), 0), without(dofs.n_local(), 0);
            for (int e = 0; e < mesh.n_local_elements(); ++e)
            {
                const ModelElementTypes::Entry &entry = in.element_types->get(mesh.subdomain(e));
                const ElementType &type = *entry.type;
                const int *local = nullptr;
                const double *sign = nullptr;
                unsigned int n = 0;
                dofs.element_dofs(e, local, sign, n);
                const MeshElement element{&mesh, &catalog, nullptr, nullptr, e};
                const auto indices = type.element_indices(dofs, element, FieldKind::Dof, entry.domain);
                for (std::size_t field = 0; field < type.spec().fields.size(); ++field)
                {
                    const FieldSpec &spec = type.field(field);
                    if (spec.kind != FieldKind::Dof || !ElementType::active(spec, entry.domain))
                        continue;
                    auto &target = type.has_stiffness(field, entry.domain) ? with : without;
                    const ElementType::FieldRange &range = indices.field[field];
                    for (unsigned int i = 0; i < range.count; ++i)
                        target[static_cast<std::size_t>(local[range[i]])] = 1;
                }
            }
            dofs.reduce_or(with);
            dofs.reduce_or(without);
            for (std::size_t l = 0; l < dofs.n_local(); ++l)
                if (without[l] && !with[l] && dofs.owned(l))
                    rows.emplace(dofs.global(l), ConstraintRow{});
            return rows;
        }

        /** @brief All the rows of all the ranks, on every rank. */
        ConstraintSet gather(const RowMap &rows, MPI_Comm comm)
        {
            // Serialized: dof, value, n masters, then (master, coefficient) pairs; doubles hold the indices exactly.
            std::vector<double> mine;
            for (const auto &[dof, row] : rows)
            {
                mine.push_back(static_cast<double>(dof));
                mine.push_back(row.value);
                mine.push_back(static_cast<double>(row.masters.size()));
                for (const auto &[m, c] : row.masters)
                {
                    mine.push_back(static_cast<double>(m));
                    mine.push_back(c);
                }
            }
            int n_ranks = 1;
            MPI_Comm_size(comm, &n_ranks);
            int count = static_cast<int>(mine.size());
            std::vector<int> counts(static_cast<std::size_t>(n_ranks)), displs(static_cast<std::size_t>(n_ranks), 0);
            MPI_Allgather(&count, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
            for (int r = 1; r < n_ranks; ++r)
                displs[static_cast<std::size_t>(r)] = displs[static_cast<std::size_t>(r) - 1] + counts[static_cast<std::size_t>(r) - 1];
            std::vector<double> all(static_cast<std::size_t>(displs.back() + counts.back()));
            MPI_Allgatherv(mine.data(), count, MPI_DOUBLE, all.data(), counts.data(), displs.data(), MPI_DOUBLE, comm);
            ConstraintSet set;
            for (std::size_t i = 0; i < all.size();)
            {
                const auto dof = static_cast<PetscInt>(all[i]);
                ConstraintRow row;
                row.value = all[i + 1];
                const auto n = static_cast<std::size_t>(all[i + 2]);
                i += 3;
                for (std::size_t k = 0; k < n; ++k, i += 2)
                    row.masters.emplace_back(static_cast<PetscInt>(all[i]), all[i + 1]);
                set.emplace(dof, std::move(row));
            }
            return set;
        }
    } // namespace

    CouplingPairs coupling_pairs(const ModelArrays &arrays, std::size_t n_points)
    {
        std::map<int, std::vector<std::int64_t>> sets;
        for (std::size_t p = 0; p < n_points; ++p)
            for (const auto &[name, array] : arrays.couplings)
            {
                (void)name;
                const int cp = array(p);
                if (cp > 0)
                    sets[cp].push_back(static_cast<std::int64_t>(p));
            }
        CouplingPairs pairs;
        for (const auto &[cp, points] : sets)
        {
            const std::int64_t master = *std::min_element(points.begin(), points.end());
            auto &list = pairs[cp];
            for (const std::int64_t p : points)
                if (p != master)
                    list.emplace_back(master, p);
        }
        return pairs;
    }

    ConstraintSet build_constraints(const ConstraintInputs &in)
    {
        logging::info("Computing nodal DOF constraints");
        const fem::DofSpace &dofs = *in.dofs;
        const MPI_Comm comm = in.mesh->comm();

        const RowMap couplings = coupling_constraints(in);
        const RowMap points = point_constraints(in);
        const RowMap pins = static_constraints(in);

        RowMap rows = couplings;
        std::size_t overridden = 0;
        for (const auto &[dof, row] : points)
        {
            if (couplings.count(dof))
            {
                ++overridden;
                logging::warn("Point constraint overrides coupling constraint (dof={})", dof);
            }
            rows[dof] = row;
        }
        for (const auto &[dof, row] : pins)
            if (!points.count(dof))
                rows[dof] = row; // A Dirichlet value wins over the pin of a DOF without stiffness.
        (void)overridden;

        // The summary: counts per source and per variable, all ranks.
        std::vector<std::size_t> per_variable(dofs.n_vars(), 0);
        {
            // The variable of an owned global DOF: the variable blocks of this rank are consecutive.
            std::vector<PetscInt> first(dofs.n_vars() + 1, dofs.first_owned());
            std::vector<PetscInt> counts(dofs.n_vars(), 0);
            for (std::size_t l = 0; l < dofs.n_local(); ++l)
                if (dofs.owned(l))
                    ++counts[dofs.variable_of(l)];
            for (unsigned int v = 0; v < dofs.n_vars(); ++v)
                first[v + 1] = first[v] + counts[v];
            for (const auto &[dof, row] : rows)
            {
                (void)row;
                const auto it = std::upper_bound(first.begin(), first.end(), dof);
                const auto v = static_cast<std::size_t>(std::distance(first.begin(), it)) - 1;
                if (v < per_variable.size())
                    ++per_variable[v];
            }
        }
        std::array<unsigned long long, 4> totals = {pins.size(), points.size(), couplings.size(), rows.size()};
        MPI_Allreduce(MPI_IN_PLACE, totals.data(), 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
        std::vector<unsigned long long> per_variable_all(per_variable.begin(), per_variable.end());
        MPI_Allreduce(MPI_IN_PLACE, per_variable_all.data(), static_cast<int>(per_variable_all.size()), MPI_UNSIGNED_LONG_LONG,
                      MPI_SUM, comm);
        std::vector<report::ConstraintVariableStat> stats;
        for (unsigned int v = 0; v < dofs.n_vars(); ++v)
            if (per_variable_all[v] > 0)
                stats.push_back({dofs.variable(v).label, per_variable_all[v], per_variable_all[v]});
        std::sort(stats.begin(), stats.end(), [](const auto &a, const auto &b)
                  { return a.total_count != b.total_count ? a.total_count > b.total_count : a.variable_name < b.variable_name; });
        logging::table_root(report::format_constraints_summary_table(totals[0], totals[1], totals[2], totals[3], totals[3],
                                                                     stats, 0, 0));
        return gather(rows, comm);
    }
} // namespace pytec3d
