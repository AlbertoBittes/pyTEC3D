/**
 * @file assembler.h
 * @brief Global FE matrix and RHS assembler for the DOF solution system.
 *
 * AssemblerLinear reads per-element matrices and vectors previously written by
 * ElementBuilder, scales and combines them according to the chosen time-integration
 * scheme (generalised-theta method), applies constraint rows, and scatter-adds
 * the results into the global libMesh system matrix and RHS.
 */
#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <memory>

#include <libmesh/system.h>
#include <pybind11/pybind11.h>

#include "fe_parent_object.h"
#include "apec_m_visibility.h"

class ElementBuilder;
class ApecMEquationSystems;
class ApecMLinearSystem;

/**
 * @brief Assembles the global stiffness matrix and right-hand-side vector.
 *
 * Implements the generalised-theta time-integration scheme:
 *   - Reads pre-computed element matrices (Ke, Ce) and base RHS vectors from
 *     binary @c .emat cache files.
 *   - Scales the capacitance matrix and adds transient contributions to the RHS.
 *   - Constrains element matrices via the DofMap before scatter-adding into the
 *     global sparse system.
 */
class APEC_M_HIDDEN AssemblerLinear
    : public libMesh::System::Assembly,
      public FEParentObjectBase<ApecMEquationSystems>
{
public:
    using FEParentObjectBase<ApecMEquationSystems>::FEParentObjectBase;

    /** @brief Constructs the assembler and associates it with the given equation systems. */
    AssemblerLinear(ApecMEquationSystems *equation_systems = nullptr);

    /** @brief Assembles the global system matrix and RHS.  Called by libMesh::LinearImplicitSystem::solve(). */
    void assemble() override;

    /**
     * @brief Returns true if the global system matrix must be refactored this step.
     *
     * Triggers refactorisation when boundary conditions change, or when the
     * time-step size (dt) or integration parameter (theta) differ from the
     * previous step.
     */
    bool check_update_system_matrix();

    /**
     * @brief Returns true when element matrices must be recomputed.
     *
     * Currently always returns false because element matrices are
     * recomputed only when the @c .emat cache files are absent.
     */
    bool check_update_element_matrices();

    /** @brief Computes and caches element matrices and the base RHS to @c .emat files. */
    void write_element_matrices_and_rhs();

    /** @brief Computes the nodal-solution (B, H) fields and writes them to the Nodal Solution system. */
    void write_nodal_solution();

    /** @brief Computes the element-solution (Be, He, Jt) fields and writes them to the Element Solution system. */
    void write_element_solution();
};

#endif // ASSEMBLER_H
