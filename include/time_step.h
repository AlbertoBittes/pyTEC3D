/**
 * @file time_step.h
 * @brief Time-step state of the equation systems and the time scheme it gives.
 *
 * The DOF system solves M a + C v + K u = F at the end of the step, u, v, a the solution and its first and second
 * time derivatives. The scheme writes v and a as linear functions of the unknown u and of the state at the
 * start of the step (u_n, v_n, a_n):
 *
 *   v = c_v u + (hu_v u_n + hv_v v_n + ha_v a_n),   a = c_a u + (hu_a u_n + hv_a v_n + ha_a a_n)
 *
 * so the step is (K + c_v C + c_a M) u = F - C h_v - M h_a, and v, a follow from u after the solve.
 *  - First-order formulations (time_order 1), theta method: v = (u - u_n) / (theta dt) - (1 - theta) / theta v_n.
 *  - Second-order formulations (time_order 2), Newmark (beta, gamma):
 *      a = (u - u_n - dt v_n) / (beta dt^2) - (1 / (2 beta) - 1) a_n,   v = v_n + dt ((1 - gamma) a_n + gamma a).
 *    gamma = 1/2, beta = 1/4 is the unconditionally stable average acceleration scheme.
 *  - A static step (no transient effects, or time_order 0): v = a = 0.
 * The time order of the system is the highest of its element types.
 */
#ifndef pyTEC3D_TIME_STEP_H
#define pyTEC3D_TIME_STEP_H

#include <libmesh/libmesh_common.h>

struct TimeStep
{
    libMesh::Real time = 0.0;
    libMesh::Real theta = 1.0;  ///< Theta method (first order).
    libMesh::Real beta = 0.25;  ///< Newmark beta (second order).
    libMesh::Real gamma = 0.5;  ///< Newmark gamma (second order).
    bool transient_effects = false;
    unsigned int time_step_idx = 0;
    libMesh::Real dt = 0.0;
    bool updated = true;
};

/** @brief The coefficients of v and a in u and in the state at the start of the step (see the file comment). */
struct TimeScheme
{
    libMesh::Real c_v = 0.0, hu_v = 0.0, hv_v = 0.0, ha_v = 0.0;
    libMesh::Real c_a = 0.0, hu_a = 0.0, hv_a = 0.0, ha_a = 0.0;

    /** @brief True when the step is not static: C and M enter the system. */
    bool transient() const { return c_v != 0.0 || c_a != 0.0; }

    static TimeScheme make(const TimeStep &step, unsigned int time_order)
    {
        TimeScheme s;
        const libMesh::Real dt = step.dt;
        if (!step.transient_effects || time_order == 0)
            return s;
        if (time_order == 1)
        {
            const libMesh::Real theta = step.theta;
            s.c_v = 1.0 / (theta * dt);
            s.hu_v = -1.0 / (theta * dt);
            s.hv_v = -(1.0 - theta) / theta;
            return s;
        }
        const libMesh::Real beta = step.beta;
        const libMesh::Real gamma = step.gamma;
        s.c_a = 1.0 / (beta * dt * dt);
        s.hu_a = -1.0 / (beta * dt * dt);
        s.hv_a = -1.0 / (beta * dt);
        s.ha_a = -(1.0 / (2.0 * beta) - 1.0);
        s.c_v = gamma / (beta * dt);
        s.hu_v = -gamma / (beta * dt);
        s.hv_v = 1.0 - gamma / beta;
        s.ha_v = dt * (1.0 - gamma / (2.0 * beta));
        return s;
    }
};

#endif // pyTEC3D_TIME_STEP_H
