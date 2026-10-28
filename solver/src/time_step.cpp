/**
 * @file time_step.cpp
 * @brief TimeScheme::make(): the coefficients of the theta method and of the Newmark scheme (time_step.h).
 */
#include "time_step.h"

TimeScheme TimeScheme::make(const TimeStep &step, unsigned int time_order)
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
