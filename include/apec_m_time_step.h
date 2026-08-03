/**
 * @file apec_m_time_step.h
 * @brief Shared time-step state for APEC-M equation systems.
 */
#ifndef APEC_M_TIME_STEP_H
#define APEC_M_TIME_STEP_H

#include <libmesh/libmesh_common.h>

struct ApecMTimeStep
{
    libMesh::Real time = 0.0;
    libMesh::Real theta = 1.0;
    bool transient_effects = false;
    unsigned int time_step_idx = 0;
    libMesh::Real dt = 0.0;
    bool updated = true;
};

#endif // APEC_M_TIME_STEP_H