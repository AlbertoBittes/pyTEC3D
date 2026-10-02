/**
 * @file pytec3d_time_step.h
 * @brief Shared time-step state for pyTEC3D equation systems.
 */
#ifndef pyTEC3D_TIME_STEP_H
#define pyTEC3D_TIME_STEP_H

#include <libmesh/libmesh_common.h>

struct TimeStep
{
    libMesh::Real time = 0.0;
    libMesh::Real theta = 1.0;
    bool transient_effects = false;
    unsigned int time_step_idx = 0;
    libMesh::Real dt = 0.0;
    bool updated = true;
};

#endif // pyTEC3D_TIME_STEP_H