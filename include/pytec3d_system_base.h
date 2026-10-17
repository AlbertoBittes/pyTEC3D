/**
 * @file pytec3d_system_base.h
 * @brief pyTEC3DSystemBase: the interface and the names shared by the solution systems (pytec3d_system.h).
 *
 * Kept apart so that code needing only the system types or vector names does not include the system hierarchy.
 */
#pragma once

#include <map>
#include <string>
#include <unordered_map>

#include <libmesh/id_types.h>

class pyTEC3DEquationSystems;
struct TimeStep;

/// DOF variable label -> (DOF id -> the local element that injects its nodal load).
using NodalLoadInjectorElemByVarMap =
    std::unordered_map<std::string,
                       std::unordered_map<libMesh::dof_id_type, libMesh::dof_id_type>>;

class pyTEC3DSystemBase
{
public:
    enum SystemType
    {
        DOF_SOLUTION_SYSTEM = 0,
        NODAL_SOLUTION,
        ELEMENT_SOLUTION,
        NUM_SYSTEMS
    };

    virtual ~pyTEC3DSystemBase() = default;
    virtual bool is_nodal() const = 0;
    virtual void assemble_system() = 0;
    virtual pyTEC3DEquationSystems &get_pytec3d_equation_systems() = 0;
    virtual const pyTEC3DEquationSystems &get_pytec3d_equation_systems() const = 0;
    virtual const TimeStep *current_time_step() const = 0;
    virtual const TimeStep *previous_time_step() const = 0;
    virtual void solve() = 0;
    virtual void init_variables() = 0;
    virtual void post_init() = 0;

    static const std::map<SystemType, std::string> SYSTEM_INFO_MAP;

    /// DOF system vectors of the time scheme (time_step.h): first and second time derivatives of the solution at
    /// the end of the step, and at its start.
    static constexpr const char *kVelocity = "solution_time_derivative";
    static constexpr const char *kAcceleration = "solution_second_time_derivative";
    static constexpr const char *kOldVelocity = "old_solution_time_derivative";
    static constexpr const char *kOldAcceleration = "old_solution_second_time_derivative";
    /// Mean rate of the solution over the step, (u - u_n) / dt: the dt() of the nodal and element fields of
    /// first-order models (with theta < 1 the scheme velocity oscillates from step to step).
    static constexpr const char *kMeanVelocity = "solution_mean_time_derivative";
    static const std::string &system_name_from_type(SystemType system_type);
    static SystemType system_type_from_name(const std::string &system_name);
};
