#ifndef APEC_M_EQUATION_SYSTEMS_H
#define APEC_M_EQUATION_SYSTEMS_H

#include <memory>
#include <tuple>
#include <vector>
#include <libmesh/equation_systems.h>
#include <libmesh/numeric_vector.h>
#include <pybind11/pybind11.h>

#include "fe_parent_object.h"
#include "solver_configuration.h"

using namespace libMesh;

class ElementTypesFactory;
class ApecMSystemBase;
class ApecMSolver;

class ApecMEquationSystems
    : public EquationSystems,
      public FEParentObjectBase<ApecMSolver>

{
public:
    void solve_equation_systems(const pybind11::dict &time_step);
    void configure_systems();
    void initialize_systems();
    void post_initialize_systems();
    void initialize_physics();
    void zero_all_systems();

    // Variadic template constructor: forwards all args to EquationSystems, context view is last and optional
    template <typename... Args>
    ApecMEquationSystems(Args &&...args)
        : EquationSystems(std::forward<Args>(args)...)
    {
    }
};

#endif // APEC_M_EQUATION_SYSTEMS_H
