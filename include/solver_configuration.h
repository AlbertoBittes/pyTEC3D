#ifndef SOLVER_CONFIGURATION_H
#define SOLVER_CONFIGURATION_H

#include "libmesh/petsc_linear_solver.h"
#include "libmesh/solver_configuration.h"
#include "../include/fe_data_object.h"
#include "../include/load_step.h"

#include <string>
#include <unordered_map>
#include <variant>
#include <memory>
#include <stdexcept>

using namespace libMesh;

class PetscSolverConfiguration : public SolverConfiguration, public SettingsFEDataObject
{
public:
    // Default constructor
    PetscSolverConfiguration(std::shared_ptr<LoadStep> parent = nullptr,
                             std::shared_ptr<SettingsMap> map = nullptr)
        : SettingsFEDataObject(parent, map),
          _petsc_linear_solver(nullptr) {}

    // Constructor with PetscLinearSolver
    explicit PetscSolverConfiguration(PetscLinearSolver<Number> &petsc_linear_solver,
                                      std::shared_ptr<LoadStep> parent = nullptr,
                                      std::shared_ptr<SettingsMap> map = nullptr)
        : SettingsFEDataObject(parent, map),
          _petsc_linear_solver(&petsc_linear_solver) {}

    // Getters and setters for solver settings
    void configure_solver() override;

private:
    PetscLinearSolver<Number> *_petsc_linear_solver; // Pointer to the PetscLinearSolver
};

#endif // SOLVER_CONFIGURATION_H