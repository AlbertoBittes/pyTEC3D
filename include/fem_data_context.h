#ifndef FEM_DATA_CONTEXT_H
#define FEM_DATA_CONTEXT_H

#include <vector>
#include <memory>
#include <type_traits>
#include <typeinfo>
#include <stdexcept>
#include "fe_data_object.h"
#include "solver_configuration.h"
#include "node_bc_object.h"
#include "time_step.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "loads_object.h"

class FEMDataContext
{
public:
    FEMDataContext() = default;
    virtual ~FEMDataContext() = 0;

    // Add any FEDataObjectAbstract or subclass
    template <typename T>
    void set_object(std::shared_ptr<T> obj)
    {
        static_assert(std::is_base_of<FEDataObjectAbstract, T>::value, "T must derive from FEDataObjectAbstract");
        // Remove any existing object of this type
        remove_object<T>();
        data_objects.push_back(obj);
    }

    // Get a FEDataObjectAbstract or subclass by type
    template <typename T>
    std::shared_ptr<T> get_object() const
    {
        static_assert(std::is_base_of<FEDataObjectAbstract, T>::value, "T must derive from FEDataObjectAbstract");
        for (const auto &obj : data_objects)
        {
            auto casted = std::dynamic_pointer_cast<T>(obj);
            if (casted)
                return casted;
        }
        return nullptr;
    }

    // Remove object of a given type
    template <typename T>
    void remove_object()
    {
        data_objects.erase(
            std::remove_if(data_objects.begin(), data_objects.end(),
                           [](const std::shared_ptr<FEDataObjectAbstract> &obj)
                           {
                               return std::dynamic_pointer_cast<T>(obj) != nullptr;
                           }),
            data_objects.end());
    }

    // Access to all objects
    const std::vector<std::shared_ptr<FEDataObjectAbstract>> &get_all_objects() const
    {
        return data_objects;
    }

protected:
    std::vector<std::shared_ptr<FEDataObjectAbstract>> data_objects;
};

class SystemBuilderContext : public FEMDataContext
{
public:
    // Solver configuration
    void set_solver_config(std::shared_ptr<PetscSolverConfiguration> obj) { set_object(obj); }
    std::shared_ptr<PetscSolverConfiguration> get_solver_config() const { return get_object<PetscSolverConfiguration>(); }

    // Materials
    void set_materials(std::shared_ptr<MaterialsObject> obj) { set_object(obj); }
    std::shared_ptr<MaterialsObject> get_materials() const { return get_object<MaterialsObject>(); }
};

class SystemContext : public FEMDataContext
{
public:
    void set_node_bc(std::shared_ptr<NodeBCObject> obj) { set_object(obj); }
    std::shared_ptr<NodeBCObject> get_node_bc() const { return get_object<NodeBCObject>(); }

    void set_time_step(std::shared_ptr<TimeStep> obj) { set_object(obj); }
    std::shared_ptr<TimeStep> get_time_step() const { return get_object<TimeStep>(); }
};

class ElementBuilderContext : public FEMDataContext
{
public:
    void set_materials(std::shared_ptr<MaterialsObject> obj) { set_object(obj); }
    std::shared_ptr<MaterialsObject> get_materials() const { return get_object<MaterialsObject>(); }

    void set_real_constants(std::shared_ptr<RealConstantsObject> obj) { set_object(obj); }
    std::shared_ptr<RealConstantsObject> get_real_constants() const { return get_object<RealConstantsObject>(); }

    void set_loads(std::shared_ptr<LoadsObject> obj) { set_object(obj); }
    std::shared_ptr<LoadsObject> get_loads() const { return get_object<LoadsObject>(); }
};

#endif // FEM_DATA_CONTEXT_H