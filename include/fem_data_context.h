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
#include "session.h"

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

class ElementPhysicsContext : public FEMDataContext
{
public:
    ElementPhysicsContext(std::shared_ptr<ElementTypeFactory> element_types)
    {
        set_element_types(element_types);
    }
    void set_element_types(std::shared_ptr<ElementTypeFactory> obj) { set_object(obj); }
    std::shared_ptr<ElementTypeFactory> get_element_types() const
    {
        return get_object<ElementTypeFactory>();
    }
};

class ElementBuilderContext : public ElementPhysicsContext
{
public:
    ElementBuilderContext(
        std::shared_ptr<ElementTypeFactory> element_types,
        std::shared_ptr<MaterialsObject> materials,
        std::shared_ptr<RealConstantsObject> real_constants,
        std::shared_ptr<LoadsObject> loads,
        std::shared_ptr<Session> session)
        : ElementPhysicsContext(element_types)
    {
        set_materials(materials);
        set_real_constants(real_constants);
        set_loads(loads);
        set_session(session);
    }

    void set_materials(std::shared_ptr<MaterialsObject> obj) { set_object(obj); }
    std::shared_ptr<MaterialsObject> get_materials() const { return get_object<MaterialsObject>(); }

    void set_real_constants(std::shared_ptr<RealConstantsObject> obj) { set_object(obj); }
    std::shared_ptr<RealConstantsObject> get_real_constants() const { return get_object<RealConstantsObject>(); }

    void set_loads(std::shared_ptr<LoadsObject> obj) { set_object(obj); }
    std::shared_ptr<LoadsObject> get_loads() const { return get_object<LoadsObject>(); }

    void set_session(std::shared_ptr<Session> obj) { set_object(obj); }
    std::shared_ptr<Session> get_session() const
    {
        return get_object<Session>();
    }
};

class AssemblerContext : public FEMDataContext
{
public:
    AssemblerContext(std::shared_ptr<TimeStep> time_step, std::shared_ptr<Session> session)
    {
        set_time_step(time_step);
        set_session(session);
    }

    void set_time_step(std::shared_ptr<TimeStep> obj) { set_object(obj); }
    std::shared_ptr<TimeStep> get_time_step() const { return get_object<TimeStep>(); }

    void set_session(std::shared_ptr<Session> obj) { set_object(obj); }
    std::shared_ptr<Session> get_session() const
    {
        return get_object<Session>();
    }
};

class BoundaryConditionsContext : public FEMDataContext
{
public:
    BoundaryConditionsContext(std::shared_ptr<NodeBCObject> node_bc_object, std::shared_ptr<TimeStep> time_step)
    {
        set_node_bc_object(node_bc_object);
        set_time_step(time_step);
    }

    void set_node_bc_object(std::shared_ptr<NodeBCObject> obj) { set_object(obj); }
    std::shared_ptr<NodeBCObject> get_node_bc_object() const { return get_object<NodeBCObject>(); }

    void set_time_step(std::shared_ptr<TimeStep> obj) { set_object(obj); }
    std::shared_ptr<TimeStep> get_time_step() const { return get_object<TimeStep>(); }

    // void set_session(std::shared_ptr<Session> obj) { set_object(obj); }
    // std::shared_ptr<Session> get_session() const
    // {
    //     return get_object<Session>();
    // }
};

#endif // FEM_DATA_CONTEXT_H
// #ifndef FEM_DATA_CONTEXT_H
// #define FEM_DATA_CONTEXT_H

// #include <vector>
// #include <memory>
// #include <type_traits>
// #include <typeinfo>
// #include <typeindex>
// #include <stdexcept>
// #include <unordered_map>
// #include "fe_data_object.h"
// #include "solver_configuration.h"
// #include "node_bc_object.h"
// #include "time_step.h"
// #include "materials_object.h"
// #include "real_constants_object.h"
// #include "loads_object.h"
// #include "session.h"

// class FEMDataContext
// {
// public:
//     FEMDataContext() = default;
//     virtual ~FEMDataContext() = default;

//     // Set object by type
//     template <typename T>
//     void set_object(std::shared_ptr<T> obj)
//     {
//         static_assert(std::is_base_of<FEDataObjectAbstract, T>::value, "T must derive from FEDataObjectAbstract");
//         data_objects[typeid(T)] = obj;
//     }

//     // Get object by type
//     template <typename T>
//     std::shared_ptr<T> get_object() const
//     {
//         static_assert(std::is_base_of<FEDataObjectAbstract, T>::value, "T must derive from FEDataObjectAbstract");
//         auto it = data_objects.find(typeid(T));
//         if (it == data_objects.end())
//             throw std::runtime_error("Requested object not found in FEMDataContext");
//         return std::dynamic_pointer_cast<T>(it->second);
//     }

//     // Remove object by type
//     template <typename T>
//     void remove_object()
//     {
//         data_objects.erase(typeid(T));
//     }

//     // Access to all objects (for introspection, not type-safe)
//     const std::unordered_map<std::type_index, std::shared_ptr<FEDataObjectAbstract>> &get_all_objects() const
//     {
//         return data_objects;
//     }

// protected:
//     std::unordered_map<std::type_index, std::shared_ptr<FEDataObjectAbstract>> data_objects;
// };

// // Example of a "view" context for typed access
// template <typename... Types>
// struct TypedContextView
// {
//     std::tuple<Types &...> refs;

//     template <typename Repo1, typename... Repos>
//     static TypedContextView from(const Repo1 &r1, const Repos &...rs)
//     {
//         return {getFromRepos<Types>(r1, rs...)...};
//     }

//     template <typename T>
//     T &get() { return std::get<T &>(refs); }

// private:
//     template <typename T, typename Repo>
//     static T &tryGet(const Repo &repo)
//     {
//         return *repo.template get_object<T>();
//     }

//     template <typename T, typename Repo1, typename... Rest>
//     static T &getFromRepos(const Repo1 &r1, const Rest &...rest)
//     {
//         try
//         {
//             return tryGet<T>(r1);
//         }
//         catch (...)
//         {
//             if constexpr (sizeof...(rest) == 0)
//                 throw std::runtime_error("Data not found in any repo");
//             else
//                 return getFromRepos<T>(rest...);
//         }
//     }
// };

// #endif // FEM_DATA_CONTEXT_H
