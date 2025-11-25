#ifndef FEM_DATA_REPOSITORY_H
#define FEM_DATA_REPOSITORY_H

#include <vector>
#include <memory>
#include <type_traits>
#include <typeinfo>
#include <typeindex>
#include <stdexcept>
#include <unordered_map>
#include "fe_data_object.h"

// A FEM data repository holds many FEDataObjects. E.g. a LoadStep carries information about the Loads, the
// Materials, the SolverSettings for that time step...
class FEMDataRepository
{
public:
    FEMDataRepository() = default;
    virtual ~FEMDataRepository() = default;

    // Set object by type
    template <typename T>
    void set_object(std::shared_ptr<T> obj)
    {
        static_assert(std::is_base_of<FEDataObjectAbstract, T>::value, "T must derive from FEDataObjectAbstract");
        data_objects[typeid(T)] = std::dynamic_pointer_cast<FEDataObjectAbstract>(obj);
    }

    // Get object by type
    template <typename T>
    std::shared_ptr<T> get_object() const
    {
        static_assert(std::is_base_of<FEDataObjectAbstract, T>::value, "T must derive from FEDataObjectAbstract");
        auto it = data_objects.find(typeid(T));
        if (it == data_objects.end())
            return nullptr; // <--- Return nullptr instead of error
        return std::dynamic_pointer_cast<T>(it->second);
    }

    // Remove object by type
    template <typename T>
    void remove_object()
    {
        data_objects.erase(typeid(T));
    }

    // Access to all objects (for introspection, not type-safe)
    const std::unordered_map<std::type_index, std::shared_ptr<FEDataObjectAbstract>> &get_all_objects() const
    {
        return data_objects;
    }

    // void add_all_from(std::shared_ptr<FEMDataRepository> other)
    // {
    //     if (!other)
    //         return;
    //     for (const auto &pair : other->get_all_objects())
    //     {
    //         data_objects[pair.first] = pair.second;
    //     }
    // }

    // void move_all_from(const std::shared_ptr<FEMDataRepository> &other)
    // {
    //     if (!other)
    //         return;

    //     for (auto &pair : other->data_objects)
    //     {
    //         data_objects[pair.first] = std::move(pair.second);
    //     }
    //     other->data_objects.clear();
    // }

    // void move_all_from(FEMDataRepository &other)
    // {
    //     for (auto &pair : other.data_objects)
    //     {
    //         data_objects[pair.first] = std::move(pair.second);
    //     }
    //     other.data_objects.clear();
    // }

protected:
    std::unordered_map<std::type_index, std::shared_ptr<FEDataObjectAbstract>> data_objects;
};

// The typed context views are the main way of interfacing data between classes in APEC-M
// E.g. the Assembler only needs details about the TimeStep integration scheme, it doesnt need
// the Loads, the Materials, etc.
// The user can create "slices" of data from the repositories, to be assignet to the services that
// can properly use the data
// Example of a "view" context for typed access
template <typename... Types>
struct TypedContextView
{
    std::tuple<std::shared_ptr<Types>...> ptrs;

    template <typename... Repos>
    static TypedContextView from(Repos &&...repos)
    {
        return TypedContextView{std::make_tuple(get_from_repos<Types>(std::forward<Repos>(repos)...)...)};
    }

    template <typename T>
    std::shared_ptr<T> get() const { return std::get<std::shared_ptr<T>>(ptrs); }

    template <typename... QueryTypes>
    TypedContextView<QueryTypes...> project() const
    {
        return TypedContextView<QueryTypes...>{std::make_tuple(get<QueryTypes>()...)};
    }

    template <typename... QueryTypes>
    std::shared_ptr<TypedContextView<QueryTypes...>> project_view() const
    {
        auto view = std::make_shared<TypedContextView<QueryTypes...>>();
        view->ptrs = std::make_tuple(get<QueryTypes>()...);

        const bool has_nullptr = std::apply([](const auto &...ptrs)
                                            { return ((ptrs == nullptr) || ...); },
                                            view->ptrs);

        if (has_nullptr)
            return nullptr;

        return view;
    }

private:
    // Overload for reference
    template <typename T, typename Repo>
    static std::shared_ptr<T> try_get(const Repo &repo)
    {
        return repo.template get_object<T>();
    }
    // Overload for shared_ptr
    template <typename T, typename Repo>
    static std::shared_ptr<T> try_get(const std::shared_ptr<Repo> &repo)
    {
        return repo ? repo->template get_object<T>() : nullptr;
    }

    template <typename T, typename Repo1, typename... Rest>
    static std::shared_ptr<T> get_from_repos(Repo1 &&r1, Rest &&...rest)
    {
        auto ptr = try_get<T>(std::forward<Repo1>(r1));
        if (ptr)
            return ptr;
        if constexpr (sizeof...(rest) == 0)
            libmesh_error_msg("Data not found in any repo");
        else
            return get_from_repos<T>(std::forward<Rest>(rest)...);
    }
};

#endif // FEM_DATA_REPOSITORY_H
