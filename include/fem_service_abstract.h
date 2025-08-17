#ifndef FEM_SERVICE_ABSTRACT_H
#define FEM_SERVICE_ABSTRACT_H

#include <memory>
#include <tuple>
#include "fem_data_repository.h"

// Forward declare TypedContextView
template <typename... Types>
struct TypedContextView;

// A FEMService uses the data contained in a TypedContextView and performs its tasks.
// E.g. The ElementBuilder, given the Loads, Materials, RealConstants, ElementTypes,
// Builds the element matrices to be added to the global stiffness matrix
template <typename... ContextTypes>
class FEMServiceAbstract
{
public:
    using ContextView = TypedContextView<ContextTypes...>;

    virtual ~FEMServiceAbstract() = default;

    void set_context_view(std::shared_ptr<ContextView> view) { ctx_view = view; }
    std::shared_ptr<ContextView> get_context_view() const { return ctx_view; }

    // Convenience: get a specific context by type
    template <typename T>
    std::shared_ptr<T> context() { return ctx_view->template get<T>(); }

protected:
    FEMServiceAbstract(std::shared_ptr<ContextView> view = nullptr) : ctx_view(view) {}

    std::shared_ptr<ContextView> ctx_view;
};

template <typename... ContextTypes>
class FEMServiceWithOldContext : public FEMServiceAbstract<ContextTypes...>
{
public:
    using Base = FEMServiceAbstract<ContextTypes...>;
    using ContextView = typename Base::ContextView;

    FEMServiceWithOldContext(std::shared_ptr<ContextView> view = nullptr,
                             std::shared_ptr<ContextView> old_view = nullptr)
        : Base(view), old_ctx_view(old_view) {}

    void set_old_context_view(std::shared_ptr<ContextView> view) { old_ctx_view = view; }
    std::shared_ptr<ContextView> get_old_context_view() const { return old_ctx_view; }

    // Convenience: get a specific context by type from old context
    template <typename T>
    std::shared_ptr<T> old_context() { return old_ctx_view ? old_ctx_view->template get<T>() : nullptr; }

    // Update both: current becomes old, new_view becomes current
    void update_context_views(std::shared_ptr<ContextView> new_view)
    {
        std::swap(this->ctx_view, old_ctx_view);
        this->ctx_view = new_view;
    }

protected:
    std::shared_ptr<ContextView> old_ctx_view;
};
#endif // FEM_SERVICE_ABSTRACT_H