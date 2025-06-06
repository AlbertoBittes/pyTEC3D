#ifndef FEM_SERVICE_ABSTRACT_H
#define FEM_SERVICE_ABSTRACT_H

#include <memory>
#include <tuple>
#include "fem_data_repository.h"

// Forward declare TypedContextView
template <typename... Types>
struct TypedContextView;

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

#endif // FEM_SERVICE_ABSTRACT_H