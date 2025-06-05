#ifndef FEM_SERVICE_ABSTRACT_H
#define FEM_SERVICE_ABSTRACT_H

#include <memory>
#include <vector>

// Forward declarations (add more as needed)

template <typename ContextType>
class FEMServiceAbstract
{
    static_assert(std::is_base_of<FEMDataContext, ContextType>::value, "ContextType must derive from FEMDataContext");

public:
    virtual ~FEMServiceAbstract() = 0;

    // Example: set/get context
    void set_context(std::shared_ptr<ContextType> context) { ctx = context; }
    std::shared_ptr<ContextType> get_context() const { return ctx; }

    // Add your service interface methods here
    // virtual void run() = 0;

protected:
    FEMServiceAbstract(std::shared_ptr<ContextType> context = nullptr)
        : ctx(context) {}

    std::shared_ptr<ContextType> ctx;
};

#endif // FEM_SERVICE_ABSTRACT_H