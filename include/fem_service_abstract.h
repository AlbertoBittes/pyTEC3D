#ifndef FEM_SERVICE_ABSTRACT_H
#define FEM_SERVICE_ABSTRACT_H

#include <memory>
#include <vector>
#include "fem_data_context.h"

// Forward declarations (add more as needed)
class FEMDataContext;

class FEMServiceAbstract
{
public:
    virtual ~FEMServiceAbstract() = 0;

    // Example: set/get context
    void set_context(std::shared_ptr<FEMDataContext> context) { ctx = context; }
    std::shared_ptr<FEMDataContext> get_context() const { return ctx; }

    // Add your service interface methods here
    // virtual void run() = 0;

protected:
    FEMServiceAbstract(std::shared_ptr<FEMDataContext> context = nullptr)
        : ctx(context) {}

    std::shared_ptr<FEMDataContext> ctx;
};

#endif // FEM_SERVICE_ABSTRACT_H