#ifndef TRANSIENT_ANALYSIS_H
#define TRANSIENT_ANALYSIS_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "static_analysis.h"

// class ApecMLinearSystem;
// class ApecMExplicitSystem;

class TransientAnalysis : public StaticAnalysis
{
public:
    explicit TransientAnalysis(APEC_M *program = nullptr,
                               std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);

    bool solve_transient() override;

protected:
    TransientAnalysis(AnalysisBase::Type type,
                      APEC_M *program,
                      std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);
};

#endif // TRANSIENT_ANALYSIS_H
