#ifndef TRANSIENT_ANALYSIS_H
#define TRANSIENT_ANALYSIS_H

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "static_analysis.h"

// class ApecMLinearSystem;
// class ApecMExplicitSystem;

class TransientAnalysis : public virtual StaticAnalysis
{
public:
    explicit TransientAnalysis(APEC_M *program = nullptr,
                               std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);

    std::vector<std::filesystem::path> list_load_step_files() const;
    std::size_t load_step_file_count() const;

    bool solve_transient() override;

private:
    unsigned int cumulative_substep_index = 0;

protected:
    TransientAnalysis(AnalysisBase::Type type,
                      APEC_M *program,
                      std::string analysis_name,
                      std::shared_ptr<AnalysisBase::ContextView> context_view = nullptr);
};

#endif // TRANSIENT_ANALYSIS_H
