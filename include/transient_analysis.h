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
    using StaticAnalysis::StaticAnalysis;

    bool solve_transient() override;
};

#endif // TRANSIENT_ANALYSIS_H
