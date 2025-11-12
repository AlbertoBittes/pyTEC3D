#ifndef APEC_M_ANALYSIS_H
#define APEC_M_ANALYSIS_H

#include <memory>
#include <string>
#include <utility>

#include "analysis_base.h"

class StaticAnalysis : public AnalysisBase
{
public:
    explicit StaticAnalysis(ApecMEquationSystems *eq_systems = nullptr);

    bool solve(const std::string &identifier) override;
};

class APECMAnalysis : public AnalysisBase
{
public:
    explicit APECMAnalysis(ApecMEquationSystems *eq_systems = nullptr);

    bool solve(const std::string &identifier) override;
};

class DOFInterpolationAnalysis : public AnalysisBase
{
public:
    explicit DOFInterpolationAnalysis(ApecMEquationSystems *eq_systems = nullptr);

    bool solve(const std::string &identifier) override;
};

#endif // APEC_M_ANALYSIS_H
