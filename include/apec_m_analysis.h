#ifndef APEC_M_ANALYSIS_H
#define APEC_M_ANALYSIS_H

#include <string>
#include <utility>

class ApecMAnalysis
{
public:
    ApecMAnalysis() = default;

    explicit ApecMAnalysis(std::string analysis_name)
        : name(std::move(analysis_name)) {}

    virtual ~ApecMAnalysis() = default;

    const std::string &get_name() const { return name; }
    void set_name(const std::string &analysis_name) { name = analysis_name; }

private:
    std::string name;
};

#endif // APEC_M_ANALYSIS_H
