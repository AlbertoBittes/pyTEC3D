#include "apec_m_settings.h"

#include <algorithm>
#include <cctype>

ApecMSettingsObject::ApecMSettingsObject()
    : BinaryDataFEDataObject()
{
    parse_columns();
}

ApecMSettingsObject::ApecMSettingsObject(std::ifstream &file)
    : BinaryDataFEDataObject(file)
{
    parse_columns();
}

void ApecMSettingsObject::read_data(std::ifstream &file)
{
    BinaryDataFEDataObject::read_data(file);
    parse_columns();
}

void ApecMSettingsObject::parse_columns()
{
    auto normalize = [](std::string value) {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return value;
    };

    if (has_column("type"))
    {
        auto type_string = normalize(get_string("type"));
        if (type_string == "evolution")
            mode = AnalysisMode::EVOLUTION;
        else
            mode = AnalysisMode::EQUILIBRIUM;
    }
    else
        mode = AnalysisMode::EQUILIBRIUM;

    if (has_column("target_central_axis"))
    {
        const auto &axis = get_real_column("target_central_axis");
        if (axis.size() >= 3)
            target_central_axis = libMesh::Point(axis[0], axis[1], axis[2]);
        else
            libMesh::out << "WARNING: 'target_central_axis' requires three coordinates." << std::endl;
    }

    if (has_column("R0"))
        R0 = get_real("R0");

    if (has_column("max_iterations"))
    {
        int iteration_count = get_int("max_iterations");
        if (iteration_count < 1)
            iteration_count = 1;
        max_iterations = static_cast<unsigned int>(iteration_count);
    }

    if (has_column("min_distance_to_axis"))
        min_distance = get_real("min_distance_to_axis");

    if (has_column("cpgain"))
        cp_gain = get_real("cpgain");

    if (has_column("tolerance"))
        tolerance_value = get_real("tolerance");

    if (has_column("bradx"))
        bradx_value = get_real("bradx");
    if (has_column("brady"))
        brady_value = get_real("brady");
    if (has_column("bverx"))
        bverx_value = get_real("bverx");
    if (has_column("bvery"))
        bvery_value = get_real("bvery");

    if (has_column("alpha"))
        alpha_value = get_real("alpha");

    plasma_current_profile.clear();
    beta_p_profile.clear();

    if (has_column("I_plasma"))
    {
        const auto &values = get_real_column("I_plasma");
        plasma_current_profile.assign(values.begin(), values.end());
    }

    if (has_column("beta_p"))
    {
        const auto &values = get_real_column("beta_p");
        beta_p_profile.assign(values.begin(), values.end());
    }

    if (!plasma_current_profile.empty() && !beta_p_profile.empty() && plasma_current_profile.size() != beta_p_profile.size())
    {
        const auto limit = std::min(plasma_current_profile.size(), beta_p_profile.size());
        plasma_current_profile.resize(limit);
        beta_p_profile.resize(limit);
        libMesh::out << "WARNING: 'I_plasma' and 'beta_p' column sizes differ. Truncating to " << limit << " entries." << std::endl;
    }

    if (plasma_current_profile.empty())
        plasma_current_profile.push_back(0.0);

    if (beta_p_profile.empty())
        beta_p_profile.resize(plasma_current_profile.size(), 0.0);
}
