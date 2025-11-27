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
    target_central_axis = Point();
    R0 = 0.0;
    max_picard_iterations_ = 1;
    min_distance = 0.3;
    cp_gain = 0.0;
    tolerance_value = 1e-6;
    bradx_value = 0.0;
    brady_value = 0.0;
    bverx_value = 0.0;
    bvery_value = 0.0;
    alpha_value = 1.0;
    plasma_current = 0.0;
    beta_p = 0.0;

    const auto get_axis_component = [this](const std::string &key, unsigned int idx)
    {
        if (this->has_column(key))
        {
            const auto &vals = this->get_real_column(key);
            if (vals.size() > idx)
                target_central_axis(idx) = vals[idx];
        }
    };

    if (this->has_column("target_central_axis_x"))
        target_central_axis(0) = this->get_real("target_central_axis_x");
    else
        get_axis_component("target_central_axis", 0);

    if (this->has_column("target_central_axis_y"))
        target_central_axis(1) = this->get_real("target_central_axis_y");
    else
        get_axis_component("target_central_axis", 1);

    if (this->has_column("target_central_axis_z"))
        target_central_axis(2) = this->get_real("target_central_axis_z");
    else
        get_axis_component("target_central_axis", 2);

    if (this->has_column("R0"))
        R0 = this->get_real("R0");

    if (this->has_column("max_picard_iterations"))
        max_picard_iterations_ = static_cast<unsigned int>(this->get_int("max_picard_iterations"));

    if (this->has_column("min_distance"))
        min_distance = this->get_real("min_distance");

    if (this->has_column("cp_gain"))
        cp_gain = this->get_real("cp_gain");

    if (this->has_column("tolerance_value"))
        tolerance_value = this->get_real("tolerance_value");

    if (this->has_column("bradx_value"))
        bradx_value = this->get_real("bradx_value");

    if (this->has_column("brady_value"))
        brady_value = this->get_real("brady_value");

    if (this->has_column("bverx_value"))
        bverx_value = this->get_real("bverx_value");

    if (this->has_column("bvery_value"))
        bvery_value = this->get_real("bvery_value");

    if (this->has_column("alpha_value"))
        alpha_value = this->get_real("alpha_value");

    if (this->has_column("plasma_current"))
        plasma_current = this->get_real("plasma_current");

    if (this->has_column("beta_p"))
        beta_p = this->get_real("beta_p");
}
