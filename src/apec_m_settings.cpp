#include "apec_m_settings.h"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <vector>

#include <libmesh/libmesh.h>

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
    target_central_axis_element_id = 0;
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

    const auto *data_map = this->get_data_map();
    if (data_map && !data_map->columns.empty())
    {
        const std::vector<std::string> required_columns = {
            "R0",
            "max_picard_iterations",
            "tolerance_value",
            "alpha_value",
            "plasma_current",
            "beta_p"};
        std::vector<std::string> missing;
        missing.reserve(required_columns.size());
        for (const auto &name : required_columns)
        {
            if (!this->has_column(name))
                missing.push_back(name);
        }
        if (!missing.empty())
        {
            std::ostringstream oss;
            oss << "ApecMSettingsObject: missing required columns: ";
            for (std::size_t i = 0; i < missing.size(); ++i)
            {
                if (i > 0)
                    oss << ", ";
                oss << missing[i];
            }
            libmesh_error_msg(oss.str());
        }
    }

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

    if (this->has_column("target_central_axis_element_id"))
        target_central_axis_element_id = static_cast<dof_id_type>(this->get_int("target_central_axis_element_id"));

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

    if (this->has_column("print_full_output"))
        print_full_output_ = static_cast<bool>(this->get_int("print_full_output"));

    if (this->has_column("main_mesh"))
        main_mesh_ = this->get_string("main_mesh");

    if (this->has_column("static_mesh"))
        static_mesh_ = this->get_string("static_mesh");
}
