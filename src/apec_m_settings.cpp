#include "apec_m_settings.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <sstream>
#include <type_traits>
#include <variant>
#include <vector>

#include <libmesh/libmesh.h>

namespace
{
    std::vector<std::vector<dof_id_type>> parse_toroidal_slice_map_columnwise(
        const NDArray &array,
        const std::string &column_name)
    {
        if (!std::holds_alternative<std::vector<int>>(array.values))
            libmesh_error_msg("ApecMSettingsObject: '" + column_name + "' must store integer IDs.");

        auto normalized_shape = array.shape;
        while (normalized_shape.size() > 2 && normalized_shape.front() == 1u)
            normalized_shape.erase(normalized_shape.begin());

        if (normalized_shape.size() != 2)
            libmesh_error_msg("ApecMSettingsObject: '" + column_name + "' must be a 2D map.");

        const std::size_t n_rows = static_cast<std::size_t>(normalized_shape[0]);
        const std::size_t n_cols = static_cast<std::size_t>(normalized_shape[1]);
        if (n_rows == 0 || n_cols == 0)
            return {};

        if (n_cols < 2)
            libmesh_error_msg("ApecMSettingsObject: '" + column_name + "' must include IDs and at least one toroidal slice.");

        const auto &row_major_values = array.view<int>();
        const std::size_t expected_size = n_rows * n_cols;
        if (row_major_values.size() != expected_size)
            libmesh_error_msg("ApecMSettingsObject: '" + column_name + "' has inconsistent shape and data size.");

        // Input arrays are row-oriented [2D_id, slice_1, ...]. Cache as column vectors:
        // [2D_ids, slice_1_ids, ...].
        std::vector<std::vector<dof_id_type>> columnwise(n_cols, std::vector<dof_id_type>(n_rows));
        for (std::size_t row = 0; row < n_rows; ++row)
            for (std::size_t col = 0; col < n_cols; ++col)
                columnwise[col][row] = static_cast<dof_id_type>(row_major_values[row * n_cols + col]);

        return columnwise;
    }
} // namespace

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
    alpha_mp = 1.0001;
    alpha_np = 0.83197;
    plasma_current = 0.0;
    beta_p = 0.0;
    toroidal_slices_nodes_.clear();
    toroidal_slices_elems_.clear();

    const auto *data_map = this->get_data_map();
    if (data_map && !data_map->columns.empty())
    {
        const std::vector<std::string> required_columns = {
            "R0",
            "max_picard_iterations",
            "tolerance_value",
            "alpha_mp",
            "alpha_np",
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

    if (this->has_column("alpha_mp"))
        alpha_mp = this->get_real("alpha_mp");

    if (this->has_column("alpha_np"))
        alpha_np = this->get_real("alpha_np");

    if (this->has_column("plasma_current"))
        plasma_current = this->get_real("plasma_current");

    if (this->has_column("beta_p"))
        beta_p = this->get_real("beta_p");

    if (this->has_column("print_full_output"))
        print_full_output_ = static_cast<bool>(this->get_int("print_full_output"));

    if (this->has_column("evolution_mesh"))
        evolution_mesh_ = this->get_string("evolution_mesh");

    if (this->has_column("equilibrium_mesh"))
        equilibrium_mesh_ = this->get_string("equilibrium_mesh");

    if (this->has_column("static_mesh"))
        static_mesh_ = this->get_string("static_mesh");

    if (this->has_column("AF0"))
        AF0 = this->get_real("AF0");

    if (data_map)
    {
        auto it_nodes = data_map->columns.find("toroidal_slices_nodes");
        if (it_nodes != data_map->columns.end())
            toroidal_slices_nodes_ = parse_toroidal_slice_map_columnwise(it_nodes->second, "toroidal_slices_nodes");

        auto it_elems = data_map->columns.find("toroidal_slices_elems");
        if (it_elems != data_map->columns.end())
            toroidal_slices_elems_ = parse_toroidal_slice_map_columnwise(it_elems->second, "toroidal_slices_elems");
    }
}
