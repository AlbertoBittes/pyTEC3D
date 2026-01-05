#include "apec_m_settings.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <sstream>
#include <type_traits>
#include <variant>
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

    if (this->has_column("evolution_mesh"))
        evolution_mesh_ = this->get_string("evolution_mesh");

    if (this->has_column("equilibrium_mesh"))
        equilibrium_mesh_ = this->get_string("equilibrium_mesh");

    if (this->has_column("n_slices"))
        n_slices_ = static_cast<unsigned int>(this->get_int("n_slices"));

    const auto parse_toroidal_map = [this](const std::string &key,
                                           std::vector<std::vector<dof_id_type>> &out)
    {
        out.clear();
        if (!this->has_column(key))
            return;

        const auto &array = this->get_column(key);
        if (array.shape.size() != 2 && array.shape.size() != 3)
            libmesh_error_msg("ApecMSettingsObject: '" + key + "' must be a 2D array.");

        std::size_t entry_dim = 0;
        std::size_t rows = 0;
        std::size_t cols = 0;
        if (array.shape.size() == 3)
        {
            entry_dim = array.shape[0];
            rows = array.shape[1];
            cols = array.shape[2];
        }
        else
        {
            rows = array.shape[0];
            cols = array.shape[1];
        }

        if (rows == 0 || cols == 0)
            return;
        if (array.shape.size() == 3 && entry_dim == 0)
            libmesh_error_msg("ApecMSettingsObject: '" + key + "' has no entries.");

        out.assign(rows, {});
        const std::size_t entry_stride = rows * cols;
        const std::size_t base_offset = 0;

        std::visit(
            [&](const auto &vec)
            {
                using VecType = std::decay_t<decltype(vec)>;
                using ElemType = typename VecType::value_type;
                if constexpr (std::is_same_v<ElemType, int>)
                {
                    for (std::size_t r = 0; r < rows; ++r)
                    {
                        auto &row = out[r];
                        row.reserve(cols);
                        for (std::size_t c = 0; c < cols; ++c)
                        {
                            const std::size_t idx = base_offset + r * cols + c;
                            const int value = vec[idx];
                            if (c == 0)
                            {
                                if (value < 0)
                                    libmesh_error_msg("ApecMSettingsObject: invalid input id in '" + key + "'.");
                                row.push_back(static_cast<dof_id_type>(value));
                            }
                            else if (value >= 0)
                            {
                                row.push_back(static_cast<dof_id_type>(value));
                            }
                        }
                    }
                }
                else if constexpr (std::is_same_v<ElemType, Real>)
                {
                    for (std::size_t r = 0; r < rows; ++r)
                    {
                        auto &row = out[r];
                        row.reserve(cols);
                        for (std::size_t c = 0; c < cols; ++c)
                        {
                            const std::size_t idx = base_offset + r * cols + c;
                            const Real value = vec[idx];
                            if (!std::isfinite(static_cast<double>(value)))
                                continue;
                            const auto rounded = static_cast<long long>(std::llround(value));
                            if (c == 0)
                            {
                                if (rounded < 0)
                                    libmesh_error_msg("ApecMSettingsObject: invalid input id in '" + key + "'.");
                                row.push_back(static_cast<dof_id_type>(rounded));
                            }
                            else if (rounded >= 0)
                            {
                                row.push_back(static_cast<dof_id_type>(rounded));
                            }
                        }
                    }
                }
                else
                {
                    libmesh_error_msg("ApecMSettingsObject: '" + key + "' must store numeric data.");
                }
            },
            array.values);
    };

    parse_toroidal_map("toroidal_projection_map_nodes", toroidal_projection_map_nodes_);
    parse_toroidal_map("toroidal_projection_map_elems", toroidal_projection_map_elems_);
}
