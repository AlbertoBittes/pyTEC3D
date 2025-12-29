#include "dof_couplings_factory.h"
#include "variables_map.h"

#include <algorithm>
#include <type_traits>

namespace
{
    const std::vector<int> &extract_int_array(const Globals::NDArray &array)
    {
        if (!std::holds_alternative<std::vector<int>>(array.values))
            libmesh_error_msg("Dof couplings: expected integer array.");
        return std::get<std::vector<int>>(array.values);
    }

    const std::vector<std::string> &extract_string_array(const Globals::NDArray &array)
    {
        if (!std::holds_alternative<std::vector<std::string>>(array.values))
            libmesh_error_msg("Dof couplings: expected string array.");
        return std::get<std::vector<std::string>>(array.values);
    }

    unsigned int to_unsigned(int value, const std::string &label)
    {
        if (value < 0)
            libmesh_error_msg("Dof couplings: negative value encountered in '" + label + "'.");
        return static_cast<unsigned int>(value);
    }

    std::size_t stride_from_shape(const Globals::NDArray &array, std::size_t entry_count, const std::string &label)
    {
        if (!array.shape.empty())
        {
            if (array.shape[0] != entry_count)
                libmesh_error_msg("Dof couplings: shape mismatch for '" + label + "'.");

            std::size_t stride = 1;
            for (std::size_t dim = 1; dim < array.shape.size(); ++dim)
                stride *= array.shape[dim];
            return stride;
        }

        const std::size_t value_count = BinaryDataFEDataObject::flattened_size(array);
        if (entry_count == 0)
            return 0;
        if (value_count % entry_count != 0)
            libmesh_error_msg("Dof couplings: flattened size mismatch for '" + label + "'.");
        return value_count / entry_count;
    }

    std::vector<std::string> to_global_labels(const std::string &label)
    {
        const auto matches_global = VariablesMap::get_keys_by_predicate(
            VariablesMap::GLOBAL_VARIABLE_INFO_MAP,
            [&label](const VariableInfo &info)
            {
                return info.global_label == label;
            });

        if (!matches_global.empty())
            return {label};

        auto matches_api = VariablesMap::get_keys_by_predicate(
            VariablesMap::GLOBAL_VARIABLE_INFO_MAP,
            [&label](const VariableInfo &info)
            {
                return info.api_label == label;
            });

        if (matches_api.empty())
            libmesh_error_msg("Dof couplings: unknown variable label '" + label + "'.");

        std::sort(matches_api.begin(), matches_api.end(),
                  [](GlobalVariableName lhs, GlobalVariableName rhs)
                  {
                      return VariablesMap::get_info_by_key(lhs).component_index <
                             VariablesMap::get_info_by_key(rhs).component_index;
                  });

        std::vector<std::string> labels;
        labels.reserve(matches_api.size());
        for (const auto &key : matches_api)
            labels.push_back(VariablesMap::get_info_by_key(key).global_label);
        return labels;
    }

    void append_unique(std::vector<std::string> &values, const std::string &label)
    {
        if (std::find(values.begin(), values.end(), label) == values.end())
            values.push_back(label);
    }
} // namespace

DofCouplingsObject::DofCouplingsObject(std::ifstream &file)
    : BinaryDataFEDataObject(file)
{
    parse_couplings_from_map();
}

void DofCouplingsObject::read_data(std::ifstream &file)
{
    BinaryDataFEDataObject::read_data(file);
    parse_couplings_from_map();
}

void DofCouplingsObject::parse_couplings_from_map()
{
    couplings.clear();

    if (empty() || this->get_data_map()->columns.empty())
        return;

    if (!has_column("array_name") || !has_column("cp_id") || !has_column("pairs") ||
        !has_column("roots") || !has_column("isolated") || !has_column("coupled_variables"))
        libmesh_error_msg("Dof couplings: required keys missing from settings data.");

    const auto &array_name_array = get_column("array_name");
    const auto &cp_id_array = get_column("cp_id");
    const auto &pairs_array = get_column("pairs");
    const auto &roots_array = get_column("roots");
    const auto &isolated_array = get_column("isolated");
    const auto &variables_array = get_column("coupled_variables");

    const auto &array_names = extract_string_array(array_name_array);
    const auto &cp_ids = extract_int_array(cp_id_array);
    const auto &pair_values = extract_int_array(pairs_array);
    const auto &root_values = extract_int_array(roots_array);
    const auto &isolated_values = extract_int_array(isolated_array);
    const auto &variable_values = extract_string_array(variables_array);

    const std::size_t count = entry_count();
    if (array_names.size() < count || cp_ids.size() < count)
        libmesh_error_msg("Dof couplings: scalar column size mismatch.");

    const std::size_t pair_stride = stride_from_shape(pairs_array, count, "pairs");
    if (pair_stride % 2 != 0)
        libmesh_error_msg("Dof couplings: pairs stride must be divisible by 2.");
    const std::size_t roots_stride = stride_from_shape(roots_array, count, "roots");
    const std::size_t isolated_stride = stride_from_shape(isolated_array, count, "isolated");
    const std::size_t variables_stride = stride_from_shape(variables_array, count, "coupled_variables");

    if (pair_values.size() < count * pair_stride)
        libmesh_error_msg("Dof couplings: pair values size mismatch.");
    if (root_values.size() < count * roots_stride)
        libmesh_error_msg("Dof couplings: root values size mismatch.");
    if (isolated_values.size() < count * isolated_stride)
        libmesh_error_msg("Dof couplings: isolated values size mismatch.");
    if (variable_values.size() < count * variables_stride)
        libmesh_error_msg("Dof couplings: coupled variables size mismatch.");

    couplings.reserve(count);

    for (std::size_t i = 0; i < count; ++i)
    {
        DofCouplingEntry entry;
        entry.array_name = array_names[i];
        entry.cp_id = cp_ids[i];

        const std::size_t pair_offset = i * pair_stride;
        const std::size_t pairs_count = pair_stride / 2;
        entry.pairs.reserve(pairs_count);
        for (std::size_t p = 0; p < pairs_count; ++p)
        {
            const std::size_t base = pair_offset + p * 2;
            entry.pairs.emplace_back(
                to_unsigned(pair_values[base], "pairs"),
                to_unsigned(pair_values[base + 1], "pairs"));
        }

        const std::size_t roots_offset = i * roots_stride;
        entry.roots.reserve(roots_stride);
        for (std::size_t r = 0; r < roots_stride; ++r)
            entry.roots.push_back(to_unsigned(root_values[roots_offset + r], "roots"));

        const std::size_t isolated_offset = i * isolated_stride;
        entry.isolated.reserve(isolated_stride);
        for (std::size_t r = 0; r < isolated_stride; ++r)
            entry.isolated.push_back(
                to_unsigned(isolated_values[isolated_offset + r], "isolated"));

        const std::size_t variables_offset = i * variables_stride;
        entry.coupled_variables.reserve(variables_stride);
        for (std::size_t v = 0; v < variables_stride; ++v)
        {
            const auto &raw_label = variable_values[variables_offset + v];
            for (const auto &global_label : to_global_labels(raw_label))
                append_unique(entry.coupled_variables, global_label);
        }

        couplings.push_back(std::move(entry));
    }
}

const std::vector<DofCouplingEntry> &DofCouplingsObject::get_couplings() const
{
    return couplings;
}

void DofCouplingsObject::print_data_info() const
{
    if (couplings.empty())
    {
        libMesh::out << "  DOF couplings: none configured." << std::endl;
        return;
    }

    libMesh::out << "  DOF couplings: " << couplings.size() << std::endl;
    for (const auto &entry : couplings)
    {
        libMesh::out << "    " << entry.array_name << " cp_id=" << entry.cp_id
                     << " pairs=" << entry.pairs.size()
                     << " roots=" << entry.roots.size()
                     << " isolated=" << entry.isolated.size()
                     << " variables=" << entry.coupled_variables.size()
                     << std::endl;
    }
}
