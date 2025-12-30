#include "fe_data_object.h"
#include "nd_array.h"

#include <algorithm>
#include <cmath>
#include <type_traits>

void BinaryDataFEDataObject::read_data(std::ifstream &file)
{
    if (!file.is_open())
        libmesh_error_msg("Failed to open file stream.");

    data_map->columns.clear();
    data_map->entry_count = 0;

    std::uint32_t entry_count = 0;
    file.read(reinterpret_cast<char *>(&entry_count), sizeof(entry_count));
    if (entry_count == 0)
    {
        data_map->entry_count = 0;
        return;
    }

    std::uint32_t num_keys = 0;
    file.read(reinterpret_cast<char *>(&num_keys), sizeof(num_keys));

    if (num_keys == 0)
    {
        data_map->entry_count = 0;
        return;
    }

    std::size_t deduced_entry_count = entry_count > 1 ? entry_count : 0;

    for (std::uint32_t i = 0; i < num_keys; ++i)
    {
        std::uint32_t key_size = 0;
        file.read(reinterpret_cast<char *>(&key_size), sizeof(key_size));

        std::string key(key_size, '\0');
        file.read(&key[0], key_size);

        int32_t payload_type = 0;
        file.read(reinterpret_cast<char *>(&payload_type), sizeof(payload_type));
        if (payload_type != 3)
            libmesh_error_msg("BinaryDataFEDataObject: unexpected payload type while reading key '" + key + "'.");

        std::uint32_t dtype_len = 0;
        file.read(reinterpret_cast<char *>(&dtype_len), sizeof(dtype_len));
        std::string dtype(dtype_len, '\0');
        file.read(&dtype[0], dtype_len);

        std::uint32_t ndim = 0;
        file.read(reinterpret_cast<char *>(&ndim), sizeof(ndim));

        NDArray array_value;
        array_value.shape.resize(ndim);

        std::size_t expected_elements = 1;
        for (std::uint32_t dim = 0; dim < ndim; ++dim)
        {
            std::uint32_t extent = 0;
            file.read(reinterpret_cast<char *>(&extent), sizeof(extent));
            array_value.shape[dim] = extent;
            expected_elements *= static_cast<std::size_t>(extent);
        }

        int32_t element_kind = 0;
        file.read(reinterpret_cast<char *>(&element_kind), sizeof(element_kind));

        switch (element_kind)
        {
        case 0:
        {
            std::vector<std::int64_t> buffer(expected_elements);
            if (expected_elements)
                file.read(reinterpret_cast<char *>(buffer.data()), buffer.size() * sizeof(std::int64_t));

            std::vector<int> values(buffer.size());
            for (std::size_t idx = 0; idx < buffer.size(); ++idx)
                values[idx] = static_cast<int>(buffer[idx]);

            array_value.values = std::move(values);
            break;
        }
        case 1:
        {
            std::vector<double> buffer(expected_elements);
            if (expected_elements)
                file.read(reinterpret_cast<char *>(buffer.data()), buffer.size() * sizeof(double));

            std::vector<Real> values(buffer.size());
            for (std::size_t idx = 0; idx < buffer.size(); ++idx)
                values[idx] = static_cast<Real>(buffer[idx]);

            array_value.values = std::move(values);
            break;
        }
        case 2:
        {
            std::uint32_t count = 0;
            file.read(reinterpret_cast<char *>(&count), sizeof(count));

            if (count != expected_elements)
                libmesh_error_msg("BinaryDataFEDataObject: string array size mismatch for key '" + key + "'.");

            std::vector<std::string> values(count);
            for (std::uint32_t idx = 0; idx < count; ++idx)
            {
                std::uint32_t str_len = 0;
                file.read(reinterpret_cast<char *>(&str_len), sizeof(str_len));

                std::string value(str_len, '\0');
                if (str_len)
                    file.read(&value[0], str_len);
                values[idx] = std::move(value);
            }

            array_value.values = std::move(values);
            break;
        }
        default:
            libmesh_error_msg("BinaryDataFEDataObject: unsupported element kind encountered for key '" + key + "'.");
        }

        const std::size_t value_count = std::visit([](const auto &vec)
                                                   { return NDArray::vector_size(vec); }, array_value.values);
        if (value_count != expected_elements)
            libmesh_error_msg("BinaryDataFEDataObject: flattened size mismatch while reading key '" + key + "'.");

        const std::size_t entry_span = NDArray::infer_entry_span(array_value.shape, value_count);
        if (entry_count > 1)
        {
            const bool matches_shape = array_value.shape.empty()
                                           ? (value_count == entry_count)
                                           : std::any_of(array_value.shape.begin(),
                                                         array_value.shape.end(),
                                                         [&](std::uint32_t dim)
                                                         { return dim == entry_count; });
            if (!matches_shape)
                libmesh_error_msg("BinaryDataFEDataObject: provided entry count does not match any dimension of key '" + key + "'.");
        }
        else
            deduced_entry_count = std::max(deduced_entry_count, entry_span);

        data_map->columns.emplace(std::move(key), std::move(array_value));
    }

    data_map->entry_count = (entry_count > 1) ? entry_count : deduced_entry_count;
}

const NDArray &BinaryDataFEDataObject::get_column(const std::string &key) const
{
    auto it = data_map->columns.find(key);
    if (it == data_map->columns.end())
        libmesh_error_msg("BinaryDataFEDataObject: column '" + key + "' not found.");
    return it->second;
}

bool BinaryDataFEDataObject::has_column(const std::string &key) const
{
    return data_map->columns.find(key) != data_map->columns.end();
}

std::size_t BinaryDataFEDataObject::entry_count() const
{
    return data_map->entry_count;
}

std::size_t BinaryDataFEDataObject::flattened_size(const NDArray &array)
{
    return std::visit([](const auto &vec)
                      { return NDArray::vector_size(vec); }, array.values);
}

void BinaryDataFEDataObject::ensure_index_in_bounds(const std::string &key, std::size_t index, std::size_t size)
{
    if (index >= size)
        libmesh_error_msg("BinaryDataFEDataObject: index " + std::to_string(index) + " out of bounds for column '" + key + "'.");
}

int BinaryDataFEDataObject::get_int(const std::string &key, std::size_t index) const
{
    const auto &column = get_column(key);
    const std::size_t size = flattened_size(column);
    ensure_index_in_bounds(key, index, size);

    return std::visit(
        [&](const auto &vec) -> int
        {
            using VecType = std::decay_t<decltype(vec)>;
            using ElemType = typename VecType::value_type;

            if constexpr (std::is_same_v<ElemType, int>)
                return vec[index];
            else if constexpr (std::is_same_v<ElemType, Real>)
                return static_cast<int>(std::llround(vec[index]));
            else if constexpr (std::is_same_v<ElemType, std::string>)
                return std::stoi(vec[index]);
        },
        column.values);
}

Real BinaryDataFEDataObject::get_real(const std::string &key, std::size_t index) const
{
    const auto &column = get_column(key);
    const std::size_t size = flattened_size(column);
    ensure_index_in_bounds(key, index, size);

    return std::visit(
        [&](const auto &vec) -> Real
        {
            using VecType = std::decay_t<decltype(vec)>;
            using ElemType = typename VecType::value_type;

            if constexpr (std::is_same_v<ElemType, int>)
                return static_cast<Real>(vec[index]);
            else if constexpr (std::is_same_v<ElemType, Real>)
                return vec[index];
            else if constexpr (std::is_same_v<ElemType, std::string>)
                return static_cast<Real>(std::stod(vec[index]));
        },
        column.values);
}

std::string BinaryDataFEDataObject::get_string(const std::string &key, std::size_t index) const
{
    const auto &column = get_column(key);
    const std::size_t size = flattened_size(column);
    ensure_index_in_bounds(key, index, size);

    return std::visit(
        [&](const auto &vec) -> std::string
        {
            using VecType = std::decay_t<decltype(vec)>;
            using ElemType = typename VecType::value_type;

            if constexpr (std::is_same_v<ElemType, int>)
                return std::to_string(vec[index]);
            else if constexpr (std::is_same_v<ElemType, Real>)
                return std::to_string(vec[index]);
            else if constexpr (std::is_same_v<ElemType, std::string>)
                return vec[index];
        },
        column.values);
}

const std::vector<int> &BinaryDataFEDataObject::get_int_column(const std::string &key) const
{
    const auto &column = get_column(key);
    if (!std::holds_alternative<std::vector<int>>(column.values))
        libmesh_error_msg("BinaryDataFEDataObject: column '" + key + "' does not store integers.");
    return std::get<std::vector<int>>(column.values);
}

const std::vector<Real> &BinaryDataFEDataObject::get_real_column(const std::string &key) const
{
    const auto &column = get_column(key);
    if (!std::holds_alternative<std::vector<Real>>(column.values))
        libmesh_error_msg("BinaryDataFEDataObject: column '" + key + "' does not store real values.");
    return std::get<std::vector<Real>>(column.values);
}

const std::vector<std::string> &BinaryDataFEDataObject::get_string_column(const std::string &key) const
{
    const auto &column = get_column(key);
    if (!std::holds_alternative<std::vector<std::string>>(column.values))
        libmesh_error_msg("BinaryDataFEDataObject: column '" + key + "' does not store strings.");
    return std::get<std::vector<std::string>>(column.values);
}

void BinaryDataFEDataObject::print_data_info() const
{
    const auto *map = this->get_data_map();
    if (!map || map->columns.empty())
    {
        libMesh::out << "  Binary data: no columns loaded." << std::endl;
        return;
    }

    libMesh::out << "  Binary data entries: " << entry_count() << std::endl;
    libMesh::out << "  Columns:" << std::endl;
    for (const auto &column_pair : map->columns)
    {
        libMesh::out << "    " << column_pair.first << " shape=" << NDArray::format_shape(column_pair.second.shape) << " values=";
        std::visit(
            [&](const auto &vec)
            {
                NDArray::print_vector_contents(vec);
            },
            column_pair.second.values);
        libMesh::out << std::endl;
    }
}
