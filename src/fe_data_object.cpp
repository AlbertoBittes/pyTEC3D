#include "fe_data_object.h"
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkStringArray.h>

void DataframeFEDataObject::read_data(std::ifstream &file)
{
    if (!file.is_open())
    {
        libmesh_error_msg("Failed to open file stream.");
    }

    data_map->clear();

    // Read the number of rows and columns
    dof_id_type num_rows, num_columns;
    file.read(reinterpret_cast<char *>(&num_rows), sizeof(num_rows));
    file.read(reinterpret_cast<char *>(&num_columns), sizeof(num_columns));

    // Read column metadata (name and type)
    std::vector<std::string> column_names(num_columns);
    std::vector<int32_t> column_types(num_columns); // Use int32_t for type indicator

    for (dof_id_type i = 0; i < num_columns; ++i)
    {
        // Read the column name
        dof_id_type name_length;
        file.read(reinterpret_cast<char *>(&name_length), sizeof(name_length));

        std::string column_name(name_length, '\0');
        file.read(&column_name[0], name_length);

        // Read the column type
        int32_t column_type;
        file.read(reinterpret_cast<char *>(&column_type), sizeof(column_type));

        column_names[i] = column_name;
        column_types[i] = column_type;
    }

    // Read the data row by row
    for (dof_id_type row = 0; row < num_rows; ++row)
    {
        std::vector<std::variant<int, Real, std::string>> row_data;

        for (dof_id_type col = 0; col < num_columns; ++col)
        {
            switch (static_cast<int>(column_types[col]))
            {
            case 0: // int
            {
                int32_t value;
                file.read(reinterpret_cast<char *>(&value), sizeof(value));
                row_data.push_back(static_cast<int>(value)); // If your variant expects int
                break;
            }
            case 1: // double
            {
                Real value;
                file.read(reinterpret_cast<char *>(&value), sizeof(value));
                row_data.push_back(value);
                break;
            }
            case 2: // string
            {
                dof_id_type string_length;
                file.read(reinterpret_cast<char *>(&string_length), sizeof(string_length));

                std::string value(string_length, '\0');
                file.read(&value[0], string_length);
                row_data.push_back(value);
                break;
            }
            default:
                libmesh_error_msg("Unsupported column type encountered in binary file.");
            }
        }

        data_map->push_back(std::move(row_data));
    }
}

template <typename T>
void DataframeMapFEDataObject<T>::populate_dataframe_map()
{
    dataframe_map->clear();

    if (data_map->empty())
        return;

    for (const auto &row : *data_map)
    {
        if (row.empty())
            libmesh_error_msg("Encountered an empty row in the data map.");

        // First column is always int
        int id = std::get<int>(row[0]);

        // Remaining columns are of type T
        std::vector<T> values;
        for (size_t col = 1; col < row.size(); ++col)
        {
            if constexpr (std::is_same<T, Real>::value)
            {
                if (std::holds_alternative<Real>(row[col]))
                    values.push_back(std::get<Real>(row[col]));
                else
                    libmesh_error_msg("Expected Real in data map row.");
            }
            else if constexpr (std::is_same<T, int>::value)
            {
                if (std::holds_alternative<int>(row[col]))
                    values.push_back(std::get<int>(row[col]));
                else
                    libmesh_error_msg("Expected int in data map row.");
            }
            else if constexpr (std::is_same<T, std::string>::value)
            {
                if (std::holds_alternative<std::string>(row[col]))
                    values.push_back(std::get<std::string>(row[col]));
                else
                    libmesh_error_msg("Expected string in data map row.");
            }
            else
            {
                static_assert(sizeof(T) == 0, "Unsupported type for DataframeMapFEDataObject");
            }
        }

        // Store the vector<T> in the map with int key
        (*dataframe_map)[id] = std::move(values);
    }
}

void SettingsFEDataObject::read_data(std::ifstream &file)
{
    this->data_map->clear();

    if (!file.is_open())
    {
        libmesh_error_msg("Failed to open file stream.");
    }

    // Read the number of settings
    dof_id_type num_settings;
    file.read(reinterpret_cast<char *>(&num_settings), sizeof(num_settings));

    // Read each setting
    for (dof_id_type i = 0; i < num_settings; ++i)
    {
        // Read the size of the key (string label)
        dof_id_type key_size;
        file.read(reinterpret_cast<char *>(&key_size), sizeof(key_size));

        // Read the key (string label)
        std::string key(key_size, '\0');
        file.read(&key[0], key_size);

        // Read the type of the value
        int32_t value_type;
        file.read(reinterpret_cast<char *>(&value_type), sizeof(value_type));

        // Read the value based on its type
        if (value_type == 0)
        { // int
            int32_t int_value;
            file.read(reinterpret_cast<char *>(&int_value), sizeof(int_value));
            (*data_map)[key] = int_value;
        }
        else if (value_type == 1)
        { // double
            double double_value;
            file.read(reinterpret_cast<char *>(&double_value), sizeof(double_value));
            (*data_map)[key] = double_value;
        }
        else if (value_type == 2)
        { // string
            dof_id_type string_size;
            file.read(reinterpret_cast<char *>(&string_size), sizeof(string_size));
            std::string string_value(string_size, '\0');
            file.read(&string_value[0], string_size);
            (*data_map)[key] = string_value;
        }
        else
        {
            libmesh_error_msg("Unsupported value type encountered in binary file.");
        }
    }
}

template <typename TMap>
void vtkFieldDataFEDataObject<TMap>::read_from_vtk_field_data(vtkFieldData *field_data, const std::vector<std::string> &field_names)
{
    if (!field_data)
        libmesh_error_msg("vtkFieldDataFEDataObject::read_from_vtk_field_data: null vtkFieldData.");

    vtkIdType num_rows = 0;
    if (!field_names.empty())
    {
        auto *arr = field_data->GetAbstractArray(field_names.front().c_str());
        if (arr)
            num_rows = arr->GetNumberOfTuples();
    }

    auto &data_map = *(this->data_map);
    data_map.clear();

    // Use unsigned int keys
    for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
        data_map[irow] = std::tuple<>();

    for (const auto &base_name : field_names)
    {
        vtkAbstractArray *val_array = field_data->GetAbstractArray(base_name.c_str());
        vtkIntArray *shape_array = vtkIntArray::SafeDownCast(field_data->GetAbstractArray((base_name + "_shape").c_str()));
        vtkIntArray *offset_array = vtkIntArray::SafeDownCast(field_data->GetAbstractArray((base_name + "_offset").c_str()));

        if (!val_array)
        {
            libmesh_warning("Field '" << base_name << "' not found — skipping.");
            continue;
        }

        // --------------------------
        // SCALAR FIELDS
        // --------------------------
        if (!shape_array || !offset_array)
        {
            if (auto *int_arr = vtkIntArray::SafeDownCast(val_array))
            {
                for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
                    append_to_tuple(data_map[irow], int_arr->GetValue(irow));
            }
            else if (auto *dbl_arr = vtkDoubleArray::SafeDownCast(val_array))
            {
                for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
                    append_to_tuple(data_map[irow], dbl_arr->GetValue(irow));
            }
            else if (auto *str_arr = vtkStringArray::SafeDownCast(val_array))
            {
                for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
                    append_to_tuple(data_map[irow], std::string(str_arr->GetValue(irow)));
            }
            else
            {
                libmesh_error_msg("Unsupported scalar array type for field '" << base_name << "'.");
            }
            continue;
        }

        // --------------------------
        // ITERABLE FIELDS
        // --------------------------
        read_iterable_field(data_map, val_array, shape_array, offset_array, num_rows, base_name);
    }
}

// Append a new element to a tuple
template <typename TMap>
template <typename Tuple, typename T>
void vtkFieldDataFEDataObject<TMap>::append_to_tuple(Tuple &t, T &&value)
{
    t = std::tuple_cat(t, std::make_tuple(std::forward<T>(value)));
}

// Helper: read iterable fields
template <typename TMap>
void vtkFieldDataFEDataObject<TMap>::read_iterable_field(
    std::unordered_map<unsigned int, std::tuple<int>, std::hash<unsigned int>, std::equal_to<vtkGenericDataArrayLookupHelper<vtkGenericDataArray<vtkAOSDataArrayTemplate<unsigned int>, unsigned int>::SelfType>::ValueType>, std::allocator<std::pair<const unsigned int, std::tuple<int>>>> &data_map,
    vtkAbstractArray *val_array,
    vtkIntArray *shape_array,
    vtkIntArray *offset_array,
    vtkIdType num_rows,
    const std::string &base_name)
{
    const int *shape_data = shape_array->GetPointer(0);
    const int *offset_data = offset_array->GetPointer(0);

    // Integers
    if (auto *int_arr = vtkIntArray::SafeDownCast(val_array))
    {
        const int *vals = int_arr->GetPointer(0);
        vtkIdType cursor = 0;
        for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
        {
            int shape_start = offset_data[irow];
            int ndim = shape_data[shape_start];

            std::vector<int> dims(ndim);
            vtkIdType total_elems = 1;
            for (int j = 0; j < ndim; ++j)
            {
                dims[j] = shape_data[shape_start + 1 + j];
                total_elems *= dims[j];
            }

            std::vector<int> flat(vals + cursor, vals + cursor + total_elems);
            cursor += total_elems;

            append_to_tuple(data_map[irow], reshape(flat, dims));
        }
    }
    // Doubles
    else if (auto *dbl_arr = vtkDoubleArray::SafeDownCast(val_array))
    {
        const double *vals = dbl_arr->GetPointer(0);
        vtkIdType cursor = 0;
        for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
        {
            int shape_start = offset_data[irow];
            int ndim = shape_data[shape_start];

            std::vector<int> dims(ndim);
            vtkIdType total_elems = 1;
            for (int j = 0; j < ndim; ++j)
            {
                dims[j] = shape_data[shape_start + 1 + j];
                total_elems *= dims[j];
            }

            std::vector<double> flat(vals + cursor, vals + cursor + total_elems);
            cursor += total_elems;

            append_to_tuple(data_map[irow], reshape(flat, dims));
        }
    }
    // Strings
    else if (auto *str_arr = vtkStringArray::SafeDownCast(val_array))
    {
        vtkIdType cursor = 0;
        for (unsigned int irow = 0; irow < static_cast<unsigned int>(num_rows); ++irow)
        {
            int shape_start = offset_data[irow];
            int ndim = shape_data[shape_start];

            std::vector<int> dims(ndim);
            vtkIdType total_elems = 1;
            for (int j = 0; j < ndim; ++j)
            {
                dims[j] = shape_data[shape_start + 1 + j];
                total_elems *= dims[j];
            }

            std::vector<std::string> flat;
            flat.reserve(total_elems);
            for (vtkIdType k = 0; k < total_elems; ++k)
                flat.push_back(str_arr->GetValue(cursor++));
            append_to_tuple(data_map[irow], reshape(flat, dims));
        }
    }
    else
    {
        libmesh_error_msg("Unsupported iterable array type for field '" + base_name + "'.");
    }
}

// Reshape flattened vector to nested vector
template <typename TMap>
template <typename T>
std::vector<std::vector<T>> vtkFieldDataFEDataObject<TMap>::reshape(
    const std::vector<T> &flat,
    const std::vector<int> &dims,
    size_t level,
    size_t offset)
{
    if (dims.empty())
        return {};

    if (level == dims.size() - 1)
    {
        std::vector<std::vector<T>> leaf;
        size_t chunk = flat.size() / dims[level];
        for (int i = 0; i < dims[level]; ++i)
            leaf.emplace_back(flat.begin() + offset + i * chunk,
                              flat.begin() + offset + (i + 1) * chunk);
        return leaf;
    }

    size_t stride = flat.size() / dims[level];
    std::vector<std::vector<T>> result;
    for (int i = 0; i < dims[level]; ++i)
    {
        auto sub = reshape(flat, dims, level + 1, offset + i * stride);
        result.insert(result.end(), sub.begin(), sub.end());
    }
    return result;
}

// Explicit instantiation if needed:
// template class vtkFieldDataFEDataObject<YourMapType>;
// template class vtkFieldDataFEDataObject<std::unordered_map<unsigned int, std::tuple<int, double, std::string>>>;
