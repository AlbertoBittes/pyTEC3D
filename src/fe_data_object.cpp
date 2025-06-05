#include "fe_data_object.h"

void DataframeFEDataObject::read_data(std::ifstream &file)
{
    if (!file.is_open())
    {
        libmesh_error_msg("Failed to open file stream.");
    }

    // Clear the existing data_map
    data_map->clear();

    // Read the number of rows and columns
    uint32_t num_rows, num_columns;
    file.read(reinterpret_cast<char *>(&num_rows), sizeof(num_rows));
    file.read(reinterpret_cast<char *>(&num_columns), sizeof(num_columns));

    // Read column metadata (name and type)
    std::vector<std::string> column_names(num_columns);
    std::vector<int> column_types(num_columns); // 0=int, 1=double, 2=string

    for (uint32_t i = 0; i < num_columns; ++i)
    {
        // Read the column name
        uint32_t name_length;
        file.read(reinterpret_cast<char *>(&name_length), sizeof(name_length));

        std::string column_name(name_length, '\0');
        file.read(&column_name[0], name_length);

        // Read the column type
        int column_type;
        file.read(reinterpret_cast<char *>(&column_type), sizeof(column_type));

        column_names[i] = column_name;
        column_types[i] = column_type;
    }

    // Read the data row by row
    for (uint32_t row = 0; row < num_rows; ++row)
    {
        std::vector<std::variant<int, Real, std::string>> row_data;

        for (uint32_t col = 0; col < num_columns; ++col)
        {
            switch (column_types[col])
            {
            case 0: // int
            {
                int value;
                file.read(reinterpret_cast<char *>(&value), sizeof(value));
                row_data.push_back(value);
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
                uint32_t string_length;
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

        // Add the row to the data_map
        data_map->push_back(std::move(row_data));
    }
}

void DataframeMatLikeFEDataObject::populate_mat_like_map()
{
    // Clear the existing mat_like_map
    mat_like_map->clear();

    // Ensure the data_map is not empty
    if (data_map->empty())
    {
        return;
        // libmesh_error_msg("Data map is empty. Cannot populate mat_like_map.");
    }

    // Iterate over each row in the data_map
    for (const auto &row : *data_map)
    {
        if (row.empty())
        {
            libmesh_error_msg("Encountered an empty row in the data map.");
        }

        // Extract the first column as the ID (cast to dof_id_type)
        dof_id_type id = std::get<int>(row[0]);

        // Extract the remaining columns as Real values
        std::vector<Real> values;
        for (size_t col = 1; col < row.size(); ++col)
        {
            if (std::holds_alternative<Real>(row[col]))
            {
                values.push_back(std::get<Real>(row[col]));
            }
        }

        // Populate the mat_like_map
        (*mat_like_map)[id] = std::move(values);
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
    int32_t num_settings;
    file.read(reinterpret_cast<char *>(&num_settings), sizeof(num_settings));

    // Read each setting
    for (int32_t i = 0; i < num_settings; ++i)
    {
        // Read the size of the key (string label)
        int32_t key_size;
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
            int32_t string_size;
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