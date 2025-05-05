#ifndef FE_DATA_OBJECT_H
#define FE_DATA_OBJECT_H

#include <memory>
#include <string>
#include <vector>
#include <variant>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "libmesh/libmesh.h"
#include "../include/load_step.h" // Include the definition of LoadStep

typedef std::unordered_map<dof_id_type, std::vector<Real>> MaterialsMap;
typedef std::unordered_map<dof_id_type, std::vector<Real>> RealConstantsMap;
typedef std::unordered_map<std::string, std::variant<int, double, std::string>> SettingsMap;

template <typename T>
class FEDataObject
{
public:
    // Virtual destructor for proper cleanup of derived classes
    virtual ~FEDataObject() = default;

    // Pure virtual method to read data from a binary file stream
    virtual void read_data(std::ifstream &file) = 0;

    // Getter for the parent LoadStep
    virtual std::shared_ptr<LoadStep> get_parent_load_step() const
    {
        return parent_load_step;
    }

    // Setter for the parent LoadStep
    virtual void set_parent_load_step(std::shared_ptr<LoadStep> parent)
    {
        parent_load_step = std::move(parent);
    }

    // Getter for the data map
    virtual std::shared_ptr<T> get_data_map() const
    {
        return data_map;
    }

    // Setter for the data map
    virtual void set_data_map(std::shared_ptr<T> map)
    {
        data_map = std::move(map);
    }

    template <typename T>
    inline T read_value(std::ifstream &file)
    {
        T value;
        file.read(reinterpret_cast<char *>(&value), sizeof(T));
        return value;
    }

protected:
    // Protected constructor to prevent direct instantiation
    FEDataObject(std::shared_ptr<LoadStep> parent = nullptr,
                 std::shared_ptr<BinaryReader> reader = nullptr,
                 std::shared_ptr<T> map = nullptr)
        : parent_load_step(std::move(parent)),
          data_map(map ? std::move(map) : std::make_shared<T>()) {}

    // Attributes
    std::shared_ptr<LoadStep> parent_load_step; // Stores the parent LoadStep
    std::shared_ptr<T> data_map;                // Stores the data map (type-specific)
};

class CSVFEDataObject : public FEDataObject<MaterialsMap>
{
public:
    // Constructor
    CSVFEDataObject(std::shared_ptr<LoadStep> parent = nullptr,
                    std::shared_ptr<MaterialsMap> map = nullptr)
        : FEDataObject(parent, map) {}

    // Implement the read_data method
    void read_data(std::ifstream &file) override
    {
        auto data_map = std::make_shared<MaterialsMap>();

        if (!file.is_open())
        {
            libmesh_error_msg("Failed to open file: ");
        }

        // Read the number of rows and columns
        uint32_t num_rows = read_value<uint32_t>(file);
        uint32_t num_columns = read_value<uint32_t>(file);

        // Read the data row by row
        for (uint32_t row = 0; row < num_rows; ++row)
        {
            // Read the first column as the ID
            uint32_t id = read_value<uint32_t>(file);

            // Cast the ID to dof_id_type
            dof_id_type dof_id = static_cast<dof_id_type>(id);

            // Read the remaining columns as Real values
            std::vector<Real> values(num_columns - 1);
            for (uint32_t col = 0; col < num_columns - 1; ++col)
            {
                Real value = read_value<Real>(file);
                values[col] = value;
            }

            // Store the values in the map
            (*data_map)[dof_id] = std::move(values);
        }

        file.close();
    }
};

class SettingsFEDataObject : public FEDataObject<SettingsMap>
{
public:
    // Constructor
    SettingsFEDataObject(std::shared_ptr<LoadStep> parent = nullptr,
                         std::shared_ptr<SettingsMap> map = nullptr)
        : FEDataObject(parent, map) {}

    template <typename T>
    T get_setting(const std::string &key) const
    {
        auto it = data_map->find(key);
        if (it == data_map->end())
        {
            libmesh_error_msg("Setting not found: " + key);
        }
        return std::get<T>(it->second);
    }

    template <typename T>
    void set_setting(const std::string &key, T value)
    {
        (*data_map)[key] = value;
    }
    // Implement the read_data method
    void read_data(std::ifstream &file) override
    {
        auto data_map = std::make_shared<SettingsMap>();

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
};

#endif // FE_DATA_OBJECT_H