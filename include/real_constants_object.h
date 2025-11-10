#ifndef REAL_CONSTANTS_OBJECT_H
#define REAL_CONSTANTS_OBJECT_H

#include <algorithm>
#include <array>
#include <iomanip>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <vector>
#include "fe_data_object.h" // Include the base class

// Holds the RealConstants numbers (R in APDL) and their values to be used by ElementBuilder during
// element matrices' construction
class RealConstantsObject : public BinaryDataFEDataObject
{
public:
    explicit RealConstantsObject(std::ifstream &file)
        : BinaryDataFEDataObject(file)
    {
        build_constants_map();
    }

    // void read_data(std::ifstream &file) override
    // {
    //     BinaryDataFEDataObject::read_data(file);
    //     build_constants_map();
    // }

    // Check if a constant exists in the map
    bool has_constant(dof_id_type constant_id) const { return constants_map.find(constant_id) != constants_map.end(); }

    // Get the values of a constant by ID
    const std::vector<Real> &get_constant(dof_id_type constant_id) const
    {
        auto it = constants_map.find(constant_id);
        if (it == constants_map.end())
        {
            libmesh_error_msg("Constant not found: " + std::to_string(constant_id));
        }
        return it->second;
    }

    const std::unordered_map<dof_id_type, std::vector<Real>> &get_constants_map() const { return constants_map; }

    void print_data_info() const override
    {
        if (constants_map.empty())
        {
            libMesh::out << "  Real constants: none defined." << std::endl;
            return;
        }

        const auto format_value = [](Real value)
        {
            std::ostringstream oss;
            oss.setf(std::ios::scientific);
            oss << std::setprecision(6) << value;
            return oss.str();
        };

        libMesh::out << "  Real constants sets: " << constants_map.size() << std::endl;
        for (const auto &entry : constants_map)
        {
            libMesh::out << "    Constant ID " << entry.first << ":" << std::endl;
            const auto &values = entry.second;
            const std::size_t printable = std::min(values.size(), PROPERTY_COLUMNS.size());
            for (std::size_t i = 0; i < printable; ++i)
            {
                libMesh::out << "      " << PROPERTY_COLUMNS[i] << " = " << format_value(values[i]) << std::endl;
            }
            if (values.size() > PROPERTY_COLUMNS.size())
                libMesh::out << "      +" << (values.size() - PROPERTY_COLUMNS.size()) << " extra values" << std::endl;
        }
    }

private:
    static constexpr std::array<const char *, 22> PROPERTY_COLUMNS = {
        "CARE", "TURN", "LENG", "DIRZ", "FILL", "SC",
        "NC", "RAD", "TX", "TY", "TZ", "R",
        "SYM", "VELOX", "VELOY", "VELOZ", "OMEGAX", "OMEGAY",
        "OMEGAZ", "XLOC", "YLOC", "ZLOC"};

    void build_constants_map()
    {
        constants_map.clear();
        if (empty())
            return;

        if (!has_column("ansys_real_constant"))
            libmesh_error_msg("RealConstantsObject: missing 'ansys_real_constant' column.");

        const auto &ids = get_int_column("ansys_real_constant");
        const std::size_t count = entry_count();

        if (ids.size() < count)
            libmesh_error_msg("RealConstantsObject: inconsistent 'ansys_real_constant' column length.");

        for (std::size_t row = 0; row < count; ++row)
        {
            std::vector<Real> values;
            values.reserve(PROPERTY_COLUMNS.size());
            for (const auto *name : PROPERTY_COLUMNS)
            {
                if (!has_column(name))
                    libmesh_error_msg(std::string("RealConstantsObject: missing column '") + name + "'.");
                const auto &col = get_real_column(name);
                if (row >= col.size())
                    libmesh_error_msg(std::string("RealConstantsObject: column '") + name + "' has insufficient entries.");
                values.push_back(col[row]);
            }
            constants_map[ids[row]] = std::move(values);
        }
    }

    std::unordered_map<dof_id_type, std::vector<Real>> constants_map;
};

#endif // REAL_CONSTANTS_OBJECT_H
