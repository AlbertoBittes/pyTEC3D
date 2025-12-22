#ifndef BINARY_LABELS_H
#define BINARY_LABELS_H

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>

#include <libmesh/libmesh.h>

enum class BinaryLabelId : std::uint32_t
{
    ELEMENT_TYPES = 0,
    TIME_PARAMETERS = 1,
    SOLVER_SETTINGS = 2,
    MATERIALS = 3,
    REAL_CONSTANTS = 4,
    POINT_DATA = 5,
    CELL_DATA = 6,
    PERIODIC_BOUNDARIES = 7,
    DOF_COUPLINGS = 8,
    APEC_SETTINGS = 9,
    DOF_INTERPOLATION_SETTINGS = 10,
    UNKNOWN = std::numeric_limits<std::uint32_t>::max()
};

struct BinaryLabelTraits
{
    static constexpr std::array<std::pair<BinaryLabelId, std::string_view>, 11> LABEL_TABLE = {{
        {BinaryLabelId::ELEMENT_TYPES, "element_types"},
        {BinaryLabelId::TIME_PARAMETERS, "time_parameters"},
        {BinaryLabelId::SOLVER_SETTINGS, "solver_settings"},
        {BinaryLabelId::MATERIALS, "materials"},
        {BinaryLabelId::REAL_CONSTANTS, "real_constants"},
        {BinaryLabelId::POINT_DATA, "node"},
        {BinaryLabelId::CELL_DATA, "elem"},
        {BinaryLabelId::PERIODIC_BOUNDARIES, "periodic_boundaries"},
        {BinaryLabelId::DOF_COUPLINGS, "dof_couplings"},
        {BinaryLabelId::APEC_SETTINGS, "apec_settings"},
        {BinaryLabelId::DOF_INTERPOLATION_SETTINGS, "dof_interpolation_settings"},
    }};

    static constexpr std::string_view to_string(BinaryLabelId id)
    {
        for (const auto &[label, name] : LABEL_TABLE)
            if (label == id)
                return name;
        return "unknown";
    }

    static BinaryLabelId from_raw(std::uint32_t raw)
    {
        for (const auto &[label, _] : LABEL_TABLE)
            if (static_cast<std::uint32_t>(label) == raw)
                return label;

        libmesh_error_msg("BinaryLabelTraits: unknown binary label id " + std::to_string(raw));
        return BinaryLabelId::UNKNOWN;
    }

    static bool is_vtk_label(BinaryLabelId id)
    {
        return id == BinaryLabelId::POINT_DATA || id == BinaryLabelId::CELL_DATA;
    }

    static bool is_point_data(BinaryLabelId id) { return id == BinaryLabelId::POINT_DATA; }
    static bool is_cell_data(BinaryLabelId id) { return id == BinaryLabelId::CELL_DATA; }
};

#endif // BINARY_LABELS_H
