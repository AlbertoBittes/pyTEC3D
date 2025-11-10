#include "nd_array.h"

#include <sstream>

namespace Globals
{

NDArray::NDArray() = default;

NDArray::NDArray(std::vector<std::uint32_t> s,
                 std::variant<std::vector<int>, std::vector<libMesh::Real>, std::vector<std::string>> v)
    : shape(std::move(s)),
      values(std::move(v))
{
}

std::size_t NDArray::flattened_index(const std::vector<std::uint32_t> &indices) const
{
    if (indices.size() != shape.size())
        libmesh_error_msg("NDArray: index dimensionality mismatch");

    std::size_t offset = 0;
    std::size_t stride = 1;
    for (std::size_t dim = shape.size(); dim-- > 0;)
    {
        if (indices[dim] >= shape[dim])
            libmesh_error_msg("NDArray: index out of bounds");
        offset += indices[dim] * stride;
        stride *= shape[dim];
    }
    return offset;
}

std::string NDArray::format_shape(const std::vector<std::uint32_t> &shape)
{
    if (shape.empty())
        return "[]";

    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < shape.size(); ++i)
    {
        if (i)
            oss << ", ";
        oss << shape[i];
    }
    oss << "]";
    return oss.str();
}

std::size_t NDArray::infer_entry_span(const std::vector<std::uint32_t> &shape, std::size_t value_count)
{
    if (value_count == 0)
        return 0;

    if (shape.empty())
        return value_count;

    for (auto it = shape.rbegin(); it != shape.rend(); ++it)
    {
        if (*it == 0)
            return 0;
        if (*it > 1)
            return *it;
    }

    return value_count;
}

} // namespace Globals
