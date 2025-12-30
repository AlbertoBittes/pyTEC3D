#ifndef ND_ARRAY_H
#define ND_ARRAY_H

#include <cstdint>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <libmesh/libmesh.h>

class NDArray
{
public:
    NDArray();
    NDArray(std::vector<std::uint32_t> s,
            std::variant<std::vector<int>, std::vector<libMesh::Real>, std::vector<std::string>> v);

    std::size_t flattened_index(const std::vector<std::uint32_t> &indices) const;

    template <typename... Indices>
    auto operator()(Indices... idxs) const -> std::variant<int, libMesh::Real, std::string>
    {
        std::vector<std::uint32_t> indices{static_cast<std::uint32_t>(idxs)...};
        auto flat = flattened_index(indices);

        if (std::holds_alternative<std::vector<int>>(values))
            return std::get<std::vector<int>>(values)[flat];
        if (std::holds_alternative<std::vector<libMesh::Real>>(values))
            return std::get<std::vector<libMesh::Real>>(values)[flat];
        if (std::holds_alternative<std::vector<std::string>>(values))
            return std::get<std::vector<std::string>>(values)[flat];

        libmesh_error_msg("NDArray: unsupported value type");
        return std::variant<int, libMesh::Real, std::string>{};
    }

    template <typename T>
    const std::vector<T> &view() const
    {
        if (!std::holds_alternative<std::vector<T>>(values))
            libmesh_error_msg("NDArray: requested view type does not match stored type.");
        return std::get<std::vector<T>>(values);
    }

    static std::string format_shape(const std::vector<std::uint32_t> &shape);
    template <typename Vec>
    static void print_vector_contents(const Vec &vec)
    {
        constexpr std::size_t max_preview = 8;
        libMesh::out << "[";
        const std::size_t limit = std::min<std::size_t>(vec.size(), max_preview);
        for (std::size_t i = 0; i < limit; ++i)
        {
            if (i)
                libMesh::out << ", ";
            if constexpr (std::is_same_v<typename Vec::value_type, std::string>)
                libMesh::out << '"' << vec[i] << '"';
            else
                libMesh::out << vec[i];
        }
        if (vec.size() > max_preview)
            libMesh::out << ", ... (" << vec.size() << " total)";
        libMesh::out << "]";
    }

    static std::size_t infer_entry_span(const std::vector<std::uint32_t> &shape, std::size_t value_count);

    template <typename Vec>
    static std::size_t vector_size(const Vec &vec)
    {
        return vec.size();
    }

    std::vector<std::uint32_t> shape;
    std::variant<std::vector<int>, std::vector<libMesh::Real>, std::vector<std::string>> values;
};

#endif // ND_ARRAY_H
