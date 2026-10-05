/**
 * @file element_type_registry.h
 * @brief JSON-defined base element types, physical fields and domains.
 *
 * The registry is built once when the solver starts. Users then instantiate concrete
 * element types (see element_type_instance.h) by naming a base type plus options.
 */
#pragma once

#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace pytec3d
{
    enum class FieldKind : unsigned char
    {
        Dof,     ///< Unknown solved for by the system.
        Nodal,   ///< Post-processed nodal field.
        Element, ///< Post-processed element (cell) field.
    };

    /** @brief One physical field of a base element type. */
    struct FieldSpec
    {
        std::string name;
        FieldKind kind = FieldKind::Dof;
        bool is_vector = false;
        std::vector<unsigned int> active_components; ///< Vector fields only; 0=x, 1=y, 2=z.
        std::vector<std::string> domains;            ///< Empty means active in every domain.
        std::string finite_elements;                 ///< libMesh FE family name overriding the default; empty keeps the default.
        std::string nodal_loads;                     ///< Name of the per-node VTK array loading this DOF field (e.g. "AMPS").
        std::string expression;                      ///< Nodal/element fields: sum of terms evaluated from the DOF fields.
    };

    /** @brief One weak-form kernel: integral of test^T * coeff * trial, e.g. {"curl(A)", "nu[][]", "curl(A)"}. Without a trial it is a source: integral of test^T * coeff with a vector coefficient. */
    struct KernelSpec
    {
        std::string name;
        std::string trial;
        std::string coeff; ///< Optional leading '-' followed by a name from BaseElementType::coefficients.
        std::string test;
        std::vector<std::string> domains; ///< Empty means active in every domain.
    };

    /** @brief A base element type as defined in element_types.json. */
    struct BaseElementType
    {
        std::string name;
        unsigned int dimension = 0;
        unsigned int time_order = 0;
        std::vector<std::string> coefficients; ///< Coefficient names usable in kernels and field expressions.
        std::vector<KernelSpec> kernels;
        std::vector<FieldSpec> fields; ///< JSON order: dof, nodal, element.
    };

    class ElementTypeRegistry
    {
    public:
        /** @brief Builds the registry from the JSON compiled into the module. */
        static ElementTypeRegistry from_embedded();

        /** @brief Builds the registry from raw JSON text (aborts on malformed input). */
        static ElementTypeRegistry from_strings(const std::string &element_types_json,
                                                const std::string &fields_json,
                                                const std::string &domains_json);

        /** @brief Raw element_types.json content. */
        const nlohmann::ordered_json &element_types_json() const { return element_types_json_; }

        bool has_type(const std::string &name) const;
        const BaseElementType &get_type(const std::string &name) const;
        std::vector<std::string> type_names() const;

        const std::vector<std::string> &domains() const { return domains_; }
        unsigned int domain_index(const std::string &name) const;
        const std::string &domain_name(unsigned int index) const;

    private:
        ElementTypeRegistry() = default;

        nlohmann::ordered_json element_types_json_;
        std::map<std::string, bool> field_is_vector_;
        std::vector<std::string> domains_;
        std::map<std::string, BaseElementType> types_;
    };
} // namespace pytec3d
