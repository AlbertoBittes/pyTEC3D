/**
 * @file element_type_parser.h
 * @brief Builds ElementType objects from element_types.json, fields.json and domains.json.
 *
 * All the text of the JSON files is interpreted here, and nowhere else: operand expressions such as
 * "grad(dt(V))", coefficient names such as "-sigma[][]", field definitions and domain names. The result is
 * an ElementType made of plain data, with the kernels instantiated for the requested options.
 *
 * element_types.json layout:
 *
 *   "parameters":  [ "mat.mur[]", "rc.turns", "me.Js[]", "mn.AMPS()", ... ]
 *       Data the kernels can use: mat. material row, rc. real constants, me. cell-data array of the grid,
 *       mn. point-data array. The type is the suffix: none scalar, [] vector, [][] tensor, () as in the grid.
 *   "expressions": { "nu[][]": "global(diag(inv(mu0 * mur[])))", ... }
 *       Derived parameters, in dependency order (see expression.h); usable in kernels and later expressions.
 *   "element_types": { "<name>": {
 *       "dimension": 2 | 3,
 *       "nonlinear_laws":  { "BH": {"type": "bh_curve", "linear": "nu[][]"} },     // used as "BH()"
 *       "dof_fields":      { "A": {"active_components": ["z"], "finite_elements": "...", "nodal_loads": "AMPS()", "domain": [...]} },
 *       "dof_kernels":     { "Kaa": {"trial": "curl(A)", "coeff": "BH()", "test": "curl(A)", "domain": [...]},
 *                            "Fs":  {"coeff": "Js_global[]", "test": "A"} },        // no trial: source
 *       "nodal_fields":    { "B": {} },
 *       "nodal_kernels":   { "B": [ {"trial": "curl(A)", "coeff": "...", "domain": [...]}, ... ] },
 *       "element_fields":  { "Jt": {"active_components": ["z"]} },
 *       "element_kernels": { "Jt": [ ... ] } } }
 */
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "element_type.h"

namespace pytec3d
{
    class ElementTypeParser
    {
    public:
        /** @brief Reads the JSON compiled into the module. */
        static ElementTypeParser from_embedded();

        /** @brief Reads raw JSON text (aborts on malformed input). */
        static ElementTypeParser from_strings(const std::string &element_types_json,
                                              const std::string &fields_json,
                                              const std::string &domains_json);

        bool has_type(const std::string &name) const;
        std::vector<std::string> type_names() const;

        const std::vector<std::string> &domains() const { return domains_; }
        unsigned int domain_index(const std::string &name) const;
        const std::string &domain_name(unsigned int index) const;

        /** @brief Instantiates the type @p name for @p options; aborts on any inconsistency. */
        std::unique_ptr<ElementType> build(const std::string &name, const ElementTypeOptions &options) const;

    private:
        ElementTypeParser() = default;

        nlohmann::ordered_json element_types_json_;
        std::vector<ParameterDefinition> declarations_;                    ///< Top-level "parameters".
        std::vector<std::pair<std::string, std::string>> expressions_; ///< Top-level "expressions", in order.
        std::map<std::string, bool> field_is_vector_;
        std::vector<std::string> domains_;
    };
} // namespace pytec3d
