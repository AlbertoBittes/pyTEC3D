/**
 * @file element_type_parser.h
 * @brief Builds ElementType objects from element_types.json, with the registries of registries.h.
 *
 * The names of fields, domains, element types, parameters and laws are enums generated at build time
 * from fields.json, domains.json and element_types.json (src/tools/generate_registries.py). What is read
 * here, at run time, is the structure of each element type: the operand expressions such as
 * "grad(dt(V))", the coefficient names such as "-sigma[][]" and the field options. The result is an
 * ElementType made of plain data, with the kernels instantiated for the requested options.
 *
 * element_types.json layout:
 *
 *   "parameters":  [ "mat.mur[]", "rc.turns", "me.Js[]", "mn.AMPS()", ... ]     (registry: ParameterId)
 *       Data the kernels can use: mat. material property, rc. real constant, me. cell-data array of the
 *       grid, mn. point-data array. The type is the suffix: none scalar, [] vector, [][] tensor, () as
 *       declared in fields.json. A material property or real constant is a key of the model data JSON.
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
 *
 * A kernel with a "domain" list only applies in those domains, and only there are the parameters of its
 * coefficient needed: a material property or real constant that no active kernel uses can stay undefined.
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "element_type.h"
#include "registries.h"

namespace pytec3d
{
    class ElementTypeParser
    {
    public:
        /** @brief Reads the JSON compiled into the module. */
        static ElementTypeParser from_embedded();

        /** @brief Reads raw element_types.json text (aborts on malformed input). */
        static ElementTypeParser from_string(const std::string &element_types_json);

        /** @brief Instantiates the type @p id for @p options; aborts on any inconsistency. */
        std::unique_ptr<ElementType> build(registry::ElementTypeId id, const ElementTypeOptions &options) const;

    private:
        ElementTypeParser() = default;

        nlohmann::ordered_json element_types_json_;
    };
} // namespace pytec3d
