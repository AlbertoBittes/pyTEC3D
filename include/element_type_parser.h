/**
 * @file element_type_parser.h
 * @brief Builds ElementType objects from element_types.json, with the registries of registries.h.
 *
 * The names of fields, domains, element types, parameters and laws are enums generated at build time
 * from fields.json, domains.json and element_types.json (src/tools/generate_registries.py). What is read
 * here, at run time, is the structure of each element type: the operand expressions such as
 * "grad(dt(V))", the coefficient names such as "-sigma" and the field options. The result is an
 * ElementType made of plain data, with the kernels instantiated for the requested options.
 *
 * element_types.json layout, "element_types": { "<name>": {
 *       "dimension": 2 | 3,
 *       "constants":   { "mu0": "4 * pi * 1e-7", "ex[]": "RealVectorValue(1, 0, 0)" },
 *       "parameters":  { "mat.mur[][]": {}, "rc.turns": {"domain": ["STRANDED_COIL"]}, "me.Js": {}, "mn.AMPS": {} },
 *       "expressions": { "nu[][]": "(mu0 * mur).inverse()", "kvv": {"value": "...", "domain": [...]} },
 *     The three are read by generate_registries.py, not here (see its documentation). Constants and
 *     expressions are C++; parameters are data the kernels can use: mat. material property, rc. real
 *     constant (keys of the model data JSON), me. cell-data array, mn. point-data array. A "domain" list
 *     limits the domains where a parameter exists in the element type. Kernels name parameters without
 *     suffix.
 *       "nonlinear_laws":  { "BH": {"type": "bh_curve", "linear": "nu"} },         // used as "BH"
 *       "dof_fields":      { "A": {"active_components": ["z"], "finite_elements": "...", "nodal_loads": "AMPS", "domain": [...]} },
 *       "dof_kernels":     { "Kaa": {"trial": "curl(A)", "coeff": "BH", "test": "curl(A)", "domain": [...]},
 *                            "Fs":  {"coeff": "Js", "test": "A"} },                 // no trial: source
 *       "nodal_fields":    { "B": {} },
 *       "nodal_kernels":   { "B": [ {"trial": "curl(A)", "coeff": "...", "domain": [...]}, ... ] },
 *       "element_fields":  { "Jt": {"active_components": ["z"]} },
 *       "element_kernels": { "Jt": [ ... ] } } }
 *
 * A kernel with a "domain" list only applies in those domains, and only there are the parameters of its
 * coefficient needed: a material property or real constant that no active kernel uses can stay undefined.
 * A kernel must not be active in a domain where its parameter does not exist.
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
