#!/usr/bin/env python3
"""Generates registries.h: compile-time enums for the JSON registries.

Inputs (src/resources): domains.json, fields.json, element_types.json.
Output: one header with an enum, a name table and a lookup function for each registry:

  Domain          domains.json "domains"
  FieldId         fields.json (a field is a DOF unknown, a post-processed field or a grid array)
  ElementTypeId   element_types.json "element_types"
  MaterialProperty / RealConstant   the "mat." / "rc." entries of "parameters"
  ParameterId     every parameter: the declared ones, then the expressions, in file order
  LawId           the "nonlinear_laws" names of the element types

Changing any of the JSON files means recompiling the solver module. Only the model data JSON is dynamic.
"""
import argparse
import json
import keyword
import re
import sys
from collections import OrderedDict

IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
RESERVED_PARAMETER_NAMES = {"pi", "mu0", "dimension", "axisymmetric", "ex", "ey", "ez"}
CPP_KEYWORDS = {
    "alignas", "alignof", "and", "asm", "auto", "bool", "break", "case", "catch", "char", "class", "const",
    "constexpr", "continue", "default", "delete", "do", "double", "else", "enum", "explicit", "export",
    "extern", "false", "float", "for", "friend", "goto", "if", "inline", "int", "long", "mutable", "namespace",
    "new", "noexcept", "not", "nullptr", "operator", "or", "private", "protected", "public", "register",
    "return", "short", "signed", "sizeof", "static", "struct", "switch", "template", "this", "throw", "true",
    "try", "typedef", "typename", "union", "unsigned", "using", "virtual", "void", "volatile", "while", "xor",
    "COUNT",
}
TYPE_BY_SUFFIX = {"": "Scalar", "[]": "Vector", "[][]": "Tensor", "()": None}


def fail(message):
    sys.exit(f"generate_registries.py: {message}")


def load(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=OrderedDict)


def check_enumerator(name, context):
    if not IDENTIFIER.match(name) or keyword.iskeyword(name) or name in CPP_KEYWORDS:
        fail(f"{context}: '{name}' cannot be used as a C++ enumerator.")


def unique(names, context):
    seen = set()
    for name in names:
        if name in seen:
            fail(f"{context}: duplicate entry '{name}'.")
        seen.add(name)
    return list(names)


def split_key(text, context):
    """'nu[][]' -> ('nu', '[][]')."""
    text = "".join(text.split())
    match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)(\[\]\[\]|\[\]|\(\))?$", text)
    if not match:
        fail(f"{context}: invalid parameter '{text}'.")
    return match.group(1), match.group(2) or ""


def read_registries(resource_dir):
    domains = load(f"{resource_dir}/domains.json").get("domains")
    if not isinstance(domains, list):
        fail("domains.json must contain a 'domains' array.")
    domains = unique(domains, "domains.json")

    fields = load(f"{resource_dir}/fields.json")
    for name, kind in fields.items():
        if kind not in ("scalar", "vector"):
            fail(f"fields.json: field '{name}' must be \"scalar\" or \"vector\".")

    root = load(f"{resource_dir}/element_types.json")
    element_types = root.get("element_types")
    if not isinstance(element_types, dict) or not element_types:
        fail("element_types.json must contain an 'element_types' object.")
    for name, definition in element_types.items():
        if definition.get("dimension") not in (2, 3):
            fail(f"element type '{name}': dimension must be 2 or 3.")

    # Parameters: declared first, then the expressions.
    parameters = []  # dicts: key, name, type, source, property, array, expression
    for declaration in root.get("parameters", []):
        context = f"parameter '{declaration}'"
        prefix, dot, rest = "".join(declaration.split()).partition(".")
        if not dot:
            fail(f"{context} must start with mat., rc., me. or mn.")
        name, suffix = split_key(rest, context)
        entry = {"key": name + suffix, "name": name, "type": TYPE_BY_SUFFIX[suffix], "expression": ""}
        if prefix in ("mat", "rc"):
            if suffix == "()":
                fail(f"{context}: a {'material' if prefix == 'mat' else 'real-constant'} parameter has a fixed type.")
            entry["source"] = "Material" if prefix == "mat" else "Constant"
        elif prefix in ("me", "mn"):
            if name not in fields:
                fail(f"{context}: '{name}' is not a field of fields.json.")
            entry["source"] = "Cell" if prefix == "me" else "Point"
            entry["array"] = name
            vector = fields[name] == "vector"
            if suffix == "()":
                entry["type"] = "Vector" if vector else "Scalar"
            elif (entry["type"] == "Vector") != vector or entry["type"] == "Tensor":
                fail(f"{context}: field '{name}' is a {fields[name]}.")
            if prefix == "mn" and vector:
                fail(f"{context}: point data is a scalar array per node.")
        else:
            fail(f"{context}: unknown source '{prefix}' (expected mat, rc, me or mn).")
        parameters.append(entry)

    for key, text in root.get("expressions", {}).items():
        name, suffix = split_key(key, f"expression '{key}'")
        if suffix == "()":
            fail(f"expression '{key}': an expression is a scalar, a vector[] or a tensor[][].")
        if not isinstance(text, str):
            fail(f"expression '{key}' must be a string.")
        parameters.append({"key": name + suffix, "name": name, "type": TYPE_BY_SUFFIX[suffix],
                           "source": "Expression", "expression": text})

    unique([p["key"] for p in parameters], "parameters and expressions")
    for entry in parameters:
        if entry["name"] in RESERVED_PARAMETER_NAMES:
            fail(f"parameter '{entry['key']}': '{entry['name']}' is reserved.")
        check_enumerator(entry["name"], f"parameter '{entry['key']}'")
    unique([p["name"] for p in parameters], "parameter names (a name can be used with one type suffix only)")

    material_properties = [p for p in parameters if p["source"] == "Material"]
    real_constants = [p for p in parameters if p["source"] == "Constant"]
    for index, entry in enumerate(material_properties):
        entry["property"] = index
    for index, entry in enumerate(real_constants):
        entry["property"] = index

    # Nonlinear laws: names are global, a name has one type.
    laws = OrderedDict()
    for type_name, definition in element_types.items():
        for law, law_definition in definition.get("nonlinear_laws", {}).items():
            kind = law_definition.get("type")
            if kind != "bh_curve":
                fail(f"element type '{type_name}', law '{law}': the only law type is 'bh_curve'.")
            laws[law] = kind
            check_enumerator(law, f"law '{law}'")

    unique([p["name"] for p in parameters] + list(laws), "parameter and law names")

    for name in list(fields) + list(domains) + list(element_types):
        check_enumerator(name, f"registry name '{name}'")

    return {
        "domains": domains,
        "fields": fields,
        "element_types": OrderedDict((name, d["dimension"]) for name, d in element_types.items()),
        "parameters": parameters,
        "material_properties": material_properties,
        "real_constants": real_constants,
        "laws": laws,
    }


def cpp_string(text):
    return json.dumps(text)


def emit_enum(out, enum, names, doc):
    out.append(f"/// {doc}")
    out.append(f"enum class {enum} : unsigned char")
    out.append("{")
    for name in names:
        out.append(f"    {name},")
    out.append("    COUNT")
    out.append("};")
    out.append(f"inline constexpr std::size_t kNum{enum} = {len(names)};")
    out.append("")
    out.append(f"inline constexpr std::array<std::string_view, {len(names)}> k{enum}Names = {{")
    for name in names:
        out.append(f"    {cpp_string(name)},")
    out.append("};")
    out.append(f"inline constexpr std::string_view name({enum} id) {{ return k{enum}Names[static_cast<std::size_t>(id)]; }}")
    out.append(f"inline std::optional<{enum}> parse_{enum}(std::string_view text)")
    out.append("{")
    out.append(f"    for (std::size_t i = 0; i < k{enum}Names.size(); ++i)")
    out.append(f"        if (k{enum}Names[i] == text)")
    out.append(f"            return static_cast<{enum}>(i);")
    out.append("    return std::nullopt;")
    out.append("}")
    out.append("")


def render(registries):
    out = [
        "// Generated by src/tools/generate_registries.py from src/resources/*.json. Do not edit.",
        "#pragma once",
        "",
        "#include <array>",
        "#include <cstddef>",
        "#include <optional>",
        "#include <string_view>",
        "",
        '#include "parameter_value.h"',
        "",
        "namespace pytec3d::registry",
        "{",
        "",
    ]

    emit_enum(out, "Domain", registries["domains"], "Physical regions (domains.json).")
    emit_enum(out, "FieldId", list(registries["fields"]), "Fields and grid arrays (fields.json).")
    out.append("inline constexpr std::array<bool, kNumFieldId> kFieldIsVector = {")
    for kind in registries["fields"].values():
        out.append(f"    {'true' if kind == 'vector' else 'false'},")
    out.append("};")
    out.append("inline constexpr bool is_vector(FieldId id) { return kFieldIsVector[static_cast<std::size_t>(id)]; }")
    out.append("")

    names = list(registries["element_types"])
    emit_enum(out, "ElementTypeId", names, "Element types (element_types.json).")
    out.append("inline constexpr std::array<unsigned int, kNumElementTypeId> kElementTypeDimensions = {")
    for dimension in registries["element_types"].values():
        out.append(f"    {dimension},")
    out.append("};")
    out.append("inline constexpr unsigned int dimension(ElementTypeId id) { return kElementTypeDimensions[static_cast<std::size_t>(id)]; }")
    out.append("")

    # The material and real-constant properties are a model data key each.
    for enum, key, doc in (("MaterialProperty", "material_properties", 'Material properties ("mat." parameters).'),
                           ("RealConstant", "real_constants", 'Real constants ("rc." parameters).')):
        entries = registries[key]
        emit_enum(out, enum, [e["name"] for e in entries], doc)
        out.append(f"inline constexpr std::array<ParameterType, {len(entries)}> k{enum}Types = {{")
        for entry in entries:
            out.append(f"    ParameterType::{entry['type']},")
        out.append("};")
        out.append(f"inline constexpr ParameterType type({enum} id) {{ return k{enum}Types[static_cast<std::size_t>(id)]; }}")
        out.append("")

    out.append("enum class ParameterSource : unsigned char")
    out.append("{")
    out.append("    Material,   ///< Property of the material (MaterialProperty).")
    out.append("    Constant,   ///< Real constant (RealConstant).")
    out.append("    Cell,       ///< Cell-data array of the grid, named like a field.")
    out.append("    Point,      ///< Point-data array of the grid; only usable as the nodal load of a DOF field.")
    out.append("    Expression, ///< Derived from other parameters.")
    out.append("};")
    out.append("")

    parameters = registries["parameters"]
    emit_enum(out, "ParameterId", [p["name"] for p in parameters], "Parameters: declared first, then expressions.")
    out.append("struct ParameterInfo")
    out.append("{")
    out.append("    std::string_view key; ///< As written, with its type suffix: nu[][].")
    out.append("    ParameterType type;")
    out.append("    ParameterSource source;")
    out.append("    unsigned int property;       ///< MaterialProperty or RealConstant, by source.")
    out.append("    FieldId array;               ///< Grid array of Cell and Point parameters.")
    out.append("    std::string_view expression; ///< Text of an Expression parameter.")
    out.append("};")
    out.append(f"inline constexpr std::array<ParameterInfo, {len(parameters)}> kParameters = {{{{")
    for entry in parameters:
        array = f"FieldId::{entry['array']}" if "array" in entry else "FieldId::COUNT"
        out.append(f"    {{{cpp_string(entry['key'])}, ParameterType::{entry['type']}, ParameterSource::{entry['source']}, "
                   f"{entry.get('property', 0)}u, {array}, {cpp_string(entry['expression'])}}},")
    out.append("}};")
    out.append("inline constexpr const ParameterInfo &info(ParameterId id) { return kParameters[static_cast<std::size_t>(id)]; }")
    out.append("/// Parameter written as in the JSON files, with its type suffix.")
    out.append("inline std::optional<ParameterId> parse_parameter(std::string_view key)")
    out.append("{")
    out.append("    for (std::size_t i = 0; i < kParameters.size(); ++i)")
    out.append("        if (kParameters[i].key == key)")
    out.append("            return static_cast<ParameterId>(i);")
    out.append("    return std::nullopt;")
    out.append("}")
    out.append("")

    out.append("enum class LawKind : unsigned char")
    out.append("{")
    out.append("    BHCurve, ///< H(B) of the BH curve of the element material.")
    out.append("};")
    emit_enum(out, "LawId", list(registries["laws"]), 'Nonlinear laws ("nonlinear_laws"), used in kernels as name().')
    out.append(f"inline constexpr std::array<LawKind, {len(registries['laws'])}> kLawKinds = {{")
    for _ in registries["laws"]:
        out.append("    LawKind::BHCurve,")
    out.append("};")
    out.append("inline constexpr LawKind kind(LawId id) { return kLawKinds[static_cast<std::size_t>(id)]; }")
    out.append("")
    out.append("} // namespace pytec3d::registry")
    out.append("")
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resources", required=True, help="Directory with the JSON registries.")
    parser.add_argument("--output", required=True, help="Header to write.")
    args = parser.parse_args()

    text = render(read_registries(args.resources))
    try:
        with open(args.output, encoding="utf-8") as handle:
            if handle.read() == text:
                return  # Leave the timestamp alone: nothing recompiles.
    except FileNotFoundError:
        pass
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(text)


if __name__ == "__main__":
    main()
