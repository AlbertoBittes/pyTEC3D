#!/usr/bin/env python3
"""Generates registries.h: compile-time enums for the JSON registries.

Inputs (src/resources): domains.json, fields.json, element_types.json.
Output: one header with an enum, a name table and a lookup function for each registry:

  Domain          domains.json "domains"
  FieldId         fields.json (a field is a DOF unknown, a post-processed field or a grid array)
  ElementTypeId   element_types.json "element_types"
  MaterialProperty / RealConstant   the "mat." / "rc." parameters of the element types
  ParameterId     the parameters of all the element types: the declared ones, then the expressions
  LawId           the "nonlinear_laws" names of the element types

Each element type declares its own data:

  "constants":   { "mu0": "4 * pi * 1e-7", "ex[]": "RealVectorValue(1, 0, 0)" }
  "parameters":  { "mat.mur[][]": {}, "rc.turns": {"domain": ["STRANDED_COIL"]}, "me.Js": {}, "mn.AMPS": {} }
  "expressions": { "nu[][]": "(mu0 * mur).inverse()", "sigma[][]": {"value": "rho.inverse()", "domain": [...]} }

The type of a constant, a mat./rc. parameter or an expression is declared once, by the suffix of its key
(none Real, [] RealVectorValue, [][] RealTensorValue); me./mn. parameters take the type of their field.
Everywhere else a parameter is written by name. A "domain" list limits the domains where a parameter exists
in the element type. A name has the same meaning in every element type that has it, except for the value of
an expression.

Constants and expressions are C++ (libMesh Real, RealVectorValue, RealTensorValue; the std math functions
and libMesh::outer_product without their namespace). Each expression becomes an expression of
registries.h (evaluate_expression), so the compiler checks it.

Changing any of the JSON files means recompiling the solver module. Only the model data JSON is dynamic.
"""
import argparse
import json
import keyword
import re
import sys
from collections import OrderedDict

IDENTIFIER = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
# Names an expression can use without a namespace.
STD_FUNCTIONS = ("abs", "sqrt", "cbrt", "pow", "exp", "log", "sin", "cos", "tan", "atan2", "hypot", "min", "max")
LIBMESH_NAMES = ("Real", "RealVectorValue", "RealTensorValue", "pi", "outer_product")
RESERVED_NAMES = {"axisymmetric", "source", "std", "libMesh", *STD_FUNCTIONS, *LIBMESH_NAMES}
CPP_KEYWORDS = {
    "alignas", "alignof", "and", "asm", "auto", "bool", "break", "case", "catch", "char", "class", "const",
    "constexpr", "continue", "default", "delete", "do", "double", "else", "enum", "explicit", "export",
    "extern", "false", "float", "for", "friend", "goto", "if", "inline", "int", "long", "mutable", "namespace",
    "new", "noexcept", "not", "nullptr", "operator", "or", "private", "protected", "public", "register",
    "return", "short", "signed", "sizeof", "static", "struct", "switch", "template", "this", "throw", "true",
    "try", "typedef", "typename", "union", "unsigned", "using", "virtual", "void", "volatile", "while", "xor",
    "COUNT",
}
TYPE_BY_SUFFIX = {"": "Scalar", "[]": "Vector", "[][]": "Tensor"}
CPP_TYPE = {"Scalar": "Real", "Vector": "RealVectorValue", "Tensor": "RealTensorValue"}
MAX_DOMAINS = 64  # A set of domains is a 64-bit mask.


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
    match = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)(\[\]\[\]|\[\])?$", text)
    if not match:
        fail(f"{context}: invalid name '{text}'.")
    return match.group(1), match.group(2) or ""


def read_domain_mask(options, allowed_keys, domains, context):
    """Bit mask of the "domain" list of an entry; all the domains when there is none."""
    if not isinstance(options, dict):
        fail(f"{context} must be an object.")
    for key in options:
        if key not in allowed_keys:
            fail(f"{context}: unknown key '{key}'.")
    if "domain" not in options:
        return (1 << len(domains)) - 1
    names = options["domain"]
    if not isinstance(names, list) or not names or not all(isinstance(n, str) for n in names):
        fail(f"{context}: 'domain' must be a non-empty array of domain names.")
    mask = 0
    for name in names:
        if name not in domains:
            fail(f"{context}: unknown domain '{name}' (see domains.json).")
        mask |= 1 << domains.index(name)
    return mask


def read_declaration(declaration, fields, context):
    """'mat.mur[][]' -> {name, type, source[, array]}."""
    prefix, dot, rest = "".join(declaration.split()).partition(".")
    if not dot:
        fail(f"{context} must start with mat., rc., me. or mn.")
    if prefix in ("mat", "rc"):
        name, suffix = split_key(rest, context)
        return {"name": name, "type": TYPE_BY_SUFFIX[suffix], "source": "Material" if prefix == "mat" else "Constant"}
    if prefix in ("me", "mn"):
        if not IDENTIFIER.match(rest):
            fail(f"{context}: write the name of the field only, the type is the one of fields.json.")
        if rest not in fields:
            fail(f"{context}: '{rest}' is not a field of fields.json.")
        vector = fields[rest] == "vector"
        if prefix == "mn" and vector:
            fail(f"{context}: point data is a scalar array per node.")
        return {"name": rest, "type": "Vector" if vector else "Scalar", "source": "Cell" if prefix == "me" else "Point",
                "array": rest}
    fail(f"{context}: unknown source '{prefix}' (expected mat, rc, me or mn).")


def cpp_expression(entry, index, entries, context):
    """C++ text of an expression. The parameters it names, which must come before it, are read from `source`."""
    where = f"{context}, expression '{entry['name']}'"
    known = {other["name"]: (position, other) for position, other in enumerate(entries)}

    def translate(match):
        token = match.group(0)
        if token == "axisymmetric":
            return "source.axisymmetric()"
        if token not in known:
            return token
        position, other = known[token]
        if position >= index:
            fail(f"{where}: '{token}' is not defined before it.")
        if other["source"] == "Point":
            fail(f"{where}: '{token}' is point data and cannot be used in expressions.")
        if other["domains"] & entry["domains"] != entry["domains"]:
            fail(f"{where}: '{token}' does not exist in every domain of the expression.")
        return f"source.get_{other['type'].lower()}(ParameterId::{token})"

    # Not a member (.det) and not part of a number (1e-7).
    return re.sub(r"(?<![\w.])[A-Za-z_]\w*", translate, entry["text"])


def read_element_type(type_name, definition, domains, fields):
    context = f"element type '{type_name}'"

    constants = []
    for key, value in definition.get("constants", {}).items():
        name, suffix = split_key(key, f"{context}, constant '{key}'")
        if isinstance(value, bool) or not isinstance(value, (int, float, str)):
            fail(f"{context}, constant '{key}' must be a number or a C++ expression.")
        constants.append({"name": name, "type": TYPE_BY_SUFFIX[suffix], "text": value if isinstance(value, str) else repr(value)})

    entries = []
    declarations = definition.get("parameters", {})
    if not isinstance(declarations, dict):
        fail(f"{context}: 'parameters' must be an object.")
    for declaration, options in declarations.items():
        entry = read_declaration(declaration, fields, f"{context}, parameter '{declaration}'")
        entry["domains"] = read_domain_mask(options, {"domain"}, domains, f"{context}, parameter '{declaration}'")
        entries.append(entry)

    expressions = definition.get("expressions", {})
    if not isinstance(expressions, dict):
        fail(f"{context}: 'expressions' must be an object.")
    for key, value in expressions.items():
        name, suffix = split_key(key, f"{context}, expression '{key}'")
        options = value if isinstance(value, dict) else {"value": value}
        if not isinstance(options.get("value"), str):
            fail(f"{context}, expression '{key}' must be a string or an object with a string 'value'.")
        entries.append({"name": name, "type": TYPE_BY_SUFFIX[suffix], "source": "Expression", "text": options["value"],
                        "domains": read_domain_mask(options, {"value", "domain"}, domains, f"{context}, expression '{key}'")})

    for item in constants + entries:
        if item["name"] in RESERVED_NAMES:
            fail(f"{context}: '{item['name']}' is reserved.")
        check_enumerator(item["name"], context)
    unique([item["name"] for item in constants + entries], f"{context}: constants, parameters and expressions")
    for index, entry in enumerate(entries):
        if entry["source"] == "Expression":
            entry["code"] = cpp_expression(entry, index, entries, context)

    return {"dimension": definition.get("dimension"), "constants": constants, "entries": entries}


def read_registries(resource_dir):
    domains = load(f"{resource_dir}/domains.json").get("domains")
    if not isinstance(domains, list):
        fail("domains.json must contain a 'domains' array.")
    domains = unique(domains, "domains.json")
    if len(domains) > MAX_DOMAINS:
        fail(f"domains.json: at most {MAX_DOMAINS} domains.")

    fields = load(f"{resource_dir}/fields.json")
    for name, kind in fields.items():
        if kind not in ("scalar", "vector"):
            fail(f"fields.json: field '{name}' must be \"scalar\" or \"vector\".")

    root = load(f"{resource_dir}/element_types.json")
    definitions = root.get("element_types")
    if not isinstance(definitions, dict) or not definitions:
        fail("element_types.json must contain an 'element_types' object.")

    types = OrderedDict()
    for name, definition in definitions.items():
        if definition.get("dimension") not in (2, 3):
            fail(f"element type '{name}': dimension must be 2 or 3.")
        types[name] = read_element_type(name, definition, domains, fields)

    # A parameter name has one type and one source in every element type.
    by_name = OrderedDict()
    for type_name, element_type in types.items():
        for entry in element_type["entries"]:
            shared = {key: entry.get(key) for key in ("name", "type", "source", "array")}
            known = by_name.setdefault(entry["name"], shared)
            if known != shared:
                fail(f"element type '{type_name}': '{entry['name']}' is declared differently in another element type.")
    parameters = [p for p in by_name.values() if p["source"] != "Expression"]
    parameters += [p for p in by_name.values() if p["source"] == "Expression"]

    material_properties = [p for p in parameters if p["source"] == "Material"]
    real_constants = [p for p in parameters if p["source"] == "Constant"]
    for index, entry in enumerate(material_properties):
        entry["property"] = index
    for index, entry in enumerate(real_constants):
        entry["property"] = index

    # The domains where each parameter exists in each element type.
    for element_type in types.values():
        mask_of = {entry["name"]: entry["domains"] for entry in element_type["entries"]}
        element_type["masks"] = [mask_of.get(p["name"], 0) for p in parameters]

    # Nonlinear laws: names are global, a name has one type.
    laws = OrderedDict()
    for type_name, definition in definitions.items():
        for law, law_definition in definition.get("nonlinear_laws", {}).items():
            kind = law_definition.get("type")
            if kind != "bh_curve":
                fail(f"element type '{type_name}', law '{law}': the only law type is 'bh_curve'.")
            laws[law] = kind
            check_enumerator(law, f"law '{law}'")

    unique([p["name"] for p in parameters] + list(laws), "parameter and law names")

    for name in list(fields) + list(domains) + list(types):
        check_enumerator(name, f"registry name '{name}'")

    return {
        "domains": domains,
        "fields": fields,
        "types": types,
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


def emit_expressions(out, types):
    out.append("/**")
    out.append(" * @brief Value of the Expression parameter @p id of the element type @p type.")
    out.append(" *")
    out.append(" * Source gives axisymmetric() and the value of the parameters, get_scalar(), get_vector() and get_tensor()")
    out.append(" * by ParameterId. A parameter is read where the expression names it, so one in a branch that is not taken")
    out.append(" * is not read.")
    out.append(" */")
    out.append("template <typename Source>")
    out.append("inline ParameterValue evaluate_expression(ElementTypeId type, ParameterId id, [[maybe_unused]] Source &source)")
    out.append("{")
    for name in LIBMESH_NAMES:
        out.append(f"    using libMesh::{name};")
    for name in STD_FUNCTIONS:
        out.append(f"    using std::{name};")
    out.append("    switch (type)")
    out.append("    {")
    for type_name, element_type in types.items():
        expressions = [e for e in element_type["entries"] if e["source"] == "Expression"]
        if not expressions:
            continue
        out.append(f"    case ElementTypeId::{type_name}:")
        out.append("    {")
        for constant in element_type["constants"]:
            out.append(f"        [[maybe_unused]] const {CPP_TYPE[constant['type']]} {constant['name']} = {constant['text']};")
        out.append("        switch (id)")
        out.append("        {")
        for entry in expressions:
            out.append(f"        case ParameterId::{entry['name']}:")
            out.append(f"            // {' '.join(entry['text'].split())}")
            out.append(f"            return expression_support::expression_value<ParameterType::{entry['type']}>({entry['code']});")
        out.append("        default:")
        out.append("            break;")
        out.append("        }")
        out.append("        break;")
        out.append("    }")
    out.append("    default:")
    out.append("        break;")
    out.append("    }")
    out.append("    return ParameterValue();")
    out.append("}")
    out.append("")


def render(registries):
    out = [
        "// Generated by src/tools/generate_registries.py from src/resources/*.json. Do not edit.",
        "#pragma once",
        "",
        "#include <array>",
        "#include <cmath>",
        "#include <cstddef>",
        "#include <cstdint>",
        "#include <optional>",
        "#include <string_view>",
        "",
        '#include "expression_support.h"',
        '#include "parameter_value.h"',
        "",
        "namespace pytec3d::registry",
        "{",
        "",
    ]

    emit_enum(out, "Domain", registries["domains"], "Physical regions (domains.json).")
    out.append("/// A set of domains: bit i is the Domain with index i.")
    out.append("using DomainMask = std::uint64_t;")
    out.append("inline constexpr DomainMask bit(Domain domain) { return DomainMask(1) << static_cast<unsigned int>(domain); }")
    out.append("")

    emit_enum(out, "FieldId", list(registries["fields"]), "Fields and grid arrays (fields.json).")
    out.append("inline constexpr std::array<bool, kNumFieldId> kFieldIsVector = {")
    for kind in registries["fields"].values():
        out.append(f"    {'true' if kind == 'vector' else 'false'},")
    out.append("};")
    out.append("inline constexpr bool is_vector(FieldId id) { return kFieldIsVector[static_cast<std::size_t>(id)]; }")
    out.append("")

    types = registries["types"]
    emit_enum(out, "ElementTypeId", list(types), "Element types (element_types.json).")
    out.append("inline constexpr std::array<unsigned int, kNumElementTypeId> kElementTypeDimensions = {")
    for element_type in types.values():
        out.append(f"    {element_type['dimension']},")
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
    emit_enum(out, "ParameterId", [p["name"] for p in parameters], "Parameters of all the element types: declared first, then expressions.")
    out.append("struct ParameterInfo")
    out.append("{")
    out.append("    ParameterType type;")
    out.append("    ParameterSource source;")
    out.append("    unsigned int property; ///< MaterialProperty or RealConstant, by source.")
    out.append("    FieldId array;         ///< Grid array of Cell and Point parameters.")
    out.append("};")
    out.append(f"inline constexpr std::array<ParameterInfo, {len(parameters)}> kParameters = {{{{")
    for entry in parameters:
        array = f"FieldId::{entry['array']}" if entry.get("array") else "FieldId::COUNT"
        out.append(f"    {{ParameterType::{entry['type']}, ParameterSource::{entry['source']}, {entry.get('property', 0)}u, {array}}},")
    out.append("}};")
    out.append("inline constexpr const ParameterInfo &info(ParameterId id) { return kParameters[static_cast<std::size_t>(id)]; }")
    out.append("")

    out.append("/// The domains where each parameter exists in each element type; none if the type does not have it.")
    out.append(f"inline constexpr std::array<std::array<DomainMask, kNumParameterId>, kNumElementTypeId> kParameterDomains = {{{{")
    for type_name, element_type in types.items():
        out.append(f"    {{{{ // {type_name}")
        for parameter, mask in zip(parameters, element_type["masks"]):
            out.append(f"        {mask}ull, // {parameter['name']}")
        out.append("    }},")
    out.append("}};")
    out.append("inline constexpr DomainMask parameter_domains(ElementTypeId type, ParameterId id)")
    out.append("{")
    out.append("    return kParameterDomains[static_cast<std::size_t>(type)][static_cast<std::size_t>(id)];")
    out.append("}")
    out.append("inline constexpr bool defined(ElementTypeId type, ParameterId id) { return parameter_domains(type, id) != 0; }")
    out.append("inline constexpr bool defined(ElementTypeId type, ParameterId id, Domain domain)")
    out.append("{")
    out.append("    return domain != Domain::COUNT && (parameter_domains(type, id) & bit(domain)) != 0;")
    out.append("}")
    out.append("")

    emit_expressions(out, types)

    out.append("enum class LawKind : unsigned char")
    out.append("{")
    out.append("    BHCurve, ///< H(B) of the BH curve of the element material.")
    out.append("};")
    emit_enum(out, "LawId", list(registries["laws"]), 'Nonlinear laws ("nonlinear_laws"), used in kernels like a parameter.')
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
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
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
