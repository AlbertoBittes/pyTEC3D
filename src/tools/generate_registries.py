#!/usr/bin/env python3
"""Compiles the JSON registries into C++ headers: enums, compiled expressions and element type tables.

Inputs (src/pytec3d/resources): domains.json, fields.json, element_types.json, nonlinear_laws.json.
Outputs:
  registries.h     enums with name tables and lookups (Domain, FieldId, VariableId, ElementTypeId,
                   MaterialProperty, RealConstant, ParameterId, CurveId, LawId, LawOutputId), the
                   parameter tables and evaluate_expression().
  registries.cpp   the definitions of the lookups that are not constexpr (parse_*(), default_value()).
  formulations.h   one formulation::ElementTypeSpec per element type (formulation_spec.h) and evaluate_law().

Everything is checked here, so a formulation that builds is consistent: names, types of every product,
domains of the parameters, operators against the finite elements. Changing any JSON file means rebuilding the
solver module; only the model data JSON is read at run time.

element_types.json, "element_types": { "<name>": {

  "inherits": "<element type>"
      Starts from another element type: each section given here replaces the parent's entries key by key
      (a dof_kernels "K" list, an element_kernels "Jt" list, a parameter...); "dimension" replaces it.
  "dimension": 2 | 3,
  "time_order": 0 | 1 | 2      highest time derivative of the formulation (default 1): at most that many dt().
  "out_of_plane_length": "thickness"
      2D only: the scalar parameter giving the out-of-plane length L of the elements (times 2 pi r when
      axisymmetric); L = 1 without it. It scales the integration weights and the out-of-plane gradient.

  "constants":   { "mu0": "4 * pi * 1e-7", "ex[]": "RealVectorValue(1, 0, 0)" }
  "parameters":  { "mat.mur[][]": {}, "rc.turns": {"domain": ["STRANDED_COIL"]}, "me.Js": {}, "mn.AMPS": {},
                   "mat.M[]": {"default": [0, 0, 0]} }
  "expressions": { "nu[][]": "(mu0 * mur).inverse()", "sigma[][]": {"value": "rho.inverse()", "domain": [...]} }
      The suffix of a key gives the type: none Real, [] RealVectorValue, [][] RealTensorValue; me./mn.
      parameters take the type of their field. mat. material property, rc. real constant (keys of the model
      data), me. cell-data array, mn. point-data array (only usable as nodal loads). A "domain" list limits the
      domains where a parameter exists; a "default" (mat. and rc.) is the value used when the model data leaves
      it undefined: a number, its components, or C++ of the constants. Constants and expressions are C++ on the
      libMesh types (the std math functions, outer_product, I the identity tensor, axisymmetric); each
      expression is compiled into evaluate_expression(). Everywhere else a parameter is written by name.

  "dof_fields":     { "A": {"finite_elements": "NEDELEC_ONE", "active_components": ["z"],
                            "nodal_loads": "AMPS", "domain": [...]} }
  "nodal_fields":   { "B": {"active_components": [...], "domain": [...]} }
  "element_fields": { "Jt": {"active_components": ["z"], "domain": [...]} }
      Fields of fields.json. A vector field has one active component or all those of the dimension (default).
      DOF fields are LAGRANGE (LAGRANGE_VEC for a vector) unless finite_elements says NEDELEC_ONE; nodal fields
      L2_LAGRANGE(_VEC); element fields MONOMIAL CONSTANT. Lagrange families take the order of the options.
      Variables: every field is one libMesh variable. A vector with all the components of the element dimension
      (2D or 3D) is a vector variable named like the field (LAGRANGE_VEC, L2_LAGRANGE_VEC, NEDELEC_ONE,
      MONOMIAL_VEC); a single component is a scalar variable <field>_<c> (e.g. A_z of a 2D element). nodal_loads (scalar DOF fields) is an
      mn. parameter added to the right-hand side on the nodes. A nodal and an element field can share a name.
      The active components of a scalar field are the directions it varies along: by default those of the
      element (its grad is the in-plane gradient on a 2D element); ["z"] on a 2D element makes it the out-of-plane
      component of a 3D field, e.g. a potential over the out-of-plane length L: its grad is -phi / L e_z.

  "dof_kernels": { "K": [kernel, ...], "F": [...] }
      The key says where the kernels go:
        K  blocks test^T * coeff * trial: in K, or in C with a dt() trial, in M with a dt(dt()) trial.
        F  right-hand side: test^T * coeff.
      The coefficients are linear parameters. Where a law of nonlinear_laws.json replaces one (on the elements whose model
      data selects the law), the term is left out of K / F and the law's terms stand for it (see nonlinear_laws.json).
  "nodal_kernels":   { "B": [kernel, ...] }   L2 projection onto the field: M^-1 * sum of test^T * coeff * trial(u),
                                              the test being the field's shape functions ("test" may name a nodal
                                              field with the same finite elements).
  "element_kernels": { "Jt": [kernel, ...] }  mean over the quadrature points of sum of coeff * trial(u).
      A kernel: {"trial": "grad(dt(V))", "coeff": "-sigma", "test": "A", "domain": [...]}. An operand is a DOF
      field with at most one of grad, curl, div and, for the trial, up to time_order dt(): in nodal and element
      kernels dt() and dt(dt()) mean the first and second time derivatives of the solution instead of the solution. A coefficient is an optional sign
      and a parameter, expression or law output; none means 1. A kernel without trial in a field is the
      coefficient itself (e.g. Jt = Js).

} }

nonlinear_laws.json, "nonlinear_laws": { "<name>": {

  "source": "mat" | "rc" | "me",       the ids of the model data "nonlinear_laws" section: MATERIAL, CONSTANT or TYPE ids.
                                       An entry {"type": "<name>", ...} under an id selects the law on those cells.
  "input": "B",                        a nodal or element field, computed at the quadrature points by its kernels
                                       (element kernels if the element type has an element field of that name, else
                                       the nodal ones), neither averaged nor projected.
  "curve": {"x": "B", "y": "H"},       optional: the entry of the law carries the arrays "<x>" and "<y>" of a
                                       curve y(x), extended by its first and last segments.
  "expressions": { "b": "B.norm()", ... },  C++, in order.
  "coefficients": { "nu[][]": {"secant": "...", "tangent": "..."}, "nu_gauge": {"secant": "..."} } } }
      The coefficients (parameters or expressions of the element types) the law replaces, with their type,
      where the law is selected: secant S(x) and optionally
      tangent d(S(x) x)/dx, x the input. Expressions see the input by its name, the curve (curve.value(x),
      curve.slope(x)), the constants and parameters of the element type.

  A law applies to the element types that have a parameter it replaces. Every kernel with such a coefficient
  P gets variants used where the law is selected, in place of the kernel:
      K block t^T P L(u), d dt() on the trial:  residual  t^T S L(d^d u / dt^d),
                                                Jacobian  c_d t^T T L, T the tangent if L is the input's own trial
                                                (B = curl(A)) and d = 0, the secant otherwise;
      F vector t^T P:                           residual  - t^T S (no Jacobian: its input is not the unknown);
      nodal / element kernel:                   S in place of P; the field is then evaluated at every output.
  c_0 = 1 and c_1, c_2 are the coefficients of v and a on u in the time scheme. The residual of the system is
  K u + C v + M a - F + R, its Jacobian K + c_1 C + c_2 M + J.
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
RESERVED_NAMES = {"axisymmetric", "source", "input", "curve", "I", "std", "libMesh", *STD_FUNCTIONS, *LIBMESH_NAMES}
CPP_KEYWORDS = {
    "alignas", "alignof", "and", "asm", "auto", "bool", "break", "case", "catch", "char", "class", "const",
    "constexpr", "continue", "default", "delete", "do", "double", "else", "enum", "explicit", "export",
    "extern", "false", "float", "for", "friend", "goto", "if", "inline", "int", "long", "mutable", "namespace",
    "new", "noexcept", "not", "nullptr", "operator", "or", "private", "protected", "public", "register",
    "return", "short", "signed", "sizeof", "static", "struct", "switch", "template", "this", "throw", "true",
    "try", "typedef", "typename", "union", "unsigned", "using", "virtual", "void", "volatile", "while", "xor",
    "COUNT",
}
S, V, T = "Scalar", "Vector", "Tensor"
TYPE_BY_SUFFIX = {"": S, "[]": V, "[][]": T}
CPP_TYPE = {S: "Real", V: "RealVectorValue", T: "RealTensorValue"}
MAX_DOMAINS = 64  # A set of domains is a 64-bit mask.
COMPONENTS = ("x", "y", "z")
FIELD_KINDS = {"dof": "Dof", "nodal": "Nodal", "element": "Element"}
OPERATORS = {"grad": "Grad", "curl": "Curl", "div": "Div"}
# Where a dof kernel goes: stiffness (or time-derivative) matrices, right-hand side.
DOF_TARGETS = ("K", "F")
# The model data rows a law is selected by: materials, real constants, element types (TYPE).
LAW_SOURCES = {"mat": "Material", "rc": "Constant", "me": "Cell"}
MAX_LAWS = 32  # A set of laws is a 32-bit mask.
DOF_FAMILIES = ("LAGRANGE", "LAGRANGE_VEC", "NEDELEC_ONE")
# Families with the same shape functions share one FE: the key is the FE family used.
FE_OF_FAMILY = {"LAGRANGE": "LAGRANGE", "L2_LAGRANGE": "LAGRANGE",
                "LAGRANGE_VEC": "LAGRANGE_VEC", "L2_LAGRANGE_VEC": "LAGRANGE_VEC", "NEDELEC_ONE": "NEDELEC_ONE"}

# Products of the kernels: (test, coefficient, trial) -> Real for a block; (test, coefficient) -> Real for a
# source; (coefficient, trial) -> value for a term without test.
BLOCK_PRODUCTS = {(S, S, S), (V, S, V), (V, T, V), (S, V, V), (V, V, S), (T, S, T)}
SOURCE_PRODUCTS = {(S, S), (V, V), (T, T)}
VALUE_PRODUCTS = {(S, S): S, (S, V): V, (S, T): T, (T, V): V, (V, S): V}


def fail(message):
    sys.exit(f"generate_registries.py: {message}")


def load(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle, object_pairs_hook=OrderedDict)


def check_keys(options, allowed, context):
    if not isinstance(options, dict):
        fail(f"{context} must be an object.")
    for key in options:
        if key not in allowed:
            fail(f"{context}: unknown key '{key}' (allowed: {', '.join(sorted(allowed))}).")


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


def read_string(options, key, context, required=True):
    if key not in options:
        if required:
            fail(f"{context} needs '{key}'.")
        return None
    if not isinstance(options[key], str):
        fail(f"{context}: '{key}' must be a string.")
    return options[key]


def read_domain_mask(options, domains, context):
    """Bit mask of the "domain" list of an entry; all the domains when there is none."""
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


def domain_names(mask, domains):
    return ", ".join(name for index, name in enumerate(domains) if mask >> index & 1)


# ---- parameters and expressions --------------------------------------------------------------------------


def read_declaration(declaration, fields, context):
    """'mat.mur[][]' -> {name, type, source[, array]}."""
    prefix, dot, rest = "".join(declaration.split()).partition(".")
    if not dot:
        fail(f"{context} must start with mat., rc., me. or mn.")
    if prefix in ("mat", "rc"):
        name, suffix = split_key(rest, context)
        if name == "law":
            fail(f"{context}: 'law' is the model data key that selects the nonlinear laws.")
        return {"name": name, "type": TYPE_BY_SUFFIX[suffix], "source": "Material" if prefix == "mat" else "Constant"}
    if prefix in ("me", "mn"):
        if not IDENTIFIER.match(rest):
            fail(f"{context}: write the name of the field only, the type is the one of fields.json.")
        if rest not in fields:
            fail(f"{context}: '{rest}' is not a field of fields.json.")
        vector = fields[rest] == "vector"
        if prefix == "mn" and vector:
            fail(f"{context}: point data is a scalar array per node.")
        return {"name": rest, "type": V if vector else S, "source": "Cell" if prefix == "me" else "Point", "array": rest}
    fail(f"{context}: unknown source '{prefix}' (expected mat, rc, me or mn).")


def cpp_literal(value, type_, context):
    """C++ text of a default value: a number, [3 numbers], [[3 numbers] x 3], or a C++ expression of the constants."""
    if isinstance(value, str):
        return value
    number = lambda v: isinstance(v, (int, float)) and not isinstance(v, bool)
    if type_ == S and number(value):
        return repr(float(value))
    if type_ == V and isinstance(value, list) and len(value) == 3 and all(map(number, value)):
        return "RealVectorValue(" + ", ".join(repr(float(v)) for v in value) + ")"
    if type_ == T and isinstance(value, list) and len(value) == 3 and all(
            isinstance(row, list) and len(row) == 3 and all(map(number, row)) for row in value):
        return "RealTensorValue(" + ", ".join(repr(float(v)) for row in value for v in row) + ")"
    fail(f"{context}: the default of a {type_.lower()} is a C++ expression or {'a number' if type_ == S else 'its components'}.")


def translate_parameters(text, known, where, check=None):
    """C++ text with the parameters it names read from `source`; `check(name, entry)` validates each use."""

    def translate(match):
        token = match.group(0)
        if token == "axisymmetric":
            return "source.axisymmetric()"
        if token not in known:
            return token
        entry = known[token]
        if entry["source"] == "Point":
            fail(f"{where}: '{token}' is point data and cannot be used in expressions.")
        if check:
            check(token, entry)
        return f"source.get_{entry['type'].lower()}(ParameterId::{token})"

    # Not a member (.det) and not part of a number (1e-7).
    return re.sub(r"(?<![\w.])[A-Za-z_]\w*", translate, text)


def cpp_expression(entry, index, entries, context):
    """C++ text of an expression. The parameters it names must come before it and exist where it does."""
    where = f"{context}, expression '{entry['name']}'"
    position = {other["name"]: i for i, other in enumerate(entries)}

    def check(token, other):
        if position[token] >= index:
            fail(f"{where}: '{token}' is not defined before it.")
        if other["domains"] & entry["domains"] != entry["domains"]:
            fail(f"{where}: '{token}' does not exist in every domain of the expression.")

    return translate_parameters(entry["text"], {e["name"]: e for e in entries}, where, check)


# ---- laws ------------------------------------------------------------------------------------------------


def read_laws(resource_dir):
    """nonlinear_laws.json: name -> {source, input, curve (name, x, y) or None, locals, replaces}."""
    definitions = load(f"{resource_dir}/nonlinear_laws.json").get("nonlinear_laws")
    if not isinstance(definitions, dict):
        fail("nonlinear_laws.json must contain a 'nonlinear_laws' object.")
    laws = OrderedDict()
    for name, options in definitions.items():
        context = f"nonlinear_laws.json, law '{name}'"
        check_enumerator(name, context)
        check_keys(options, {"source", "input", "curve", "expressions", "coefficients"}, context)
        source = read_string(options, "source", context)
        if source not in LAW_SOURCES:
            fail(f"{context}: source must be one of {', '.join(LAW_SOURCES)} (the model data rows that select it).")
        curve = None
        if "curve" in options:
            check_keys(options["curve"], {"x", "y"}, f"{context}, curve")
            x, y = (read_string(options["curve"], key, f"{context}, curve") for key in ("x", "y"))
            if x == y or "type" in (x, y):
                fail(f"{context}, curve: x and y must be two names other than 'type'.")
            curve = (name, x, y)
        locals_ = options.get("expressions", {})
        check_keys(locals_, set(locals_), f"{context}, expressions")
        for local, text in locals_.items():
            check_enumerator(local, f"{context}, expression '{local}'")
            if local in RESERVED_NAMES:
                fail(f"{context}: expression '{local}' is reserved.")
            if not isinstance(text, str):
                fail(f"{context}, expression '{local}' must be a C++ expression.")
        replaces = OrderedDict()
        if not isinstance(options.get("coefficients"), dict) or not options["coefficients"]:
            fail(f"{context}: 'coefficients' names the coefficients the law replaces, with their secant (and tangent).")
        for key, value in options["coefficients"].items():
            where = f"{context}, coefficient '{key}'"
            parameter, suffix = split_key(key, where)
            check_keys(value, {"secant", "tangent"}, where)
            replaces[parameter] = {"type": TYPE_BY_SUFFIX[suffix], "secant": read_string(value, "secant", where),
                                   "tangent": read_string(value, "tangent", where, required=False)}
        laws[name] = {"name": name, "source": source, "input": read_string(options, "input", context), "curve": curve,
                      "locals": locals_, "replaces": replaces}
    if len(laws) > MAX_LAWS:
        fail(f"nonlinear_laws.json: at most {MAX_LAWS} laws.")
    return laws


# ---- element types ---------------------------------------------------------------------------------------


class ElementType:
    """One element type of element_types.json, checked and reduced to the tables of formulation_spec.h."""

    SECTIONS = {"inherits", "dimension", "time_order", "out_of_plane_length", "constants", "parameters", "expressions",
                "dof_fields", "dof_kernels", "nodal_fields", "nodal_kernels", "element_fields", "element_kernels"}

    def __init__(self, name, definition, domains, fields_json, laws, law_index):
        self.name = name
        self.context = f"element type '{name}'"
        self.definition = definition
        self.domains = domains
        self.fields_json = fields_json
        self.all_domains = (1 << len(domains)) - 1
        check_keys(definition, self.SECTIONS, self.context)
        self.dimension = definition.get("dimension")
        if self.dimension not in (2, 3):
            fail(f"{self.context}: dimension must be 2 or 3.")
        self.time_order = definition.get("time_order", 1)
        if self.time_order not in (0, 1, 2):
            fail(f"{self.context}: time_order must be 0 (static), 1 or 2.")

        self.read_parameters(definition)
        self.length = definition.get("out_of_plane_length")
        if self.length is not None:
            where = f"{self.context}, out_of_plane_length"
            if self.dimension != 2:
                fail(f"{where}: only 2D element types have an out-of-plane length.")
            entry = self.parameters.get(self.length) if isinstance(self.length, str) else None
            if entry is None or entry["type"] != S or entry["source"] == "Point":
                fail(f"{where}: must name a scalar parameter or expression of the element type.")
            if entry["domains"] != self.all_domains:
                fail(f"{where}: '{self.length}' must exist in every domain.")
        self.fields = []
        self.fes = []
        self.nablas = []
        self.kernels = []
        self.terms = []
        self.evaluations = []
        self.laws = []
        self.field_terms = {}  # (kind, name) -> its term definitions, for the inputs of the laws
        for kind in FIELD_KINDS:
            self.read_fields(kind, definition.get(f"{kind}_fields", {}))
        self.bind_laws(laws, law_index)
        self.read_dof_kernels(definition.get("dof_kernels", {}))
        for kind, evaluation in (("nodal", "L2Projection"), ("element", "Average")):
            self.read_field_kernels(kind, evaluation, definition.get(f"{kind}_kernels", {}))
        self.read_law_inputs()

    # -- parameters --

    def read_parameters(self, definition):
        context = self.context
        self.constants = []
        for key, value in definition.get("constants", {}).items():
            name, suffix = split_key(key, f"{context}, constant '{key}'")
            if isinstance(value, bool) or not isinstance(value, (int, float, str)):
                fail(f"{context}, constant '{key}' must be a number or a C++ expression.")
            self.constants.append({"name": name, "type": TYPE_BY_SUFFIX[suffix],
                                   "text": value if isinstance(value, str) else repr(value)})

        entries = []
        declarations = definition.get("parameters", {})
        check_keys(declarations, set(declarations), f"{context}, parameters")
        for declaration, options in declarations.items():
            where = f"{context}, parameter '{declaration}'"
            entry = read_declaration(declaration, self.fields_json, where)
            check_keys(options, {"domain", "default"}, where)
            entry["domains"] = read_domain_mask(options, self.domains, where)
            if "default" in options:
                if entry["source"] not in ("Material", "Constant"):
                    fail(f"{where}: only mat. and rc. parameters have a default.")
                entry["default"] = cpp_literal(options["default"], entry["type"], where)
            entries.append(entry)

        expressions = definition.get("expressions", {})
        check_keys(expressions, set(expressions), f"{context}, expressions")
        for key, value in expressions.items():
            name, suffix = split_key(key, f"{context}, expression '{key}'")
            options = value if isinstance(value, dict) else {"value": value}
            check_keys(options, {"value", "domain"}, f"{context}, expression '{key}'")
            if not isinstance(options.get("value"), str):
                fail(f"{context}, expression '{key}' must be a string or an object with a string 'value'.")
            entries.append({"name": name, "type": TYPE_BY_SUFFIX[suffix], "source": "Expression", "text": options["value"],
                            "domains": read_domain_mask(options, self.domains, f"{context}, expression '{key}'")})

        for item in self.constants + entries:
            if item["name"] in RESERVED_NAMES:
                fail(f"{context}: '{item['name']}' is reserved.")
            check_enumerator(item["name"], context)
        unique([item["name"] for item in self.constants + entries], f"{context}: constants, parameters and expressions")
        for index, entry in enumerate(entries):
            if entry["source"] == "Expression":
                entry["code"] = cpp_expression(entry, index, entries, context)
        self.entries = entries
        self.parameters = {entry["name"]: entry for entry in entries}

    # -- fields and finite elements --

    def read_fields(self, kind, fields):
        check_keys(fields, set(fields), f"{self.context}, {kind}_fields")
        allowed = {"dof": {"active_components", "finite_elements", "nodal_loads", "domain"},
                   "nodal": {"active_components", "domain"},
                   "element": {"active_components", "domain"}}[kind]
        for name, options in fields.items():
            context = f"{self.context}, {kind} field '{name}'"
            check_keys(options, allowed, context)
            if name not in self.fields_json:
                fail(f"{context} is not declared in fields.json.")
            physical_vector = self.fields_json[name] == "vector"

            requested = options.get("active_components")
            if requested is not None:
                if not isinstance(requested, list) or not requested or any(c not in COMPONENTS for c in requested):
                    fail(f"{context}: active_components is a list of x, y, z.")
                requested = sorted(set(COMPONENTS.index(c) for c in requested))
            # The active components of a scalar are the directions it varies along: those of the element, or only
            # z on a 2D element (the out-of-plane component of a 3D field, e.g. a potential over the length L).
            out_of_plane = False
            if not physical_vector:
                if requested is not None and requested != list(range(self.dimension)):
                    if self.dimension != 2 or requested != [2]:
                        fail(f"{context}: a scalar field varies along all the directions of the element, or only z on a "
                             "2D element.")
                    out_of_plane = True
                components = [0]
            elif requested is None or requested == list(range(self.dimension)):
                components = list(range(self.dimension))
            elif len(requested) == 1:
                components = requested
            else:
                fail(f"{context}: a vector field has one active component or all those of the dimension.")
            single_component = physical_vector and len(components) == 1
            multi = len(components) > 1


            field = {"name": name, "kind": kind, "value": V if physical_vector else S, "components": components,
                     "single_component": single_component, "out_of_plane": out_of_plane, "fe": None, "family": None,
                     "order": "Option", "variable": None, "nodal_load": None,
                     "domains": read_domain_mask(options, self.domains, context)}

            if kind in ("dof", "nodal"):
                family = options.get("finite_elements", "LAGRANGE")
                if family not in DOF_FAMILIES:
                    fail(f"{context}: finite_elements must be one of {', '.join(DOF_FAMILIES)}.")
                if family == "LAGRANGE" and multi:
                    family = "LAGRANGE_VEC"
                if family == "LAGRANGE_VEC" and not multi:
                    fail(f"{context}: LAGRANGE_VEC needs a vector field with all its components.")
                if family == "NEDELEC_ONE" and (self.dimension != 3 or not multi):
                    fail(f"{context}: NEDELEC_ONE needs a vector field with all its components on a 3D element.")
                if kind == "nodal":
                    family = "L2_" + family
                field["family"] = family
                field["order"] = "First" if family == "NEDELEC_ONE" else "Option"
                field["fe"] = self.fe_index(FE_OF_FAMILY[family])
            else:
                field["family"] = "MONOMIAL_VEC" if multi else "MONOMIAL"
                field["order"] = "Constant"
            # Every field is one libMesh variable: a vector field with all the components of the element dimension
            # (2D or 3D) a vector variable named like the field, a single component the scalar variable <field>_<c>.
            field["variable"] = f"{name}_{COMPONENTS[components[0]]}" if single_component else name

            loads = options.get("nodal_loads")
            if loads is not None:
                if physical_vector:
                    fail(f"{context}: nodal_loads only apply to scalar dof fields.")
                parameter = self.parameters.get(loads)
                if not parameter or parameter["source"] != "Point":
                    fail(f"{context}: nodal_loads '{loads}' is not a point-data (mn.) parameter of the element type.")
                field["nodal_load"] = parameter["array"]

            field["index"] = len(self.fields)
            self.fields.append(field)
        unique([f["name"] for f in self.fields if f["kind"] == kind], f"{self.context}, {kind}_fields")

    def fe_index(self, family):
        for index, fe in enumerate(self.fes):
            if fe["family"] == family:
                return index
        self.fes.append({"family": family, "vector": family != "LAGRANGE",
                         "order": "First" if family == "NEDELEC_ONE" else "Option"})
        return len(self.fes) - 1

    def field(self, name, context, kind):
        for field in self.fields:
            if field["name"] == name and field["kind"] == kind:
                return field
        fail(f"{context}: '{name}' is not a {kind} field of the element type.")

    # -- operands and nablas --

    def parse_operand(self, text, context, allow_dt):
        """'grad(dt(V))' -> (dof field, operator, number of dt())."""
        rest = "".join(text.split())
        operator, dt = "Identity", 0
        while "(" in rest:
            open_at = rest.index("(")
            if not rest.endswith(")"):
                fail(f"{context}: unbalanced parentheses in '{text}'.")
            name, rest = rest[:open_at], rest[open_at + 1:-1]
            if name == "dt":
                if not allow_dt:
                    fail(f"{context}: dt() is only allowed on the trial, in '{text}'.")
                dt += 1
                if dt > self.time_order:
                    fail(f"{context}: '{text}' has {dt} dt(), more than the time_order {self.time_order} of the element type.")
            elif name in OPERATORS:
                if operator != "Identity":
                    fail(f"{context}: only one of grad, curl, div in '{text}'.")
                operator = OPERATORS[name]
            else:
                fail(f"{context}: unknown operator '{name}' in '{text}'.")
        return self.field(rest, context, "dof"), operator, dt

    def nabla(self, field, operator, context):
        """Index of the nabla of `operator` on the FE of `field` (created once), and its output type."""
        fe = self.fes[field["fe"]]
        component = field["components"][0] if field["single_component"] else -1
        out_of_plane = field.get("out_of_plane", False)
        name = f"{operator.lower()}({field['name']})"
        if not fe["vector"]:
            embedded = field["single_component"]
            if operator == "Identity":
                output = V if embedded else S
            elif operator == "Grad" and not embedded or operator == "Curl" and embedded and self.dimension == 2 and component == 2:
                output = V
            else:
                fail(f"{context}: {name} is not defined for this field on a {self.dimension}D element "
                     "(scalar FE: grad of a scalar; curl of the out-of-plane component of a 2D element).")
        else:
            nedelec = fe["family"] == "NEDELEC_ONE"
            if operator == "Identity" or operator == "Curl" and self.dimension == 3:
                output = V
            elif operator == "Div" and not nedelec:
                output = S
            elif operator == "Grad" and not nedelec:
                output = T
            else:
                fail(f"{context}: {name} is not defined for {fe['family']} on a {self.dimension}D element.")
        key = (field["fe"], operator, component, out_of_plane)
        for index, nabla in enumerate(self.nablas):
            if nabla["key"] == key:
                return index, output
        self.nablas.append({"key": key, "fe": field["fe"], "op": operator, "component": component, "out_of_plane": out_of_plane,
                            "output": output})
        return len(self.nablas) - 1, output

    # -- laws --

    def bind_laws(self, laws, law_index):
        """The laws of nonlinear_laws.json that replace parameters of the element type, compiled on its parameters."""
        self.law_outputs = OrderedDict()
        self.replacements = {}  # parameter -> [{law, law_id, secant, tangent}]
        self.law_by_name = {}
        for name, law in laws.items():
            present = [(parameter, r) for parameter, r in law["replaces"].items() if parameter in self.parameters]
            if not present:
                continue
            context = f"{self.context}, law '{name}'"
            parameters = self.parameters
            used_names = set(parameters) | {c["name"] for c in self.constants} | {law["input"]}
            locals_ = []
            for local, text in law["locals"].items():
                if local in used_names:
                    fail(f"{context}: local '{local}' hides another name of the element type.")
                used_names.add(local)
                locals_.append({"name": local, "code": translate_parameters(text, parameters, context)})
            outputs = []
            for parameter, replaced in present:
                entry = parameters[parameter]
                if entry["source"] == "Point":
                    fail(f"{context}: '{parameter}' is point data, which cannot be nonlinear.")
                if entry["type"] != replaced["type"]:
                    fail(f"{context}: '{parameter}' is a {entry['type'].lower()} in the element type, a "
                         f"{replaced['type'].lower()} in the law.")
                names = {}
                for role in ("secant", "tangent"):
                    if replaced[role] is None:
                        continue
                    output = {"name": f"{name}_{parameter}_{role}", "type": replaced["type"], "law": name,
                              "parameter": parameter, "index": len(outputs),
                              "code": translate_parameters(replaced[role], parameters, context)}
                    outputs.append(output)
                    self.law_outputs[output["name"]] = output
                    names[role] = output["name"]
                self.replacements.setdefault(parameter, []).append(
                    {"law": name, "law_id": law_index[name], "secant": names["secant"], "tangent": names.get("tangent")})
            bound = {"name": name, "input_name": law["input"], "curve_key": law["curve"], "locals": locals_,
                     "outputs": outputs}
            self.laws.append(bound)
            self.law_by_name[name] = bound

    def law_operand(self, name):
        """The trial of the input of law `name` when the input is that trial alone (B = curl(A)); None otherwise.
        A K kernel with this trial uses the tangent of the law in J, any other its secant."""
        input_name = self.law_by_name[name]["input_name"]
        for kind in ("element", "nodal"):
            terms = self.definition.get(f"{kind}_kernels", {}).get(input_name)
            if terms is not None:
                break
        else:
            return None
        if len(terms) != 1 or "".join((terms[0].get("coeff") or "").split()) or not isinstance(terms[0].get("trial"), str):
            return None
        return "".join(terms[0]["trial"].split())

    def read_law_inputs(self):
        """The input of a law is a field evaluated at the quadrature points with its own kernels, not averaged nor projected."""
        for law in self.laws:
            context = f"{self.context}, law '{law['name']}'"
            name = law["input_name"]
            key = ("element", name) if ("element", name) in self.field_terms else ("nodal", name)
            if key not in self.field_terms:
                fail(f"{context}: input '{name}' is not a nodal or element field with kernels.")
            field = self.field(name, context, key[0])
            evaluation = {"kind": "Quadrature", "field": field, "first_term": len(self.terms), "mass_kernel": -1}
            terms = self.field_terms[key]
            for position, term in enumerate(terms):
                options = {k: v for k, v in term.items() if k != "test"}
                self.read_term(options, "value", field, f"{context}, input term {position}")
            evaluation["n_terms"] = len(terms)
            law["input"] = field
            law["evaluation"] = len(self.evaluations)
            self.evaluations.append(evaluation)

    # -- coefficients, kernels, terms --

    def coefficient(self, text, context):
        if text is None or not "".join(text.split()):
            return {"kind": "One", "name": None, "type": S, "sign": 1.0}
        text = "".join(text.split())
        sign = 1.0
        if text[0] in "+-":
            sign = -1.0 if text[0] == "-" else 1.0
            text = text[1:]
        if text in self.parameters:
            parameter = self.parameters[text]
            if parameter["source"] == "Point":
                fail(f"{context}: '{text}' is point data; it can only be the nodal_loads of a dof field.")
            return {"kind": "Parameter", "name": text, "type": parameter["type"], "sign": sign}
        fail(f"{context}: '{text}' is not a parameter or expression of the element type.")

    def kernel(self, trial, coefficient, test, output):
        key = (trial, test, coefficient["kind"], coefficient["name"], coefficient["sign"])
        for index, kernel in enumerate(self.kernels):
            if kernel["key"] == key:
                return index
        self.kernels.append({"key": key, "trial": trial, "test": test, "coefficient": coefficient, "output": output})
        return len(self.kernels) - 1

    def read_dof_kernels(self, groups):
        check_keys(groups, set(DOF_TARGETS), f"{self.context}, dof_kernels")
        # The terms of the laws, derived from the K and F kernels whose coefficient a law replaces: Jacobian blocks
        # (the tangent when the trial is the law's own input, the secant otherwise), residual blocks (secant) and
        # residual vectors (- the secant of an F kernel).
        derived = {"Tangent": [], "Residual": [], "ResidualRhs": []}
        for target, terms in groups.items():
            context = f"{self.context}, dof_kernels '{target}'"
            if not isinstance(terms, list) or not terms:
                fail(f"{context} must be a non-empty array of kernels.")
            if any(not isinstance(t, dict) for t in terms):
                fail(f"{context}: each kernel is an object.")
            # A group is one evaluation per shape: blocks (with trial), vectors (without).
            for with_trial in (True, False):
                selected = [(i, t) for i, t in enumerate(terms) if ("trial" in t) == with_trial]
                if not selected:
                    continue
                evaluation = {"kind": "Block" if with_trial else "Rhs", "field": None, "first_term": len(self.terms),
                              "n_terms": len(selected), "mass_kernel": -1}
                for position, term in selected:
                    parsed = self.read_term(term, target, None, f"{context}, kernel {position}")
                    for r in parsed["replacements"]:
                        if with_trial:
                            own = parsed["term"]["trial_dt"] == 0 and parsed["trial_text"] == self.law_operand(r["law"])
                            derived["Tangent"].append((parsed, r["tangent"] if own and r["tangent"] else r["secant"],
                                                       parsed["sign"], r["law_id"]))
                            derived["Residual"].append((parsed, r["secant"], parsed["sign"], r["law_id"]))
                        else:
                            derived["ResidualRhs"].append((parsed, r["secant"], -parsed["sign"], r["law_id"]))
                self.evaluations.append(evaluation)
        for key, kind, stage in (("Tangent", "Block", "Tangent"), ("Residual", "Block", "Residual"),
                                 ("ResidualRhs", "Rhs", "Residual")):
            if not derived[key]:
                continue
            evaluation = {"kind": kind, "field": None, "first_term": len(self.terms), "n_terms": len(derived[key]),
                          "mass_kernel": -1}
            for parsed, output, sign, law_id in derived[key]:
                self.law_variant(parsed, output, sign, stage, law_id)
            self.evaluations.append(evaluation)

    def read_field_kernels(self, kind, evaluation_kind, groups):
        check_keys(groups, set(groups), f"{self.context}, {kind}_kernels")
        for name, terms in groups.items():
            context = f"{self.context}, {kind}_kernels '{name}'"
            field = self.field(name, context, kind)
            if not isinstance(terms, list) or not terms or any(not isinstance(t, dict) for t in terms):
                fail(f"{context} must be a non-empty array of kernels.")
            self.field_terms[(kind, name)] = terms
            evaluation = {"kind": evaluation_kind, "field": field, "first_term": len(self.terms), "mass_kernel": -1}
            if evaluation_kind == "L2Projection":
                identity, _ = self.nabla(field, "Identity", context)
                evaluation["mass_kernel"] = self.kernel(identity, self.coefficient(None, context), identity, None)
            # A kernel whose coefficient a law replaces has a variant with the law's secant, used where the law is.
            for position, term in enumerate(terms):
                parsed = self.read_term(term, kind, field, f"{context}, kernel {position}")
                for r in parsed["replacements"]:
                    self.law_variant(parsed, r["secant"], parsed["sign"], "Linear", r["law_id"])
            evaluation["n_terms"] = len(self.terms) - evaluation["first_term"]
            self.evaluations.append(evaluation)
        for field in self.fields:
            if field["kind"] == kind and (kind, field["name"]) not in self.field_terms:
                fail(f"{self.context}: {kind} field '{field['name']}' has no {kind}_kernels.")

    def read_term(self, options, target, field, context):
        """A kernel placed in `target`: a dof_kernels key, "nodal", "element", or "value" (a law input).

        Appends its term and returns what its law variants need (law_variant)."""
        if "nonlinear" in options:
            fail(f"{context}: 'nonlinear' was removed: the laws of nonlinear_laws.json say which parameters they replace.")
        check_keys(options, {"trial", "coeff", "test", "domain"}, context)
        trial_text = read_string(options, "trial", context, required=False)
        test_text = read_string(options, "test", context, required=False)
        coefficient = self.coefficient(read_string(options, "coeff", context, required=False), context)
        ctype = coefficient["type"]

        trial = test = None
        trial_field = test_field = None
        trial_dt = 0
        domains = read_domain_mask(options, self.domains, context) & (field["domains"] if field else self.all_domains)
        if trial_text is not None:
            trial_field, operator, trial_dt = self.parse_operand(trial_text, f"{context} trial", True)
            trial, trial_type = self.nabla(trial_field, operator, f"{context} trial")
            domains &= trial_field["domains"]

        output = None
        stage = "Linear"
        if target in DOF_TARGETS:
            if test_text is None:
                fail(f"{context}: a dof kernel has a test.")
            test_field, operator, _ = self.parse_operand(test_text, f"{context} test", False)
            test, test_type = self.nabla(test_field, operator, f"{context} test")
            domains &= test_field["domains"]
            if trial is not None and (test_type, ctype, trial_type) not in BLOCK_PRODUCTS:
                fail(f"{context}: no product {test_type} test, {ctype} coefficient, {trial_type} trial.")
            if trial is None and (test_type, ctype) not in SOURCE_PRODUCTS:
                fail(f"{context}: no product {test_type} test, {ctype} coefficient.")
            if target == "K" and trial is None:
                fail(f"{context}: K kernels have a trial.")
            if target == "F" and trial is not None:
                fail(f"{context}: F kernels have no trial.")
        else:
            if trial_dt and target == "value":
                fail(f"{context}: the input of a law has no dt() trial.")
            output = VALUE_PRODUCTS.get((ctype, trial_type)) if trial is not None else ctype
            if output is None:
                fail(f"{context}: no product {ctype} coefficient, {trial_type} trial.")
            if output != field["value"]:
                fail(f"{context}: the kernel gives a {output.lower()}, field '{field['name']}' is a {field['value'].lower()}.")
            if target == "nodal":
                test, test_type = self.nabla(field, "Identity", context)
                if test_text is not None:
                    basis = self.field(test_text, f"{context} test", "nodal")
                    if basis["fe"] != field["fe"] or basis["components"] != field["components"]:
                        fail(f"{context}: the test '{test_text}' must have the finite elements of '{field['name']}'.")
                if trial is not None and (test_type, ctype, trial_type) not in BLOCK_PRODUCTS:
                    fail(f"{context}: no product {test_type} field, {ctype} coefficient, {trial_type} trial.")
                if trial is None and (test_type, ctype) not in SOURCE_PRODUCTS:
                    fail(f"{context}: no product {test_type} field, {ctype} coefficient.")
                output = None
            elif test_text is not None:
                fail(f"{context}: an element field kernel has no test.")

        if domains == 0:
            fail(f"{context}: the kernel applies in no domain (check the domains of its fields).")
        self.check_coefficient_domains(coefficient, domains, context)

        replacements = self.replacements.get(coefficient["name"], []) if coefficient["kind"] == "Parameter" else []
        if replacements and target == "value":
            fail(f"{context}: the input of a law cannot use '{coefficient['name']}', which a law replaces.")
        term = {
            "kernel": self.kernel(trial, coefficient, test, output),
            "trial_field": trial_field["index"] if trial_field else -1,
            "test_field": test_field["index"] if test_field else -1,
            "trial_dt": trial_dt, "stage": stage, "domains": domains,
            "law": -1, "replaced_by": sum(1 << r["law_id"] for r in replacements),
        }
        self.terms.append(term)
        return {"trial": trial, "test": test, "output": output, "sign": coefficient["sign"], "term": term,
                "replacements": replacements, "trial_text": "".join(trial_text.split()) if trial_text else None}

    def law_variant(self, parsed, output_name, sign, stage, law_id):
        """The term `parsed` with the law output `output_name` as coefficient, used where law `law_id` is."""
        output = self.law_outputs[output_name]
        coefficient = {"kind": "LawOutput", "name": output_name, "type": output["type"], "sign": sign}
        self.terms.append(dict(parsed["term"], kernel=self.kernel(parsed["trial"], coefficient, parsed["test"], parsed["output"]),
                               stage=stage, law=law_id, replaced_by=0))

    def check_coefficient_domains(self, coefficient, domains, context):
        if coefficient["kind"] == "One":
            return
        name = coefficient["name"]
        missing = domains & ~self.parameters[name]["domains"]
        if missing:
            fail(f"{context}: '{name}' is not defined in {domain_names(missing, self.domains)}, where the kernel applies.")


def resolve_inheritance(definitions):
    """An element type with "inherits" starts from its parent: each of its sections replaces the parent's key by key."""
    resolved = OrderedDict()

    def resolve(name, chain):
        if name in resolved:
            return resolved[name]
        if name in chain:
            fail(f"element type '{name}': circular inherits.")
        definition = definitions.get(name)
        if not isinstance(definition, dict):
            fail(f"element type '{chain[-1]}' inherits from '{name}', which is not an element type.")
        parent_name = definition.get("inherits")
        if parent_name is None:
            result = OrderedDict(definition)
        else:
            if not isinstance(parent_name, str):
                fail(f"element type '{name}': inherits is the name of an element type.")
            result = OrderedDict((k, OrderedDict(v) if isinstance(v, dict) else v)
                                 for k, v in resolve(parent_name, chain + [name]).items())
            for key, value in definition.items():
                if key == "inherits":
                    continue
                if isinstance(value, dict) and isinstance(result.get(key), dict):
                    result[key].update(value)
                else:
                    result[key] = value
        result.pop("inherits", None)
        resolved[name] = result
        return result

    for name in definitions:
        resolve(name, [])
    return OrderedDict((name, resolved[name]) for name in definitions)


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

    laws = read_laws(resource_dir)
    law_index = {name: i for i, name in enumerate(laws)}
    types = OrderedDict((name, ElementType(name, definition, domains, fields, laws, law_index))
                        for name, definition in resolve_inheritance(definitions).items())

    # A parameter name has one type and one source in every element type.
    by_name = OrderedDict()
    for type_name, element_type in types.items():
        for entry in element_type.entries:
            shared = {key: entry.get(key) for key in ("name", "type", "source", "array")}
            known = by_name.setdefault(entry["name"], shared)
            if known != shared:
                fail(f"element type '{type_name}': '{entry['name']}' is declared differently in another element type.")
    parameters = [p for p in by_name.values() if p["source"] != "Expression"]
    parameters += [p for p in by_name.values() if p["source"] == "Expression"]
    parameter_index = {p["name"]: i for i, p in enumerate(parameters)}

    material_properties = [p for p in parameters if p["source"] == "Material"]
    real_constants = [p for p in parameters if p["source"] == "Constant"]
    for index, entry in enumerate(material_properties):
        entry["property"] = index
    for index, entry in enumerate(real_constants):
        entry["property"] = index
    for element_type in types.values():
        mask_of = {entry["name"]: entry["domains"] for entry in element_type.entries}
        element_type.masks = [mask_of.get(p["name"], 0) for p in parameters]

    # Curves: one per law that has one, given with the law in the "nonlinear_laws" section of the model data.
    curves = OrderedDict((law["curve"], {"index": i}) for i, law in enumerate(x for x in laws.values() if x["curve"]))
    # Law outputs: global names, one meaning each.
    law_outputs = OrderedDict()
    for type_name, element_type in types.items():
        for law in element_type.laws:
            for output in law["outputs"]:
                shared = {"law": law["name"], "type": output["type"]}
                if law_outputs.setdefault(output["name"], shared) != shared:
                    fail(f"element type '{type_name}', law '{law['name']}': output '{output['name']}' belongs to another "
                         "law or has another type in another element type.")
    output_index = {name: i for i, name in enumerate(law_outputs)}

    # libMesh variables: every label the element types use, in the order of fields.json.
    variables = OrderedDict()
    nodal_loads = OrderedDict()
    for field_name, kind in fields.items():
        labels = [(field_name, field_name, 0, "Scalar")] if kind == "scalar" else \
            [(field_name, field_name, 0, "Vector")] + [(f"{field_name}_{c}", field_name, i, "Component") for i, c in enumerate(COMPONENTS)]
        for label, api_label, component, variable_kind in labels:
            variables.setdefault(label, {"label": label, "api": api_label, "field": field_name, "component": component,
                                         "kind": variable_kind})
    for type_name, element_type in types.items():
        for field in element_type.fields:
            if field["nodal_load"]:
                label = field["variable"]
                if nodal_loads.setdefault(label, field["nodal_load"]) != field["nodal_load"]:
                    fail(f"element type '{type_name}': '{label}' is loaded by another array in another element type.")

    names = [p["name"] for p in parameters] + list(law_outputs) + list(laws)
    unique(names, "parameter, law output and law names")
    for name in list(fields) + list(domains) + list(types) + list(variables):
        check_enumerator(name, f"registry name '{name}'")

    return {
        "domains": domains,
        "fields": fields,
        "types": types,
        "parameters": parameters,
        "parameter_index": parameter_index,
        "material_properties": material_properties,
        "real_constants": real_constants,
        "laws": laws,
        "law_index": law_index,
        "law_outputs": law_outputs,
        "output_index": output_index,
        "curves": curves,
        "variables": variables,
        "nodal_loads": nodal_loads,
    }


# ---- registries.h ----------------------------------------------------------------------------------------


def cpp_string(text):
    return json.dumps(text)


def emit_enum(out, src, enum, names, doc):
    out.append(f"/// {doc}")
    out.append(f"enum class {enum} : unsigned short")
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
    out.append(f"/// The {enum} named @p text; nullopt when there is none.")
    out.append(f"std::optional<{enum}> parse_{enum}(std::string_view text);")
    out.append("")
    src.append(f"std::optional<{enum}> parse_{enum}(std::string_view text)")
    src.append("{")
    src.append(f"    for (std::size_t i = 0; i < k{enum}Names.size(); ++i)")
    src.append(f"        if (k{enum}Names[i] == text)")
    src.append(f"            return static_cast<{enum}>(i);")
    src.append("    return std::nullopt;")
    src.append("}")
    src.append("")


def emit_scope(out, indent):
    for name in LIBMESH_NAMES:
        out.append(f"{indent}using libMesh::{name};")
    for name in STD_FUNCTIONS:
        out.append(f"{indent}using std::{name};")
    out.append(f"{indent}[[maybe_unused]] const RealTensorValue I(1, 0, 0, 0, 1, 0, 0, 0, 1);")


def emit_constants(out, element_type, indent):
    for constant in element_type.constants:
        out.append(f"{indent}[[maybe_unused]] const {CPP_TYPE[constant['type']]} {constant['name']} = {constant['text']};")


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
    emit_scope(out, "    ")
    out.append("    switch (type)")
    out.append("    {")
    for type_name, element_type in types.items():
        expressions = [e for e in element_type.entries if e["source"] == "Expression"]
        if not expressions:
            continue
        out.append(f"    case ElementTypeId::{type_name}:")
        out.append("    {")
        emit_constants(out, element_type, "        ")
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


def emit_defaults(out, src, types):
    out.append("/**")
    out.append(" * @brief Default value of the mat./rc. parameter @p id of the element type @p type, used when the model data")
    out.append(" * leaves it undefined; nullopt when it has none.")
    out.append(" */")
    out.append("std::optional<ParameterValue> default_value(ElementTypeId type, ParameterId id);")
    out.append("")
    out = src  # the definition goes to registries.cpp
    out.append("std::optional<ParameterValue> default_value(ElementTypeId type, ParameterId id)")
    out.append("{")
    emit_scope(out, "    ")
    out.append("    switch (type)")
    out.append("    {")
    for type_name, element_type in types.items():
        defaults = [e for e in element_type.entries if "default" in e]
        if not defaults:
            continue
        out.append(f"    case ElementTypeId::{type_name}:")
        out.append("    {")
        emit_constants(out, element_type, "        ")
        out.append("        switch (id)")
        out.append("        {")
        for entry in defaults:
            out.append(f"        case ParameterId::{entry['name']}:")
            out.append(f"            return expression_support::expression_value<ParameterType::{entry['type']}>({entry['default']});")
        out.append("        default:")
        out.append("            break;")
        out.append("        }")
        out.append("        break;")
        out.append("    }")
    out.append("    default:")
    out.append("        break;")
    out.append("    }")
    out.append("    return std::nullopt;")
    out.append("}")
    out.append("")


def render_registries(registries):
    """The text of registries.h and of registries.cpp (the definitions of the lookups that are not constexpr)."""
    src = [
        "// Generated by src/tools/generate_registries.py from src/pytec3d/resources/*.json. Do not edit.",
        '#include "registries.h"',
        "",
        "#include <cmath>",
        "",
        "namespace pytec3d::registry",
        "{",
        "",
    ]
    out = [
        "// Generated by src/tools/generate_registries.py from src/pytec3d/resources/*.json. Do not edit.",
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

    emit_enum(out, src, "Domain", registries["domains"], "Physical regions (domains.json).")
    out.append("/// A set of domains: bit i is the Domain with index i.")
    out.append("using DomainMask = std::uint64_t;")
    out.append("inline constexpr DomainMask bit(Domain domain) { return DomainMask(1) << static_cast<unsigned int>(domain); }")
    out.append("")

    emit_enum(out, src, "FieldId", list(registries["fields"]), "Fields and grid arrays (fields.json).")
    out.append("inline constexpr std::array<bool, kNumFieldId> kFieldIsVector = {")
    for kind in registries["fields"].values():
        out.append(f"    {'true' if kind == 'vector' else 'false'},")
    out.append("};")
    out.append("inline constexpr bool is_vector(FieldId id) { return kFieldIsVector[static_cast<std::size_t>(id)]; }")
    out.append("")

    variables = registries["variables"]
    emit_enum(out, src, "VariableId", list(variables), "libMesh variable labels the fields can have, in system order.")
    out.append("enum class VariableKind : unsigned char")
    out.append("{")
    out.append("    Scalar,    ///< A scalar field.")
    out.append("    Vector,    ///< A whole vector field: Lagrange (a DOF per component on each node) or NEDELEC_ONE.")
    out.append("    Component, ///< One component of a vector field.")
    out.append("};")
    out.append("struct VariableInfo")
    out.append("{")
    out.append("    std::string_view label;     ///< Of the libMesh variable.")
    out.append("    std::string_view api_label; ///< Of the grid arrays (solution and _ constraint arrays).")
    out.append("    FieldId field;")
    out.append("    unsigned char component;    ///< Physical component of a Component variable.")
    out.append("    VariableKind kind;")
    out.append("    std::string_view nodal_load; ///< Point array loading it; empty if none.")
    out.append("};")
    out.append(f"inline constexpr std::array<VariableInfo, {len(variables)}> kVariables = {{{{")
    for variable in variables.values():
        load = registries["nodal_loads"].get(variable["label"], "")
        out.append(f"    {{{cpp_string(variable['label'])}, {cpp_string(variable['api'])}, FieldId::{variable['field']}, "
                   f"{variable['component']}, VariableKind::{variable['kind']}, {cpp_string(load)}}},")
    out.append("}};")
    out.append("inline constexpr const VariableInfo &info(VariableId id) { return kVariables[static_cast<std::size_t>(id)]; }")
    out.append("")

    types = registries["types"]
    emit_enum(out, src, "ElementTypeId", list(types), "Element types (element_types.json).")
    out.append("inline constexpr std::array<unsigned int, kNumElementTypeId> kElementTypeDimensions = {")
    for element_type in types.values():
        out.append(f"    {element_type.dimension},")
    out.append("};")
    out.append("inline constexpr unsigned int dimension(ElementTypeId id) { return kElementTypeDimensions[static_cast<std::size_t>(id)]; }")
    out.append("")

    # The material and real-constant properties are a model data key each.
    for enum, key, doc in (("MaterialProperty", "material_properties", 'Material properties ("mat." parameters).'),
                           ("RealConstant", "real_constants", 'Real constants ("rc." parameters).')):
        entries = registries[key]
        emit_enum(out, src, enum, [e["name"] for e in entries], doc)
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
    emit_enum(out, src, "ParameterId", [p["name"] for p in parameters], "Parameters of all the element types: declared first, then expressions.")
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
    out.append("inline constexpr std::array<std::array<DomainMask, kNumParameterId>, kNumElementTypeId> kParameterDomains = {{")
    for type_name, element_type in types.items():
        out.append(f"    {{{{ // {type_name}")
        for parameter, mask in zip(parameters, element_type.masks):
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
    emit_defaults(out, src, types)


    curves = registries["curves"]
    emit_enum(out, src, "CurveId", [f"{law}_{x}_{y}" for law, x, y in curves],
              "Curves of the laws: the arrays x and y of a law entry of the model data \"laws\" section.")
    out.append("struct CurveSource")
    out.append("{")
    out.append("    std::string_view law;")
    out.append("    std::string_view x;")
    out.append("    std::string_view y;")
    out.append("};")
    out.append(f"inline constexpr std::array<CurveSource, {len(curves)}> kCurveSources = {{{{")
    for law, x, y in curves:
        out.append(f"    {{{cpp_string(law)}, {cpp_string(x)}, {cpp_string(y)}}},")
    out.append("}};")
    out.append("")

    laws = registries["laws"]
    emit_enum(out, src, "LawId", list(laws), "Nonlinear laws (nonlinear_laws.json).")
    out.append("/// What the ids of the model data \"laws\" section are, for a law: MATERIAL, CONSTANT or TYPE ids.")
    out.append("enum class LawSource : unsigned char")
    out.append("{")
    out.append("    Material, ///< materials (mat.)")
    out.append("    Constant, ///< constants (rc.)")
    out.append("    Cell,     ///< element_types, by TYPE id (me.)")
    out.append("};")
    out.append("struct LawInfo")
    out.append("{")
    out.append("    LawSource source;")
    out.append("    CurveId curve; ///< COUNT for a law without curve.")
    out.append("    std::string_view input; ///< The field the law is a function of (its curve's x).")
    out.append("};")
    out.append(f"inline constexpr std::array<LawInfo, {len(laws)}> kLaws = {{{{")
    for law in laws.values():
        curve = f"CurveId::{'_'.join(law['curve'])}" if law["curve"] else "CurveId::COUNT"
        out.append(f"    {{LawSource::{LAW_SOURCES[law['source']]}, {curve}, {cpp_string(law['input'])}}},")
    out.append("}};")
    out.append("inline constexpr const LawInfo &info(LawId id) { return kLaws[static_cast<std::size_t>(id)]; }")
    out.append("/// A set of laws: bit i is the LawId with index i.")
    out.append("using LawMask = std::uint32_t;")
    out.append("inline constexpr LawMask bit(LawId law) { return LawMask(1) << static_cast<unsigned int>(law); }")
    out.append("")
    law_outputs = registries["law_outputs"]
    emit_enum(out, src, "LawOutputId", list(law_outputs), "Outputs of the nonlinear laws, used in kernels like a parameter.")
    out.append("struct LawOutputInfo")
    out.append("{")
    out.append("    LawId law;")
    out.append("    ParameterType type;")
    out.append("};")
    out.append(f"inline constexpr std::array<LawOutputInfo, {len(law_outputs)}> kLawOutputs = {{{{")
    for output in law_outputs.values():
        out.append(f"    {{LawId::{output['law']}, ParameterType::{output['type']}}},")
    out.append("}};")
    out.append("inline constexpr const LawOutputInfo &info(LawOutputId id) { return kLawOutputs[static_cast<std::size_t>(id)]; }")
    out.append("")
    out.append("} // namespace pytec3d::registry")
    out.append("")
    src.append("} // namespace pytec3d::registry")
    src.append("")
    return "\n".join(out), "\n".join(src)


# ---- formulations.h --------------------------------------------------------------------------------------


def span(name, count):
    return f"{{{name}, {count}}}" if count else "{nullptr, 0}"


def emit_laws(out, types, registries):
    out.append("/**")
    out.append(" * @brief Outputs of the law @p law of the element type @p type at one quadrature point.")
    out.append(" *")
    out.append(" * @p input is the value of the law's input field at the quadrature point; @p curve_data the curve of the row selecting the law (null without curve); @p source")
    out.append(" * gives the parameters (as for evaluate_expression); @p outputs receives one value per output of the law,")
    out.append(" * in the order of LawSpec::outputs.")
    out.append(" */")
    out.append("template <typename Source>")
    out.append("inline void evaluate_law(registry::ElementTypeId type,")
    out.append("                         registry::LawId law,")
    out.append("                         const ParameterValue &input,")
    out.append("                         [[maybe_unused]] const CurveData *curve_data,")
    out.append("                         [[maybe_unused]] Source &source,")
    out.append("                         ParameterValue *outputs)")
    out.append("{")
    out.append("    using registry::ParameterId;")
    emit_scope(out, "    ")
    out.append("    switch (type)")
    out.append("    {")
    for type_name, element_type in types.items():
        if not element_type.laws:
            continue
        out.append(f"    case registry::ElementTypeId::{type_name}:")
        out.append("        switch (law)")
        out.append("        {")
        for law in element_type.laws:
            input_field = law["input"]
            out.append(f"        case registry::LawId::{law['name']}:")
            out.append("        {")
            emit_constants(out, element_type, "            ")
            if law["curve_key"]:
                out.append("            [[maybe_unused]] const CurveData &curve = *curve_data;")
            accessor = "input.vector()" if input_field["value"] == V else "input.scalar()"
            out.append(f"            const {CPP_TYPE[input_field['value']]} &{input_field['name']} = {accessor};")
            for local in law["locals"]:
                out.append(f"            const auto {local['name']} = {local['code']};")
            for output in law["outputs"]:
                out.append(f"            outputs[{output['index']}] = expression_support::expression_value<ParameterType::{output['type']}>({output['code']});")
            out.append("            return;")
            out.append("        }")
        out.append("        default:")
        out.append("            break;")
        out.append("        }")
        out.append("        break;")
    out.append("    default:")
    out.append("        break;")
    out.append("    }")
    out.append("}")
    out.append("")


def render_formulations(registries):
    types = registries["types"]
    out = [
        "// Generated by src/tools/generate_registries.py from src/pytec3d/resources/element_types.json and nonlinear_laws.json. Do not edit.",
        "#pragma once",
        "",
        "#include <array>",
        "#include <cmath>",
        "",
        '#include "curve.h"',
        '#include "expression_support.h"',
        '#include "formulation_spec.h"',
        '#include "registries.h"',
        "",
        "namespace pytec3d::formulation",
        "{",
        "namespace generated",
        "{",
    ]
    domains = registries["domains"]
    for type_name, t in types.items():
        p = f"k{type_name}"
        out.append(f"// ---- {type_name} " + "-" * max(0, 100 - len(type_name)))
        out.append(f"inline constexpr FieldSpec {p}Fields[] = {{")
        for f in t.fields:
            components = list(f["components"]) + [0] * (3 - len(f["components"]))
            family = f"libMesh::{f['family']}" if f["family"] else "libMesh::INVALID_FE"
            order = f"OrderRule::{f['order']}"
            load = f"registry::FieldId::{f['nodal_load']}" if f["nodal_load"] else "registry::FieldId::COUNT"
            out.append(f"    {{registry::FieldId::{f['name']}, FieldKind::{FIELD_KINDS[f['kind']]}, ValueType::{f['value']}, "
                       f"{len(f['components'])}, {{{', '.join(map(str, components))}}}, {f['fe'] if f['fe'] is not None else -1}, "
                       f"{family}, {order}, registry::VariableId::{f['variable']}, "
                       f"{f['domains']}ull, {load}}}, // {domain_names(f['domains'], domains)}")
        out.append("};")
        if t.fes:
            out.append(f"inline constexpr FESpec {p}FEs[] = {{")
            for fe in t.fes:
                out.append(f"    {{libMesh::{fe['family']}, {'true' if fe['vector'] else 'false'}, OrderRule::{fe['order']}}},")
            out.append("};")
        if t.nablas:
            out.append(f"inline constexpr NablaSpec {p}Nablas[] = {{")
            for n in t.nablas:
                out.append(f"    {{{n['fe']}, Operator::{n['op']}, {n['component']}, {'true' if n['out_of_plane'] else 'false'}, "
                           f"ValueType::{n['output']}}},")
            out.append("};")
        if t.kernels:
            out.append(f"inline constexpr KernelSpec {p}Kernels[] = {{")
            for k in t.kernels:
                c = k["coefficient"]
                if c["kind"] == "Parameter":
                    index = registries["parameter_index"][c["name"]]
                elif c["kind"] == "LawOutput":
                    index = registries["output_index"][c["name"]]
                else:
                    index = 0
                trial = k["trial"] if k["trial"] is not None else -1
                test = k["test"] if k["test"] is not None else -1
                output = k["output"] or S
                out.append(f"    {{{trial}, {test}, {{CoefficientKind::{c['kind']}, {index}, ValueType::{c['type']}, {c['sign']}}}, "
                           f"ValueType::{output}}}, // {c['name'] or '1'}")
            out.append("};")
        if t.terms:
            out.append(f"inline constexpr TermSpec {p}Terms[] = {{")
            for term in t.terms:
                out.append(f"    {{{term['kernel']}, {term['trial_field']}, {term['test_field']}, {term['trial_dt']}, "
                           f"Stage::{term['stage']}, {term['law']}, {term['replaced_by']}u, {term['domains']}ull}},")
            out.append("};")
        if t.evaluations:
            out.append(f"inline constexpr EvaluationSpec {p}Evaluations[] = {{")
            for e in t.evaluations:
                field = e["field"]["index"] if e["field"] else -1
                out.append(f"    {{EvaluationKind::{e['kind']}, {field}, {e['mass_kernel']}, {e['first_term']}, {e['n_terms']}}},")
            out.append("};")
        for law in t.laws:
            outputs = law["outputs"]
            out.append(f"inline constexpr registry::LawOutputId {p}{law['name']}Outputs[] = {{"
                       + ", ".join(f"registry::LawOutputId::{o['name']}" for o in outputs) + "};")
        if t.laws:
            out.append(f"inline constexpr LawSpec {p}Laws[] = {{")
            for law in t.laws:
                out.append(f"    {{registry::LawId::{law['name']}, {law['evaluation']}, "
                           f"{span(f'{p}{law['name']}Outputs', len(law['outputs']))}}},")
            out.append("};")
        out.append("")

    out.append("} // namespace generated")
    out.append("")
    out.append("inline constexpr std::array<ElementTypeSpec, registry::kNumElementTypeId> kElementTypes = {{")
    for type_name, t in types.items():
        p = f"generated::k{type_name}"
        length = registries["parameter_index"][t.length] if t.length else -1
        out.append(f"    {{registry::ElementTypeId::{type_name}, {t.dimension}, {t.time_order}, {length}, "
                   f"{span(p + 'Fields', len(t.fields))}, {span(p + 'FEs', len(t.fes))}, {span(p + 'Nablas', len(t.nablas))}, "
                   f"{span(p + 'Kernels', len(t.kernels))}, {span(p + 'Terms', len(t.terms))}, "
                   f"{span(p + 'Evaluations', len(t.evaluations))}, {span(p + 'Laws', len(t.laws))}}},")
    out.append("}};")
    out.append("inline constexpr const ElementTypeSpec &spec(registry::ElementTypeId id) { return kElementTypes[static_cast<std::size_t>(id)]; }")
    out.append("")
    emit_laws(out, types, registries)
    out.append("} // namespace pytec3d::formulation")
    out.append("")
    return "\n".join(out)


def write_if_changed(path, text):
    try:
        with open(path, encoding="utf-8") as handle:
            if handle.read() == text:
                return  # Leave the timestamp alone: nothing recompiles.
    except FileNotFoundError:
        pass
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--resources", required=True, help="Directory with the JSON registries.")
    parser.add_argument("--output", required=True, help="Directory for registries.h, registries.cpp and formulations.h.")
    args = parser.parse_args()

    registries = read_registries(args.resources)
    header, source = render_registries(registries)
    write_if_changed(f"{args.output}/registries.h", header)
    write_if_changed(f"{args.output}/registries.cpp", source)
    write_if_changed(f"{args.output}/formulations.h", render_formulations(registries))


if __name__ == "__main__":
    main()
