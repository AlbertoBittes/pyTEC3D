/**
 * @file element_type_parser.cpp
 * @brief ElementTypeParser: JSON parsing, validation and ElementType instantiation.
 */
#include "element_type_parser.h"

#include <algorithm>
#include <cctype>
#include <functional>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <unordered_map>

#include <libmesh/fe_interface.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/string_to_enum.h>

#include "embedded_resources.h"
#include "logging.h"
#include "parameter_columns.h"

namespace pytec3d
{
    namespace
    {
        using json = nlohmann::ordered_json;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementTypeParser: " + message);
        }

        json parse_json(const std::string &text, const char *what)
        {
            try
            {
                return json::parse(text);
            }
            catch (const json::exception &ex)
            {
                fail(std::string("cannot parse ") + what + ": " + ex.what());
            }
        }

        std::vector<std::string> read_string_list(const json &value, const std::string &context)
        {
            if (!value.is_array())
                fail(context + " must be an array of strings.");
            std::vector<std::string> out;
            for (const auto &item : value)
            {
                if (!item.is_string())
                    fail(context + " must be an array of strings.");
                out.push_back(item.get<std::string>());
            }
            return out;
        }

        unsigned int component_index(const std::string &label, const std::string &context)
        {
            if (label == "x")
                return 0;
            if (label == "y")
                return 1;
            if (label == "z")
                return 2;
            fail(context + ": invalid component '" + label + "' (expected x, y or z).");
        }

        std::string strip_spaces(const std::string &text)
        {
            std::string out;
            for (char c : text)
                if (!std::isspace(static_cast<unsigned char>(c)))
                    out.push_back(c);
            return out;
        }

        bool is_identifier(const std::string &text)
        {
            if (text.empty() || !(std::isalpha(static_cast<unsigned char>(text[0])) || text[0] == '_'))
                return false;
            return std::all_of(text.begin(), text.end(), [](char c)
                               { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; });
        }

        std::vector<unsigned int> read_domains(const json &definition, const std::string &context, const std::vector<std::string> &domains)
        {
            std::vector<unsigned int> out;
            if (!definition.is_object() || !definition.contains("domain"))
                return out;
            for (const std::string &name : read_string_list(definition["domain"], context + " domain"))
            {
                const auto it = std::find(domains.begin(), domains.end(), name);
                if (it == domains.end())
                    fail(context + ": unknown domain '" + name + "'.");
                out.push_back(static_cast<unsigned int>(it - domains.begin()));
            }
            return out;
        }

        // ---- expressions -------------------------------------------------------------------------------------

        struct OperandText
        {
            NablaOperator op = NablaOperator::Identity;
            bool time_derivative = false;
            std::string field;
        };

        // "A", "curl(A)", "dt(A)", "grad(dt(V))": at most one spatial operator and one dt().
        OperandText parse_operand(const std::string &expression, const std::string &context)
        {
            std::string text = strip_spaces(expression);
            OperandText out;
            bool has_operator = false;

            for (std::size_t open = text.find('('); open != std::string::npos; open = text.find('('))
            {
                if (text.back() != ')')
                    fail(context + ": unbalanced parentheses in '" + expression + "'.");
                const std::string name = text.substr(0, open);
                const std::string inner = text.substr(open + 1, text.size() - open - 2);

                if (name == "dt")
                {
                    if (out.time_derivative)
                        fail(context + ": repeated dt() in '" + expression + "'.");
                    out.time_derivative = true;
                }
                else
                {
                    if (has_operator)
                        fail(context + ": only one of grad, curl, div is allowed in '" + expression + "'.");
                    if (name == "grad")
                        out.op = NablaOperator::Grad;
                    else if (name == "curl")
                        out.op = NablaOperator::Curl;
                    else if (name == "div")
                        out.op = NablaOperator::Div;
                    else
                        fail(context + ": unknown operator '" + name + "' in '" + expression + "'.");
                    has_operator = true;
                }
                text = inner;
            }

            if (!is_identifier(text))
                fail(context + ": '" + expression + "' does not reduce to a field name.");
            out.field = text;
            return out;
        }

        // ---- parameters --------------------------------------------------------------------------------------

        struct ParameterKey
        {
            std::string name; ///< Without suffix.
            ParameterType type = ParameterType::Scalar;
            std::string text; ///< As written: "nu[][]", "Js[]", "amps()", "turns".
        };

        // "name", "name[]", "name[][]" or "name()".
        ParameterKey split_key(const std::string &expression, const std::string &context)
        {
            const std::string text = strip_spaces(expression);
            std::size_t end = 0;
            while (end < text.size() && (std::isalnum(static_cast<unsigned char>(text[end])) || text[end] == '_'))
                ++end;

            ParameterKey key;
            key.name = text.substr(0, end);
            key.text = text;
            const std::string suffix = text.substr(end);
            if (suffix.empty())
                key.type = ParameterType::Scalar;
            else if (suffix == "[]")
                key.type = ParameterType::Vector;
            else if (suffix == "[][]")
                key.type = ParameterType::Tensor;
            else if (suffix == "()")
                key.type = ParameterType::Unknown;
            else
                fail(context + ": invalid parameter '" + expression + "'.");
            if (!is_identifier(key.name))
                fail(context + ": invalid parameter name '" + expression + "'.");
            return key;
        }

        std::string lowercase(std::string text)
        {
            std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return text;
        }

        struct CatalogEntry
        {
            const char *name; ///< Lowercase.
            ParameterType type;
            std::vector<unsigned int> columns;
        };

        const std::vector<CatalogEntry> &material_catalog()
        {
            static const std::vector<CatalogEntry> catalog = {
                {"mur", ParameterType::Vector, {MAT_MUR_X, MAT_MUR_Y, MAT_MUR_Z}},
                {"relative_permeability", ParameterType::Vector, {MAT_MUR_X, MAT_MUR_Y, MAT_MUR_Z}},
                {"rho", ParameterType::Vector, {MAT_RHO_X, MAT_RHO_Y, MAT_RHO_Z}},
                {"resistivity", ParameterType::Vector, {MAT_RHO_X, MAT_RHO_Y, MAT_RHO_Z}},
                {"m", ParameterType::Vector, {MAT_M_X, MAT_M_Y, MAT_M_Z}},
                {"magnetization", ParameterType::Vector, {MAT_M_X, MAT_M_Y, MAT_M_Z}},
            };
            return catalog;
        }

        const std::vector<CatalogEntry> &constant_catalog()
        {
            static const std::vector<CatalogEntry> catalog = {
                {"thk", ParameterType::Scalar, {CONST_THK}},
                {"thickness", ParameterType::Scalar, {CONST_THK}},
                {"cs", ParameterType::Scalar, {CONST_CS}},
                {"cross_section", ParameterType::Scalar, {CONST_CS}},
                {"turns", ParameterType::Scalar, {CONST_TURNS}},
                {"vol", ParameterType::Scalar, {CONST_VOL}},
                {"volume", ParameterType::Scalar, {CONST_VOL}},
                {"rad", ParameterType::Scalar, {CONST_RAD}},
                {"radius", ParameterType::Scalar, {CONST_RAD}},
                {"r", ParameterType::Scalar, {CONST_R}},
                {"resistance", ParameterType::Scalar, {CONST_R}},
                {"s", ParameterType::Scalar, {CONST_S}},
                {"symmetry", ParameterType::Scalar, {CONST_S}},
                {"t", ParameterType::Vector, {CONST_T_X, CONST_T_Y, CONST_T_Z}},
                {"current_direction", ParameterType::Vector, {CONST_T_X, CONST_T_Y, CONST_T_Z}},
            };
            return catalog;
        }

        bool is_builtin_name(const std::string &name)
        {
            for (const char *builtin : {"pi", "mu0", "dimension", "axisymmetric", "ex", "ey", "ez"})
                if (name == builtin)
                    return true;
            return false;
        }

        // "mat.mur[]", "rc.turns", "me.Js[]", "mn.AMPS()": the grid and model data that parameters read.
        ParameterDefinition parse_declaration(const std::string &declaration)
        {
            const std::string context = "parameter '" + declaration + "'";
            const std::string text = strip_spaces(declaration);
            const std::size_t dot = text.find('.');
            if (dot == std::string::npos)
                fail(context + " must start with mat., rc., me. or mn.");

            const std::string prefix = text.substr(0, dot);
            const ParameterKey key = split_key(text.substr(dot + 1), context);

            ParameterDefinition definition;
            definition.name = key.text;
            definition.type = key.type;

            if (prefix == "mat" || prefix == "rc")
            {
                if (key.type == ParameterType::Unknown)
                    fail(context + ": a " + (prefix == "mat" ? "material" : "real-constant") + " parameter has a fixed type.");
                const auto &catalog = prefix == "mat" ? material_catalog() : constant_catalog();
                const std::string wanted = lowercase(key.name);
                const auto entry = std::find_if(catalog.begin(), catalog.end(), [&wanted](const CatalogEntry &e) { return wanted == e.name; });
                if (entry == catalog.end())
                    fail(context + ": unknown " + (prefix == "mat" ? "material" : "real constant") + " '" + key.name + "'.");
                if (entry->type != key.type)
                    fail(context + ": '" + key.name + "' is a " + type_name(entry->type) + ", write it " +
                         (entry->type == ParameterType::Vector ? "with []" : "without a suffix") + ".");
                definition.source = prefix == "mat" ? ParameterSource::Material : ParameterSource::Constant;
                definition.columns = entry->columns;
            }
            else if (prefix == "me" || prefix == "mn")
            {
                definition.source = prefix == "me" ? ParameterSource::Cell : ParameterSource::Point;
                definition.array = key.name;
                if (definition.source == ParameterSource::Point && key.type != ParameterType::Unknown && key.type != ParameterType::Scalar)
                    fail(context + ": point data is a scalar array per node.");
            }
            else
                fail(context + ": unknown source '" + prefix + "' (expected mat, rc, me or mn).");
            return definition;
        }

        // What a kernel can name: parameters, laws and, for nodal loads, point data.
        struct ParameterTable
        {
            const std::vector<ParameterDefinition> &parameters;
            const std::unordered_map<std::string, unsigned int> &index;
            const std::vector<NonlinearLaw> &laws;
            const std::unordered_map<std::string, unsigned int> &law_index;
            const std::map<std::string, std::string> &point_arrays; ///< "AMPS()" -> array name.
        };

        // "-nu[][]", "BH()", "turns": an optional sign and a parameter or law; empty is 1.
        ParameterRef parse_parameter(const std::string &expression, const std::string &context, const ParameterTable &table)
        {
            ParameterRef out;
            std::string text = strip_spaces(expression);
            if (text.empty())
                return out;

            if (text[0] == '-' || text[0] == '+')
            {
                out.sign = text[0] == '-' ? -1.0 : 1.0;
                text.erase(0, 1);
            }
            out.name = text;

            if (text.size() > 2 && text.compare(text.size() - 2, 2, "()") == 0)
            {
                const auto law = table.law_index.find(text.substr(0, text.size() - 2));
                if (law != table.law_index.end())
                {
                    out.kind = ParameterKind::Law;
                    out.index = law->second;
                    out.type = ParameterType::Tensor;
                    return out;
                }
            }

            const auto it = table.index.find(text);
            if (it != table.index.end())
            {
                out.kind = ParameterKind::Parameter;
                out.index = it->second;
                out.type = table.parameters[it->second].type;
                return out;
            }
            if (table.point_arrays.count(text))
                fail(context + ": '" + text + "' is point data; it can only be the nodal_loads of a dof field.");
            fail(context + ": unknown parameter '" + text + "'.");
        }

        // ---- fields ------------------------------------------------------------------------------------------

        libMesh::FEFamily default_family(FieldKind kind, bool vector)
        {
            switch (kind)
            {
            case FieldKind::Dof:
                return vector ? libMesh::LAGRANGE_VEC : libMesh::LAGRANGE;
            case FieldKind::Nodal:
                return vector ? libMesh::L2_LAGRANGE_VEC : libMesh::L2_LAGRANGE;
            case FieldKind::Element:
                return vector ? libMesh::MONOMIAL_VEC : libMesh::MONOMIAL;
            }
            fail("unknown field kind.");
        }

        bool is_lagrange_type(libMesh::FEFamily family)
        {
            return family == libMesh::LAGRANGE || family == libMesh::LAGRANGE_VEC || family == libMesh::L2_LAGRANGE ||
                   family == libMesh::L2_LAGRANGE_VEC;
        }

        libMesh::FEFamily scalar_family(libMesh::FEFamily vector_family)
        {
            switch (vector_family)
            {
            case libMesh::LAGRANGE_VEC:
                return libMesh::LAGRANGE;
            case libMesh::L2_LAGRANGE_VEC:
                return libMesh::L2_LAGRANGE;
            case libMesh::MONOMIAL_VEC:
                return libMesh::MONOMIAL;
            default:
                fail("FE family has no component-wise scalar family.");
            }
        }

        // One libMesh scalar variable per component, except for an H(curl) field that is a single variable.
        std::vector<LibMeshVariable> libmesh_variables(const Field &field)
        {
            static const char *const suffixes[3] = {"_x", "_y", "_z"};
            std::vector<LibMeshVariable> out;

            if (field.vector_output && field.family == libMesh::NEDELEC_ONE)
                out.push_back({field.name + "z", field.family, field.order});
            else if (field.vector_output)
                for (unsigned int component : field.active_components)
                    out.push_back({field.name + suffixes[component], scalar_family(field.family), field.order});
            else if (field.physical_vector)
                out.push_back({field.name + suffixes[field.active_components.at(0)], field.family, field.order});
            else
                out.push_back({field.name, field.family, field.order});
            return out;
        }

        Field parse_field(const std::string &type_name,
                          const std::string &name,
                          const json &definition,
                          FieldKind kind,
                          unsigned int dimension,
                          const ElementTypeOptions &options,
                          const std::map<std::string, bool> &field_is_vector,
                          const std::vector<std::string> &domains,
                          const std::map<std::string, std::string> &point_arrays)
        {
            const std::string context = "element type '" + type_name + "', field '" + name + "'";

            const auto vector_it = field_is_vector.find(name);
            if (vector_it == field_is_vector.end())
                fail(context + " is not declared in fields.json.");
            if (!definition.is_object())
                fail(context + " must be an object.");

            Field field;
            field.name = name;
            field.kind = kind;
            field.physical_vector = vector_it->second;

            std::vector<unsigned int> requested;
            if (definition.contains("active_components"))
                for (const std::string &label : read_string_list(definition["active_components"], context + " active_components"))
                    requested.push_back(component_index(label, context));
            std::sort(requested.begin(), requested.end());
            if (std::adjacent_find(requested.begin(), requested.end()) != requested.end())
                fail(context + ": duplicate active_components.");
            if (!field.physical_vector && !requested.empty())
                fail(context + ": a scalar field cannot define active_components.");
            if (!requested.empty() && requested.back() >= 3)
                fail(context + ": invalid active_components.");

            // A vector field has every component of the dimension; a single component makes it a scalar.
            if (!field.physical_vector)
                field.active_components = {0};
            else if (requested.size() == 1)
                field.active_components = requested;
            else
            {
                field.active_components.resize(dimension);
                std::iota(field.active_components.begin(), field.active_components.end(), 0u);
            }
            const bool multi_component = field.active_components.size() > 1;

            field.family = default_family(kind, multi_component);
            if (definition.contains("finite_elements"))
            {
                if (!definition["finite_elements"].is_string())
                    fail(context + ": 'finite_elements' must be a string.");
                field.family = libMesh::Utility::string_to_enum<libMesh::FEFamily>(definition["finite_elements"].get<std::string>());
                // "LAGRANGE" on a multi-component field means component-wise Lagrange.
                if (multi_component && field.family == libMesh::LAGRANGE)
                    field.family = libMesh::LAGRANGE_VEC;
            }
            field.vector_output = libMesh::FEInterface::field_type(field.family) == libMesh::TYPE_VECTOR;
            if (field.vector_output != multi_component)
                fail(context + ": FE family does not match the number of active components.");

            if (is_lagrange_type(field.family))
                field.order = static_cast<libMesh::Order>(options.order);
            else if (kind == FieldKind::Element)
                field.order = libMesh::CONSTANT;
            else
                field.order = libMesh::FIRST;

            if (definition.contains("nodal_loads"))
            {
                if (kind != FieldKind::Dof || !definition["nodal_loads"].is_string())
                    fail(context + ": 'nodal_loads' is a string and only applies to dof fields.");
                const std::string key = strip_spaces(definition["nodal_loads"].get<std::string>());
                const auto point = point_arrays.find(key);
                if (point == point_arrays.end())
                    fail(context + ": nodal_loads '" + key + "' is not a point-data (mn.) parameter.");
                field.nodal_loads = point->second;
            }
            field.domains = read_domains(definition, context, domains);
            field.libmesh_variables = libmesh_variables(field);
            return field;
        }

        std::vector<Field> parse_fields(const std::string &type_name,
                                        const json &type_json,
                                        const char *section,
                                        FieldKind kind,
                                        unsigned int dimension,
                                        const ElementTypeOptions &options,
                                        const std::map<std::string, bool> &field_is_vector,
                                        const std::vector<std::string> &domains,
                                        const std::map<std::string, std::string> &point_arrays)
        {
            std::vector<Field> out;
            if (!type_json.contains(section))
                return out;
            const json &fields = type_json[section];
            if (!fields.is_object())
                fail("element type '" + type_name + "': '" + section + "' must be an object.");
            out.reserve(fields.size());
            for (auto it = fields.begin(); it != fields.end(); ++it)
            {
                Field field = parse_field(type_name, it.key(), it.value(), kind, dimension, options, field_is_vector, domains, point_arrays);
                field.index = static_cast<unsigned int>(out.size());
                out.push_back(std::move(field));
            }
            return out;
        }

        // ---- kernels -----------------------------------------------------------------------------------------

        const Field *find_by_name(const std::vector<Field> &fields, const std::string &name)
        {
            for (const Field &field : fields)
                if (field.name == name)
                    return &field;
            return nullptr;
        }

        // Creates (or finds) the operator of a field; supplied by the element type being built.
        using NablaFactory = std::function<const Nabla *(const Field &, NablaOperator)>;

        OperandRef make_operand(const std::string &expression,
                                const std::vector<Field> &dof_fields,
                                const std::string &context,
                                bool allow_dt,
                                const NablaFactory &nabla_of)
        {
            const OperandText text = parse_operand(expression, context);
            const Field *field = find_by_name(dof_fields, text.field);
            if (!field)
                fail(context + ": '" + text.field + "' is not a dof field of the element type.");
            if (text.time_derivative && !allow_dt)
                fail(context + ": dt() is only allowed on the trial side.");
            return {nabla_of(*field, text.op), text.time_derivative};
        }

        std::string read_string(const json &definition, const char *key, const std::string &context, bool required)
        {
            if (!definition.contains(key))
            {
                if (required)
                    fail(context + " needs a string '" + key + "'.");
                return {};
            }
            if (!definition[key].is_string())
                fail(context + ": '" + key + "' must be a string.");
            return definition[key].get<std::string>();
        }

        std::vector<DofKernel> parse_dof_kernels(const std::string &type_name,
                                                 const json &type_json,
                                                 const std::vector<Field> &dof_fields,
                                                 const std::vector<std::string> &domains,
                                                 const ParameterTable &table,
                                                 const NablaFactory &nabla_of)
        {
            std::vector<DofKernel> out;
            if (!type_json.contains("dof_kernels"))
                return out;
            const json &kernels = type_json["dof_kernels"];
            if (!kernels.is_object())
                fail("element type '" + type_name + "': 'dof_kernels' must be an object.");

            for (auto it = kernels.begin(); it != kernels.end(); ++it)
            {
                const std::string context = "element type '" + type_name + "', dof kernel '" + it.key() + "'";
                const json &definition = it.value();
                if (!definition.is_object())
                    fail(context + " must be an object.");

                DofKernel kernel;
                kernel.name = it.key();
                kernel.test = make_operand(read_string(definition, "test", context, true), dof_fields, context + " test", false, nabla_of);
                kernel.has_trial = definition.contains("trial");
                if (kernel.has_trial)
                    kernel.trial = make_operand(read_string(definition, "trial", context, true), dof_fields, context + " trial", true, nabla_of);
                kernel.parameter = parse_parameter(read_string(definition, "coeff", context, false), context, table);
                kernel.domains = read_domains(definition, context, domains);

                const ParameterRef &parameter = kernel.parameter;
                if (kernel.has_trial && parameter.type == ParameterType::Vector)
                    fail(context + ": a vector parameter needs a kernel without trial.");
                if (!kernel.has_trial && parameter.type == ParameterType::Tensor)
                    fail(context + ": a kernel without trial needs a scalar or vector parameter.");
                if (parameter.kind == ParameterKind::Law && (!kernel.has_trial || kernel.trial.time_derivative))
                    fail(context + ": a law needs a trial without dt().");
                out.push_back(std::move(kernel));
            }
            return out;
        }

        std::vector<FieldKernel> parse_field_kernels(const std::string &type_name,
                                                     const json &type_json,
                                                     const char *section,
                                                     const std::vector<Field> &fields,
                                                     const std::vector<Field> &dof_fields,
                                                     const std::vector<std::string> &domains,
                                                     const ParameterTable &table,
                                                     const NablaFactory &nabla_of)
        {
            std::vector<FieldKernel> out;
            if (type_json.contains(section))
            {
                const json &kernels = type_json[section];
                if (!kernels.is_object())
                    fail("element type '" + type_name + "': '" + section + "' must be an object.");

                for (auto it = kernels.begin(); it != kernels.end(); ++it)
                {
                    const std::string context = "element type '" + type_name + "', " + section + " '" + it.key() + "'";
                    const Field *field = find_by_name(fields, it.key());
                    if (!field)
                        fail(context + ": no such field.");
                    if (!it.value().is_array())
                        fail(context + " must be an array of terms.");

                    unsigned int term = 0;
                    for (const json &definition : it.value())
                    {
                        const std::string term_context = context + ", term " + std::to_string(term);
                        if (!definition.is_object())
                            fail(term_context + " must be an object.");

                        FieldKernel kernel;
                        kernel.name = it.key() + "[" + std::to_string(term++) + "]";
                        kernel.field = field;
                        kernel.trial = make_operand(read_string(definition, "trial", term_context, true), dof_fields, term_context + " trial", true, nabla_of);
                        kernel.parameter = parse_parameter(read_string(definition, "coeff", term_context, false), term_context, table);
                        kernel.domains = read_domains(definition, term_context, domains);
                        if (kernel.parameter.type == ParameterType::Vector)
                            fail(term_context + ": a field parameter is a scalar or a tensor.");
                        if (kernel.parameter.kind == ParameterKind::Law && kernel.trial.time_derivative)
                            fail(term_context + ": a law needs a trial without dt().");
                        out.push_back(std::move(kernel));
                    }
                }
            }

            for (const Field &field : fields)
            {
                const bool has_term = std::any_of(out.begin(), out.end(), [&field](const FieldKernel &k)
                                                  { return k.field == &field; });
                if (!has_term)
                    fail("element type '" + type_name + "': field '" + field.name + "' has no kernel in '" + section + "'.");
            }
            return out;
        }
    } // namespace

    ElementTypeParser ElementTypeParser::from_embedded()
    {
        return from_strings(resources::ELEMENT_TYPES_JSON, resources::FIELDS_JSON, resources::DOMAINS_JSON);
    }

    ElementTypeParser ElementTypeParser::from_strings(const std::string &element_types_json,
                                                      const std::string &fields_json,
                                                      const std::string &domains_json)
    {
        ElementTypeParser parser;

        const json domains = parse_json(domains_json, "domains.json");
        if (!domains.is_object() || !domains.contains("domains"))
            fail("domains.json must contain a 'domains' array.");
        parser.domains_ = read_string_list(domains["domains"], "domains.json 'domains'");

        const json fields = parse_json(fields_json, "fields.json");
        if (!fields.is_object())
            fail("fields.json must be an object.");
        for (auto it = fields.begin(); it != fields.end(); ++it)
        {
            if (!it.value().is_object() || !it.value().contains("vector") || !it.value()["vector"].is_boolean())
                fail("fields.json: field '" + it.key() + "' needs a boolean 'vector' entry.");
            parser.field_is_vector_[it.key()] = it.value()["vector"].get<bool>();
        }

        parser.element_types_json_ = parse_json(element_types_json, "element_types.json");
        const json &root = parser.element_types_json_;
        if (!root.is_object() || !root.contains("element_types") || !root["element_types"].is_object())
            fail("element_types.json must contain an 'element_types' object.");

        std::set<std::string> names;
        const auto claim = [&names](const std::string &key, const std::string &context)
        {
            if (is_builtin_name(key) || !names.insert(key).second)
                fail(context + ": '" + key + "' is already defined or reserved.");
        };

        if (root.contains("parameters"))
            for (const std::string &declaration : read_string_list(root["parameters"], "element_types.json 'parameters'"))
            {
                ParameterDefinition definition = parse_declaration(declaration);
                claim(definition.name, "parameter '" + declaration + "'");
                parser.declarations_.push_back(std::move(definition));
            }

        if (root.contains("expressions"))
        {
            if (!root["expressions"].is_object())
                fail("element_types.json: 'expressions' must be an object.");
            for (auto it = root["expressions"].begin(); it != root["expressions"].end(); ++it)
            {
                const std::string context = "expression '" + it.key() + "'";
                const ParameterKey key = split_key(it.key(), context);
                if (key.type == ParameterType::Unknown)
                    fail(context + ": an expression is a scalar, a vector[] or a tensor[][].");
                if (!it.value().is_string())
                    fail(context + " must be a string.");
                claim(key.text, context);
                parser.expressions_.emplace_back(key.text, it.value().get<std::string>());
            }
        }
        return parser;
    }

    bool ElementTypeParser::has_type(const std::string &name) const
    {
        return element_types_json_["element_types"].contains(name);
    }

    std::vector<std::string> ElementTypeParser::type_names() const
    {
        std::vector<std::string> names;
        for (auto it = element_types_json_["element_types"].begin(); it != element_types_json_["element_types"].end(); ++it)
            names.push_back(it.key());
        return names;
    }

    unsigned int ElementTypeParser::domain_index(const std::string &name) const
    {
        const auto it = std::find(domains_.begin(), domains_.end(), name);
        if (it == domains_.end())
            fail("unknown domain '" + name + "'.");
        return static_cast<unsigned int>(it - domains_.begin());
    }

    const std::string &ElementTypeParser::domain_name(unsigned int index) const
    {
        if (index >= domains_.size())
            fail("domain index " + std::to_string(index) + " out of range.");
        return domains_[index];
    }

    std::unique_ptr<ElementType> ElementTypeParser::build(const std::string &name, const ElementTypeOptions &options) const
    {
        if (!has_type(name))
            fail("unknown element type '" + name + "'.");
        const json &type_json = element_types_json_["element_types"][name];
        const std::string context = "element type '" + name + "'";

        if (!type_json.contains("dimension") || !type_json["dimension"].is_number_unsigned())
            fail(context + " needs an unsigned 'dimension'.");
        const unsigned int dimension = type_json["dimension"].get<unsigned int>();
        if (dimension != 2 && dimension != 3)
            fail(context + ": dimension must be 2 or 3.");
        if (options.order < static_cast<unsigned int>(libMesh::FIRST))
            fail(context + ": order must be >= 1.");
        if (options.axisymmetric && dimension != 2)
            fail(context + ": axisymmetric requires a 2D element type.");

        std::unique_ptr<ElementType> type(new ElementType());
        type->name_ = name;
        type->dimension_ = dimension;
        type->options_ = options;
        type->domains_ = domains_;

        // Parameters: the declared ones, then the expressions that build on them in order.
        type->parameters_ = declarations_;
        std::unordered_map<std::string, unsigned int> index;
        std::map<std::string, std::string> point_arrays;
        for (unsigned int i = 0; i < type->parameters_.size(); ++i)
        {
            const ParameterDefinition &definition = type->parameters_[i];
            if (definition.source == ParameterSource::Point)
                point_arrays[definition.name] = definition.array;
            else
                index[definition.name] = i;
        }
        for (const auto &[key, text] : expressions_)
        {
            const std::string expression_context = "element type '" + name + "', expression '" + key + "'";
            const Expression::Resolver resolve = [&](const std::string &parameter) -> std::optional<unsigned int>
            {
                const auto it = index.find(parameter);
                if (it != index.end())
                    return it->second;
                if (point_arrays.count(parameter))
                    fail(expression_context + ": '" + parameter + "' is point data and cannot be used in expressions.");
                return std::nullopt;
            };

            ParameterDefinition definition;
            definition.name = key;
            definition.type = split_key(key, expression_context).type;
            definition.source = ParameterSource::Expression;
            definition.expression = std::make_shared<const Expression>(Expression::parse(text, expression_context, resolve, dimension, options.axisymmetric));
            index[key] = static_cast<unsigned int>(type->parameters_.size());
            type->parameters_.push_back(std::move(definition));
        }

        // Nonlinear laws: used in kernels as "name()".
        std::unordered_map<std::string, unsigned int> law_index;
        if (type_json.contains("nonlinear_laws"))
        {
            const json &laws = type_json["nonlinear_laws"];
            if (!laws.is_object())
                fail(context + ": 'nonlinear_laws' must be an object.");
            for (auto it = laws.begin(); it != laws.end(); ++it)
            {
                const std::string law_context = context + ", law '" + it.key() + "'";
                const json &definition = it.value();
                if (!is_identifier(it.key()) || !definition.is_object())
                    fail(law_context + " must be an object with an identifier as name.");
                if (index.count(it.key() + "()") || point_arrays.count(it.key() + "()"))
                    fail(law_context + ": a parameter has the same name.");

                NonlinearLaw law;
                law.name = it.key();
                if (read_string(definition, "type", law_context, true) != "bh_curve")
                    fail(law_context + ": the only law type is 'bh_curve'.");
                law.type = LawType::BHCurve;

                const std::string linear = strip_spaces(read_string(definition, "linear", law_context, true));
                const auto linear_it = index.find(linear);
                if (linear_it == index.end() || type->parameters_[linear_it->second].type != ParameterType::Tensor)
                    fail(law_context + ": 'linear' must be a tensor parameter.");
                law.linear = linear_it->second;

                law_index[law.name] = static_cast<unsigned int>(type->laws_.size());
                type->laws_.push_back(std::move(law));
            }
        }

        const ParameterTable table{type->parameters_, index, type->laws_, law_index, point_arrays};

        type->dof_fields_ = parse_fields(name, type_json, "dof_fields", FieldKind::Dof, dimension, options, field_is_vector_, domains_, point_arrays);
        type->nodal_fields_ = parse_fields(name, type_json, "nodal_fields", FieldKind::Nodal, dimension, options, field_is_vector_, domains_, point_arrays);
        type->element_fields_ = parse_fields(name, type_json, "element_fields", FieldKind::Element, dimension, options, field_is_vector_, domains_, point_arrays);

        // One quadrature rule and one FE per (family, order), shared by every field using it.
        unsigned int highest_order = static_cast<unsigned int>(libMesh::FIRST);
        for (const Field &field : type->dof_fields_)
            highest_order = std::max(highest_order, static_cast<unsigned int>(field.order));
        type->quadrature_ = std::make_unique<libMesh::QGauss>(dimension, static_cast<libMesh::Order>(2 * highest_order + 1));

        const auto attach_fe = [&type](Field &field)
        {
            const std::pair<int, int> key{static_cast<int>(field.family), static_cast<int>(field.order)};
            auto it = type->fes_.find(key);
            if (it == type->fes_.end())
            {
                std::unique_ptr<libMesh::FEAbstract> fe = libMesh::FEAbstract::build(type->dimension_, libMesh::FEType(field.order, field.family));
                fe->attach_quadrature_rule(type->quadrature_.get());
                it = type->fes_.emplace(key, std::move(fe)).first;
            }
            field.fe = it->second.get();
        };
        for (Field &field : type->dof_fields_)
            attach_fe(field);
        for (Field &field : type->nodal_fields_)
            attach_fe(field);

        // Fields are complete and never move again: operators and kernels point into them.
        const NablaFactory nabla_of = [&type, &options](const Field &field, NablaOperator op) -> const Nabla *
        {
            for (const Nabla &existing : type->nablas_)
                if (&existing.field() == &field && existing.op() == op)
                    return &existing;
            type->nablas_.emplace_back(op, field, type->dimension_, options.axisymmetric);
            return &type->nablas_.back();
        };
        for (const std::vector<Field> *fields : {&type->dof_fields_, &type->nodal_fields_})
            for (const Field &field : *fields)
                nabla_of(field, NablaOperator::Identity);

        type->dof_kernels_ = parse_dof_kernels(name, type_json, type->dof_fields_, domains_, table, nabla_of);
        type->nodal_kernels_ = parse_field_kernels(name, type_json, "nodal_kernels", type->nodal_fields_, type->dof_fields_, domains_, table, nabla_of);
        type->element_kernels_ = parse_field_kernels(name, type_json, "element_kernels", type->element_fields_, type->dof_fields_, domains_, table, nabla_of);

        type->geometry_fe_ = libMesh::FEBase::build(dimension, libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE));
        type->geometry_fe_->attach_quadrature_rule(type->quadrature_.get());
        type->geometry_fe_->get_JxW();
        type->geometry_fe_->get_xyz();

        // The FEs to reinit in each domain, each once.
        type->fes_by_domain_.resize(domains_.size());
        for (unsigned int domain = 0; domain < domains_.size(); ++domain)
        {
            std::vector<libMesh::FEAbstract *> &active = type->fes_by_domain_[domain];
            for (const std::vector<Field> *fields : {&type->dof_fields_, &type->nodal_fields_})
                for (const Field &field : *fields)
                    if (field.active_in(domain) && std::find(active.begin(), active.end(), field.fe) == active.end())
                        active.push_back(field.fe);
        }
        return type;
    }
} // namespace pytec3d
