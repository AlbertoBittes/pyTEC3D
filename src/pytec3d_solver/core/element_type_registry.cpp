/**
 * @file element_type_registry.cpp
 * @brief ElementTypeRegistry: parsing and validation of the JSON element-type definitions.
 */
#include "element_type_registry.h"

#include <algorithm>

#include "embedded_resources.h"
#include "logging.h"
#include "nabla.h"

namespace pytec3d
{
    namespace
    {
        using json = nlohmann::ordered_json;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("ElementTypeRegistry: " + message);
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

        FieldSpec parse_field(const std::string &type_name,
                              const std::string &name,
                              const json &definition,
                              FieldKind kind,
                              unsigned int dimension,
                              const std::map<std::string, bool> &field_is_vector,
                              const std::vector<std::string> &domains)
        {
            const std::string context = "element type '" + type_name + "', field '" + name + "'";

            const auto vec_it = field_is_vector.find(name);
            if (vec_it == field_is_vector.end())
                fail(context + " is not declared in fields.json.");

            FieldSpec field;
            field.name = name;
            field.kind = kind;
            field.is_vector = vec_it->second;

            // Nodal/element fields are an expression string, or an object that adds metadata.
            if (definition.is_string())
                field.expression = definition.get<std::string>();
            if (definition.is_object())
            {
                if (definition.contains("expression"))
                {
                    if (!definition["expression"].is_string())
                        fail(context + ": 'expression' must be a string.");
                    field.expression = definition["expression"].get<std::string>();
                }
                if (definition.contains("active_components"))
                {
                    for (const auto &label : read_string_list(definition["active_components"], context + " active_components"))
                        field.active_components.push_back(component_index(label, context));
                }
                if (definition.contains("finite_elements"))
                {
                    if (!definition["finite_elements"].is_string())
                        fail(context + ": 'finite_elements' must be a string.");
                    field.finite_elements = definition["finite_elements"].get<std::string>();
                }
                for (const char *key : {"nodal_loads"})
                {
                    if (!definition.contains(key))
                        continue;
                    if (!definition[key].is_string())
                        fail(context + ": '" + key + "' must be a string.");
                    field.nodal_loads = definition[key].get<std::string>();
                }
                if (definition.contains("domain"))
                {
                    field.domains = read_string_list(definition["domain"], context + " domain");
                    for (const auto &domain : field.domains)
                    {
                        if (std::find(domains.begin(), domains.end(), domain) == domains.end())
                            fail(context + ": unknown domain '" + domain + "'.");
                    }
                }
            }

            if (!field.is_vector)
            {
                if (!field.active_components.empty())
                    fail(context + ": scalar field cannot define active_components.");
                return field;
            }

            if (field.active_components.empty())
            {
                for (unsigned int c = 0; c < dimension; ++c)
                    field.active_components.push_back(c);
            }
            std::sort(field.active_components.begin(), field.active_components.end());
            if (std::adjacent_find(field.active_components.begin(), field.active_components.end()) != field.active_components.end())
                fail(context + ": duplicate active_components.");
            return field;
        }

        void parse_section(BaseElementType &type,
                           const json &type_json,
                           const char *section,
                           FieldKind kind,
                           const std::map<std::string, bool> &field_is_vector,
                           const std::vector<std::string> &domains)
        {
            if (!type_json.contains(section))
                return;
            const json &fields = type_json[section];
            if (!fields.is_object())
                fail("element type '" + type.name + "': '" + section + "' must be an object.");
            for (auto it = fields.begin(); it != fields.end(); ++it)
                type.fields.push_back(parse_field(type.name, it.key(), it.value(), kind, type.dimension, field_is_vector, domains));
        }

        const FieldSpec *find_dof_field(const BaseElementType &type, const std::string &name)
        {
            for (const auto &field : type.fields)
                if (field.kind == FieldKind::Dof && field.name == name)
                    return &field;
            return nullptr;
        }

        void check_operand(const BaseElementType &type, const OperatorExpression &operand, const std::string &context)
        {
            if (!find_dof_field(type, operand.variable))
                fail(context + ": '" + operand.variable + "' is not a dof field of the element type.");
        }

        bool has_coefficient(const BaseElementType &type, const std::string &name)
        {
            return std::find(type.coefficients.begin(), type.coefficients.end(), name) != type.coefficients.end();
        }

        void parse_kernels(BaseElementType &type, const json &type_json, const std::vector<std::string> &domains)
        {
            const std::string type_context = "element type '" + type.name + "'";

            if (type_json.contains("time_order"))
            {
                if (!type_json["time_order"].is_number_unsigned())
                    fail(type_context + ": 'time_order' must be an unsigned integer.");
                type.time_order = type_json["time_order"].get<unsigned int>();
            }
            if (type_json.contains("coefficients"))
                type.coefficients = read_string_list(type_json["coefficients"], type_context + " coefficients");

            if (type_json.contains("kernels"))
            {
                if (!type_json["kernels"].is_object())
                    fail(type_context + ": 'kernels' must be an object.");
                for (auto it = type_json["kernels"].begin(); it != type_json["kernels"].end(); ++it)
                {
                    const std::string context = type_context + ", kernel '" + it.key() + "'";
                    const json &definition = it.value();
                    for (const char *key : {"coeff", "test"})
                        if (!definition.is_object() || !definition.contains(key) || !definition[key].is_string())
                            fail(context + " needs a string '" + key + "'.");
                    if (definition.contains("trial") && !definition["trial"].is_string())
                        fail(context + ": 'trial' must be a string.");

                    KernelSpec kernel;
                    kernel.name = it.key();
                    if (definition.contains("trial"))
                        kernel.trial = definition["trial"].get<std::string>();
                    kernel.coeff = definition["coeff"].get<std::string>();
                    kernel.test = definition["test"].get<std::string>();
                    if (definition.contains("domain"))
                    {
                        kernel.domains = read_string_list(definition["domain"], context + " domain");
                        for (const auto &domain : kernel.domains)
                            if (std::find(domains.begin(), domains.end(), domain) == domains.end())
                                fail(context + ": unknown domain '" + domain + "'.");
                    }

                    const std::string coeff_name = !kernel.coeff.empty() && kernel.coeff[0] == '-' ? kernel.coeff.substr(1) : kernel.coeff;
                    if (!has_coefficient(type, coeff_name))
                        fail(context + ": coefficient '" + coeff_name + "' is not listed in 'coefficients'.");
                    if (!kernel.trial.empty())
                        check_operand(type, OperatorExpression::parse(kernel.trial), context + " trial");
                    const OperatorExpression test = OperatorExpression::parse(kernel.test);
                    check_operand(type, test, context + " test");
                    if (test.time_derivative)
                        fail(context + ": dt() is only allowed on the trial side.");
                    type.kernels.push_back(std::move(kernel));
                }
            }

            for (const auto &field : type.fields)
            {
                if (field.kind == FieldKind::Dof)
                    continue;
                const std::string context = type_context + ", field '" + field.name + "'";
                if (field.expression.empty())
                    fail(context + " needs an expression.");
                for (const auto &term : ExpressionTerm::parse_sum(field.expression, type.coefficients))
                    check_operand(type, term.operand, context);
            }
        }
    } // namespace

    ElementTypeRegistry ElementTypeRegistry::from_embedded()
    {
        return from_strings(resources::ELEMENT_TYPES_JSON, resources::FIELDS_JSON, resources::DOMAINS_JSON);
    }

    ElementTypeRegistry ElementTypeRegistry::from_strings(const std::string &element_types_json,
                                                          const std::string &fields_json,
                                                          const std::string &domains_json)
    {
        ElementTypeRegistry registry;

        const json domains = parse_json(domains_json, "domains.json");
        if (!domains.is_object() || !domains.contains("domains"))
            fail("domains.json must contain a 'domains' array.");
        registry.domains_ = read_string_list(domains["domains"], "domains.json 'domains'");

        const json fields = parse_json(fields_json, "fields.json");
        if (!fields.is_object())
            fail("fields.json must be an object.");
        for (auto it = fields.begin(); it != fields.end(); ++it)
        {
            if (!it.value().is_object() || !it.value().contains("vector") || !it.value()["vector"].is_boolean())
                fail("fields.json: field '" + it.key() + "' needs a boolean 'vector' entry.");
            registry.field_is_vector_[it.key()] = it.value()["vector"].get<bool>();
        }

        registry.element_types_json_ = parse_json(element_types_json, "element_types.json");
        if (!registry.element_types_json_.is_object() || !registry.element_types_json_.contains("element_types"))
            fail("element_types.json must contain an 'element_types' object.");

        const json &types = registry.element_types_json_["element_types"];
        for (auto it = types.begin(); it != types.end(); ++it)
        {
            BaseElementType type;
            type.name = it.key();
            const json &type_json = it.value();

            if (!type_json.contains("dimension") || !type_json["dimension"].is_number_unsigned())
                fail("element type '" + type.name + "' needs an unsigned 'dimension'.");
            type.dimension = type_json["dimension"].get<unsigned int>();
            if (type.dimension != 2 && type.dimension != 3)
                fail("element type '" + type.name + "': dimension must be 2 or 3.");

            parse_section(type, type_json, "dof_fields", FieldKind::Dof, registry.field_is_vector_, registry.domains_);
            parse_section(type, type_json, "nodal_fields", FieldKind::Nodal, registry.field_is_vector_, registry.domains_);
            parse_section(type, type_json, "element_fields", FieldKind::Element, registry.field_is_vector_, registry.domains_);
            parse_kernels(type, type_json, registry.domains_);

            registry.types_.emplace(type.name, std::move(type));
        }
        return registry;
    }

    bool ElementTypeRegistry::has_type(const std::string &name) const
    {
        return types_.find(name) != types_.end();
    }

    const BaseElementType &ElementTypeRegistry::get_type(const std::string &name) const
    {
        const auto it = types_.find(name);
        if (it == types_.end())
            fail("unknown element type '" + name + "'.");
        return it->second;
    }

    std::vector<std::string> ElementTypeRegistry::type_names() const
    {
        std::vector<std::string> names;
        names.reserve(types_.size());
        for (const auto &entry : types_)
            names.push_back(entry.first);
        return names;
    }

    unsigned int ElementTypeRegistry::domain_index(const std::string &name) const
    {
        const auto it = std::find(domains_.begin(), domains_.end(), name);
        if (it == domains_.end())
            fail("unknown domain '" + name + "'.");
        return static_cast<unsigned int>(it - domains_.begin());
    }

    const std::string &ElementTypeRegistry::domain_name(unsigned int index) const
    {
        if (index >= domains_.size())
            fail("domain index " + std::to_string(index) + " out of range.");
        return domains_[index];
    }
} // namespace pytec3d
