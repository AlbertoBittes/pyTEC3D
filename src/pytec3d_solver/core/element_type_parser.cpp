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

#include <libmesh/fe_interface.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/string_to_enum.h>

#include "embedded_resources.h"
#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using json = nlohmann::ordered_json;
        using registry::Domain;
        using registry::FieldId;
        using registry::LawId;

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

        std::vector<Domain> read_domains(const json &definition, const std::string &context)
        {
            std::vector<Domain> out;
            if (!definition.is_object() || !definition.contains("domain"))
                return out;
            for (const std::string &name : read_string_list(definition["domain"], context + " domain"))
            {
                const auto domain = registry::parse_Domain(name);
                if (!domain)
                    fail(context + ": unknown domain '" + name + "' (see domains.json).");
                out.push_back(*domain);
            }
            return out;
        }

        // ---- operands ----------------------------------------------------------------------------------------

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

        // What a kernel can name besides the parameters of the registry: the laws of the element type.
        struct ParameterTable
        {
            registry::ElementTypeId type;
            const std::vector<std::optional<NonlinearLaw>> &laws;
        };

        // "-nu", "BH", "turns": an optional sign and the name of a parameter or law; empty is 1.
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

            if (const auto law = registry::parse_LawId(text))
            {
                if (!table.laws[static_cast<std::size_t>(*law)])
                    fail(context + ": the element type does not define the law '" + text + "'.");
                out.kind = ParameterKind::Law;
                out.index = static_cast<unsigned int>(*law);
                out.type = ParameterType::Tensor;
                return out;
            }

            const auto id = registry::parse_ParameterId(text);
            if (!id)
                fail(context + ": unknown parameter '" + text + "'.");
            if (!registry::defined(table.type, *id))
                fail(context + ": the element type does not have the parameter '" + text + "'.");
            const registry::ParameterInfo &info = registry::info(*id);
            if (info.source == ParameterSource::Point)
                fail(context + ": '" + text + "' is point data; it can only be the nodal_loads of a dof field.");
            out.kind = ParameterKind::Parameter;
            out.index = static_cast<unsigned int>(*id);
            out.type = info.type;
            return out;
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
            const std::string name(field.name());
            std::vector<LibMeshVariable> out;

            if (field.vector_output && field.family == libMesh::NEDELEC_ONE)
                out.push_back({name + "z", field.family, field.order});
            else if (field.vector_output)
                for (unsigned int component : field.active_components)
                    out.push_back({name + suffixes[component], scalar_family(field.family), field.order});
            else if (field.physical_vector)
                out.push_back({name + suffixes[field.active_components.at(0)], field.family, field.order});
            else
                out.push_back({name, field.family, field.order});
            return out;
        }

        Field parse_field(const std::string &type_name,
                          const std::string &name,
                          const json &definition,
                          FieldKind kind,
                          unsigned int dimension,
                          const ElementTypeOptions &options)
        {
            const std::string context = "element type '" + type_name + "', field '" + name + "'";

            const auto id = registry::parse_FieldId(name);
            if (!id)
                fail(context + " is not declared in fields.json.");
            if (!definition.is_object())
                fail(context + " must be an object.");

            Field field;
            field.id = *id;
            field.kind = kind;
            field.physical_vector = registry::is_vector(*id);

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
                const auto parameter = registry::parse_ParameterId(key);
                if (!parameter || registry::info(*parameter).source != ParameterSource::Point)
                    fail(context + ": nodal_loads '" + key + "' is not a point-data (mn.) parameter.");
                if (!registry::defined(*registry::parse_ElementTypeId(type_name), *parameter))
                    fail(context + ": the element type does not have the parameter '" + key + "'.");
                field.nodal_loads = registry::info(*parameter).array;
            }
            field.domains = read_domains(definition, context);
            field.libmesh_variables = libmesh_variables(field);
            return field;
        }

        std::vector<Field> parse_fields(const std::string &type_name,
                                        const json &type_json,
                                        const char *section,
                                        FieldKind kind,
                                        unsigned int dimension,
                                        const ElementTypeOptions &options)
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
                Field field = parse_field(type_name, it.key(), it.value(), kind, dimension, options);
                field.index = static_cast<unsigned int>(out.size());
                out.push_back(std::move(field));
            }
            return out;
        }

        // ---- kernels -----------------------------------------------------------------------------------------

        const Field *find_by_name(const std::vector<Field> &fields, const std::string &name)
        {
            const auto id = registry::parse_FieldId(name);
            if (!id)
                return nullptr;
            for (const Field &field : fields)
                if (field.id == *id)
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
                kernel.domains = read_domains(definition, context);

                const ParameterRef &parameter = kernel.parameter;
                if (kernel.has_trial && parameter.type == ParameterType::Vector &&
                    kernel.test.nabla->scalar_valued() == kernel.trial.nabla->scalar_valued())
                    fail(context + ": a vector parameter couples a scalar operand with a vector operand.");
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
                        kernel.domains = read_domains(definition, term_context);
                        if (kernel.parameter.type == ParameterType::Vector && !kernel.trial.nabla->scalar_valued())
                            fail(term_context + ": a vector parameter needs a scalar trial.");
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
                    fail("element type '" + type_name + "': field '" + std::string(field.name()) + "' has no kernel in '" + section + "'.");
            }
            return out;
        }

        // A kernel needs its parameter (or the linear parameter of its law) in every domain where it is active.
        template <typename Kernel>
        void check_parameter_domains(const ElementType &type, const std::vector<Kernel> &kernels)
        {
            for (const Kernel &kernel : kernels)
            {
                const ParameterRef &reference = kernel.parameter;
                if (reference.kind == ParameterKind::One)
                    continue;
                const registry::ParameterId id = reference.kind == ParameterKind::Law ? type.law(reference.law()).linear : reference.parameter();
                for (std::size_t index = 0; index < registry::kNumDomain; ++index)
                {
                    const Domain domain = static_cast<Domain>(index);
                    if (type.active(kernel, domain) && !registry::defined(type.id(), id, domain))
                        fail("element type '" + std::string(type.name()) + "', kernel '" + kernel.name + "' is active in domain " +
                             std::string(registry::name(domain)) + " where the parameter '" + std::string(registry::name(id)) +
                             "' is not defined.");
                }
            }
        }
    } // namespace

    ElementTypeParser ElementTypeParser::from_embedded()
    {
        return from_string(resources::ELEMENT_TYPES_JSON);
    }

    ElementTypeParser ElementTypeParser::from_string(const std::string &element_types_json)
    {
        ElementTypeParser parser;
        parser.element_types_json_ = parse_json(element_types_json, "element_types.json");
        const json &root = parser.element_types_json_;
        if (!root.is_object() || !root.contains("element_types") || !root["element_types"].is_object())
            fail("element_types.json must contain an 'element_types' object.");
        return parser;
    }

    std::unique_ptr<ElementType> ElementTypeParser::build(registry::ElementTypeId id, const ElementTypeOptions &options) const
    {
        const std::string name(registry::name(id));
        const json &type_json = element_types_json_["element_types"][name];
        const std::string context = "element type '" + name + "'";

        const unsigned int dimension = registry::dimension(id);
        if (options.order < static_cast<unsigned int>(libMesh::FIRST))
            fail(context + ": order must be >= 1.");
        if (options.axisymmetric && dimension != 2)
            fail(context + ": axisymmetric requires a 2D element type.");

        std::unique_ptr<ElementType> type(new ElementType());
        type->id_ = id;
        type->dimension_ = dimension;
        type->options_ = options;

        // Nonlinear laws: used in kernels by name, like a parameter.
        type->laws_.assign(registry::kNumLawId, std::nullopt);
        if (type_json.contains("nonlinear_laws"))
        {
            const json &laws = type_json["nonlinear_laws"];
            if (!laws.is_object())
                fail(context + ": 'nonlinear_laws' must be an object.");
            for (auto it = laws.begin(); it != laws.end(); ++it)
            {
                const std::string law_context = context + ", law '" + it.key() + "'";
                const json &definition = it.value();
                if (!definition.is_object())
                    fail(law_context + " must be an object.");
                const auto law_id = registry::parse_LawId(it.key());
                if (!law_id)
                    fail(law_context + " is not in the registry.");

                NonlinearLaw law;
                law.id = *law_id;
                const std::string linear = strip_spaces(read_string(definition, "linear", law_context, true));
                const auto linear_id = registry::parse_ParameterId(linear);
                if (!linear_id || registry::info(*linear_id).source == ParameterSource::Point ||
                    registry::info(*linear_id).type != ParameterType::Tensor)
                    fail(law_context + ": 'linear' must be a tensor parameter.");
                law.linear = *linear_id;
                type->laws_[static_cast<std::size_t>(*law_id)] = law;
            }
        }

        const ParameterTable table{id, type->laws_};

        type->dof_fields_ = parse_fields(name, type_json, "dof_fields", FieldKind::Dof, dimension, options);
        type->nodal_fields_ = parse_fields(name, type_json, "nodal_fields", FieldKind::Nodal, dimension, options);
        type->element_fields_ = parse_fields(name, type_json, "element_fields", FieldKind::Element, dimension, options);

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

        type->dof_kernels_ = parse_dof_kernels(name, type_json, type->dof_fields_, table, nabla_of);
        type->nodal_kernels_ = parse_field_kernels(name, type_json, "nodal_kernels", type->nodal_fields_, type->dof_fields_, table, nabla_of);
        type->element_kernels_ = parse_field_kernels(name, type_json, "element_kernels", type->element_fields_, type->dof_fields_, table, nabla_of);
        check_parameter_domains(*type, type->dof_kernels_);
        check_parameter_domains(*type, type->nodal_kernels_);
        check_parameter_domains(*type, type->element_kernels_);

        type->geometry_fe_ = libMesh::FEBase::build(dimension, libMesh::FEType(libMesh::FIRST, libMesh::LAGRANGE));
        type->geometry_fe_->attach_quadrature_rule(type->quadrature_.get());
        type->geometry_fe_->get_JxW();
        type->geometry_fe_->get_xyz();

        // The FEs to reinit in each domain, each once.
        type->fes_by_domain_.resize(registry::kNumDomain);
        for (std::size_t index = 0; index < registry::kNumDomain; ++index)
        {
            const Domain domain = static_cast<Domain>(index);
            std::vector<libMesh::FEAbstract *> &active = type->fes_by_domain_[index];
            for (const std::vector<Field> *fields : {&type->dof_fields_, &type->nodal_fields_})
                for (const Field &field : *fields)
                    if (field.active_in(domain) && std::find(active.begin(), active.end(), field.fe) == active.end())
                        active.push_back(field.fe);
        }
        return type;
    }
} // namespace pytec3d
