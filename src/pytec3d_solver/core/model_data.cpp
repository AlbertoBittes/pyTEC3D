/**
 * @file model_data.cpp
 * @brief ModelData: reading of the model data JSON.
 */
#include "model_data.h"

#include <cmath>
#include <set>

#include <nlohmann/json.hpp>

#include "logging.h"

namespace pytec3d
{
    namespace
    {
        using json = nlohmann::json;

        [[noreturn]] void fail(const std::string &message)
        {
            logging::error_and_abort("Model data: " + message);
        }

        const json *section(const json &root, const char *name)
        {
            const auto it = root.find(name);
            if (it == root.end() || it->is_null())
                return nullptr;
            if (!it->is_object())
                fail("'" + std::string(name) + "' must be an object.");
            return &*it;
        }

        unsigned int parse_id(const std::string &key, const std::string &context)
        {
            std::size_t used = 0;
            long value = -1;
            try
            {
                value = std::stol(key, &used);
            }
            catch (const std::exception &)
            {
                used = 0;
            }
            if (used != key.size() || value <= 0)
                fail(context + ": the id '" + key + "' must be a positive integer.");
            return static_cast<unsigned int>(value);
        }

        double finite_number(const json &value, const std::string &context)
        {
            if (!value.is_number())
                fail(context + " must be a number.");
            const double number = value.get<double>();
            if (!std::isfinite(number))
                fail(context + " must be finite.");
            return number;
        }

        ParameterValue read_value(const json &value, ParameterType type, const std::string &context)
        {
            ParameterValue out;
            out.type = type;
            switch (type)
            {
            case ParameterType::Scalar:
                out.scalar = finite_number(value, context);
                return out;
            case ParameterType::Vector:
                if (value.is_number())
                {
                    const double number = finite_number(value, context);
                    for (unsigned int c = 0; c < 3; ++c)
                        out.vector(c) = number;
                    return out;
                }
                if (!value.is_array() || value.size() != 3)
                    fail(context + " must be a number or an array of 3 numbers.");
                for (unsigned int c = 0; c < 3; ++c)
                    out.vector(c) = finite_number(value[c], context + "[" + std::to_string(c) + "]");
                return out;
            case ParameterType::Tensor:
                if (!value.is_array() || value.size() != 3)
                    fail(context + " must be a 3x3 array.");
                for (unsigned int i = 0; i < 3; ++i)
                {
                    if (!value[i].is_array() || value[i].size() != 3)
                        fail(context + " must be a 3x3 array.");
                    for (unsigned int j = 0; j < 3; ++j)
                        out.tensor(i, j) = finite_number(value[i][j], context + "[" + std::to_string(i) + "][" + std::to_string(j) + "]");
                }
                return out;
            }
            fail(context + ": unknown value type.");
        }

        // A name of a registry, or its index.
        template <typename Id, typename Parse>
        Id read_registry_entry(const json &value, std::size_t count, Parse parse, const std::string &context, const char *what)
        {
            if (value.is_string())
            {
                const auto id = parse(value.get<std::string>());
                if (!id)
                    fail(context + ": '" + value.get<std::string>() + "' is not a " + what + ".");
                return *id;
            }
            if (value.is_number_integer() && value.get<long>() >= 0 && static_cast<std::size_t>(value.get<long>()) < count)
                return static_cast<Id>(value.get<long>());
            fail(context + ": expected a " + what + " name or an index below " + std::to_string(count) + ".");
        }

        bool read_flag(const json &value, const std::string &context)
        {
            if (value.is_boolean())
                return value.get<bool>();
            if (value.is_number_integer() && (value.get<int>() == 0 || value.get<int>() == 1))
                return value.get<int>() == 1;
            fail(context + " must be 0 or 1.");
        }

        template <typename Row, typename Parse, typename Type>
        void read_property_rows(const json &root,
                                const char *name,
                                const char *what,
                                Parse parse,
                                Type type_of,
                                std::map<unsigned int, Row> &rows)
        {
            const json *rows_json = section(root, name);
            if (!rows_json)
                return;
            for (auto it = rows_json->begin(); it != rows_json->end(); ++it)
            {
                const std::string context = std::string(name) + "[" + it.key() + "]";
                const unsigned int id = parse_id(it.key(), name);
                if (!it.value().is_object())
                    fail(context + " must be an object of " + what + ".");
                Row row;
                for (auto entry = it.value().begin(); entry != it.value().end(); ++entry)
                {
                    const auto property = parse(entry.key());
                    if (!property)
                        fail(context + ": '" + entry.key() + "' is not a " + what + " (see the \"parameters\" of element_types.json).");
                    row.set(*property, read_value(entry.value(), type_of(*property), context + "." + entry.key()));
                }
                rows.emplace(id, std::move(row));
            }
        }
    } // namespace

    ModelData ModelData::parse(const std::string &json_text)
    {
        json root;
        try
        {
            root = json::parse(json_text);
        }
        catch (const json::exception &ex)
        {
            fail(std::string("cannot parse the JSON: ") + ex.what());
        }
        if (!root.is_object())
            fail("the JSON must be an object.");

        ModelData data;

        if (const json *types = section(root, "element_types"))
        {
            static const std::set<std::string> known = {"element_type", "option1", "option2", "option3", "option4", "option5"};
            for (auto it = types->begin(); it != types->end(); ++it)
            {
                const std::string context = "element_types[" + it.key() + "]";
                const unsigned int id = parse_id(it.key(), "element_types");
                const json &definition = it.value();
                if (!definition.is_object())
                    fail(context + " must be an object.");
                for (auto entry = definition.begin(); entry != definition.end(); ++entry)
                    if (!known.count(entry.key()))
                        fail(context + ": unknown key '" + entry.key() + "'.");
                for (const char *required : {"element_type", "option1"})
                    if (!definition.contains(required))
                        fail(context + " needs '" + required + "'.");

                ElementTypeSpec spec;
                spec.type = read_registry_entry<registry::ElementTypeId>(
                    definition["element_type"], registry::kNumElementTypeId, [](const std::string &n)
                    { return registry::parse_ElementTypeId(n); },
                    context + ".element_type", "element type");
                spec.domain = read_registry_entry<registry::Domain>(
                    definition["option1"], registry::kNumDomain, [](const std::string &n)
                    { return registry::parse_Domain(n); },
                    context + ".option1", "domain");
                if (definition.contains("option2"))
                {
                    if (!definition["option2"].is_number_integer() || definition["option2"].get<int>() < 1)
                        fail(context + ".option2 (FE order) must be an integer >= 1.");
                    spec.order = definition["option2"].get<unsigned int>();
                }
                if (definition.contains("option3"))
                    spec.axisymmetric = read_flag(definition["option3"], context + ".option3 (axisymmetric)");
                data.element_types_.emplace(id, spec);
            }
        }

        read_property_rows(
            root, "materials", "material property",
            [](const std::string &n)
            { return registry::parse_MaterialProperty(n); },
            [](registry::MaterialProperty p)
            { return registry::type(p); }, data.materials_);
        read_property_rows(
            root, "constants", "real constant",
            [](const std::string &n)
            { return registry::parse_RealConstant(n); },
            [](registry::RealConstant p)
            { return registry::type(p); }, data.constants_);

        if (const json *couplings = section(root, "dof_couplings"))
            for (auto it = couplings->begin(); it != couplings->end(); ++it)
            {
                const std::string context = "dof_couplings[" + it.key() + "]";
                const int cp_id = static_cast<int>(parse_id(it.key(), "dof_couplings"));
                if (!it.value().is_array())
                    fail(context + " must be an array of variable names.");
                for (const json &variable : it.value())
                {
                    if (!variable.is_string())
                        fail(context + " must be an array of variable names.");
                    data.dof_couplings_.push_back({cp_id, variable.get<std::string>()});
                }
            }

        if (const json *curves = section(root, "bh_curves"))
            for (auto it = curves->begin(); it != curves->end(); ++it)
            {
                const std::string context = "bh_curves[" + it.key() + "]";
                const unsigned int id = parse_id(it.key(), "bh_curves");
                const json &curve = it.value();
                if (!curve.is_object() || !curve.contains("H") || !curve.contains("B") || !curve["H"].is_array() || !curve["B"].is_array())
                    fail(context + " must be an object with the arrays 'H' and 'B'.");
                if (curve["H"].size() != curve["B"].size() || curve["H"].size() < 2)
                    fail(context + ": 'H' and 'B' must have the same length, at least 2.");

                std::vector<libMesh::Real> h_values, b_values;
                for (std::size_t i = 0; i < curve["H"].size(); ++i)
                {
                    h_values.push_back(finite_number(curve["H"][i], context + ".H[" + std::to_string(i) + "]"));
                    b_values.push_back(finite_number(curve["B"][i], context + ".B[" + std::to_string(i) + "]"));
                }
                data.bh_curves_.emplace(id, BHCurve(std::move(h_values), std::move(b_values)));
            }

        return data;
    }

    const MaterialRow *ModelData::material(unsigned int id) const
    {
        const auto it = materials_.find(id);
        return it == materials_.end() ? nullptr : &it->second;
    }

    const ConstantRow *ModelData::constants(unsigned int id) const
    {
        const auto it = constants_.find(id);
        return it == constants_.end() ? nullptr : &it->second;
    }
} // namespace pytec3d
