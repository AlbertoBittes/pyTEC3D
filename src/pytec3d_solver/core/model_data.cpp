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
            ParameterValue out = ParameterValue::zero(type);
            switch (type)
            {
            case ParameterType::Scalar:
                out.scalar() = finite_number(value, context);
                return out;
            case ParameterType::Vector:
                if (value.is_number())
                {
                    const double number = finite_number(value, context);
                    for (unsigned int c = 0; c < 3; ++c)
                        out.vector()(c) = number;
                    return out;
                }
                if (!value.is_array() || value.size() != 3)
                    fail(context + " must be a number or an array of 3 numbers.");
                for (unsigned int c = 0; c < 3; ++c)
                    out.vector()(c) = finite_number(value[c], context + "[" + std::to_string(c) + "]");
                return out;
            case ParameterType::Tensor:
                if (!value.is_array() || value.size() != 3)
                    fail(context + " must be a 3x3 array.");
                for (unsigned int i = 0; i < 3; ++i)
                {
                    if (!value[i].is_array() || value[i].size() != 3)
                        fail(context + " must be a 3x3 array.");
                    for (unsigned int j = 0; j < 3; ++j)
                        out.tensor()(i, j) = finite_number(value[i][j], context + "[" + std::to_string(i) + "][" + std::to_string(j) + "]");
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

        [[noreturn]] void moved_law_key(const std::string &context)
        {
            fail(context + ": the nonlinear laws are given in the \"nonlinear_laws\" section, keyed by the same id: "
                           "\"nonlinear_laws\": {\"<id>\": {\"type\": \"bh_curve\", \"B\": [...], \"H\": [...]}}.");
        }

        std::vector<libMesh::Real> read_curve_array(const json &entry, const std::string &key, const std::string &context)
        {
            if (!entry.contains(key) || !entry[key].is_array())
                fail(context + " needs the array '" + key + "'.");
            std::vector<libMesh::Real> values;
            for (std::size_t i = 0; i < entry[key].size(); ++i)
                values.push_back(finite_number(entry[key][i], context + "." + key + "[" + std::to_string(i) + "]"));
            return values;
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
                    if (entry.key() == "law")
                        moved_law_key(context);
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
            static const std::set<std::string> known = {"element_type", "option1", "option2", "option3"};
            for (auto it = types->begin(); it != types->end(); ++it)
            {
                const std::string context = "element_types[" + it.key() + "]";
                const unsigned int id = parse_id(it.key(), "element_types");
                const json &definition = it.value();
                if (!definition.is_object())
                    fail(context + " must be an object.");
                for (auto entry = definition.begin(); entry != definition.end(); ++entry)
                {
                    if (entry.key() == "option4" || entry.key() == "law")
                        moved_law_key(context);
                    if (!known.count(entry.key()))
                        fail(context + ": unknown key '" + entry.key() + "'.");
                }
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

        if (root.contains("bh_curves"))
            fail("bh_curves: the curves are given with their law in the \"nonlinear_laws\" section: "
                 "\"nonlinear_laws\": {\"<MATERIAL id>\": {\"type\": \"bh_curve\", \"B\": [...], \"H\": [...]}}.");

        if (const json *laws = section(root, "nonlinear_laws"))
            for (auto it = laws->begin(); it != laws->end(); ++it)
            {
                const unsigned int id = parse_id(it.key(), "nonlinear_laws");
                const bool many = it.value().is_array();
                const json entries = many ? it.value() : json::array({it.value()});
                for (std::size_t n = 0; n < entries.size(); ++n)
                {
                    const std::string context = "nonlinear_laws[" + it.key() + "]" + (many ? "[" + std::to_string(n) + "]" : "");
                    const json &entry = entries[n];
                    if (!entry.is_object() || !entry.contains("type"))
                        fail(context + " must be an object with the \"type\" of the law (a law of nonlinear_laws.json).");
                    const auto law = read_registry_entry<registry::LawId>(
                        entry["type"], registry::kNumLawId, [](const std::string &text)
                        { return registry::parse_LawId(text); },
                        context + ".type", "law of nonlinear_laws.json");
                    const registry::LawInfo &info = registry::info(law);
                    registry::LawMask &selected = data.laws_[static_cast<std::size_t>(info.source)][id];
                    if (selected & registry::bit(law))
                        fail(context + ": the law '" + std::string(registry::name(law)) + "' is given twice for this id.");
                    selected |= registry::bit(law);

                    std::set<std::string> keys = {"type"};
                    if (info.curve != registry::CurveId::COUNT)
                    {
                        const registry::CurveSource &curve = registry::kCurveSources[static_cast<std::size_t>(info.curve)];
                        const std::string x_key(curve.x), y_key(curve.y);
                        keys.insert({x_key, y_key});
                        auto x = read_curve_array(entry, x_key, context);
                        auto y = read_curve_array(entry, y_key, context);
                        if (x.size() != y.size() || x.size() < 2)
                            fail(context + ": '" + x_key + "' and '" + y_key + "' must have the same length, at least 2.");
                        for (std::size_t i = 1; i < x.size(); ++i)
                            if (x[i] <= x[i - 1])
                                fail(context + ": '" + x_key + "' must be strictly increasing.");
                        data.curves_[static_cast<std::size_t>(info.curve)].emplace(id, CurveData(std::move(x), std::move(y)));
                    }
                    for (auto key = entry.begin(); key != entry.end(); ++key)
                        if (!keys.count(key.key()))
                            fail(context + ": unknown key '" + key.key() + "' for the law '" +
                                 std::string(registry::name(law)) + "'.");
                }
            }

        return data;
    }

    registry::LawMask ModelData::laws(unsigned int type_id, unsigned int material_id, unsigned int constants_id) const
    {
        registry::LawMask out = 0;
        for (const auto &[source, id] : {std::pair{registry::LawSource::Cell, type_id},
                                         std::pair{registry::LawSource::Material, material_id},
                                         std::pair{registry::LawSource::Constant, constants_id}})
        {
            const auto &selected = laws_[static_cast<std::size_t>(source)];
            if (const auto it = selected.find(id); it != selected.end())
                out |= it->second;
        }
        return out;
    }

    std::string ModelData::summary() const
    {
        std::set<int> coupling_ids;
        for (const DofCoupling &coupling : dof_couplings_)
            coupling_ids.insert(coupling.cp_id);
        std::string out = std::to_string(element_types_.size()) + " element type row(s), " +
                          std::to_string(materials_.size()) + " material(s), " + std::to_string(constants_.size()) +
                          " constants row(s), " + std::to_string(coupling_ids.size()) + " DOF coupling id(s)";
        static const char *id_names[] = {"MATERIAL", "CONSTANT", "TYPE"};
        std::string laws;
        for (std::size_t source = 0; source < laws_.size(); ++source)
            for (const auto &[id, mask] : laws_[source])
                for (std::size_t l = 0; l < registry::kNumLawId; ++l)
                    if (mask & registry::bit(static_cast<registry::LawId>(l)))
                        laws += (laws.empty() ? "" : ", ") + std::string(registry::name(static_cast<registry::LawId>(l))) +
                                " on " + id_names[source] + " " + std::to_string(id);
        return out + (laws.empty() ? ", no nonlinear law" : ", nonlinear laws: " + laws);
    }

    const MaterialRow *ModelData::material(unsigned int id) const
    {
        const auto it = materials_.find(id);
        return it == materials_.end() ? nullptr : &it->second;
    }

    const CurveData *ModelData::curve(registry::CurveId id, unsigned int row_id) const
    {
        const auto &curves = curves_[static_cast<std::size_t>(id)];
        const auto it = curves.find(row_id);
        return it == curves.end() ? nullptr : &it->second;
    }

    const ConstantRow *ModelData::constants(unsigned int id) const
    {
        const auto it = constants_.find(id);
        return it == constants_.end() ? nullptr : &it->second;
    }
} // namespace pytec3d
