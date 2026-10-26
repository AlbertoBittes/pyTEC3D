/**
 * @file model_data.h
 * @brief The model data JSON: element types, materials, real constants, DOF couplings and nonlinear laws.
 *
 * The caller holds the data as a JSON object (the Python model: the dict model.data) and hands it to the solver
 * with the model (model_input.h); it is read here. Only these sections are read; the keys that name something are checked against the registries
 * (registries.h):
 *
 *   "element_types":  { "<TYPE id>": { "element_type": 0 | "AV2D_nodal",       // index or name in element_types.json
 *                                      "option1": 6 | "STRANDED_COIL",         // domain: index or name in domains.json
 *                                      "option2": 1,                           // FE order (default 1)
 *                                      "option3": 0 } }                        // axisymmetric (default 0)
 *   "materials":      { "<MATERIAL id>": { "mur": [[1,0,0],[0,1,0],[0,0,1]], "rho": [[..]], "M": [0,0,0] } }  // "mat."
 *   "constants":      { "<CONSTANT id>": { "thickness": 1.0, "T": [0, 0, 1], ... } }                         // "rc."
 *   "dof_couplings":  { "<CP id>": ["A", "VOLT"] }
 *   "nonlinear_laws": { "<id>": { "type": "bh_curve", "B": [0, ...], "H": [0, ...] } }   // or an array of such entries
 *
 * A "nonlinear_laws" entry selects the law "type" of nonlinear_laws.json (name or index) on the cells of its id, which is a MATERIAL,
 * CONSTANT or TYPE id as the source of the law says (mat, rc, me). The entry carries the data of the law: the x and
 * y arrays of its curve, if it has one (x strictly increasing). The laws selected on an element are those of its
 * TYPE, MATERIAL and CONSTANT ids that its element type has; they replace the parameters they say (see nonlinear_laws.json).
 * A formulation variant (e.g. time-integrated potentials) is an element type of its own, not an option.
 *
 * A material property or real constant is a number, a 3-vector or a 3x3 tensor as its type says (a number
 * stands for a vector with equal components). Any of them can be left out: it is only an error when a
 * kernel active in the domain of an element needs it (see ElementParameterValues).
 */
#pragma once

#include <nlohmann/json_fwd.hpp>

#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "curve.h"
#include "parameter_value.h"
#include "registries.h"

namespace pytec3d
{
    /** @brief The values of one material or of one set of real constants, each of which may be undefined. */
    template <typename Id, std::size_t N>
    class PropertyRow
    {
    public:
        const ParameterValue *find(Id id) const
        {
            const auto &value = values_[static_cast<std::size_t>(id)];
            return value ? &*value : nullptr;
        }

        void set(Id id, const ParameterValue &value) { values_[static_cast<std::size_t>(id)] = value; }

    private:
        std::array<std::optional<ParameterValue>, N> values_;
    };

    using MaterialRow = PropertyRow<registry::MaterialProperty, registry::kNumMaterialProperty>;
    using ConstantRow = PropertyRow<registry::RealConstant, registry::kNumRealConstant>;

    /** @brief What a TYPE id of the mesh stands for. */
    struct ElementTypeSpec
    {
        registry::ElementTypeId type = registry::ElementTypeId::COUNT;
        registry::Domain domain = registry::Domain::COUNT;
        unsigned int order = 1;
        bool axisymmetric = false;
    };

    struct DofCoupling
    {
        int cp_id = 0;
        std::string variable;
    };

    class ModelData
    {
    public:
        /** @brief Reads the model data; aborts with a message on any inconsistency. */
        static ModelData from_json(const nlohmann::json &root);

        /** @brief from_json() of JSON text. */
        static ModelData parse(const std::string &json_text);

        /** @brief By TYPE id. */
        const std::map<unsigned int, ElementTypeSpec> &element_types() const { return element_types_; }

        /** @brief The material or constants @p id; null when the model does not define it. */
        const MaterialRow *material(unsigned int id) const;
        const ConstantRow *constants(unsigned int id) const;

        /** @brief The curve @p id of the law entry of id @p row_id (MATERIAL, CONSTANT or TYPE id as the source of the
         * law says); null when undefined. */
        const CurveData *curve(registry::CurveId id, unsigned int row_id) const;

        /** @brief The laws selected on the TYPE, MATERIAL and CONSTANT ids @p type_id, @p material_id, @p constants_id. */
        registry::LawMask laws(unsigned int type_id, unsigned int material_id, unsigned int constants_id) const;


        const std::vector<DofCoupling> &dof_couplings() const { return dof_couplings_; }

        /** @brief All the materials and constants rows, by id. */
        const std::map<unsigned int, MaterialRow> &material_rows() const { return materials_; }
        const std::map<unsigned int, ConstantRow> &constant_rows() const { return constants_; }

        /** @brief The laws selected by id, by LawSource (Material, Constant, Cell). */
        const std::array<std::map<unsigned int, registry::LawMask>, 3> &law_selection() const { return laws_; }

        /** @brief One line for the log: the number of rows of each section and the nonlinear laws by id. */
        std::string summary() const;

    private:
        std::map<unsigned int, ElementTypeSpec> element_types_;
        std::map<unsigned int, MaterialRow> materials_;
        std::map<unsigned int, ConstantRow> constants_;
        std::array<std::map<unsigned int, CurveData>, registry::kNumCurveId> curves_; ///< By CurveId, then row id.
        std::vector<DofCoupling> dof_couplings_;
        std::array<std::map<unsigned int, registry::LawMask>, 3> laws_; ///< By LawSource, then id.
    };
} // namespace pytec3d
