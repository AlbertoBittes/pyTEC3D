/**
 * @file model_data.h
 * @brief The model data JSON: element types, materials, real constants, DOF couplings and BH curves.
 *
 * The Python model keeps the data as a JSON-like dict (model.data); the solver receives it as text and reads
 * it here. Only these sections are read; the keys that name something are checked against the registries
 * (registries.h):
 *
 *   "element_types":  { "<TYPE id>": { "element_type": 0 | "AV2D_nodal",       // index or name in element_types.json
 *                                      "option1": 6 | "STRANDED_COIL",         // domain: index or name in domains.json
 *                                      "option2": 1,                           // FE order (default 1)
 *                                      "option3": 0,                           // axisymmetric (default 0)
 *                                      "law": "<me. law>" } }                  // nonlinear laws of source me (optional)
 *   "materials":      { "<MATERIAL id>": { "mur": [[1,0,0],[0,1,0],[0,0,1]], "rho": [[..]], "M": [0,0,0],   // "mat." parameters
 *                                          "law": "ferromagnetic" } }          // nonlinear laws of source mat (optional)
 *   "constants":      { "<CONSTANT id>": { "thickness": 1.0, "T": [0, 0, 1], ..., "law": ... } }   // "rc." parameters
 *   "dof_couplings":  { "<CP id>": ["A", "VOLT"] }
 *   "bh_curves":      { "<MATERIAL id>": { "B": [0, ...], "H": [0, ...] } }   // the curves of the laws (CurveId)
 *
 * "law" is a law of laws.json (name or index) or an array of them; each must have the source of the row. The
 * laws selected on an element are those of its TYPE, MATERIAL and CONSTANT rows that its element type has; they
 * replace the parameters they say (see laws.json). A curve section, and the names of its x and y arrays, are
 * those the laws declare; it is keyed by the ids of the rows of the laws' source, and every row selecting a law
 * with a curve needs one. A formulation variant (e.g. time-integrated potentials) is an element type of its
 * own, not an option.
 *
 * A material property or real constant is a number, a 3-vector or a 3x3 tensor as its type says (a number
 * stands for a vector with equal components). Any of them can be left out: it is only an error when a
 * kernel active in the domain of an element needs it (see ElementParameterValues).
 */
#pragma once

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
        registry::LawMask laws = 0; ///< Selected by the "law" key.

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
        registry::LawMask laws = 0; ///< Laws of source me selected by the "law" key.
    };

    struct DofCoupling
    {
        int cp_id = 0;
        std::string variable;
    };

    class ModelData
    {
    public:
        /** @brief Reads the JSON text of model.data; aborts with a message on any inconsistency. */
        static ModelData parse(const std::string &json_text);

        /** @brief By TYPE id. */
        const std::map<unsigned int, ElementTypeSpec> &element_types() const { return element_types_; }

        /** @brief The material or constants @p id; null when the model does not define it. */
        const MaterialRow *material(unsigned int id) const;
        const ConstantRow *constants(unsigned int id) const;

        /** @brief The curve @p id of row @p row_id (of the source of the laws reading it); null when undefined. */
        const CurveData *curve(registry::CurveId id, unsigned int row_id) const;

        /** @brief The laws selected by the TYPE, MATERIAL and CONSTANT rows @p type_id, @p material_id, @p constants_id. */
        registry::LawMask laws(unsigned int type_id, unsigned int material_id, unsigned int constants_id) const;


        const std::vector<DofCoupling> &dof_couplings() const { return dof_couplings_; }

    private:
        std::map<unsigned int, ElementTypeSpec> element_types_;
        std::map<unsigned int, MaterialRow> materials_;
        std::map<unsigned int, ConstantRow> constants_;
        std::array<std::map<unsigned int, CurveData>, registry::kNumCurveId> curves_; ///< By CurveId, then row id.
        std::vector<DofCoupling> dof_couplings_;
    };
} // namespace pytec3d
