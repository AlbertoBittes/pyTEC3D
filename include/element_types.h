#ifndef ELEMENT_TYPES_H
#define ELEMENT_TYPES_H

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <vector>
#include <libmesh/libmesh.h> // Include the libMesh header
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"
#include "libmesh/string_to_enum.h"

using namespace libMesh;

class ElementType;

struct ScalarShapeFunctions
{
    const std::vector<Real> &JxW;
    const std::vector<std::vector<Real>> &phi;
    const std::vector<std::vector<RealGradient>> &dphi;

    // Constructor
    ScalarShapeFunctions(const std::vector<Real> &JxW,
                         const std::vector<std::vector<Real>> &phi,
                         const std::vector<std::vector<RealGradient>> &dphi)
        : JxW(JxW), phi(phi), dphi(dphi) {}
};

struct VectorShapeFunctions
{
    const std::vector<Real> &JxW;
    const std::vector<std::vector<RealGradient>> &phi_vec;
    const std::vector<std::vector<RealTensor>> &dphi_vec;
    const std::vector<std::vector<RealGradient>> &curl_phi_vec;

    // Constructor
    VectorShapeFunctions(const std::vector<Real> &JxW,
                         const std::vector<std::vector<RealGradient>> &phi_vec,
                         const std::vector<std::vector<RealTensor>> &dphi_vec,
                         const std::vector<std::vector<RealGradient>> &curl_phi_vec)
        : JxW(JxW), phi_vec(phi_vec), dphi_vec(dphi_vec), curl_phi_vec(curl_phi_vec) {}
};

typedef std::unordered_map<unsigned int, std::vector<unsigned int>> ElementTypesMapCSV;
typedef std::unordered_map<unsigned int, std::shared_ptr<ElementType>> ElementTypesMap;
typedef std::unordered_map<std::string, std::set<subdomain_id_type>> VariableToElementIdsMapType;
//                         label         is_active, physical components, is scalar, API label, API components, fe family, order
typedef std::map<std::string, std::tuple<bool, unsigned int, bool, std::string, unsigned int, std::string, std::string>> VariableData;
typedef std::unordered_map<std::string, VariableData> ElementTypeData;

class ElementType
{
public:
    // Constructor
    ElementType(const std::vector<unsigned int> &keyopts)
        : keyopts(keyopts)
    {
        element_type_data["DOFVariables"] = VariableData();
        element_type_data["NodalDerivedQuantitiesVariables"] = VariableData();
        element_type_data["ElementDerivedQuantitiesVariables"] = VariableData();
        element_type_data["LoadsVariables"] = VariableData();
    }

    // Virtual methods to be implemented by derived classes
    std::vector<std::string> getActiveVariables(const std::string &var_type) const
    {
        std::vector<std::string> active_vars;
        auto it = element_type_data.find(var_type);
        if (it != element_type_data.end())
        {
            for (const auto &var : it->second)
            {
                if (std::get<0>(var.second)) // Check if the variable is active
                {
                    active_vars.push_back(var.first);
                }
            }
        }
        return active_vars;
    }

    bool is_conductive() const
    {
        // Check if "VOLT" exists in element_type_data["DOFVariables"]
        auto it = element_type_data.find("DOFVariables");
        if (it != element_type_data.end())
        {
            const auto &dof_variables = it->second;
            return dof_variables.find("VOLT") != dof_variables.end();
        }
        return false; // "DOFVariables" or "VOLT" not found
    }

    // Getters
    const std::vector<unsigned int> &get_properties() const { return keyopts; }
    // Getter for ElementTypeData
    const ElementTypeData &get_element_type_data() const
    {
        return element_type_data;
    }

    const FEBase *get_fe() const { return fe.get(); }
    const FEVectorBase *get_fe_vec() const { return fe_vec.get(); }

    const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_matrix_idx_size_map() const { return matrix_idx_size_map; }
    const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_derived_matrix_idx_size_map() const { return derived_matrix_idx_size_map; }
    const std::map<std::string, std::tuple<unsigned int, unsigned int>> &get_vector_idx_size_map() const { return vector_idx_size_map; }

    // Setters
    void set_fe(std::unique_ptr<FEBase> new_fe) { fe = std::move(new_fe); }
    // Setter for ElementTypeData
    void set_element_type_data(const ElementTypeData &data)
    {
        element_type_data = data;
    }
    // Virtual methods for populating maps
    virtual void populate_matrix_idx_size_map() = 0;
    virtual void populate_derived_matrix_idx_size_map() = 0;
    virtual void populate_vector_idx_size_map() = 0;

    // Populate all maps
    void populate_all();

    void reinit_fe(Elem *elem)
    {
        fe->reinit(elem);
    }
    void reinit_fe_vec(Elem *elem)
    {
        fe_vec->reinit(elem);
    }
    // Shape functions
    inline std::shared_ptr<ScalarShapeFunctions> get_scalar_shape_functions() const
    {
        return std::make_shared<ScalarShapeFunctions>(fe->get_JxW(), fe->get_phi(), fe->get_dphi());
    }

    inline std::shared_ptr<VectorShapeFunctions> get_vector_shape_functions() const
    {
        return std::make_shared<VectorShapeFunctions>(fe_vec->get_JxW(), fe_vec->get_phi(), fe_vec->get_dphi(), fe_vec->get_curl_phi());
    }

protected:
    std::unique_ptr<FEBase> fe;           // Finite element object
    std::unique_ptr<FEVectorBase> fe_vec; // Finite element vector object
    std::vector<unsigned int> keyopts;    // Element type keyopts
    ElementTypeData element_type_data;    // Element type data
    // Maps for matrix and vector sizes
    std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> matrix_idx_size_map;
    std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> derived_matrix_idx_size_map;
    std::map<std::string, std::tuple<unsigned int, unsigned int>> vector_idx_size_map;
};

class AVElementType : public ElementType
{
public:
    // Constructor
    AVElementType(const std::vector<unsigned int> &keyopts)
        : ElementType(keyopts) // Call the parent class constructor
    {
        // Additional initialization for AVElementType (if needed)
        element_type_data["NodalDerivedQuantitiesVariables"]["B_x"] =
            std::make_tuple(true, 3, false, "B", 0, "L2_LAGRANGE", "FIRST");
        element_type_data["NodalDerivedQuantitiesVariables"]["B_y"] =
            std::make_tuple(true, 3, false, "B", 1, "L2_LAGRANGE", "FIRST");
        element_type_data["ElementDerivedQuantitiesVariables"]["J_z"] =
            std::make_tuple(true, 3, false, "J", 2, "MONOMIAL", "CONSTANT");
        element_type_data["LoadsVariables"]["Js"] =
            std::make_tuple(true, 3, false, "Js", 3, "", "");
    }
};

class AV3DElementType : public AVElementType
{
public:
    // Constructor
    AV3DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
    {
        // Additional initialization for AVElementType (if needed)
        element_type_data["NodalDerivedQuantitiesVariables"]["B_z"] =
            std::make_tuple(true, 3, false, "B", 2, "L2_LAGRANGE", "FIRST");
        element_type_data["ElementDerivedQuantitiesVariables"]["J_x"] =
            std::make_tuple(true, 3, false, "J", 0, "MONOMIAL", "CONSTANT");
        element_type_data["ElementDerivedQuantitiesVariables"]["J_y"] =
            std::make_tuple(true, 3, false, "J", 1, "MONOMIAL", "CONSTANT");
        element_type_data["DOFVariables"]["VOLT"] =
            std::make_tuple(false, 1, true, "VOLT", 0, "LAGRANGE", "FIRST");
        element_type_data["DOFVariables"]["EMF"] =
            std::make_tuple(false, 1, true, "EMF", 0, "LAGRANGE", "FIRST");

        FEType fe_type(FIRST, LAGRANGE);
        // Build a Finite Element object of the specified type.  Since the
        // FEBase::build() member dynamically creates memory we will
        // store the object as a std::unique_ptr<FEBase>.  This can be thought
        // of as a pointer that will clean up after itself.
        fe = FEBase::build(3, fe_type);
        // A 2nd order Gauss quadrature rule for numerical integration.
        QGauss *qrule = new QGauss(3, fe_type.default_quadrature_order());
        // Tell the finite element object to use our quadrature rule.
        fe->attach_quadrature_rule(qrule);

        auto sf = this->get_scalar_shape_functions();
    }
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &keyopts) : AVElementType(keyopts) // Call the parent class constructor
    {
        element_type_data["DOFVariables"]["VOLT"] =
            std::make_tuple(false, 1, true, "VOLT", 0, "LAGRANGE", "FIRST");
        // Additional initialization for AVElementType (if needed)
        FEType fe_type(FIRST, LAGRANGE);
        // Build a Finite Element object of the specified type.  Since the
        // FEBase::build() member dynamically creates memory we will
        // store the object as a std::unique_ptr<FEBase>.  This can be thought
        // of as a pointer that will clean up after itself.
        fe = FEBase::build(2, fe_type);
        // A 2nd order Gauss quadrature rule for numerical integration.
        QGauss *qrule = new QGauss(2, fe_type.default_quadrature_order());
        // Tell the finite element object to use our quadrature rule.
        fe->attach_quadrature_rule(qrule);

        auto sf = this->get_scalar_shape_functions();
    }
};

class SOLID97Type : public AV3DElementType
{
public:
    SOLID97Type(const std::vector<unsigned int> &keyopts);

    void populate_matrix_idx_size_map() override;
    void populate_derived_matrix_idx_size_map() override;
    void populate_vector_idx_size_map() override;
};

class SOLID236Type : public AV3DElementType
{
public:
    SOLID236Type(const std::vector<unsigned int> &keyopts);

    void populate_matrix_idx_size_map() override;
    void populate_derived_matrix_idx_size_map() override;
    void populate_vector_idx_size_map() override;
};

class PLANE13Type : public AV2DElementType
{
public:
    PLANE13Type(const std::vector<unsigned int> &keyopts);

    void populate_matrix_idx_size_map() override;
    void populate_derived_matrix_idx_size_map() override;
    void populate_vector_idx_size_map() override;
};

class PLANE233Type : public AV2DElementType
{
public:
    PLANE233Type(const std::vector<unsigned int> &keyopts);

    void populate_matrix_idx_size_map() override;
    void populate_derived_matrix_idx_size_map() override;
    void populate_vector_idx_size_map() override;
};

class ElementTypesCollector
{
public:
    ElementTypesCollector() {}

    ElementTypesCollector(ElementTypesMapCSV &elementTypesMapCSV)
    {
        collectElementTypes(elementTypesMapCSV);
    }

    ElementTypesMap &getElementTypes()
    {
        return elementTypes;
    }

    ElementTypesMapCSV *readElementTypes(const std::string &filename);

    VariableToElementIdsMapType *associateVariablesWithElementIds(std::string variables_set) const;

    void collectElementTypes(ElementTypesMapCSV &elementTypesMapCSV);

    ElementTypesMap elementTypes;

    std::set<std::string> getLoadsLabels() const;
};
#endif // ELEMENT_TYPES_H