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

struct ElementTypeProperties
{
    std::string element_type;
    std::vector<int> keyopts;
};

struct ScalarShapeFunctions
{
    const std::vector<Real> JxW;
    const std::vector<std::vector<Real>> phi;
    const std::vector<std::vector<RealGradient>> dphi;
};

struct VectorShapeFunctions
{
    const std::vector<Real> &JxW;
    const std::vector<std::vector<RealGradient>> &phi_vec;
    const std::vector<std::vector<RealTensor>> &dphi_vec;
    const std::vector<std::vector<RealGradient>> &curl_phi_vec;
};

typedef std::unordered_map<unsigned int, ElementTypeProperties> ElementTypesMapCSV;
typedef std::unordered_map<unsigned int, std::shared_ptr<ElementType>> ElementTypesMap;
typedef std::unordered_map<std::string, std::set<subdomain_id_type>> VariableToElementIdsMapType;

class ElementType
{
public:
    ElementType(const ElementTypeProperties &properties, const std::string &order = "default_order", const std::string &fe_family = "default_family", unsigned int quadrature_order = 0, unsigned int dimension = 0)
        : properties(properties), order(order), fe_family(fe_family), quadrature_order(quadrature_order), dimension(dimension), az_dofs(0) {}

    virtual std::vector<std::string> getActiveVariables(std::string variables_set) const = 0;

    virtual std::vector<std::string> getVariables(std::string variables_set) const = 0;

    virtual std::vector<std::string> getElementLoads() const = 0;

    virtual bool is_conductive() const = 0;

    std::string getVariableFEFamily(std::string variable);

    std::string getVariableOrder(std::string variable);

    inline ScalarShapeFunctions get_scalar_shape_functions() const
    {
        // initialize the ShapeFunctions struct
        // and fill it with the shape functions from the finite element object
        // and the finite element vector object
        return ScalarShapeFunctions{fe->get_JxW(), fe->get_phi(), fe->get_dphi()};
    }

    inline VectorShapeFunctions get_vector_shape_functions() const
    {
        // initialize the ShapeFunctions struct
        // and fill it with the shape functions from the finite element object
        // and the finite element vector object
        return VectorShapeFunctions{fe_vec->get_JxW(), fe_vec->get_phi(), fe_vec->get_dphi(), fe_vec->get_curl_phi()};
    }

    virtual void populate_matrix_idx_size_map() = 0;
    virtual void populate_derived_matrix_idx_size_map() = 0;
    virtual void populate_vector_idx_size_map() = 0;

    // think if this can be don with n_dofs
    void populate_all();

    ElementTypeProperties properties;
    std::string order;
    std::string fe_family;
    unsigned int quadrature_order;
    unsigned int dimension;
    unsigned int az_dofs;
    std::unique_ptr<FEBase> fe;
    std::unique_ptr<FEVectorBase> fe_vec;
    std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> matrix_idx_size_map;
    std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> derived_matrix_idx_size_map;
    std::map<std::string, std::tuple<unsigned int, unsigned int>> vector_idx_size_map;
};

class SOLID97Type : public ElementType
{
public:
    SOLID97Type(const ElementTypeProperties &properties);

    std::vector<std::string> getElementLoads() const override;

    void populate_matrix_idx_size_map() override;

    void populate_derived_matrix_idx_size_map() override;

    void populate_vector_idx_size_map() override;

    std::vector<std::string> getActiveVariables(std::string variables_set) const override;

    std::vector<std::string> getVariables(std::string variables_set) const override;

    inline bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
    }
};

class SOLID236Type : public ElementType
{
public:
    SOLID236Type(const ElementTypeProperties &properties);

    std::vector<std::string> getElementLoads() const override;

    std::vector<std::string> getActiveVariables(std::string variables_set) const override;

    std::vector<std::string> getVariables(std::string variables_set) const override;

    void populate_matrix_idx_size_map() override;
    void populate_derived_matrix_idx_size_map() override;
    void populate_vector_idx_size_map() override;

    inline bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
    }
};

class PLANE13Type : public ElementType
{
public:
    PLANE13Type(const ElementTypeProperties &properties);

    std::vector<std::string> getElementLoads() const override;

    void populate_matrix_idx_size_map() override;

    void populate_derived_matrix_idx_size_map() override;

    void populate_vector_idx_size_map() override;

    std::vector<std::string> getActiveVariables(std::string variables_set) const override;

    std::vector<std::string> getVariables(std::string variables_set) const override;

    inline bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
    }
};

class PLANE233Type : public ElementType
{
public:
    PLANE233Type(const ElementTypeProperties &properties);

    std::vector<std::string> getElementLoads() const override;

    std::vector<std::string> getActiveVariables(std::string variables_set) const override;

    std::vector<std::string> getVariables(std::string variables_set) const override;

    void populate_matrix_idx_size_map() override;
    void populate_derived_matrix_idx_size_map() override;
    void populate_vector_idx_size_map() override;

    inline bool is_conductive() const override
    {
        return properties.keyopts[0] == 1;
    }
};

class ElementTypesCollector
{
public:
    ElementTypesCollector() {}

    ElementTypesCollector(ElementTypesMapCSV &elementTypesMapCSV)
    {
        collectElementTypes(elementTypesMapCSV);
    }

    ElementTypesCollector(std::string &filename)
    {
        ElementTypesMapCSV *elementTypesMapCSV = readElementTypes(filename);
        collectElementTypes(*elementTypesMapCSV);
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