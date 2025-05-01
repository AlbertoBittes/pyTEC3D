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
    const std::vector<Real> &JxW;
    const std::vector<std::vector<Real>> &phi;
    const std::vector<std::vector<RealGradient>> &dphi;
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
    // Constructor
    ElementType(const ElementTypeProperties &properties, const std::string &order = "default_order", const std::string &fe_family = "default_family", unsigned int quadrature_order = 0, unsigned int dimension = 0)
        : properties(properties), order(order), fe_family(fe_family), quadrature_order(quadrature_order), dimension(dimension), az_dofs(0) {}

    // Virtual methods to be implemented by derived classes
    virtual std::vector<std::string> getActiveVariables(std::string variables_set) const = 0;
    virtual std::vector<std::string> getVariables(std::string variables_set) const = 0;
    virtual std::vector<std::string> getElementLoads() const = 0;
    std::string ElementType::getVariableFEFamily(std::string variable);
    std::string ElementType::getVariableOrder(std::string variable);
    virtual bool is_conductive() const = 0;

    // Getters
    const ElementTypeProperties &get_properties() const { return properties; }
    const std::string &get_order() const { return order; }
    const std::string &get_fe_family() const { return fe_family; }
    unsigned int get_quadrature_order() const { return quadrature_order; }
    unsigned int get_dimension() const { return dimension; }
    unsigned int get_az_dofs() const { return az_dofs; }
    const FEBase *get_fe() const { return fe.get(); }
    const FEVectorBase *get_fe_vec() const { return fe_vec.get(); }

    const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_matrix_idx_size_map() const { return matrix_idx_size_map; }
    const std::map<std::string, std::tuple<unsigned int, unsigned int, unsigned int>> &get_derived_matrix_idx_size_map() const { return derived_matrix_idx_size_map; }
    const std::map<std::string, std::tuple<unsigned int, unsigned int>> &get_vector_idx_size_map() const { return vector_idx_size_map; }

    // Setters
    void set_order(const std::string &new_order) { order = new_order; }
    void set_fe_family(const std::string &new_fe_family) { fe_family = new_fe_family; }
    void set_quadrature_order(unsigned int new_quadrature_order) { quadrature_order = new_quadrature_order; }
    void set_dimension(unsigned int new_dimension) { dimension = new_dimension; }
    void set_az_dofs(unsigned int new_az_dofs) { az_dofs = new_az_dofs; }
    void set_fe(std::unique_ptr<FEBase> new_fe) { fe = std::move(new_fe); }
    // Shape functions
    inline std::shared_ptr<ScalarShapeFunctions> &get_scalar_shape_functions() const
    {
        std::shared_ptr<ScalarShapeFunctions> shape_functions = std::make_shared<ScalarShapeFunctions>(fe->get_JxW(), fe->get_phi(), fe->get_dphi());
        return shape_functions;
    }

    inline std::shared_ptr<VectorShapeFunctions> &get_vector_shape_functions() const
    {
        std::shared_ptr<VectorShapeFunctions> shape_functions = std::make_shared<VectorShapeFunctions>(fe_vec->get_JxW(), fe_vec->get_phi(), fe_vec->get_dphi(), fe_vec->get_curl_phi());
        return shape_functions;
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

protected:
    std::unique_ptr<FEBase> fe = nullptr;           // Finite element object
    std::unique_ptr<FEVectorBase> fe_vec = nullptr; // Finite element vector object
    ElementTypeProperties properties;               // Element type properties

    std::string order;             // Order of the element
    std::string fe_family;         // Finite element family
    unsigned int quadrature_order; // Quadrature order
    unsigned int dimension;        // Dimension of the element
    unsigned int az_dofs;          // Number of DOFs

    // Maps for matrix and vector sizes
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