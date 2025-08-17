#ifndef ELEMENT_TYPES_H
#define ELEMENT_TYPES_H

#include <unordered_map>
#include <set>
#include <string>
#include <memory>
#include <vector>
#include <libmesh/libmesh.h> // Include the libMesh header
#include <libmesh/quadrature_gauss.h>
#include <libmesh/string_to_enum.h>
#include <libmesh/fe.h>

#include "fe_data_object.h"
#include "globals.h" // Include the header for global constants

using namespace libMesh;

class ElementType;

enum ElementTypeEnum : unsigned int
{
    SOLID97 = 97,
    SOLID236 = 236,
    PLANE13 = 13,
    PLANE233 = 233,
    INFIN110 = 110,
    INFIN111 = 111,
    INFIN112 = 112
};

enum PhysicalRegionEnum : unsigned int
{
    NON_MAGNETIC,
    SOFT_MAGNETIC,
    HARD_MAGNETIC,
    STRANDED_CONDUCTOR,
    SOLID_CONDUCTOR,
    SOLID_CONDUCTOR_CIRCUIT_COUPLED,
    STRANDED_COIL,
    STRANDED_COIL_VOLTAGE_FED,
    STRANDED_COIL_CIRCUIT_COUPLED,
    UNKNOWN
};

typedef std::unordered_map<unsigned int, std::unique_ptr<ElementType>> ElementTypesMap;
typedef std::map<Globals::GlobalVariableName, std::set<subdomain_id_type>> VariableToElementIdsMapType;
typedef std::map<Globals::GlobalVariableName, Globals::VariableInfo> VariableInformation;

// This is probably one of the most important classes of the program (ETYPE in APDL).
// It holds the ElementTypes ID numbers (formulations) and all the information relevant for the formulation of the problem:
// The physical region represented by the element (KEYOPT1)
//
// Axisymmetric options (KEYOPT3)
// Last, it contains an array of libMesh FE objects, which hold all the shape functions, quadrature points, etc.
// These infomration will be used to build by the ElementBuilder e.g. to build the element matrices
// Also, each element has associated DOF variables (e.g. SOLI97, see APDL, can have A_x, A_y, A_z, VOLT...)
// given a mesh containing the element type IDs for each element, the SystemBuilder will distribute properly the variables
// in the proper subdomains of the mesh. The ElementType contains info about name, number of variables, if they are active, etc.
class ElementType
{
public:
    // Constructor
    ElementType(const std::vector<unsigned int> &keyopts)
        : keyopts(keyopts)
    {
        physical_region = keyopts.empty() ? PhysicalRegionEnum::UNKNOWN : static_cast<PhysicalRegionEnum>(keyopts[1]);
        order = static_cast<Order>(keyopts[2]);
    }

    // Virtual methods to be implemented by derived classes
    inline std::vector<Globals::GlobalVariableName> get_active_variables(Globals::VariableCategory variable_cat) const
    {
        std::vector<Globals::GlobalVariableName> active_vars;
        for (const auto &pair : variable_info)
        {
            const Globals::GlobalVariableName &var_name = pair.first;
            const Globals::VariableInfo &var_data = pair.second;

            if (var_data.category == variable_cat)
            {
                active_vars.push_back(var_name);
            }
        }
        return active_vars;
    }
    // // Sets system_var_num for all variables in variable_info incrementally from 0
    // void set_system_var_num()
    // {
    //     unsigned int idx = 0;
    //     for (auto &pair : variable_info)
    //     {
    //         if (pair.second.is_active)
    //         {
    //             pair.second.system_var_num = idx++;
    //         }
    //     }
    // }

    // Getters
    inline const std::vector<unsigned int> &get_keyopts() const { return keyopts; }
    // Getter for ElementTypeData
    inline const VariableInformation &get_variable_information() const
    {
        return variable_info;
    }
    inline VariableInformation &get_variable_information()
    {
        return variable_info;
    }
    inline unsigned int get_active_dof_vars_size() const
    {
        return active_dof_vars_size;
    }
    inline PhysicalRegionEnum get_physical_region() const
    {
        return physical_region;
    }

    // Virtual methods for populating maps
    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_vector_idx_size_map() = 0;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_vector_idx_size_map() = 0;

    std::vector<std::unique_ptr<FEAbstract>> fe_array; // Finite element object

    const std::vector<unsigned int> keyopts; // Element type keyopts

    virtual void set_variable_info() = 0;

    virtual void set_fe_array() = 0;

    void set_variables_order();

protected:
    VariableInformation variable_info;
    unsigned int active_dof_vars_size;  // Number of active DOF variables
    PhysicalRegionEnum physical_region; // Flag to indicate the physical region of the element type
    void set_active_dof_vars_size();
    Order order;
};

class AVElementType : public ElementType
{
public:
    // Constructor
    AVElementType(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override;

    virtual void set_fe_array() override;
};

class AV3DElementType : public AVElementType
{
public:
    // Constructor
    AV3DElementType(const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_vector_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_vector_idx_size_map() override;

    virtual void set_variable_info() override;

    virtual void set_fe_array() override;
};

class AV2DElementType : public AVElementType
{
public:
    // Constructor
    AV2DElementType(const std::vector<unsigned int> &keyopts);

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_dof_matrix_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_matrix_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_matrix_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_vector_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_nodal_solution_vector_idx_size_map() override;

    virtual std::unique_ptr<std::vector<std::vector<unsigned int>>> get_element_solution_vector_idx_size_map() override;

    virtual void set_variable_info() override;

    virtual void set_fe_array() override;
};

class SOLID97Type : public AV3DElementType
{
public:
    SOLID97Type(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override;

    virtual void set_fe_array() override;
};

class SOLID236Type : public AV3DElementType
{
public:
    SOLID236Type(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override {};

    virtual void set_fe_array() override {};
};

class PLANE13Type : public AV2DElementType
{
public:
    PLANE13Type(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override;

    virtual void set_fe_array() override;
};

class PLANE233Type : public AV2DElementType
{
public:
    PLANE233Type(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override {};

    virtual void set_fe_array() override {};
};

// class InfFEBaseEtype
// {
// public:
//     std::vector<std::unique_ptr<FEAbstract>> inf_fe_array;
// };

class INFIN110Type : public PLANE13Type
{
public:
    INFIN110Type(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override;

    virtual void set_fe_array() override;
};

class INFIN111Type : public SOLID97Type
{
public:
    INFIN111Type(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override;

    virtual void set_fe_array() override;
};

class INFIN112Type : public SOLID236Type
{
public:
    INFIN112Type(const std::vector<unsigned int> &keyopts);

    virtual void set_variable_info() override { libmesh_error_msg("INFIN112Type::set_variable_info() not implemented yet."); };

    virtual void set_fe_array() override { libmesh_error_msg("INFIN112Type::set_fe_array() not implemented yet."); };
};
#endif // ELEMENT_TYPES_H