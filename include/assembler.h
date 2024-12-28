#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "libmesh/equation_systems.h"
#include "../include/materials_properties.h"
#include "../include/json.hpp"
#include <vtkCellData.h>

using namespace libMesh;

class AssemblerLinear : public System::Assembly
{
public:
    AssemblerLinear(EquationSystems &eq_sys, nlohmann::json config)
        : eq_sys(eq_sys), config(config)
    {
        transient = config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
        ferromagnetic = config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
        fe_family = config["Finite Elements Settings"]["fe_family"].get<std::string>();
        formulation = config["Finite Elements Settings"]["formulation"].get<std::string>();
        order = config["Finite Elements Settings"]["approx_order"].get<std::string>();
        dim = eq_sys.get_mesh().mesh_dimension();
    }

    void assemble() override;

    void update_rhs();

    // Methods to insert and get values from the maps
    void insert_ke_matrix(int key, const DenseMatrix<Number> &matrix);
    DenseMatrix<Number> &get_ke_matrix(int key);

    void insert_ce_matrix(int key, const DenseMatrix<Number> &matrix);
    DenseMatrix<Number> &get_ce_matrix(int key);

    void insert_rhs_vector(int key, const DenseVector<Number> &vector);
    DenseVector<Number> &get_rhs_vector(int key);

    void insert_old_rhs_vector(int key, const std::vector<Real> &vector);
    std::vector<Real> &get_old_rhs_vector(int key);

    void insert_dof_indices(int key, const std::vector<dof_id_type> &dof_indices);
    std::vector<dof_id_type> &get_dof_indices(int key);

private:
    EquationSystems &eq_sys;
    bool transient;
    bool ferromagnetic;
    std::string fe_family;
    std::string formulation;
    std::string order;
    int dim;
    nlohmann::json config;
    std::map<dof_id_type, DenseMatrix<Number>> ke_matrix_map;
    std::map<dof_id_type, DenseMatrix<Number>> ce_matrix_map;
    std::map<dof_id_type, DenseVector<Number>> rhs_map;
    std::map<dof_id_type, std::vector<Real>> old_rhs_map;
    // Here I assume static dofs
    std::map<dof_id_type, std::vector<dof_id_type>> dof_indices;
};

#endif // ASSEMBLER_H