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
    AssemblerLinear(EquationSystems &eq_sys, vtkCellData *cellData, nlohmann::json config, MaterialsMap &materials)
        : eq_sys(eq_sys), cell_data(cellData), config(config), materials(materials)
    {
        transient = config["Transient Analysis Settings"]["transient_analysis"].get<bool>();
        ferromagnetic = config["Finite Elements Settings"]["nonlinear_analysis_ferrom"].get<bool>();
        fe_family = config["Finite Elements Settings"]["fe_family"].get<std::string>();
        formulation = config["Finite Elements Settings"]["formulation"].get<std::string>();
        order = config["Finite Elements Settings"]["approx_order"].get<std::string>();
        dim = eq_sys.get_mesh().mesh_dimension();
    }

    void assemble() override;

private:
    EquationSystems &eq_sys;
    vtkCellData *cell_data;
    bool transient;
    bool ferromagnetic;
    std::string fe_family;
    std::string formulation;
    std::string order;
    int dim;
    nlohmann::json config;
    MaterialsMap &materials;

    void build_solid97(DenseMatrix<Number> &Ke, DenseVector<Number> &Fe,
                       DenseSubMatrix<Number> &K_no_cond, DenseSubVector<Number> &F_no_cond,
                       const Elem *elem, const std::vector<dof_id_type> &dof_indices,
                       const DofMap &dof_map, dof_id_type el_id, System &system,
                       QGauss &qrule, const std::vector<std::vector<Real>> &phi,
                       const std::vector<std::vector<RealGradient>> &dphi, const std::vector<Real> &JxW,
                       PerfLog &perf_log, vtkDataArray *jsArray, std::vector<std::vector<double>> &mu_matrix,
                       unsigned int n_dofs, std::vector<std::vector<dof_id_type>> &dof_indices_var);

    void build_solid97_cond();

    void build_solid236();

    void build_solid236_cond();

    void AssemblerLinear::build_plane13(DenseMatrix<Number> &Ke, DenseVector<Number> &Fe,
                                        const Elem *elem, const std::vector<dof_id_type> &dof_indices,
                                        const DofMap &dof_map, dof_id_type el_id,
                                        QGauss &qrule, const std::vector<std::vector<Real>> &phi,
                                        const std::vector<std::vector<RealGradient>> &dphi, const std::vector<Real> &JxW,
                                        PerfLog &perf_log, vtkDataArray *jsArray, std::vector<std::vector<double>> &mu_matrix, unsigned int n_dofs);

    void build_plane13_cond();
};

#endif // ASSEMBLER_H