#ifndef APEC_M_SETTINGS_H
#define APEC_M_SETTINGS_H

#include <vector>

#include <libmesh/point.h>

#include "fe_data_object.h"

class ApecMSettingsObject : public BinaryDataFEDataObject
{
public:
    ApecMSettingsObject();
    explicit ApecMSettingsObject(std::ifstream &file);

    void read_data(std::ifstream &file) override;

    bool empty() const { return BinaryDataFEDataObject::empty(); }

    const Point &target_axis() const { return target_central_axis; }
    dof_id_type target_axis_element_id() const { return target_central_axis_element_id; }
    Real major_radius() const { return R0; }
    unsigned int max_picard_iterations() const { return max_picard_iterations_; }
    Real min_distance_to_axis() const { return min_distance; }
    Real cpgain() const { return cp_gain; }
    Real tolerance() const { return tolerance_value; }
    Real bradx() const { return bradx_value; }
    Real brady() const { return brady_value; }
    Real bverx() const { return bverx_value; }
    Real bvery() const { return bvery_value; }
    Real alpha_mp_() const { return alpha_mp; }
    Real alpha_np_() const { return alpha_np; }
    Real plasma_currents() const { return plasma_current; }
    Real beta_profile() const { return beta_p; }
    bool print_full_output() const { return print_full_output_; }
    const std::string &evolution_mesh() const { return evolution_mesh_; }
    const std::string &equilibrium_mesh() const { return equilibrium_mesh_; }
    const std::vector<std::vector<dof_id_type>> &toroidal_slices_nodes() const { return toroidal_slices_nodes_; }
    const std::vector<std::vector<dof_id_type>> &toroidal_slices_elems() const { return toroidal_slices_elems_; }

private:
    void parse_columns();

    Point target_central_axis = Point();
    dof_id_type target_central_axis_element_id = 0;
    Real R0 = 0.0;
    unsigned int max_picard_iterations_ = 1;
    Real min_distance = 0.3;
    Real cp_gain = 0.0;
    Real tolerance_value = 1e-6;
    Real bradx_value = 0.0;
    Real brady_value = 0.0;
    Real bverx_value = 0.0;
    Real bvery_value = 0.0;
    Real alpha_mp = 1.0001;
    Real alpha_np = 0.83197;
    Real plasma_current = 0.0;
    Real beta_p = 0.0;
    bool print_full_output_ = false;
    std::string evolution_mesh_ = "";
    std::string equilibrium_mesh_ = "";
    unsigned int n_slices_ = 3;
    std::vector<std::vector<dof_id_type>> toroidal_slices_nodes_;
    std::vector<std::vector<dof_id_type>> toroidal_slices_elems_;
};

#endif // APEC_M_SETTINGS_H
