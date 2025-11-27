#ifndef APEC_M_SETTINGS_H
#define APEC_M_SETTINGS_H

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
    Real major_radius() const { return R0; }
    unsigned int max_picard_iterations() const { return max_picard_iterations_; }
    Real min_distance_to_axis() const { return min_distance; }
    Real cpgain() const { return cp_gain; }
    Real tolerance() const { return tolerance_value; }
    Real bradx() const { return bradx_value; }
    Real brady() const { return brady_value; }
    Real bverx() const { return bverx_value; }
    Real bvery() const { return bvery_value; }
    Real alpha() const { return alpha_value; }
    Real plasma_currents() const { return plasma_current; }
    Real beta_profile() const { return beta_p; }

private:
    void parse_columns();

    Point target_central_axis = Point();
    Real R0 = 0.0;
    unsigned int max_picard_iterations_ = 1;
    Real min_distance = 0.3;
    Real cp_gain = 0.0;
    Real tolerance_value = 1e-6;
    Real bradx_value = 0.0;
    Real brady_value = 0.0;
    Real bverx_value = 0.0;
    Real bvery_value = 0.0;
    Real alpha_value = 1.0;
    Real plasma_current = 0.0;
    Real beta_p = 0.0;
};

#endif // APEC_M_SETTINGS_H
