#include "apec_m_periodic_boundary.h"

#include "globals.h"
#include "variables_map.h"
#include <map>

ApecMPeriodicBoundary::ApecMPeriodicBoundary(const boundary_id_type main_b,
                                             const boundary_id_type paired_b,
                                             const std::vector<std::vector<Real>> &periodic_mtx,
                                             TransformationType t,
                                             const std::vector<std::string> &variable_labels_in,
                                             EquationSystems *eq_sys)
    : PeriodicBoundaryBase(),
      FEParentObjectBase<EquationSystems>(eq_sys),
      rotation_matrix(periodic_mtx[0][0], periodic_mtx[0][1], periodic_mtx[0][2],
                      periodic_mtx[1][0], periodic_mtx[1][1], periodic_mtx[1][2],
                      periodic_mtx[2][0], periodic_mtx[2][1], periodic_mtx[2][2]),
      translation_vector(periodic_mtx[0][3], periodic_mtx[1][3], periodic_mtx[2][3]),
      variable_labels(variable_labels_in)
{
    myboundary = main_b;
    pairedboundary = paired_b;

    if (t == INVERSE)
    {
        std::swap(myboundary, pairedboundary);
        rotation_matrix = rotation_matrix.transpose();
        translation_vector = -rotation_matrix * translation_vector;
    }

    this->set_variables();

    auto transform = build_transformation_matrix();
    if (transform.m() > 0)
    {
        if (!(transform.m() == 1 && transform.n() == 1))
            this->set_transformation_matrix(transform);
    }
}

ApecMPeriodicBoundary::ApecMPeriodicBoundary(const ApecMPeriodicBoundary &o, TransformationType t)
    : ApecMPeriodicBoundary(o.myboundary, o.pairedboundary,
                            {{o.rotation_matrix(0, 0), o.rotation_matrix(0, 1), o.rotation_matrix(0, 2), o.translation_vector(0)},
                             {o.rotation_matrix(1, 0), o.rotation_matrix(1, 1), o.rotation_matrix(1, 2), o.translation_vector(1)},
                             {o.rotation_matrix(2, 0), o.rotation_matrix(2, 1), o.rotation_matrix(2, 2), o.translation_vector(2)},
                             {0, 0, 0, 1}},
                            t,
                            o.variable_labels,
                            o.get_parent_object())
{
    // if (!o.has_transformation_matrix())
    // {
    //     this->clear_transformation_matrix();
    // }
    // else
    // {
    //     this->set_transformation_matrix(o.get_transformation_matrix());
    // }
}

Point ApecMPeriodicBoundary::get_corresponding_pos(const Point &pt) const
{
    Point pt_prime;
    pt_prime = rotation_matrix * pt + translation_vector;
    return pt_prime;
}

void ApecMPeriodicBoundary::set_variables()
{
    variables.clear();

    EquationSystems *eq_sys = this->get_parent_object();
    libmesh_error_msg_if(!eq_sys, "EquationSystems not associated with periodic boundary.");

    for (const auto &label : variable_labels)
    {
        auto keys = VariablesMap::get_keys_by_predicate(
            [&label](const VariableInfo &info)
            {
                return info.api_label == label;
            });

        for (const auto &key : keys)
        {
            auto &var_string_name = VariablesMap::get_info_by_key(key).global_label;
            if (eq_sys->get_system(0).has_variable(var_string_name))
            {
                unsigned int var_index = eq_sys->get_system(0).variable_number(var_string_name);
                set_variable(var_index);
            }
        }
    }
}

std::unique_ptr<PeriodicBoundaryBase> ApecMPeriodicBoundary::clone(TransformationType t) const
{
    return std::make_unique<ApecMPeriodicBoundary>(*this, t);
}

const TensorValue<Real> &ApecMPeriodicBoundary::get_rotation_matrix() const
{
    return rotation_matrix;
}

const VectorValue<Real> &ApecMPeriodicBoundary::get_translation_vector() const
{
    return translation_vector;
}

const std::vector<std::string> &ApecMPeriodicBoundary::get_variable_labels() const
{
    return variable_labels;
}

void ApecMPeriodicBoundary::clear_transformation_matrix()
{
    this->_transformation_matrix.reset();
}

DenseMatrix<Real> ApecMPeriodicBoundary::build_transformation_matrix() const
{
    const auto &vars_set = this->variables;
    const std::size_t n = vars_set.size();

    DenseMatrix<Real> matrix(n, n);
    matrix.zero();

    if (n == 0)
        return matrix;

    EquationSystems *eq_sys = this->get_parent_object();
    libmesh_error_msg_if(!eq_sys, "EquationSystems not associated with periodic boundary.");
    const System &system = eq_sys->get_system(0);

    std::vector<unsigned int> ordered(vars_set.begin(), vars_set.end());
    std::vector<unsigned int> component_indices(ordered.size());
    std::map<APIVariableName, std::vector<std::size_t>> label_rows;
    std::vector<APIVariableName> label_order;

    for (std::size_t i = 0; i < ordered.size(); ++i)
    {
        unsigned int var_index = ordered[i];
        std::string var_name = system.variable_name(var_index);

        auto key = VariablesMap::get_key_by_predicate_one(
            [&var_name](const VariableInfo &info)
            {
                return info.global_label == var_name;
            });

        const auto &info = VariablesMap::get_info_by_key(key);
        component_indices[i] = info.component_index;

        if (!label_rows.count(info.api_name))
            label_order.push_back(info.api_name);
        label_rows[info.api_name].push_back(i);
    }

    TensorValue<Real> rotation_T = rotation_matrix.transpose();
    const unsigned int rotation_dim = 3;

    for (APIVariableName label : label_order)
    {
        const auto &rows = label_rows[label];
        const std::size_t dim = rows.size();
        if (dim == 0)
            continue;

        for (std::size_t i = 0; i < dim; ++i)
        {
            std::size_t row = rows[i];
            for (std::size_t j = 0; j < dim; ++j)
            {
                std::size_t col = rows[j];
                Real value = 0.0;

                if (dim == rotation_dim)
                {
                    unsigned int comp_i = component_indices[row];
                    unsigned int comp_j = component_indices[col];
                    if (comp_i < rotation_dim && comp_j < rotation_dim)
                        value = rotation_T(comp_i, comp_j);
                    else
                        value = (comp_i == comp_j) ? 1.0 : 0.0;
                }
                else
                {
                    value = (i == j) ? 1.0 : 0.0;
                }

                matrix(row, col) = value;
            }
        }
    }

    return matrix;
}
