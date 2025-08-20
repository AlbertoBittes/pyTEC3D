#include "apec_m_eq_system_object.h"
#include "apec_m_system_objects.h"

void ApecMEquationSystems::apply_ramped_effects(const Substep &time_substep)
{
    // Get the current and old solution vectors
    System &main_system = this->get_system<System>("ApecMLinearSystem");
    NumericVector<Number> &current_vector = *main_system.request_vector("node_dof_constraints");
    NumericVector<Number> &old_vector = *main_system.request_vector("old_node_dof_constraints");
    NumericVector<Number> &ramped_vector = *main_system.request_vector("ramped_node_dof_constraints");

    // Call the static method to apply ramping
    ApecMSystem<System>::ramp_vectors(time_substep, old_vector, current_vector, ramped_vector);

    System &element_loads_system = this->get_system<System>("ElementLoadsSystem");
    NumericVector<Number> &element_current_vector = *element_loads_system.current_local_solution.get();
    NumericVector<Number> &element_old_vector = *element_loads_system.request_vector("old_loads");
    NumericVector<Number> &element_ramped_vector = *element_loads_system.request_vector("ramped_loads");

    ApecMSystem<System>::ramp_vectors(time_substep, element_old_vector, element_current_vector, element_ramped_vector);

    System &nodal_loads_system = this->get_system<System>("NodalLoadsSystem");
    NumericVector<Number> &nodal_current_vector = *nodal_loads_system.current_local_solution.get();
    NumericVector<Number> &nodal_old_vector = *nodal_loads_system.request_vector("old_loads");
    NumericVector<Number> &nodal_ramped_vector = *nodal_loads_system.request_vector("ramped_loads");

    ApecMSystem<System>::ramp_vectors(time_substep, nodal_old_vector, nodal_current_vector, nodal_ramped_vector);
}