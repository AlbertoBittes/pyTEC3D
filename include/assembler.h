#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "libmesh/equation_systems.h"

using namespace libMesh;

class AssemblerLinear : public System::Assembly
{
public:
    AssemblerLinear(EquationSystems &eq_sys)
        : eq_sys(eq_sys),
          elem_builder(),
          load_step(nullptr),
          old_load_step(nullptr),
          element_types_collector(eq_sys.parameters.get<std::shared_ptr<ElementTypesCollector>>("element_types_collector"))
    {
    }

    void set_elem_builder(ElementBuilder elem_builder)
    {
        this->elem_builder = elem_builder;
    }
    void set_element_types_collector(std::shared_ptr<ElementTypesCollector> element_types_collector)
    {
        this->element_types_collector = element_types_collector;
    }

    std::shared_ptr<LoadStep> get_load_step() const
    {
        return load_step;
    }
    std::shared_ptr<LoadStep> get_old_load_step() const
    {
        return old_load_step;
    }
    ElementBuilder get_elem_builder() const
    {
        return elem_builder;
    }
    const EquationSystems &get_eq_sys() const
    {
        return eq_sys;
    }
    std::shared_ptr<ElementTypesCollector> get_element_types_collector() const
    {
        return element_types_collector;
    }

    void assemble() override;

    void build_write_element_matrices();

    bool check_update_matrix();

    bool check_compute_matrix();

    void compute_derived_quantities();

    void update_rhs();

private:
    EquationSystems &eq_sys;
    ElementBuilder elem_builder;
    std::shared_ptr<LoadStep> load_step;
    std::shared_ptr<LoadStep> old_load_step;
    std::shared_ptr<ElementTypesCollector> element_types_collector;
};

#endif // ASSEMBLER_H