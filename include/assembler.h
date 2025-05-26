#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "libmesh/equation_systems.h"
#include "../include/element_builder.h" // Include the ElementBuilder header file
#include "../include/load_step.h"
#include "../include/element_types.h" // Include the ElementTypesMap header file

using namespace libMesh;

class LoadStep;

class AssemblerLinear : public System::Assembly, public FEDataObject<ElementTypesCollector, EquationSystems>
{
public:
    AssemblerLinear(EquationSystems *parent = nullptr,
                    std::shared_ptr<ElementTypesCollector> map = nullptr)
        : System::Assembly(),
          FEDataObject<ElementTypesCollector, EquationSystems>(parent, map), // Pass arguments to base
          elem_builder(std::make_shared<ElementBuilder>(this))               // Initialize elem_builder
    {
    }

    void set_elem_builder(ElementBuilder *builder)
    {
        elem_builder.reset(builder);
    }

    ElementBuilder *get_elem_builder() const
    {
        return elem_builder.get();
    }

    void assemble() override;

    void build_write_element_matrices();

    bool check_update_matrix();

    bool check_compute_matrix();

    void compute_derived_quantities();

    void update_rhs();

private:
    std::shared_ptr<ElementBuilder> elem_builder;
};

#endif // ASSEMBLER_H