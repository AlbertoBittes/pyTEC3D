#ifndef ELEMENT_TYPES_FACTORY_H
#define ELEMENT_TYPES_FACTORY_H

#include "element_types.h"
#include "fe_data_object.h"
#include <regex>
#include <variant>

// -----------------------------
// Factory to build ElementTypes from vtkFieldData
// -----------------------------
class ElementTypesFactory
    : public vtkFieldDataFEDataObject<
          std::unordered_map<unsigned int,
                             std::vector<std::variant<
                                 int, Real, std::string,
                                 std::vector<std::vector<int>>,
                                 std::vector<std::vector<Real>>,
                                 std::vector<std::vector<std::string>>>>>>
{
public:
    ElementTypesFactory();

    // Stage 1 — populate intermediate data_map
    void populate_intermediate_map_from_vtk(vtkFieldData *field_data) override;

    // Stage 2 — populate final element_types_map
    void populate_final_map() override;

    // Broadcast intermediate map
    // void broadcast(const libMesh::Parallel::Communicator &comm) override;

    // Accessors
    inline ElementTypesMap *get_element_types_map() const { return element_types_map.get(); }
    ElementType *get_element_type(subdomain_id_type element_id);
    std::unique_ptr<VariableToElementIdsMapType> associate_variables_with_element_IDs(Globals::SystemType system_type) const;

    // Add element type manually
    void add_element_type(subdomain_id_type element_id, std::unique_ptr<ElementType> element_type);

    // Factory method
    static std::unique_ptr<ElementType> create_element_type(ElementTypeEnum type,
                                                            const std::vector<unsigned int> &keyopts);

private:
    std::shared_ptr<ElementTypesMap> element_types_map;
};

#endif // ELEMENT_TYPES_FACTORY_H