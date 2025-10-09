memory >
#include <vector>
#include <string>

    // -----------------------------------------------------------------------------
    // Data container holding only primitive STL types
    // -----------------------------------------------------------------------------
    struct ApecMPeriodicBoundaryData
{
    unsigned int main_boundary;
    unsigned int paired_boundary;
    std::vector<std::vector<Real>> periodic_matrix;
    std::vector<std::string> variables;
    std::vector<std::vector<Real>> variables_transformation_matrix;
};

// -----------------------------------------------------------------------------
// Periodic boundary object
// -----------------------------------------------------------------------------
class ApecMPeriodicBoundary
{
public:
    explicit ApecMPeriodicBoundary(const ApecMPeriodicBoundaryData &data)
        : main_boundary(data.main_boundary),
          paired_boundary(data.paired_boundary),
          periodic_matrix(data.periodic_matrix),
          variables(data.variables),
          variables_transformation_matrix(data.variables_transformation_matrix)
    {
    }

    unsigned int get_main_boundary() const { return main_boundary; }
    unsigned int get_paired_boundary() const { return paired_boundary; }
    const std::vector<std::vector<Real>> &get_periodic_matrix() const { return periodic_matrix; }
    const std::vector<std::string> &get_variables() const { return variables; }
    const std::vector<std::vector<Real>> &get_variables_transformation_matrix() const { return variables_transformation_matrix; }

private:
    unsigned int main_boundary;
    unsigned int paired_boundary;
    std::vector<std::vector<Real>> periodic_matrix;
    std::vector<std::string> variables;
    std::vector<std::vector<Real>> variables_transformation_matrix;
};

// -----------------------------------------------------------------------------
// Factory class for ApecMPeriodicBoundary
// -----------------------------------------------------------------------------
class ApecMPeriodicBoundaryFactory : public FEDataObject
{
public:
    ApecMPeriodicBoundaryFactory() = default;

    // Populate data map from a vtkFieldData (originating from a vtkTable)
    void populate_data_map_from_vtk(vtkFieldData *field_data);

    // Create ApecMPeriodicBoundary objects from data map
    void create_boundaries_from_data();

    inline const std::vector<std::unique_ptr<ApecMPeriodicBoundary>> &get_boundaries() const
    {
        return boundaries;
    }

    inline const std::vector<ApecMPeriodicBoundaryData> &get_data_map() const
    {
        return data_map;
    }

private:
    std::vector<ApecMPeriodicBoundaryData> data_map;
    std::vector<std::unique_ptr<ApecMPeriodicBoundary>> boundaries;
};

#endif // APECM_PERIODIC_BOUNDARY_FACTORY_H