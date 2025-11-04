#include "binary_reader.h"

#include <fstream>
#include <memory>

#include <libmesh/explicit_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/petsc_linear_solver.h>

#include "element_builder.h"
#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "solver_configuration.h"
#include "time_step.h"
#include "vtk_data_map_object.h"
#include "apec_m_periodic_boundary_factory.h"
#include "apec_m_system_objects.h"

std::unique_ptr<std::ios> BinaryReader::open_file_stream(const std::filesystem::path &folder,
                                                         const std::string &filename,
                                                         FileMode mode)
{
    if (!std::filesystem::exists(folder))
        std::filesystem::create_directories(folder);

    auto filepath = folder / filename;
    std::ios::openmode openmode = std::ios::binary;
    if (mode == FileMode::Read)
        openmode |= std::ios::in;
    else
        openmode |= std::ios::out | std::ios::trunc;

    auto stream = std::make_unique<std::fstream>(filepath, openmode);
    if (!stream->good())
        libmesh_error_msg("BinaryReader: failed to open file '" + filepath.string() + "'.");

    return stream;
}

std::vector<std::unique_ptr<std::ios>> BinaryReader::open_file_streams(
    const std::filesystem::path &folder,
    const std::vector<std::string> &filenames,
    FileMode mode)
{
    std::vector<std::unique_ptr<std::ios>> streams;
    streams.reserve(filenames.size());
    for (const auto &name : filenames)
        streams.push_back(open_file_stream(folder, name, mode));
    return streams;
}

void BinaryReader::close_file_streams(std::vector<std::unique_ptr<std::ios>> &streams)
{
    for (auto &stream : streams)
        stream.reset();
    streams.clear();
}

void BinaryReader::read_binary_fem_data(const std::filesystem::path &filepath,
                                        const std::shared_ptr<FEMDataRepository> &repository,
                                        const std::shared_ptr<ApecMEquationSystems> &equation_systems)
{
    std::ifstream file(filepath, std::ios::binary);
    if (!file)
        libmesh_error_msg("BinaryReader: failed to open fem data binary '" + filepath.string() + "'.");

    std::uint32_t section_count = 0;
    file.read(reinterpret_cast<char *>(&section_count), sizeof(section_count));
    if (!file)
        libmesh_error_msg("BinaryReader: failed to read section count from '" + filepath.string() + "'.");

    for (std::uint32_t idx = 0; idx < section_count; ++idx)
    {
        std::uint32_t raw_label = 0;
        file.read(reinterpret_cast<char *>(&raw_label), sizeof(raw_label));
        if (!file)
            libmesh_error_msg("BinaryReader: failed to read section label from '" + filepath.string() + "'.");

        const auto label = BinaryLabelTraits::from_raw(raw_label);
        dispatch_read(label, file, repository, equation_systems);
    }
}

void BinaryReader::dispatch_read(BinaryLabelId label,
                                 std::ifstream &file,
                                 const std::shared_ptr<FEMDataRepository> &repository,
                                 const std::shared_ptr<ApecMEquationSystems> &equation_systems)
{
    switch (label)
    {
    case BinaryLabelId::TIME_PARAMETERS:
    {
        auto time_step = std::make_shared<TimeStep>(file);
        repository->set_object(time_step);
        break;
    }
    case BinaryLabelId::MATERIALS:
    {
        auto materials = std::make_shared<MaterialsObject>(file);
        repository->set_object(materials);
        break;
    }
    case BinaryLabelId::REAL_CONSTANTS:
    {
        auto reals = std::make_shared<RealConstantsObject>(file);
        repository->set_object(reals);
        break;
    }
    case BinaryLabelId::ELEMENT_TYPES:
    {
        auto element_types = std::make_shared<ElementTypesFactory>(file);
        repository->set_object(element_types);
        break;
    }
    case BinaryLabelId::PERIODIC_BOUNDARIES:
    {
        auto factory = std::make_shared<ApecMPeriodicBoundaryFactory>(file);
        repository->set_object(factory);
        break;
    }
    case BinaryLabelId::SOLVER_SETTINGS:
    {
        auto solver_config = std::make_shared<PetscSolverConfiguration>(file);
        repository->set_object(solver_config);
        break;
    }
    case BinaryLabelId::POINT_DATA:
    case BinaryLabelId::CELL_DATA:
    {
        if (!equation_systems)
            libmesh_error_msg("BinaryReader: vtk data section encountered without an equation system.");

        auto vtk_data = std::make_shared<VTKDataMapFEDataObject>(equation_systems);
        if (BinaryLabelTraits::is_cell_data(label))
        {
            vtk_data->set_is_element_data(true);
            vtk_data->add_element_loads_system();
        }
        else
        {
            vtk_data->set_is_element_data(false);
            vtk_data->add_dof_system();
            vtk_data->add_nodal_loads_system();
        }

        vtk_data->read_data(file);
        break;
    }
    case BinaryLabelId::DOF_COUPLINGS:
    case BinaryLabelId::APEC_PARAMETERS:
    {
        auto binary_data = std::make_shared<BinaryDataFEDataObject>(file);
        repository->set_object(binary_data);
        break;
    }
    default:
        libmesh_error_msg("BinaryReader: encountered unknown section label '" + std::string(BinaryLabelTraits::to_string(label)) + "'.");
    }
}

void BinaryReader::write_element_matrices(ElementBuilder &builder, std::ofstream &outFile, unsigned int matrices_set)
{
    auto &elem_matrices_vector = builder.matrices[matrices_set];

    const unsigned int elem_id = builder.elem_data.elem ? builder.elem_data.elem->id() : 0;
    outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));

    const unsigned int num_containers = static_cast<unsigned int>(elem_matrices_vector.size());
    outFile.write(reinterpret_cast<const char *>(&num_containers), sizeof(num_containers));

    for (const auto &elem_matrix : elem_matrices_vector)
    {
        const unsigned int rows = elem_matrix.m();
        const unsigned int cols = elem_matrix.n();
        outFile.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
        outFile.write(reinterpret_cast<const char *>(&cols), sizeof(cols));

        const auto &values = elem_matrix.get_values();
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

void BinaryReader::read_element_matrices(ElementBuilder &builder, std::ifstream &inFile, unsigned int matrices_set)
{
    auto &elem_matrices_vector = builder.matrices[matrices_set];

    unsigned int elem_id = 0;
    inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));

    unsigned int num_containers = 0;
    inFile.read(reinterpret_cast<char *>(&num_containers), sizeof(num_containers));
    elem_matrices_vector.resize(num_containers);

    for (auto &elem_matrix : elem_matrices_vector)
    {
        unsigned int rows = 0, cols = 0;
        inFile.read(reinterpret_cast<char *>(&rows), sizeof(rows));
        inFile.read(reinterpret_cast<char *>(&cols), sizeof(cols));

        elem_matrix.resize(rows, cols);
        auto &values = elem_matrix.get_values();
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(Number));
    }
}

void BinaryReader::write_element_vectors(ElementBuilder &builder, std::ofstream &outFile, unsigned int vectors_set)
{
    auto &elem_vectors_vector = builder.vectors[vectors_set];

    const unsigned int elem_id = builder.elem_data.elem ? builder.elem_data.elem->id() : 0;
    outFile.write(reinterpret_cast<const char *>(&elem_id), sizeof(elem_id));

    const unsigned int num_containers = static_cast<unsigned int>(elem_vectors_vector.size());
    outFile.write(reinterpret_cast<const char *>(&num_containers), sizeof(num_containers));

    for (const auto &elem_vector : elem_vectors_vector)
    {
        const unsigned int values_size = elem_vector.size();
        outFile.write(reinterpret_cast<const char *>(&values_size), sizeof(values_size));

        const auto &values = elem_vector.get_values();
        outFile.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(Number));
    }
}

void BinaryReader::read_element_vectors(ElementBuilder &builder, std::ifstream &inFile, unsigned int vectors_set)
{
    auto &elem_vectors_vector = builder.vectors[vectors_set];

    unsigned int elem_id = 0;
    inFile.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));

    unsigned int num_containers = 0;
    inFile.read(reinterpret_cast<char *>(&num_containers), sizeof(num_containers));
    elem_vectors_vector.resize(num_containers);

    for (auto &elem_vector : elem_vectors_vector)
    {
        unsigned int values_size = 0;
        inFile.read(reinterpret_cast<char *>(&values_size), sizeof(values_size));

        elem_vector.resize(values_size);
        auto &values = elem_vector.get_values();
        inFile.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(Number));
    }
}
