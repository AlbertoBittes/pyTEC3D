#include "binary_reader.h"

#include <fstream>
#include <memory>
#include <string>
#include <filesystem>

#include <libmesh/explicit_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/petsc_linear_solver.h>
#include "apec_m_equation_systems.h"
#include "apec_m_system.h"
#include "element_builder.h"
#include "element_types_factory.h"
#include "solver_configuration.h"
#include "session.h"
#include "time_step.h"
#include "vtk_data_map_object.h"

BinaryReader::BinaryReader()
{
    set_disallowed_labels_for_analysis_mode(false);
}

std::unique_ptr<std::ios> BinaryReader::open_file_stream(const std::filesystem::path &filepath,
                                                         FileMode mode)
{
    switch (mode)
    {
    case FileMode::Read:
    {
        auto stream = std::make_unique<std::ifstream>(filepath, std::ios::binary);
        if (!stream->is_open())
            libmesh_error_msg("BinaryReader: failed to open file for reading '" + filepath.string() + "'.");
        return stream;
    }
    case FileMode::Write:
    {
        const auto folder = filepath.parent_path();
        if (!folder.empty() && !std::filesystem::exists(folder))
            std::filesystem::create_directories(folder);

        auto stream = std::make_unique<std::ofstream>(filepath, std::ios::binary | std::ios::trunc);
        if (!stream->is_open())
            libmesh_error_msg("BinaryReader: failed to open file for writing '" + filepath.string() + "'.");
        return stream;
    }
    default:
        libmesh_error_msg("BinaryReader: unsupported file mode.");
    }
    return nullptr;
}

std::vector<std::unique_ptr<std::ios>> BinaryReader::open_file_streams(
    const std::filesystem::path &folder,
    const std::vector<std::string> &filenames,
    FileMode mode)
{
    std::vector<std::unique_ptr<std::ios>> streams;
    streams.reserve(filenames.size());
    for (const auto &name : filenames)
        streams.push_back(open_file_stream(folder / name, mode));
    return streams;
}

void BinaryReader::close_file_streams(std::vector<std::unique_ptr<std::ios>> &streams)
{
    for (auto &stream : streams)
        stream.reset();
    streams.clear();
}

void BinaryReader::set_disallowed_labels(std::unordered_set<BinaryLabelId> labels)
{
    disallowed_labels = std::move(labels);
}

void BinaryReader::set_disallowed_labels_for_analysis_mode(bool analysis_mode)
{
    disallowed_labels.clear();
    if (analysis_mode)
        return;

    disallowed_labels.insert(BinaryLabelId::TIME_PARAMETERS);
    disallowed_labels.insert(BinaryLabelId::POINT_DATA);
    disallowed_labels.insert(BinaryLabelId::CELL_DATA);
}

void BinaryReader::read_binary_fem_data(const std::string &filepath,
                                        const std::shared_ptr<FEMDataRepository> &repository,
                                        const std::shared_ptr<ApecMEquationSystems> &equation_systems)
{
    const std::filesystem::path path{filepath};
    if (path.empty())
        libmesh_error_msg("BinaryReader: invalid fem data filepath '" + filepath + "'.");

    if (repository)
    {
        auto session = repository->get_object<Session>();
        if (session)
            set_disallowed_labels_for_analysis_mode(session->get_mode() == Session::Mode::Analysis);
    }

    auto stream = open_file_stream(path, FileMode::Read);
    auto *file_stream = dynamic_cast<std::ifstream *>(stream.get());
    if (!file_stream)
        libmesh_error_msg("BinaryReader: unexpected stream type when opening '" + filepath + "'.");

    read_binary_fem_data_from_stream(*file_stream, repository, equation_systems);
    stream.reset();
}

void BinaryReader::read_binary_fem_data_from_stream(std::ifstream &file,
                                                    const std::shared_ptr<FEMDataRepository> &repository,
                                                    const std::shared_ptr<ApecMEquationSystems> &equation_systems)
{

    std::uint32_t section_count = 0;
    file.read(reinterpret_cast<char *>(&section_count), sizeof(section_count));
    for (std::uint32_t idx = 0; idx < section_count; ++idx)
    {
        std::uint32_t raw_label = 0;
        file.read(reinterpret_cast<char *>(&raw_label), sizeof(raw_label));

        const auto label = BinaryLabelTraits::from_raw(raw_label);
        dispatch_read(label, file, repository, equation_systems);
    }
}

void BinaryReader::dispatch_read(BinaryLabelId label,
                                 std::ifstream &file,
                                 const std::shared_ptr<FEMDataRepository> &repository,
                                 const std::shared_ptr<ApecMEquationSystems> &equation_systems)
{
    if (disallowed_labels.count(label) != 0)
        throw BinaryLabelDisallowedError("BinaryReader: section '" + std::string(BinaryLabelTraits::to_string(label)) + "' is not allowed in the current context.");

    switch (label)
    {
    case BinaryLabelId::TIME_PARAMETERS:
    {
        auto old_time_step = repository->get_object<TimeStep>();
        auto time_step = std::make_shared<TimeStep>(file, old_time_step);
        repository->set_object(time_step);
        break;
    }
    case BinaryLabelId::MATERIALS:
    {
        auto materials = std::make_shared<Materials>(file);
        repository->set_object(materials);
        break;
    }
    case BinaryLabelId::CONSTANTS:
    {
        auto reals = std::make_shared<Constants>(file);
        repository->set_object(reals);
        break;
    }
    case BinaryLabelId::ELEMENT_TYPES:
    {
        auto element_types = std::make_shared<ElementTypesFactory>(file);
        repository->set_object(element_types);
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
            libmesh_error_msg("BinaryReader: cannot read VTK data without a valid ApecMEquationSystems reference.");

        auto n_sys = equation_systems->n_systems();
        for (unsigned int i = 0; i != n_sys; ++i)
            if (!equation_systems->get_system(i).is_initialized())
                libmesh_error_msg("BinaryReader: all equation systems must be initialized before reading VTK data.");

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
