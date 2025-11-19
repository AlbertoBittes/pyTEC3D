#include "binary_reader.h"

#include <exception>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <filesystem>

#include <libmesh/explicit_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/petsc_linear_solver.h>

#include "apec_m_eq_system_object.h"
#include "apec_m_periodic_boundary_factory.h"
#include "apec_m_settings.h"
#include "apec_m_system_objects.h"
#include "element_builder.h"
#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "solver_configuration.h"
#include "time_step.h"
#include "vtk_data_map_object.h"

namespace
{
    struct EmatRecordLocation
    {
        bool found = false;
        std::streampos offset = std::streampos(-1);
    };

    bool try_extract_rank(const std::string &filename, const std::string &prefix, unsigned int &rank)
    {
        if (filename.compare(0, prefix.size(), prefix) != 0)
            return false;

        const auto start = prefix.size();
        const auto dot_pos = filename.find('.', start);
        const auto number_str = filename.substr(start, dot_pos == std::string::npos ? std::string::npos : dot_pos - start);
        if (number_str.empty())
            return false;

        try
        {
            size_t consumed = 0;
            const auto parsed_rank = std::stoul(number_str, &consumed);
            if (consumed != number_str.size())
                return false;
            rank = static_cast<unsigned int>(parsed_rank);
            return true;
        }
        catch (const std::exception &)
        {
            return false;
        }
    }

    void skip_matrix_entry(std::ifstream &stream, unsigned int num_containers)
    {
        for (unsigned int container_idx = 0; container_idx < num_containers; ++container_idx)
        {
            unsigned int rows = 0, cols = 0;
            stream.read(reinterpret_cast<char *>(&rows), sizeof(rows));
            stream.read(reinterpret_cast<char *>(&cols), sizeof(cols));
            if (!stream)
                libmesh_error_msg("BinaryReader: corrupt emat file while skipping matrix metadata.");

            const auto value_count = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
            const auto bytes_to_skip = static_cast<std::streamoff>(value_count * sizeof(Number));
            stream.seekg(bytes_to_skip, std::ios::cur);
            if (!stream)
                libmesh_error_msg("BinaryReader: corrupt emat file while skipping matrix data.");
        }
    }

    EmatRecordLocation locate_element_record(std::ifstream &stream, unsigned int target_elem_id)
    {
        EmatRecordLocation result;
        stream.clear();
        stream.seekg(0, std::ios::beg);

        while (stream)
        {
            const auto entry_start = stream.tellg();
            unsigned int elem_id = 0;
            stream.read(reinterpret_cast<char *>(&elem_id), sizeof(elem_id));
            if (!stream)
                break;

            unsigned int num_containers = 0;
            stream.read(reinterpret_cast<char *>(&num_containers), sizeof(num_containers));
            if (!stream)
                libmesh_error_msg("BinaryReader: corrupt emat file while reading container count.");

            if (elem_id == target_elem_id)
            {
                result.found = true;
                result.offset = entry_start;
                return result;
            }

            skip_matrix_entry(stream, num_containers);
        }

        stream.clear();
        return result;
    }
} // namespace

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

void BinaryReader::read_binary_fem_data(const std::string &filepath,
                                        const std::shared_ptr<FEMDataRepository> &repository,
                                        const std::shared_ptr<ApecMEquationSystems> &equation_systems)
{
    const std::filesystem::path path{filepath};
    if (path.empty())
        libmesh_error_msg("BinaryReader: invalid fem data filepath '" + filepath + "'.");

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
    {
        auto binary_data = std::make_shared<BinaryDataFEDataObject>(file);
        repository->set_object(binary_data);
        break;
    }
    case BinaryLabelId::APEC_SETTINGS:
    {
        auto settings = std::make_shared<ApecMSettingsObject>(file);
        repository->set_object(settings);
        break;
    }
    case BinaryLabelId::DOF_INTERPOLATION_SETTINGS:
    {
        auto interpolation_settings = std::make_shared<DofInterpolationSettingsObject>(file);
        repository->set_object(interpolation_settings);
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

bool BinaryReader::read_solution_matrices_for_element(ElementBuilder &builder,
                                                      dof_id_type element_id,
                                                      const std::filesystem::path &emat_folderpath,
                                                      bool include_nodal_matrices)
{
    libmesh_error_msg_if(emat_folderpath.empty(),
                         "BinaryReader: invalid emat folder path supplied for solution reconstruction.");
    libmesh_error_msg_if(!std::filesystem::exists(emat_folderpath),
                         "BinaryReader: emat folder '" + emat_folderpath.string() + "' does not exist.");
    libmesh_error_msg_if(!std::filesystem::is_directory(emat_folderpath),
                         "BinaryReader: emat path '" + emat_folderpath.string() + "' is not a directory.");

    std::map<unsigned int, std::filesystem::path> element_files;
    std::map<unsigned int, std::filesystem::path> nodal_files;

    for (const auto &entry : std::filesystem::directory_iterator(emat_folderpath))
    {
        if (!entry.is_regular_file())
            continue;

        const auto filename = entry.path().filename().string();
        unsigned int rank = 0;
        if (try_extract_rank(filename, "element_sol_mtrx_", rank))
            element_files[rank] = entry.path();
        else if (try_extract_rank(filename, "nodal_sol_mtrx_", rank))
            nodal_files[rank] = entry.path();
    }

    libmesh_error_msg_if(element_files.empty(),
                         "BinaryReader: no element solution matrix files found in '" + emat_folderpath.string() + "'.");

    const auto elem_id_uint = static_cast<unsigned int>(element_id);

    for (const auto &file_info : element_files)
    {
        auto element_stream = open_file_stream(file_info.second, FileMode::Read);
        auto *element_file = dynamic_cast<std::ifstream *>(element_stream.get());
        libmesh_error_msg_if(!element_file,
                             "BinaryReader: unexpected stream type while opening '" + file_info.second.string() + "'.");

        auto element_location = locate_element_record(*element_file, elem_id_uint);
        if (!element_location.found)
            continue;

        if (include_nodal_matrices)
        {
            auto nodal_it = nodal_files.find(file_info.first);
            libmesh_error_msg_if(nodal_it == nodal_files.end(),
                                 "BinaryReader: missing nodal solution matrix file for rank " + std::to_string(file_info.first) + ".");

            auto nodal_stream = open_file_stream(nodal_it->second, FileMode::Read);
            auto *nodal_file = dynamic_cast<std::ifstream *>(nodal_stream.get());
            libmesh_error_msg_if(!nodal_file,
                                 "BinaryReader: unexpected stream type while opening '" + nodal_it->second.string() + "'.");

            auto nodal_location = locate_element_record(*nodal_file, elem_id_uint);
            libmesh_error_msg_if(!nodal_location.found,
                                 "BinaryReader: element " + std::to_string(elem_id_uint) + " not found in nodal solution matrices file '" + nodal_it->second.string() + "'.");

            nodal_file->clear();
            nodal_file->seekg(nodal_location.offset);
            BinaryReader::read_element_matrices(builder, *nodal_file, ElementBuilder::MatricesSet::NODAL_SOLUTION_MATRICES);
        }

        element_file->clear();
        element_file->seekg(element_location.offset);
        BinaryReader::read_element_matrices(builder, *element_file, ElementBuilder::MatricesSet::ELEMENT_SOLUTION_MATRICES);
        return true;
    }

    return false;
}
