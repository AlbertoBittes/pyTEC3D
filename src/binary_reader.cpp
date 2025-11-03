#include "binary_reader.h"

#include <fstream>
#include <memory>

#include <libmesh/explicit_system.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/petsc_linear_solver.h>

#include "element_types_factory.h"
#include "materials_object.h"
#include "real_constants_object.h"
#include "solver_configuration.h"
#include "time_step.h"
#include "vtk_data_map_object.h"
#include "apec_m_periodic_boundary_factory.h"

BinaryReader::BinaryReader(Session &session) : session(session) {}

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
                                        FEMDataRepository &repository,
                                        const std::shared_ptr<ApecMEquationSystems> &equation_systems,
                                        libMesh::PetscLinearSolver<Number> *solver,
                                        const std::shared_ptr<LoadStep> &previous_load_step)
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
        dispatch_read(label, file, repository, equation_systems, solver, previous_load_step);
    }
}

void BinaryReader::dispatch_read(BinaryLabelId label,
                                 std::ifstream &file,
                       FEMDataRepository &repository,
                       const std::shared_ptr<ApecMEquationSystems> &equation_systems,
                       libMesh::PetscLinearSolver<Number> *solver,
                       const std::shared_ptr<LoadStep> &previous_load_step)
{
    switch (label)
    {
    case BinaryLabelId::TIME_PARAMETERS:
    {
        auto time_step = std::make_shared<TimeStep>(file, previous_load_step);
        repository.set_object(time_step);
        break;
    }
    case BinaryLabelId::SOLVER_SETTINGS:
    {
        if (!solver)
            libmesh_error_msg("BinaryReader: solver settings section encountered without a valid solver.");
        auto solver_cfg = std::make_shared<PetscSolverConfiguration>(file, solver);
        repository.set_object(solver_cfg);
        break;
    }
    case BinaryLabelId::MATERIALS:
    {
        auto materials = std::make_shared<MaterialsObject>(file);
        repository.set_object(materials);
        break;
    }
    case BinaryLabelId::REAL_CONSTANTS:
    {
        auto reals = std::make_shared<RealConstantsObject>(file);
        repository.set_object(reals);
        break;
    }
    case BinaryLabelId::ELEMENT_TYPES:
    {
        auto element_types = std::make_shared<ElementTypesFactory>();
        element_types->read_data(file);
        repository.set_object(element_types);
        break;
    }
    case BinaryLabelId::PERIODIC_BOUNDARIES:
    {
        if (!equation_systems)
            libmesh_error_msg("BinaryReader: periodic boundaries section encountered without an equation system.");
        auto factory = std::make_shared<ApecMPeriodicBoundaryFactory>(nullptr, equation_systems.get());
        factory->read_data(file);
        factory->populate_periodic_boundaries_from_map();
        repository.set_object(factory);
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
        // Sections recognised but not yet mapped to concrete C++ objects
        libMesh::out << "BinaryReader: skipping unsupported section '" << BinaryLabelTraits::to_string(label) << "'." << std::endl;

        // Consume payload generically to advance stream
        std::uint32_t entry_count = 0;
        file.read(reinterpret_cast<char *>(&entry_count), sizeof(entry_count));

        std::uint32_t num_keys = 0;
        file.read(reinterpret_cast<char *>(&num_keys), sizeof(num_keys));

        for (std::uint32_t k = 0; k < num_keys; ++k)
        {
            std::uint32_t key_size = 0;
            file.read(reinterpret_cast<char *>(&key_size), sizeof(key_size));
            file.ignore(static_cast<std::streamsize>(key_size));

            int32_t payload_type = 0;
            file.read(reinterpret_cast<char *>(&payload_type), sizeof(payload_type));

            std::uint32_t dtype_len = 0;
            file.read(reinterpret_cast<char *>(&dtype_len), sizeof(dtype_len));
            file.ignore(static_cast<std::streamsize>(dtype_len));

            std::uint32_t ndim = 0;
            file.read(reinterpret_cast<char *>(&ndim), sizeof(ndim));
            std::size_t total_elements = 1;
            for (std::uint32_t d = 0; d < ndim; ++d)
            {
                std::uint32_t extent = 0;
                file.read(reinterpret_cast<char *>(&extent), sizeof(extent));
                total_elements *= static_cast<std::size_t>(extent);
            }

            int32_t element_kind = 0;
            file.read(reinterpret_cast<char *>(&element_kind), sizeof(element_kind));

            if (element_kind == 0)
                file.ignore(static_cast<std::streamsize>(total_elements * sizeof(std::int64_t)));
            else if (element_kind == 1)
                file.ignore(static_cast<std::streamsize>(total_elements * sizeof(double)));
            else if (element_kind == 2)
            {
                std::uint32_t count = 0;
                file.read(reinterpret_cast<char *>(&count), sizeof(count));
                for (std::uint32_t c = 0; c < count; ++c)
                {
                    std::uint32_t strlen = 0;
                    file.read(reinterpret_cast<char *>(&strlen), sizeof(strlen));
                    file.ignore(static_cast<std::streamsize>(strlen));
                }
            }
            else
                libmesh_error_msg("BinaryReader: unsupported payload encountered while skipping section.");
        }
        break;
    }
    default:
        libmesh_error_msg("BinaryReader: encountered unknown section label '" + std::string(BinaryLabelTraits::to_string(label)) + "'.");
    }
}
