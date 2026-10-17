/**
 * @file pytec3d_equation_systems.cpp
 * @brief pyTEC3DEquationSystems: system registration, variable initialisation, and solve orchestration.
 */
#include "pytec3d_equation_systems.h"
#include "element_builder.h"
#include "pytec3d_system.h"
#include "pytec3d_model.h"
#include "element_evaluator.h"
#include "element_parameters.h"
#include "model_element_types.h"
#include "vtkio.h"
#include "logging.h"
#include "report_tables.h"
#include <libmesh/petsc_linear_solver.h>
#include <libmesh/solver_configuration.h>
#include "solver_configuration.h"

#include <cmath>
#include <limits>
#include <utility>
#include <array>
#include <map>
#include <set>
#include <vector>
#include <unordered_map>

#include <libmesh/elem.h>
#include <libmesh/node.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

pyTEC3DEquationSystems::pyTEC3DEquationSystems(libMesh::MeshBase &mesh, pyTEC3DModel &model)
    : libMesh::EquationSystems(mesh), model_(model)
{
}

pyTEC3DEquationSystems::~pyTEC3DEquationSystems() = default;

void pyTEC3DEquationSystems::load_model_data(const std::string &model_json, const SolverSettings &solver_settings)
{
    model_data_ = pytec3d::ModelData::parse(model_json);
    solver_settings_ = solver_settings;
    element_types_ = std::make_unique<ModelElementTypes>(model_data_);

    element_local_to_global_rotations_.clear();
    node_local_to_global_rotations_.clear();
}

void pyTEC3DEquationSystems::validate_model_data() const
{
    vtkUnstructuredGrid *grid = model_.get_active_vtk_grid();
    vtkCellData *cell_data = grid ? grid->GetCellData() : nullptr;
    vtkDataArray *types = cell_data ? cell_data->GetArray("TYPE") : nullptr;
    vtkDataArray *materials = cell_data ? cell_data->GetArray("MATERIAL") : nullptr;
    vtkDataArray *constants = cell_data ? cell_data->GetArray("CONSTANT") : nullptr;
    if (!types || !materials || !constants)
        pytec3d::logging::error_and_abort("The grid needs the cell arrays TYPE, MATERIAL and CONSTANT.");

    using Key = std::array<unsigned int, 3>;
    std::set<Key> checked;
    for (vtkIdType cell = 0; cell < grid->GetNumberOfCells(); ++cell)
    {
        const Key key{static_cast<unsigned int>(types->GetComponent(cell, 0)),
                      static_cast<unsigned int>(materials->GetComponent(cell, 0)),
                      static_cast<unsigned int>(constants->GetComponent(cell, 0))};
        if (!checked.insert(key).second)
            continue;

        const ModelElementTypes::Entry &entry = element_types_->get(static_cast<libMesh::subdomain_id_type>(key[0]));
        const pytec3d::ElementType &type = *entry.type;

        pytec3d::ElementParameterValues::Inputs inputs;
        inputs.material = model_data_.material(key[1]);
        inputs.constants = model_data_.constants(key[2]);
        inputs.material_id = key[1];
        inputs.constants_id = key[2];
        inputs.domain = entry.domain;
        inputs.cell_data = cell_data;
        inputs.element_id = cell;

        // Builds what the assembly will on the first element of the combination, without a solution: a
        // law on an element with a curve needs no linear parameter.
        pytec3d::ElementParameterValues values(type, inputs);
        const libMesh::Elem *elem = this->get_mesh().query_elem_ptr(static_cast<libMesh::dof_id_type>(cell));
        if (!elem)
            continue;
        pytec3d::ElementEvaluator::Inputs evaluator_inputs;
        evaluator_inputs.type = &type;
        evaluator_inputs.domain = entry.domain;
        evaluator_inputs.laws = model_data_.laws(key[0], key[1], key[2]) & type.laws();
        evaluator_inputs.elem = elem;
        evaluator_inputs.parameters = &values;
        evaluator_inputs.model_data = &model_data_;
        evaluator_inputs.material_id = key[1];
        evaluator_inputs.constants_id = key[2];
        evaluator_inputs.type_id = key[0];
        pytec3d::ElementEvaluator(evaluator_inputs).check_parameters();
    }
}

const pytec3d::ModelData &pyTEC3DEquationSystems::get_model_data() const
{
    return model_data_;
}

const SolverSettings &pyTEC3DEquationSystems::get_solver_settings() const
{
    return solver_settings_;
}

void pyTEC3DEquationSystems::build_rotation_caches(const CoordinateRotations &rotations)
{
    element_local_to_global_rotations_.clear();
    node_local_to_global_rotations_.clear();

    vtkUnstructuredGrid *vtk_grid = model_.get_active_vtk_grid();
    vtkCellData *cell_data = vtk_grid->GetCellData();
    vtkPointData *point_data = vtk_grid->GetPointData();
    vtkDataArray *coordinate_system_cell_array =
        cell_data ? vtkDataArray::SafeDownCast(cell_data->GetArray("COORDINATE_SYSTEM")) : nullptr;
    vtkDataArray *angles_array =
        point_data ? vtkDataArray::SafeDownCast(point_data->GetArray("node_angles")) : nullptr;

    const libMesh::RealTensorValue identity_rotation(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0);

    // The local elements of each coordinate system with their centres.
    std::map<unsigned int, std::vector<libMesh::dof_id_type>> elem_ids_by_cs;
    std::map<unsigned int, std::vector<libMesh::Point>> centers_by_cs;
    const libMesh::MeshBase &mesh = this->get_mesh();

    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        const auto elem_id = elem->id();
        unsigned int cs_id = 0;
        if (coordinate_system_cell_array)
            cs_id = static_cast<unsigned int>(coordinate_system_cell_array->GetComponent(static_cast<vtkIdType>(elem_id), 0));

        if (cs_id == 0)
        {
            element_local_to_global_rotations_[elem_id] = identity_rotation;
            continue;
        }

        elem_ids_by_cs[cs_id].push_back(elem_id);
        centers_by_cs[cs_id].push_back(elem->vertex_average());
    }

    for (const auto &[cs_id, centers] : centers_by_cs)
    {
        if (!rotations.at_points)
            pytec3d::logging::error_and_abort("No coordinate systems are available for COORDINATE_SYSTEM " + std::to_string(cs_id) + ".");
        const auto matrices = rotations.at_points(cs_id, centers);
        if (matrices.size() != centers.size())
            pytec3d::logging::error_and_abort("Coordinate system rotation cache size mismatch for id=" + std::to_string(cs_id));

        const auto &elem_ids = elem_ids_by_cs[cs_id];
        for (std::size_t i = 0; i < matrices.size(); ++i)
            element_local_to_global_rotations_[elem_ids[i]] = matrices[i];
    }

    if (!angles_array)
    {
        for (const auto *node : mesh.node_ptr_range())
            node_local_to_global_rotations_[node->id()] = identity_rotation;
        return;
    }

    const vtkIdType n_tuples = angles_array->GetNumberOfTuples();
    std::vector<std::array<libMesh::Real, 3>> angles(static_cast<std::size_t>(n_tuples));
    for (vtkIdType i = 0; i < n_tuples; ++i)
        angles[static_cast<std::size_t>(i)] = {angles_array->GetComponent(i, 0),
                                               angles_array->GetComponent(i, 1),
                                               angles_array->GetComponent(i, 2)};

    if (!rotations.from_angles)
        pytec3d::logging::error_and_abort("No rotations are available for the \"node_angles\" point array.");
    const auto node_rotations = rotations.from_angles(angles);
    if (node_rotations.size() != angles.size())
        pytec3d::logging::error_and_abort("Node rotation cache size mismatch for the \"node_angles\" point array.");

    for (const auto *node : mesh.node_ptr_range())
        node_local_to_global_rotations_[node->id()] = node_rotations[static_cast<std::size_t>(node->id())];
}

void pyTEC3DEquationSystems::prepare_time_step_solve(const TimeStepInput &time_step)
{
    const libMesh::Real requested_time = time_step.time;
    const libMesh::Real requested_theta = time_step.theta;
    const bool requested_transient_effect = time_step.transient_effects;
    const libMesh::Real requested_beta = time_step.beta;
    const libMesh::Real requested_gamma = time_step.gamma;

    if (!time_step_)
    {
        time_step_ = std::make_unique<TimeStep>();
        time_step_->time = requested_time;
        time_step_->theta = requested_theta;
        time_step_->beta = requested_beta;
        time_step_->gamma = requested_gamma;
        time_step_->transient_effects = requested_transient_effect;
        time_step_->time_step_idx = 0;
        time_step_->dt = requested_time;
        time_step_->updated = true;
        return;
    }

    const libMesh::Real current_time = time_step_->time;
    const libMesh::Real time_scale = std::max(std::abs(requested_time), std::abs(current_time));
    const libMesh::Real time_tol = std::numeric_limits<libMesh::Real>::epsilon() * time_scale;
    const libMesh::Real time_delta = requested_time - current_time;

    if (time_delta > time_tol)
    {
        old_time_step_ = std::move(time_step_);
        time_step_ = std::make_unique<TimeStep>();
        time_step_->time = requested_time;
        time_step_->theta = requested_theta;
        time_step_->beta = requested_beta;
        time_step_->gamma = requested_gamma;
        time_step_->transient_effects = requested_transient_effect;
        time_step_->time_step_idx = old_time_step_->time_step_idx + 1;
        time_step_->dt = requested_time - current_time;
        time_step_->updated = true;
        return;
    }

    if (std::abs(time_delta) <= time_tol)
    {
        const libMesh::Real current_theta = time_step_->theta;
        const bool current_transient_effect = time_step_->transient_effects;
        const libMesh::Real theta_tol = std::numeric_limits<libMesh::Real>::epsilon() * std::abs(current_theta);

        const auto differs = [](libMesh::Real a, libMesh::Real b)
        { return std::abs(a - b) > std::numeric_limits<libMesh::Real>::epsilon() * std::max(std::abs(a), std::abs(b)); };
        if (std::abs(requested_theta - current_theta) > theta_tol || differs(requested_beta, time_step_->beta) ||
            differs(requested_gamma, time_step_->gamma) || requested_transient_effect != current_transient_effect)
        {
            pytec3d::logging::error_and_abort("pyTEC3DEquationSystems::prepare_time_step_solve: theta/beta/gamma/transient_effects changed while time is unchanged.");
        }

        time_step_->time = requested_time;
        time_step_->theta = requested_theta;
        time_step_->transient_effects = requested_transient_effect;
        time_step_->updated = false;
        return;
    }

    pytec3d::logging::error_and_abort("pyTEC3DEquationSystems::prepare_time_step_solve: received time smaller than current time step.");
}

void pyTEC3DEquationSystems::solve_equation_systems(const TimeStepInput &time_step)
{
    auto &system = get_dof_solution_system();
    system.init_time_step(time_step);

    const TimeStep &step = *current_time_step();
    step_report_ = StepReport{};
    step_report_.time = step.time;
    step_report_.dt = step.dt;
    step_report_.index = step.time_step_idx;
    step_report_.transient = time_scheme().transient();
    step_report_.advanced = step.updated;
    step_report_.nonlinear = use_nonlinear_system_;
    step_report_.elements = get_mesh().n_active_elem();
    step_report_.total_dofs = get_system(pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM)).n_dofs();

    system.write_element_matrices();
    this->solve();
    if (use_nonlinear_system_)
        compute_law_ranges();
}

void pyTEC3DEquationSystems::compute_law_ranges()
{
    namespace reg = pytec3d::registry;
    const libMesh::System &fields =
        get_system(pyTEC3DSystemBase::system_name_from_type(pyTEC3DSystemBase::ELEMENT_SOLUTION));
    vtkCellData *cell_data = model_.get_active_vtk_grid()->GetCellData();
    vtkDataArray *materials = cell_data->GetArray("MATERIAL");
    vtkDataArray *constants = cell_data->GetArray("CONSTANT");
    if (!materials || !constants)
        return;

    // The element field variables of each law input ("B", or "B_x", "B_y", "B_z").
    std::array<std::vector<unsigned int>, reg::kNumLawId> input_vars;
    for (std::size_t law = 0; law < reg::kNumLawId; ++law)
    {
        const std::string input(reg::kLaws[law].input);
        for (unsigned int var = 0; var < fields.n_vars(); ++var)
        {
            const std::string &name = fields.variable_name(var);
            if (name == input || (name.size() == input.size() + 2 && name.compare(0, input.size(), input) == 0 &&
                                  name[input.size()] == '_' && std::string("xyz").find(name.back()) != std::string::npos))
                input_vars[law].push_back(var);
        }
    }

    // Largest element mean of |input| per (law, id), over the local elements.
    std::map<std::pair<std::size_t, unsigned int>, libMesh::Real> largest;
    std::vector<libMesh::dof_id_type> dofs;
    for (const libMesh::Elem *elem : get_mesh().active_local_element_ptr_range())
    {
        const auto cell = static_cast<vtkIdType>(elem->id());
        const std::array<unsigned int, 3> ids{static_cast<unsigned int>(materials->GetComponent(cell, 0)),
                                              static_cast<unsigned int>(constants->GetComponent(cell, 0)),
                                              static_cast<unsigned int>(elem->subdomain_id())};
        const reg::LawMask mask = model_data_.laws(ids[2], ids[0], ids[1]) &
                                  element_types_->get(elem->subdomain_id()).type->laws();
        if (!mask)
            continue;
        for (std::size_t law = 0; law < reg::kNumLawId; ++law)
        {
            if (!(mask & (reg::LawMask(1) << law)) || input_vars[law].empty())
                continue;
            libMesh::Real squares = 0.0;
            for (const unsigned int var : input_vars[law])
            {
                fields.get_dof_map().dof_indices(elem, dofs, var);
                for (const auto dof : dofs)
                    squares += std::norm(fields.current_solution(dof));
            }
            const unsigned int id = ids[static_cast<std::size_t>(reg::kLaws[law].source)];
            libMesh::Real &value = largest[{law, id}];
            value = std::max(value, std::sqrt(squares));
        }
    }

    // Every (law, id) selected in the model, in the same order on every rank.
    std::vector<LawRange> &ranges = step_report_.law_ranges;
    ranges.clear();
    static const std::array<const char *, 3> source_names = {"MATERIAL", "CONSTANT", "TYPE"};
    for (std::size_t source = 0; source < 3; ++source)
        for (const auto &[id, mask] : model_data_.law_selection()[source])
            for (std::size_t law = 0; law < reg::kNumLawId; ++law)
            {
                if (!(mask & (reg::LawMask(1) << law)) || input_vars[law].empty())
                    continue;
                LawRange range;
                range.law = std::string(reg::kLawIdNames[law]);
                range.input = std::string(reg::kLaws[law].input);
                range.id = id;
                range.source = source_names[source];
                const auto it = largest.find({law, id});
                range.max_input = it == largest.end() ? 0.0 : it->second;
                this->comm().max(range.max_input);
                if (reg::kLaws[law].curve != reg::CurveId::COUNT)
                    if (const pytec3d::CurveData *curve = model_data_.curve(reg::kLaws[law].curve, id))
                        range.curve_max = curve->x().back();
                range.extrapolated = range.curve_max > 0.0 && range.max_input > range.curve_max;
                if (range.extrapolated)
                    pytec3d::logging::warn("The {} law of {} {}: |{}| reaches {:.4g} (element mean), beyond the last "
                                           "point of its curve ({:.4g}): the law is extrapolated.",
                                           range.law, range.source, range.id, range.input, range.max_input,
                                           range.curve_max);
                ranges.push_back(std::move(range));
            }
}

void pyTEC3DEquationSystems::configure_systems()
{
    use_nonlinear_system_ = this->has_nonlinear_elements();

    for (int i = 0; i < pyTEC3DSystemBase::NUM_SYSTEMS; ++i)
    {
        pyTEC3DSystemBase::SystemType type = static_cast<pyTEC3DSystemBase::SystemType>(i);
        const std::string &sys_name = pyTEC3DSystemBase::system_name_from_type(type);

        if (type == pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM)
        {
            if (use_nonlinear_system_)
                this->add_system<pyTEC3DNonLinearSystem>(sys_name);
            else
                this->add_system<pyTEC3DLinearSystem>(sys_name);
        }
        else if (type == pyTEC3DSystemBase::NODAL_SOLUTION)
            this->add_system<pyTEC3DNodalSolutionSystem>(sys_name);
        else if (type == pyTEC3DSystemBase::ELEMENT_SOLUTION)
            this->add_system<pyTEC3DElementSolutionSystem>(sys_name);
    }
}

void pyTEC3DEquationSystems::post_initialize_systems()
{
    for (int i = 0; i < pyTEC3DSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<pyTEC3DSystemBase::SystemType>(i);
        const auto &info = pyTEC3DSystemBase::system_name_from_type(type);

        if (type == pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM)
        {
            auto *system = dynamic_cast<pyTEC3DSystemBase *>(&get_dof_solution_system());
            if (!system)
                pytec3d::logging::error_and_abort("DOF solution system is not an pyTEC3DSystemBase.");
            system->post_init();
        }
        else if (type == pyTEC3DSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DNodalSolutionSystem>(info);
            system.post_init();
        }
        else if (type == pyTEC3DSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DElementSolutionSystem>(info);
            system.post_init();
        }
    }
}

void pyTEC3DEquationSystems::initialize_systems()
{
    pytec3d::logging::info("Setting systems variables");

    for (int i = 0; i < pyTEC3DSystemBase::NUM_SYSTEMS; ++i)
    {
        const auto type = static_cast<pyTEC3DSystemBase::SystemType>(i);
        const auto &info = pyTEC3DSystemBase::system_name_from_type(type);

        if (type == pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM)
        {
            auto *system = dynamic_cast<pyTEC3DSystemBase *>(&get_dof_solution_system());
            if (!system)
                pytec3d::logging::error_and_abort("DOF solution system is not an pyTEC3DSystemBase.");
            system->init_variables();
        }
        else if (type == pyTEC3DSystemBase::NODAL_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DNodalSolutionSystem>(info);
            system.init_variables();
        }
        else if (type == pyTEC3DSystemBase::ELEMENT_SOLUTION)
        {
            auto &system = this->get_system<pyTEC3DElementSolutionSystem>(info);
            system.init_variables();
        }
    }
}

void pyTEC3DEquationSystems::initialize_physics()
{
    this->init();
    post_initialize_systems();
    pytec3d::logging::table_root(pytec3d::report::format_equation_systems_summary_table(*this));
}

bool pyTEC3DEquationSystems::uses_nonlinear_system() const
{
    return use_nonlinear_system_;
}

pyTEC3DDOFSystemBase &pyTEC3DEquationSystems::get_dof_solution_system()
{
    const auto &name = pyTEC3DSystemBase::SYSTEM_INFO_MAP.at(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<pyTEC3DNonLinearSystem>(name);
    return this->get_system<pyTEC3DLinearSystem>(name);
}

const pyTEC3DDOFSystemBase &pyTEC3DEquationSystems::get_dof_solution_system() const
{
    const auto &name = pyTEC3DSystemBase::SYSTEM_INFO_MAP.at(pyTEC3DSystemBase::DOF_SOLUTION_SYSTEM);
    if (use_nonlinear_system_)
        return this->get_system<pyTEC3DNonLinearSystem>(name);
    return this->get_system<pyTEC3DLinearSystem>(name);
}

const TimeStep *pyTEC3DEquationSystems::current_time_step() const
{
    return time_step_.get();
}

const TimeStep *pyTEC3DEquationSystems::previous_time_step() const
{
    return old_time_step_.get();
}

bool pyTEC3DEquationSystems::has_nonlinear_elements() const
{
    if (!element_types_)
        return false;
    vtkUnstructuredGrid *grid = model_.get_active_vtk_grid();
    vtkCellData *cell_data = grid ? grid->GetCellData() : nullptr;
    vtkDataArray *types = cell_data ? cell_data->GetArray("TYPE") : nullptr;
    vtkDataArray *materials = cell_data ? cell_data->GetArray("MATERIAL") : nullptr;
    vtkDataArray *constants = cell_data ? cell_data->GetArray("CONSTANT") : nullptr;
    if (!types || !materials || !constants)
        return false;

    // An element with a law selected (by its TYPE, MATERIAL or CONSTANT row) that its element type has.
    using Key = std::array<unsigned int, 3>;
    std::set<Key> checked;
    for (vtkIdType cell = 0; cell < grid->GetNumberOfCells(); ++cell)
    {
        const Key key{static_cast<unsigned int>(types->GetComponent(cell, 0)),
                      static_cast<unsigned int>(materials->GetComponent(cell, 0)),
                      static_cast<unsigned int>(constants->GetComponent(cell, 0))};
        if (!checked.insert(key).second)
            continue;
        const auto &entry = element_types_->get(static_cast<libMesh::subdomain_id_type>(key[0]));
        if (model_data_.laws(key[0], key[1], key[2]) & entry.type->laws())
            return true;
    }
    return false;
}

unsigned int pyTEC3DEquationSystems::time_order() const
{
    return element_types_ ? element_types_->time_order() : 0;
}

TimeScheme pyTEC3DEquationSystems::time_scheme() const
{
    if (!time_step_)
        pytec3d::logging::error_and_abort("pyTEC3DEquationSystems::time_scheme: no time step.");
    return TimeScheme::make(*time_step_, this->time_order());
}

const ModelElementTypes &pyTEC3DEquationSystems::element_types() const
{
    if (!element_types_)
        pytec3d::logging::error_and_abort("Element types are not available before the model data is loaded.");
    return *element_types_;
}

const pyTEC3DEquationSystems::RotationMap &pyTEC3DEquationSystems::get_element_local_to_global_rotations() const
{
    return element_local_to_global_rotations_;
}

const pyTEC3DEquationSystems::RotationMap &pyTEC3DEquationSystems::get_node_local_to_global_rotations() const
{
    return node_local_to_global_rotations_;
}
