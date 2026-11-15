/**
 * @file pytec3d_equation_systems.cpp
 * @brief pyTEC3DEquationSystems: variables, systems, time step and solve orchestration, outputs.
 */
#include "pytec3d_equation_systems.h"

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>

#include "dof_system.h"
#include "element_evaluator.h"
#include "element_parameters.h"
#include "element_type.h"
#include "logging.h"
#include "model_element_types.h"
#include "pytec3d_model.h"
#include "registries.h"
#include "report_tables.h"

namespace
{
    [[noreturn]] void fail(const std::string &message)
    {
        pytec3d::logging::error_and_abort(message);
    }

    bool lagrange_vector(pytec3d::fem::FEFamily family)
    {
        return family == pytec3d::fem::FEFamily::LAGRANGE_VEC || family == pytec3d::fem::FEFamily::L2_LAGRANGE_VEC;
    }
} // namespace

const std::string &pyTEC3DEquationSystems::system_name(SystemType type)
{
    static const std::array<std::string, NUM_SYSTEMS> names = {"DOF Solution", "Nodal Solution", "Element Solution"};
    return names.at(static_cast<std::size_t>(type));
}

pyTEC3DEquationSystems::pyTEC3DEquationSystems(pytec3d::fem::FEMesh &mesh, pyTEC3DModel &model) : mesh_(mesh), model_(model)
{
}

pyTEC3DEquationSystems::~pyTEC3DEquationSystems() = default;

void pyTEC3DEquationSystems::load_model_data(const nlohmann::json &model_data)
{
    model_data_ = pytec3d::ModelData::from_json(model_data);
    element_types_ = std::make_unique<ModelElementTypes>(model_data_);
}

const ModelElementTypes &pyTEC3DEquationSystems::element_types() const
{
    if (!element_types_)
        fail("Element types are not available before the model data is loaded.");
    return *element_types_;
}

DofSystem &pyTEC3DEquationSystems::dof_system()
{
    return *dof_system_;
}

const DofSystem &pyTEC3DEquationSystems::dof_system() const
{
    return *dof_system_;
}

const pytec3d::fem::SystemLayout &pyTEC3DEquationSystems::layout(SystemType type) const
{
    switch (type)
    {
    case DOF_SOLUTION_SYSTEM: return *dof_space_;
    case NODAL_SOLUTION: return *nodal_layout_;
    default: return *element_layout_;
    }
}

void pyTEC3DEquationSystems::validate_model_data() const
{
    const pytec3d::ModelArrays &arrays = model_.arrays();
    using Key = std::array<unsigned int, 3>;
    std::set<Key> checked;
    // No finite element is evaluated for the parameters: an element without quadrature is enough.
    const pytec3d::MeshElement element{&mesh_, nullptr, nullptr, nullptr, -1};
    for (std::size_t cell = 0; cell < arrays.type.n_tuples(); ++cell)
    {
        const Key key{static_cast<unsigned int>(arrays.type(cell)), static_cast<unsigned int>(arrays.material(cell)),
                      static_cast<unsigned int>(arrays.constant(cell))};
        if (!checked.insert(key).second)
            continue;
        const ModelElementTypes::Entry &entry = element_types_->get(static_cast<pytec3d::la::subdomain_id_type>(key[0]));
        const pytec3d::ElementType &type = *entry.type;
        pytec3d::ElementParameterValues::Inputs inputs;
        inputs.material = model_data_.material(key[1]);
        inputs.constants = model_data_.constants(key[2]);
        inputs.material_id = key[1];
        inputs.constants_id = key[2];
        inputs.domain = entry.domain;
        inputs.arrays = &arrays;
        inputs.element_id = cell;
        // What the assembly will build on the first element of the combination, without a solution: a law on an
        // element with a curve needs no linear parameter.
        pytec3d::ElementParameterValues values(type, inputs);
        pytec3d::ElementEvaluator::Inputs evaluator_inputs;
        evaluator_inputs.type = &type;
        evaluator_inputs.domain = entry.domain;
        evaluator_inputs.laws = model_data_.laws(key[0], key[1], key[2]) & type.laws();
        evaluator_inputs.elem = &element;
        evaluator_inputs.parameters = &values;
        evaluator_inputs.model_data = &model_data_;
        evaluator_inputs.material_id = key[1];
        evaluator_inputs.constants_id = key[2];
        evaluator_inputs.type_id = key[0];
        pytec3d::ElementEvaluator(evaluator_inputs).check_parameters();
    }
}

std::vector<pytec3d::fem::Variable> pyTEC3DEquationSystems::system_variables(SystemType system) const
{
    const pytec3d::FieldKind kind = system == DOF_SOLUTION_SYSTEM ? pytec3d::FieldKind::Dof
                                    : system == NODAL_SOLUTION    ? pytec3d::FieldKind::Nodal
                                                                  : pytec3d::FieldKind::Element;
    // The TYPE ids of the cells of the mesh.
    std::set<int> mesh_types;
    const pytec3d::ModelArrays &arrays = model_.arrays();
    for (std::size_t cell = 0; cell < arrays.type.n_tuples(); ++cell)
        mesh_types.insert(arrays.type(cell));

    // One variable per label, in the order of the registry, active on the TYPE ids whose element types have the field.
    std::map<pytec3d::registry::VariableId, pytec3d::fem::Variable> variables;
    std::vector<int> ids;
    for (const auto &[id, entry] : element_types().entries())
        ids.push_back(id);
    std::sort(ids.begin(), ids.end());
    for (const int id : ids)
    {
        if (!mesh_types.count(id))
            continue;
        const ModelElementTypes::Entry &entry = element_types().get(id);
        const pytec3d::ElementType &type = *entry.type;
        for (const pytec3d::FieldSpec &field : type.spec().fields)
        {
            if (field.kind != kind || !pytec3d::ElementType::active(field, entry.domain))
                continue;
            const pytec3d::SystemVariable sv = type.variable(field);
            pytec3d::fem::Variable variable;
            variable.label = std::string(pytec3d::registry::info(sv.id).label);
            variable.id = sv.id;
            variable.family = sv.family;
            variable.order = sv.order;
            variable.n_components = field.n_components;
            for (unsigned int c = 0; c < 3; ++c)
                variable.components[c] = field.components[c];
            variable.vdim = lagrange_vector(sv.family) ? static_cast<int>(field.n_components) : 1;
            auto [it, inserted] = variables.emplace(sv.id, variable);
            if (!inserted && (it->second.family != variable.family || it->second.order != variable.order ||
                              it->second.n_components != variable.n_components))
                fail("The variable '" + variable.label + "' has different finite elements in the element types of the "
                     "model (TYPE " + std::to_string(id) + ").");
            it->second.subdomains.insert(id);
        }
    }
    std::vector<pytec3d::fem::Variable> out;
    for (auto &[id, variable] : variables)
        out.push_back(std::move(variable));
    return out;
}

void pyTEC3DEquationSystems::initialize_systems()
{
    use_nonlinear_system_ = has_nonlinear_elements();
    catalog_ = std::make_unique<pytec3d::fem::FECatalog>(mesh_.dimension(), mesh_.serendipity(), mesh_.lagrange_order());

    // The quadrature of every local element at the order of its element type, computed in batches.
    std::vector<int> orders(static_cast<std::size_t>(mesh_.n_local_elements()));
    for (int e = 0; e < mesh_.n_local_elements(); ++e)
        orders[static_cast<std::size_t>(e)] = element_types().get(mesh_.subdomain(e)).type->quadrature_order();
    {
        pytec3d::logging::ScopedTimer timer("Quadrature data computed");
        quadrature_ = std::make_unique<pytec3d::fem::QuadratureData>(mesh_, orders);
    }
    {
        pytec3d::logging::ScopedTimer timer("Degrees of freedom distributed");
        dof_space_ = std::make_unique<pytec3d::fem::DofSpace>(mesh_, *catalog_, system_variables(DOF_SOLUTION_SYSTEM));
    }
    nodal_layout_ = std::make_unique<pytec3d::fem::SystemLayout>(mesh_, *catalog_, system_variables(NODAL_SOLUTION));
    element_layout_ = std::make_unique<pytec3d::fem::SystemLayout>(mesh_, *catalog_, system_variables(ELEMENT_SOLUTION));
    dof_system_ = std::make_unique<DofSystem>(*this, use_nonlinear_system_);
    pytec3d::logging::info("DOF system: {}.", use_nonlinear_system_ ? "nonlinear (cells with nonlinear laws), Newton" : "linear");
    pytec3d::logging::table_root(pytec3d::report::format_equation_systems_summary_table(*this));
}

void pyTEC3DEquationSystems::build_rotation_caches(const pytec3d::ConstRealArray &element_rotations,
                                                   const pytec3d::ConstRealArray &node_rotations)
{
    const auto rotation = [](const pytec3d::ConstRealArray &rotations, std::size_t i)
    {
        return pytec3d::la::RealTensorValue(rotations(i, 0), rotations(i, 1), rotations(i, 2), rotations(i, 3), rotations(i, 4),
                                            rotations(i, 5), rotations(i, 6), rotations(i, 7), rotations(i, 8));
    };
    const pytec3d::la::RealTensorValue identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
    element_rotations_.assign(static_cast<std::size_t>(mesh_.n_local_elements()), identity);
    if (!element_rotations.empty())
        for (int e = 0; e < mesh_.n_local_elements(); ++e)
            element_rotations_[static_cast<std::size_t>(e)] = rotation(element_rotations, static_cast<std::size_t>(mesh_.cell(e)));
    node_rotations_.clear();
    if (!node_rotations.empty())
    {
        node_rotations_.resize(node_rotations.n_tuples());
        for (std::size_t p = 0; p < node_rotations.n_tuples(); ++p)
            node_rotations_[p] = rotation(node_rotations, p);
    }
}

void pyTEC3DEquationSystems::prepare_time_step_solve(const TimeStepInput &time_step)
{
    const auto make = [&](unsigned int index, double dt)
    {
        auto step = std::make_unique<TimeStep>();
        step->time = time_step.time;
        step->theta = time_step.theta;
        step->beta = time_step.beta;
        step->gamma = time_step.gamma;
        step->transient_effects = time_step.transient_effects;
        step->time_step_idx = index;
        step->dt = dt;
        step->updated = true;
        return step;
    };
    if (!time_step_)
    {
        time_step_ = make(0, time_step.time);
        return;
    }
    const double current = time_step_->time;
    const double tol = std::numeric_limits<double>::epsilon() * std::max(std::abs(time_step.time), std::abs(current));
    const double delta = time_step.time - current;
    if (delta > tol)
    {
        old_time_step_ = std::move(time_step_);
        time_step_ = make(old_time_step_->time_step_idx + 1, delta);
        return;
    }
    if (std::abs(delta) <= tol)
    {
        const auto differs = [](double a, double b)
        { return std::abs(a - b) > std::numeric_limits<double>::epsilon() * std::max(std::abs(a), std::abs(b)); };
        if (differs(time_step.theta, time_step_->theta) || differs(time_step.beta, time_step_->beta) ||
            differs(time_step.gamma, time_step_->gamma) || time_step.transient_effects != time_step_->transient_effects)
            fail("pyTEC3DEquationSystems::prepare_time_step_solve: theta/beta/gamma/transient_effects changed while time "
                 "is unchanged.");
        time_step_->updated = false;
        return;
    }
    fail("pyTEC3DEquationSystems::prepare_time_step_solve: received time smaller than current time step.");
}

void pyTEC3DEquationSystems::solve_equation_systems(const TimeStepInput &time_step, const SolverSettings &settings)
{
    solver_settings_ = settings;
    dof_system_->apply_solver_settings(settings);
    dof_system_->init_time_step(time_step);
    const TimeStep &step = *current_time_step();
    step_report_ = StepReport{};
    step_report_.time = step.time;
    step_report_.dt = step.dt;
    step_report_.index = step.time_step_idx;
    step_report_.transient = time_scheme().transient();
    step_report_.advanced = step.updated;
    step_report_.nonlinear = use_nonlinear_system_;
    step_report_.elements = mesh_.n_cells();
    step_report_.total_dofs = static_cast<std::size_t>(dof_space_->n_global());
    dof_system_->assemble_operators();
    dof_system_->solve();
    if (use_nonlinear_system_)
        compute_law_ranges();
}

void pyTEC3DEquationSystems::compute_law_ranges()
{
    namespace reg = pytec3d::registry;
    const pytec3d::fem::SystemLayout &fields = *element_layout_;
    const pytec3d::ModelArrays &arrays = model_.arrays();
    // The element field of each law input.
    std::array<int, reg::kNumLawId> input_var;
    input_var.fill(-1);
    for (std::size_t law = 0; law < reg::kNumLawId; ++law)
        for (unsigned int v = 0; v < fields.n_vars(); ++v)
            if (reg::info(fields.variable(v).id).api_label == reg::kLaws[law].input)
                input_var[law] = static_cast<int>(v);

    std::vector<std::vector<double>> sums(fields.n_vars()), counts(fields.n_vars());
    for (unsigned int v = 0; v < fields.n_vars(); ++v)
    {
        sums[v].assign(mesh_.n_cells() * fields.variable(v).n_components, 0.0);
        counts[v].assign(mesh_.n_cells(), 0.0);
    }
    dof_system_->accumulate_fields(pytec3d::FieldKind::Element, sums, counts);

    // Largest element mean of |input| per (law, id), over the local elements.
    std::map<std::pair<std::size_t, unsigned int>, double> largest;
    for (int e = 0; e < mesh_.n_local_elements(); ++e)
    {
        const auto cell = static_cast<std::size_t>(mesh_.cell(e));
        const std::array<unsigned int, 3> ids{static_cast<unsigned int>(arrays.material(cell)),
                                              static_cast<unsigned int>(arrays.constant(cell)),
                                              static_cast<unsigned int>(mesh_.subdomain(e))};
        const reg::LawMask mask = model_data_.laws(ids[2], ids[0], ids[1]) & element_types_->get(mesh_.subdomain(e)).type->laws();
        for (std::size_t law = 0; law < reg::kNumLawId && mask; ++law)
        {
            if (!(mask & (reg::LawMask(1) << law)) || input_var[law] < 0)
                continue;
            const auto v = static_cast<std::size_t>(input_var[law]);
            const unsigned int n = fields.variable(static_cast<unsigned int>(v)).n_components;
            double squares = 0.0;
            for (unsigned int c = 0; c < n; ++c)
                squares += sums[v][cell * n + c] * sums[v][cell * n + c];
            double &value = largest[{law, ids[static_cast<std::size_t>(reg::kLaws[law].source)]}];
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
                if (!(mask & (reg::LawMask(1) << law)) || input_var[law] < 0)
                    continue;
                LawRange range;
                range.law = std::string(reg::kLawIdNames[law]);
                range.input = std::string(reg::kLaws[law].input);
                range.id = id;
                range.source = source_names[source];
                const auto it = largest.find({law, id});
                range.max_input = it == largest.end() ? 0.0 : it->second;
                MPI_Allreduce(MPI_IN_PLACE, &range.max_input, 1, MPI_DOUBLE, MPI_MAX, mesh_.comm());
                if (reg::kLaws[law].curve != reg::CurveId::COUNT)
                    if (const pytec3d::CurveData *curve = model_data_.curve(reg::kLaws[law].curve, id))
                        range.curve_max = curve->is_analytic() ? 0.0 : curve->x_max(); // analytic: no extrapolation
                range.extrapolated = range.curve_max > 0.0 && range.max_input > range.curve_max;
                if (range.extrapolated)
                    pytec3d::logging::warn("The {} law of {} {}: |{}| reaches {:.4g} (element mean), beyond the last point of "
                                           "its curve ({:.4g}): the law is extrapolated.",
                                           range.law, range.source, range.id, range.input, range.max_input, range.curve_max);
                ranges.push_back(std::move(range));
            }
}

std::vector<std::string> pyTEC3DEquationSystems::missing_outputs(SystemType system, const pytec3d::ModelArrays &arrays) const
{
    const pytec3d::fem::SystemLayout &l = layout(system);
    const bool points = system != ELEMENT_SOLUTION;
    std::vector<std::string> missing;
    std::set<std::string> seen;
    for (unsigned int v = 0; v < l.n_vars(); ++v)
    {
        const pytec3d::fem::Variable &var = l.variable(v);
        const std::string name(pytec3d::registry::info(var.id).api_label);
        const pytec3d::RealArray *array = points ? arrays.point_field(name) : arrays.cell_field(name);
        bool ok = array != nullptr;
        for (unsigned int c = 0; c < var.n_components && ok; ++c)
            ok = var.components[c] < array->n_components();
        if (!ok && seen.insert(name).second)
            missing.push_back(name + (points ? " (point)" : " (cell)"));
    }
    return missing;
}

void pyTEC3DEquationSystems::write_output(SystemType system, const pytec3d::ModelArrays &arrays) const
{
    // The values at the points (cells) of the elements, averaged over the elements that have the variable (as libMesh
    // averaged the nodal values of the elements); 0 where none has it.
    const pytec3d::fem::SystemLayout &l = layout(system);
    const bool points = system != ELEMENT_SOLUTION;
    const std::size_t n = points ? mesh_.n_points() : mesh_.n_cells();
    std::vector<std::vector<double>> sums(l.n_vars()), counts(l.n_vars());
    for (unsigned int v = 0; v < l.n_vars(); ++v)
    {
        sums[v].assign(n * l.variable(v).n_components, 0.0);
        counts[v].assign(n, 0.0);
    }
    if (system == DOF_SOLUTION_SYSTEM)
        dof_system_->accumulate_dof_fields(sums, counts);
    else
        dof_system_->accumulate_fields(system == NODAL_SOLUTION ? pytec3d::FieldKind::Nodal : pytec3d::FieldKind::Element,
                                       sums, counts);
    for (unsigned int v = 0; v < l.n_vars(); ++v)
    {
        MPI_Allreduce(MPI_IN_PLACE, sums[v].data(), static_cast<int>(sums[v].size()), MPI_DOUBLE, MPI_SUM, mesh_.comm());
        MPI_Allreduce(MPI_IN_PLACE, counts[v].data(), static_cast<int>(counts[v].size()), MPI_DOUBLE, MPI_SUM, mesh_.comm());
        const pytec3d::fem::Variable &var = l.variable(v);
        const std::string name(pytec3d::registry::info(var.id).api_label);
        const pytec3d::RealArray *array = points ? arrays.point_field(name) : arrays.cell_field(name);
        if (!array)
            continue;
        for (std::size_t i = 0; i < n; ++i)
            for (unsigned int c = 0; c < var.n_components; ++c)
                if (var.components[c] < array->n_components())
                    (*array)(i, var.components[c]) =
                        counts[v][i] > 0.0 ? sums[v][i * var.n_components + c] / counts[v][i] : 0.0;
    }
}

bool pyTEC3DEquationSystems::has_nonlinear_elements() const
{
    if (!element_types_)
        return false;
    const pytec3d::ModelArrays &arrays = model_.arrays();
    // An element with a law selected (by its TYPE, MATERIAL or CONSTANT row) that its element type has.
    using Key = std::array<unsigned int, 3>;
    std::set<Key> checked;
    for (std::size_t cell = 0; cell < arrays.type.n_tuples(); ++cell)
    {
        const Key key{static_cast<unsigned int>(arrays.type(cell)), static_cast<unsigned int>(arrays.material(cell)),
                      static_cast<unsigned int>(arrays.constant(cell))};
        if (!checked.insert(key).second)
            continue;
        const auto &entry = element_types_->get(static_cast<pytec3d::la::subdomain_id_type>(key[0]));
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
        fail("pyTEC3DEquationSystems::time_scheme: no time step.");
    return TimeScheme::make(*time_step_, time_order());
}
