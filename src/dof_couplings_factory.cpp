#include "dof_couplings_factory.h"

DofCouplingsObject::DofCouplingsObject(std::ifstream &file)
    : BinaryDataFEDataObject(file)
{
    parse_couplings_from_map();
}

void DofCouplingsObject::read_data(std::ifstream &file)
{
    BinaryDataFEDataObject::read_data(file);
    parse_couplings_from_map();
}

void DofCouplingsObject::parse_couplings_from_map()
{
    couplings.clear();

    if (empty() || this->get_data_map()->columns.empty())
        return;

    if (!has_column("cp_id") || !has_column("variable"))
        return;

    const auto &cp_id_array = get_column("cp_id");
    const auto &variable_array = get_column("variable");

    const auto &cp_ids = cp_id_array.view<int>();
    const auto &variable_values = variable_array.view<std::string>();

    const std::size_t count = entry_count();
    if (cp_ids.size() < count || variable_values.size() < count)
        libmesh_error_msg("Dof couplings: scalar column size mismatch.");

    couplings.reserve(count);

    for (std::size_t i = 0; i < count; ++i)
    {
        DofCouplingEntry entry;
        entry.cp_id = cp_ids[i];
        entry.variable = variable_values[i];

        couplings.push_back(std::move(entry));
    }
}

const std::vector<DofCouplingEntry> &DofCouplingsObject::get_couplings() const
{
    return couplings;
}

std::vector<DofCouplingEntry> &DofCouplingsObject::get_couplings()
{
    return couplings;
}

void DofCouplingsObject::print_data_info() const
{
    if (couplings.empty())
    {
        libMesh::out << "  DOF couplings: none configured." << std::endl;
        return;
    }

    libMesh::out << "  DOF couplings: " << couplings.size() << std::endl;
    for (const auto &entry : couplings)
    {
        libMesh::out << "    cp_id=" << entry.cp_id
                     << " variable=" << entry.variable
                     << std::endl;
    }
}
