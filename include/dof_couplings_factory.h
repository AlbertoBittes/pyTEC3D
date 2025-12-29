#ifndef DOF_COUPLINGS_FACTORY_H
#define DOF_COUPLINGS_FACTORY_H

#include <string>
#include <utility>
#include <vector>

#include "fe_data_object.h"

struct DofCouplingEntry
{
    std::string array_name;
    int cp_id = 0;
    std::vector<std::pair<unsigned int, unsigned int>> pairs;
    std::vector<unsigned int> roots;
    std::vector<unsigned int> isolated;
    std::vector<std::string> coupled_variables;
};

class DofCouplingsObject : public BinaryDataFEDataObject
{
public:
    DofCouplingsObject() = default;
    explicit DofCouplingsObject(std::ifstream &file);

    void read_data(std::ifstream &file) override;
    void parse_couplings_from_map();
    const std::vector<DofCouplingEntry> &get_couplings() const;
    void print_data_info() const override;

private:
    std::vector<DofCouplingEntry> couplings;
};

#endif // DOF_COUPLINGS_FACTORY_H
