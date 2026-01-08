#ifndef DOF_COUPLINGS_FACTORY_H
#define DOF_COUPLINGS_FACTORY_H

#include <string>
#include <utility>
#include <vector>

#include "fe_data_object.h"

struct DofCouplingEntry
{
    int cp_id = 0;
    std::string variable;
};

class DofCouplingsObject : public BinaryDataFEDataObject
{
public:
    DofCouplingsObject() = default;
    explicit DofCouplingsObject(std::ifstream &file);

    void read_data(std::ifstream &file) override;
    void parse_couplings_from_map();
    const std::vector<DofCouplingEntry> &get_couplings() const;
    std::vector<DofCouplingEntry> &get_couplings();
    void print_data_info() const override;

private:
    std::vector<DofCouplingEntry> couplings;
};

#endif // DOF_COUPLINGS_FACTORY_H
