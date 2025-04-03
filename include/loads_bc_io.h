#ifndef LOADS_BC_IO_H
#define LOADS_BC_IO_H

// Include necessary headers
#include <vector>
#include <unordered_map>
#include <libmesh/libmesh.h>

using namespace libMesh;

typedef std::unordered_map<std::string, Real> NodeBoundaryConditions;
typedef std::unordered_map<dof_id_type, NodeBoundaryConditions> BoundaryConditionsMap;
typedef std::unordered_map<dof_id_type, std::vector<Real>> ElementLoad;
typedef std::unordered_map<std::string, ElementLoad> LoadArraysMap;

void read_boundary_conditions(std::string &filename, BoundaryConditionsMap &node_conditions);

void read_loads_arrays(std::string &filename, LoadArraysMap &loads_arrays);

// int main() {
//     try {
//         std::string boundary_filename = "boundary_conditions.bin";
//         std::string load_filename = "load_variables.bin";

//         FileHandler boundary_fileHandler(boundary_filename);
//         BoundaryConditionsMap node_conditions;
//         read_boundary_conditions(boundary_fileHandler, node_conditions);

//         FileHandler load_fileHandler(load_filename);
//         LoadArraysMap load_variables;
//         read_loads_arrays(load_fileHandler, load_variables);

//         // Print the boundary conditions (for demonstration purposes)
//         for (const auto &node_pair : node_conditions) {
//             std::cout << "Node " << node_pair.first << ":\n";
//             for (const auto &condition_pair : node_pair.second) {
//                 std::cout << "  " << condition_pair.first << " = " << condition_pair.second << "\n";
//             }
//         }

//         // Print the load variables (for demonstration purposes)
//         for (const auto &variable_pair : load_variables) {
//             std::cout << "Load Variable: " << variable_pair.first << "\n";
//             for (const auto &load_pair : variable_pair.second) {
//                 std::cout << "  Element " << load_pair.first << ":\n";
//                 for (const auto &value : load_pair.second) {
//                     std::cout << "    " << value << "\n";
//                 }
//             }
//         }
//     } catch (const std::exception &e) {
//         std::cerr << e.what() << std::endl;
//     }

//     return 0;
// }
#endif // LOADS_BC_IO_H