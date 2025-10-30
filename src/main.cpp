#include "apec_m.h"
#include "apec_m_cli.h"

#include <iostream>
#include <memory>
#include <string>

#include <libmesh/libmesh.h>

namespace
{
void print_usage()
{
    std::cout << "Usage: apec_m [-i | -b <batch_file>]" << std::endl;
    std::cout << "  -i                 Run in interactive mode" << std::endl;
    std::cout << "  -b <batch_file>    Run commands from batch file" << std::endl;
}
}

int main(int argc, char **argv)
{
    libMesh::LibMeshInit init(argc, argv);

    bool interactive_mode = false;
    bool batch_mode = false;
    std::string batch_file;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg(argv[i]);
        if (arg == "-i")
            interactive_mode = true;
        else if (arg == "-b")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "Error: -b requires a batch file argument." << std::endl;
                print_usage();
                return 1;
            }
            batch_mode = true;
            batch_file = argv[++i];
        }
        else if (arg == "-h" || arg == "--help")
        {
            print_usage();
            return 0;
        }
        else
        {
            std::cerr << "Unknown argument: " << arg << std::endl;
            print_usage();
            return 1;
        }
    }

    if (interactive_mode && batch_mode)
    {
        std::cerr << "Error: interactive and batch modes are mutually exclusive." << std::endl;
        print_usage();
        return 1;
    }

    if (!interactive_mode && !batch_mode)
    {
        std::cerr << "Error: no mode selected." << std::endl;
        print_usage();
        return 1;
    }

    APEC_M program(init.comm());
    program.set_cli(std::make_unique<BasicCommandLineInterface>());

    if (interactive_mode)
        program.run_cli();
    else if (batch_mode)
        program.run_batch_file(batch_file);

    return 0;
}
