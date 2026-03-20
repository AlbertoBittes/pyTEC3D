#ifndef CLI_H
#define CLI_H

#include "apec_m.h"

#include <functional>
#include <string>

class BasicCommandLineInterface : public CommandLineInterface
{
public:
    enum class Command
    {
        READ_MESH,
        CLEAR_MESH,
        SET_ACTIVE_MESH,
        LOAD_DATA,
        SOLVE,
        FINISH,
        SET_MODE,
        HELP,
        EXIT,
        INVALID
    };

    BasicCommandLineInterface() = default;

    void run_interactive(APEC_M &program) override;
    void run_batch(APEC_M &program, const std::string &batch_file) override;

private:
    Command parse_command(const std::string &input) const;
    std::string extract_argument(const std::string &input) const;
    bool execute_command(APEC_M &program, const std::string &input, bool interactive);
};

#endif // CLI_H
