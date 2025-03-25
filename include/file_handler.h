#include <fstream>
#include <libmesh/libmesh.h> // Include the appropriate header for libMesh

using namespace libMesh;

class FileHandler
{
public:
    FileHandler(const std::string &filename)
    {
        inFile.open(filename, std::ios::binary);
        if (!inFile)
        {
            libmesh_error_msg("Cannot open file for reading");
        }
    }

    ~FileHandler()
    {
        if (inFile.is_open())
        {
            inFile.close();
        }
    }

    std::ifstream &getStream()
    {
        return inFile;
    }

    void closeFile()
    {
        if (inFile.is_open())
        {
            inFile.close();
        }
    }

private:
    std::ifstream inFile;
};