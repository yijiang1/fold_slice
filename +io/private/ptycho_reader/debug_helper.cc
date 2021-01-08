#include <fstream>
#include <memory>

namespace debug {
    std::unique_ptr<std::ofstream> out; //!< pointer to debug log file
}
