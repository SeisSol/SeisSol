#include "Stream.hpp"
namespace seissol::parallel::runtime {

std::mutex StreamRuntime::mutexCPU = std::mutex();

} // namespace seissol::parallel::runtime
