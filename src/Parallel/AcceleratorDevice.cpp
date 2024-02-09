#include "Parallel/AcceleratorDevice.h"
#include "Parallel/MPI.h"
#include "utils/logger.h"
#include <string>
#include <sstream>

namespace seissol {
void AcceleratorDevice::bindSyclDevice(int deviceId) {
}

void AcceleratorDevice::bindNativeDevice(int deviceId) {
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  {
    std::ostringstream info;
    info << "Device API: " << device.api->getApiName();
    infoMessages.push_back(info.str());
  }
  {
    std::ostringstream info;
    info << "Device (rank=0): " << device.api->getDeviceName(deviceId);
    infoMessages.push_back(info.str());
  }
  device.api->setDevice(deviceId);
}

void AcceleratorDevice::printInfo() {
  const auto rank = seissol::MPI::mpi.rank();
  for (const auto& warn : warnMessages) {
    logWarning(rank) << warn.c_str();
  }
  for (const auto& info : infoMessages) {
    logInfo(rank) << info.c_str();
  }
}
} // namespace seissol
