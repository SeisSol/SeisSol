#include "Parallel/AcceleratorDevice.h"
#include "Parallel/MPI.h"
#include "utils/logger.h"
#include <string>
#include <sstream>

namespace seissol {

void AcceleratorDevice::bindNativeDevice(int deviceId) {
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->setDevice(deviceId);
}

} // namespace seissol
