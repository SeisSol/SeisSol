#include "Parallel/AcceleratorDevice.h"
#include "Parallel/MPI.h"
#include "utils/logger.h"
#include <string>
#include <sstream>

namespace seissol {
void AcceleratorDevice::bindSyclDevice(int deviceId) {
  try {
#ifdef __DPCPP_COMPILER
    syclDevice = sycl::device(sycl::gpu_selector_v);
#elif defined(SYCL_PLATFORM_NVIDIA)
    syclDevice = sycl::make_device<sycl::backend::cuda>(deviceId);
#elif defined(SYCL_PLATFORM_AMD)
    syclDevice = sycl::make_device<sycl::backend::hip>(deviceId);
#else
    syclDevice = sycl::device(sycl::cpu_selector());
#endif // __DPCPP_COMPILER
  } catch (sycl::exception const& err) {
    logWarning() << "SYCL error: " << err.what();
#ifdef __DPCPP_COMPILER
    syclDevice = sycl::device(sycl::cpu_selector_v);
#else
    syclDevice = sycl::device(sycl::cpu_selector());
#endif
  }

  std::ostringstream info;
  info << "SYCL device (GPU: " << std::boolalpha << syclDevice.is_gpu() << "): " << syclDevice.get_info<sycl::info::device::name>();

  const auto rank = seissol::MPI::mpi.rank();
  logInfo(rank) << info.str();

  sycl::property_list property{sycl::property::queue::in_order()};
  syclDefaultQueue = sycl::queue(syclDevice, property);
}

void AcceleratorDevice::bindNativeDevice(int deviceId) {
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  const auto rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Device API:" << device.api->getApiName().c_str();
  logInfo(rank) << "Device (rank=0):" << device.api->getDeviceName(deviceId).c_str();
  device.api->setDevice(deviceId);
}
} // namespace seissol
