#include "Parallel/AcceleratorDevice.h"
#include "Parallel/MPI.h"
#include "utils/logger.h"
#include <string>
#include <sstream>

namespace seissol {
void AcceleratorDevice::bindSyclDevice(int deviceId) {
  std::ostringstream info;

  try {
#ifdef __DPCPP_COMPILER
    syclDevice = sycl::device(sycl::gpu_selector_v);
#else
#if defined(SYCL_PLATFORM_NVIDIA)
    syclDevice = sycl::make_device<sycl::backend::cuda>(deviceId);
#elif defined(SYCL_PLATFORM_AMD)
    syclDevice = sycl::make_device<sycl::backend::hip>(deviceId);
#else
    syclDevice = sycl::device(sycl::cpu_selector());
#endif
#endif // __DPCPP_COMPILER
  } catch (sycl::exception const& err) {
    info << "[SYCL] " << err.what() << "; ";
#ifdef __DPCPP_COMPILER
    syclDevice = sycl::device(sycl::cpu_selector_v);
#else
    syclDevice = sycl::device(sycl::cpu_selector());
#endif
  }

  info << "[SYCL] GPU device: " << std::boolalpha << syclDevice.is_gpu() << "; ";
  info << "[SYCL] Using Device: " << syclDevice.get_info<sycl::info::device::name>();

  const auto rank = seissol::MPI::mpi.rank();
  logInfo(rank) << info.str();

  sycl::property_list property{sycl::property::queue::in_order()};
  syclDefaultQueue = sycl::queue(syclDevice, property);
}

void AcceleratorDevice::bindNativeDevice(int deviceId) {
  device::DeviceInstance& device = device::DeviceInstance::getInstance();

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp critical
    { device.api->setDevice(deviceId); }
  }
#else
  device.api->setDevice(deviceId);
#endif
}
} // namespace seissol
