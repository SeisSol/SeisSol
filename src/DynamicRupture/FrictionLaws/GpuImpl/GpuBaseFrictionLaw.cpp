#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "Parallel/MPI.h"
#include "utils/logger.h"
#include <device.h>
#include <sstream>

namespace seissol::dr::friction_law::gpu {
GpuBaseFrictionLaw::GpuBaseFrictionLaw(dr::DRParameters* drParameters)
    : FrictionSolver(drParameters) {}

GpuBaseFrictionLaw::~GpuBaseFrictionLaw() {
  free(faultStresses, queue);
  free(tractionResults, queue);
  free(stateVariableBuffer, queue);
  free(strengthBuffer, queue);
  free(devTimeWeights, queue);
  free(devDeltaT, queue);
  free(resampleMatrix, queue);
}

void GpuBaseFrictionLaw::setDeviceId(int currDeviceId) {
  std::ostringstream info;

  try {
#if defined(SYCL_PLATFORM_NVIDIA)
    device = sycl::make_device<sycl::backend::cuda>(currDeviceId);
#elif defined(SYCL_PLATFORM_AMD)
    device = sycl::make_device<sycl::backend::hip>(currDeviceId);
#else
    device = sycl::device(sycl::cpu_selector());
#endif
    device::DeviceInstance& nativeDevice = device::DeviceInstance::getInstance();
    nativeDevice.api->setDevice(currDeviceId);
  } catch (sycl::exception const& err) {
    info << "[SYCL] " << err.what() << "; ";
    device = sycl::device(sycl::cpu_selector());
  }

  info << "[SYCL] GPU device: " << std::boolalpha << device.is_gpu() << "; ";
  info << "[SYCL] Using Device: " << device.get_info<sycl::info::device::name>();

  const auto rank = seissol::MPI::mpi.rank();
  logInfo(rank) << info.str();

  sycl::property_list property{sycl::property::queue::in_order()};
  queue = sycl::queue(device, property);
}

void GpuBaseFrictionLaw::allocateAuxiliaryMemory() {
  faultStresses = static_cast<FaultStresses*>(
      sycl::malloc_shared(maxClusterSize * sizeof(FaultStresses), queue));

  tractionResults = static_cast<TractionResults*>(
      sycl::malloc_shared(maxClusterSize * sizeof(TractionResults), queue));

  {
    const size_t requiredNumBytes = misc::numPaddedPoints * maxClusterSize * sizeof(real);
    using StateVariableType = decltype(stateVariableBuffer);
    stateVariableBuffer =
        reinterpret_cast<StateVariableType>(sycl::malloc_shared(requiredNumBytes, queue));

    using StrengthBufferType = decltype(stateVariableBuffer);
    strengthBuffer =
        reinterpret_cast<StrengthBufferType>(sycl::malloc_shared(requiredNumBytes, queue));
  }

  {
    const size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(double);
    devTimeWeights = static_cast<double*>(sycl::malloc_shared(requiredNumBytes, queue));
  }

  {
    const size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(real);
    devDeltaT = static_cast<real*>(sycl::malloc_shared(requiredNumBytes, queue));
  }
}

void GpuBaseFrictionLaw::copyStaticDataToDevice() {
  constexpr auto dim0 = misc::dimSize<init::resample, 0>();
  constexpr auto dim1 = misc::dimSize<init::resample, 1>();
  const size_t requiredNumBytes = dim0 * dim1 * sizeof(real);

  resampleMatrix = static_cast<real*>(sycl::malloc_shared(requiredNumBytes, queue));
  queue.memcpy(resampleMatrix, &init::resample::Values[0], requiredNumBytes).wait();
}
} // namespace seissol::dr::friction_law::gpu
