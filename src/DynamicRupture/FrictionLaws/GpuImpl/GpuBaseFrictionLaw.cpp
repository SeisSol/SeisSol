#include "DynamicRupture/FrictionLaws/GpuImpl/GpuBaseFrictionLaw.h"
#include "utils/logger.h"
#include <omp.h>
#include <algorithm>
#include <sstream>

namespace seissol::dr::friction_law::gpu {
GpuBaseFrictionLaw::GpuBaseFrictionLaw(dr::DRParameters& drParameters)
    : FrictionSolver(drParameters) {
  checkOffloading();
}

GpuBaseFrictionLaw::~GpuBaseFrictionLaw() {
  omp_target_free(faultStresses, deviceId);
  omp_target_free(tractionResults, deviceId);
  omp_target_free(stateVariableBuffer, deviceId);
  omp_target_free(strengthBuffer, deviceId);
  omp_target_free(resampleMatrix, deviceId);
}

void GpuBaseFrictionLaw::allocateAuxiliaryMemory() {
  faultStresses = reinterpret_cast<FaultStresses*>(
      omp_target_alloc(maxClusterSize * sizeof(FaultStresses), deviceId));
  tractionResults = reinterpret_cast<TractionResults*>(
      omp_target_alloc(maxClusterSize * sizeof(TractionResults), deviceId));

  size_t requiredNumBytes = misc::numPaddedPoints * maxClusterSize * sizeof(real);
  stateVariableBuffer =
      reinterpret_cast<decltype(stateVariableBuffer)>(omp_target_alloc(requiredNumBytes, deviceId));
  strengthBuffer =
      reinterpret_cast<decltype(strengthBuffer)>(omp_target_alloc(requiredNumBytes, deviceId));
}

void GpuBaseFrictionLaw::copyStaticDataToDevice() {
  constexpr auto dim0 = misc::dimSize<init::resample, 0>();
  constexpr auto dim1 = misc::dimSize<init::resample, 1>();
  size_t requiredNumBytes = dim0 * dim1 * sizeof(real);

  resampleMatrix = reinterpret_cast<real*>(omp_target_alloc(requiredNumBytes, deviceId));
  auto hostId = omp_get_initial_device();

  omp_target_memcpy(
      resampleMatrix, &init::resample::Values[0], requiredNumBytes, 0, 0, deviceId, hostId);
}

void GpuBaseFrictionLaw::checkOffloading() {
  bool canOffload = false;
  #pragma omp target map(tofrom : canOffload)
  {
    if (!omp_is_initial_device()) {
      canOffload = true;
    }
  }
  std::ostringstream info;
  info << "Device offloading: " << std::boolalpha << canOffload;
  logInfo() << info.str();
}
} // namespace seissol::dr::friction_law::gpu
