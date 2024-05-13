#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.hpp"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/precision.hpp"
#include "Parallel/AcceleratorDevice.h"
#include <cstddef>
#include <hipSYCL/sycl/usm.hpp>
#include <init.h>

namespace seissol::dr::friction_law::gpu {
FrictionSolverDetails::FrictionSolverDetails(
    seissol::initializer::parameters::DRParameters* drParameters)
    : FrictionSolverInterface(drParameters) {}

FrictionSolverDetails::~FrictionSolverDetails() {
  if (maxClusterSize == 0)
    return;

  free(faultStresses, queue);
  free(tractionResults, queue);
  free(stateVariableBuffer, queue);
  free(strengthBuffer, queue);
  free(devTimeWeights, queue);
  free(devSpaceWeights, queue);
  free(resampleMatrix, queue);
  queue.wait_and_throw();
}

void FrictionSolverDetails::initSyclQueue() {
  auto& instance = seissol::AcceleratorDevice::getInstance();
  queue = instance.getSyclDefaultQueue();
}

void FrictionSolverDetails::allocateAuxiliaryMemory() {
  if (maxClusterSize == 0)
    return;

  faultStresses = static_cast<FaultStresses*>(
      sycl::malloc_device(maxClusterSize * sizeof(FaultStresses), queue));

  tractionResults = static_cast<TractionResults*>(
      sycl::malloc_device(maxClusterSize * sizeof(TractionResults), queue));

  {
    const size_t requiredNumBytes = misc::numPaddedPoints * maxClusterSize * sizeof(real);
    using StateVariableType = decltype(stateVariableBuffer);
    stateVariableBuffer =
        reinterpret_cast<StateVariableType>(sycl::malloc_device(requiredNumBytes, queue));

    using StrengthBufferType = decltype(stateVariableBuffer);
    strengthBuffer =
        reinterpret_cast<StrengthBufferType>(sycl::malloc_device(requiredNumBytes, queue));
  }

  {
    const size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(double);
    devTimeWeights = static_cast<double*>(sycl::malloc_device(requiredNumBytes, queue));
  }

  {
    const size_t requiredNumBytes = misc::numPaddedPoints * sizeof(real);
    devSpaceWeights = static_cast<real*>(sycl::malloc_device(requiredNumBytes, queue));
  }
}

void FrictionSolverDetails::copyStaticDataToDevice() {
  if (maxClusterSize == 0)
    return;

  {
    constexpr auto dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto dim1 = misc::dimSize<init::resample, 1>();
    const size_t requiredNumBytes = dim0 * dim1 * sizeof(real);

    resampleMatrix = static_cast<real*>(sycl::malloc_device(requiredNumBytes, queue));
    queue.memcpy(resampleMatrix, &init::resample::Values[0], requiredNumBytes);
  }

  {
    const size_t requiredNumBytes = misc::numPaddedPoints * sizeof(real);
    queue.memcpy(devSpaceWeights, &spaceWeights[0], requiredNumBytes);
  }

  queue.wait_and_throw();
}
} // namespace seissol::dr::friction_law::gpu
