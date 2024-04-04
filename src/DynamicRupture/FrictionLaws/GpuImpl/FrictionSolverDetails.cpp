#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "Parallel/AcceleratorDevice.h"
#include <device.h>
#include <omp.h>

#include "utils/logger.h"

namespace seissol::dr::friction_law::gpu {
FrictionSolverDetails::FrictionSolverDetails(dr::DRParameters* drParameters)
    : FrictionSolverInterface(drParameters) {}

FrictionSolverDetails::~FrictionSolverDetails() {
  delete[] queue;
  if (maxClusterSize == 0)
    return;

  omp_free(faultStresses);
  omp_free(tractionResults);
  omp_free(stateVariableBuffer);
  omp_free(strengthBuffer);
  omp_free(devTimeWeights);
  omp_free(devSpaceWeights);
  omp_free(resampleMatrix);
}

void FrictionSolverDetails::initSyclQueue() {
  
}

void FrictionSolverDetails::allocateAuxiliaryMemory() {
  if (maxClusterSize == 0)
    return;
  
  queue = new int[chunkcount];

  chunksize = (maxClusterSize + chunkcount - 1) / chunkcount;

  faultStresses = static_cast<FaultStresses*>(
      omp_aligned_alloc(ALIGNMENT, maxClusterSize * sizeof(FaultStresses)));

  tractionResults = static_cast<TractionResults*>(
      omp_aligned_alloc(ALIGNMENT, maxClusterSize * sizeof(TractionResults)));

  {
    const size_t requiredNumBytes = misc::numPaddedPoints * maxClusterSize * sizeof(real);
    using StateVariableType = decltype(stateVariableBuffer);
    stateVariableBuffer =
        reinterpret_cast<StateVariableType>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));

    using StrengthBufferType = decltype(stateVariableBuffer);
    strengthBuffer =
        reinterpret_cast<StrengthBufferType>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
  }

  {
    const size_t requiredNumBytes = CONVERGENCE_ORDER * sizeof(double);
    devTimeWeights = static_cast<double*>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
  }

  {
    const size_t requiredNumBytes = misc::numPaddedPoints * sizeof(real);
    devSpaceWeights = static_cast<real*>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
  }
}

void FrictionSolverDetails::copyStaticDataToDevice() {
  if (maxClusterSize == 0)
    return;

  {
    constexpr auto dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto dim1 = misc::dimSize<init::resample, 1>();
    const size_t requiredNumBytes = dim0 * dim1 * sizeof(real);

    resampleMatrix = static_cast<real*>(omp_aligned_alloc(ALIGNMENT, requiredNumBytes));
    memcpy(resampleMatrix, &init::resample::Values[0], requiredNumBytes);
  }

  {
    const size_t requiredNumBytes = misc::numPaddedPoints * sizeof(real);
    memcpy(devSpaceWeights, &spaceWeights[0], requiredNumBytes);
  }
}
} // namespace seissol::dr::friction_law::gpu
