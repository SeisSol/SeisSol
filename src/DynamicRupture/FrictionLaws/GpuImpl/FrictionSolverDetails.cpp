// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "Common/Constants.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Parallel/AcceleratorDevice.h"
#include <cstddef>
#include <init.h>

#ifndef __DPCPP_COMPILER
#include <hipSYCL/sycl/usm.hpp>
#endif

namespace seissol::dr::friction_law::gpu {
FrictionSolverDetails::FrictionSolverDetails(
    seissol::initializer::parameters::DRParameters* drParameters)
    : FrictionSolverInterface(drParameters) {}

FrictionSolverDetails::~FrictionSolverDetails() {
  if (maxClusterSize == 0) {
    return;
  }

  free(devTimeWeights, queue);
  free(devSpaceWeights, queue);
  free(resampleMatrix, queue);
  free(data, queue);
  queue.wait_and_throw();
}

void FrictionSolverDetails::initSyclQueue() {
  auto& instance = seissol::AcceleratorDevice::getInstance();
  queue = instance.getInorderSyclQueue();
}

void FrictionSolverDetails::allocateAuxiliaryMemory() {
  if (maxClusterSize == 0) {
    return;
  }

  {
    const size_t requiredNumBytes = ConvergenceOrder * sizeof(double);
    devTimeWeights = static_cast<double*>(sycl::malloc_device(requiredNumBytes, queue));
  }

  {
    const size_t requiredNumBytes = misc::NumPaddedPoints * sizeof(real);
    devSpaceWeights = static_cast<real*>(sycl::malloc_device(requiredNumBytes, queue));
  }

  {
    data = static_cast<FrictionLawData*>(sycl::malloc_device(sizeof(FrictionLawData), queue));
  }
}

void FrictionSolverDetails::copyStaticDataToDevice() {
  if (maxClusterSize == 0) {
    return;
  }

  {
    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    const size_t requiredNumBytes = Dim0 * Dim1 * sizeof(real);

    resampleMatrix = static_cast<real*>(sycl::malloc_device(requiredNumBytes, queue));
    queue.memcpy(resampleMatrix, &init::resample::Values[0], requiredNumBytes);
  }

  {
    const size_t requiredNumBytes = misc::NumPaddedPoints * sizeof(real);
    queue.memcpy(devSpaceWeights, &spaceWeights[0], requiredNumBytes);
  }

  queue.wait_and_throw();
}
} // namespace seissol::dr::friction_law::gpu
