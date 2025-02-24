// SPDX-FileCopyrightText: 2022-2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "Common/Constants.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include <cstddef>
#include <init.h>

#include "Memory/MemoryAllocator.h"

#include "DynamicRupture/FrictionLaws/TPCommon.h"

namespace seissol::dr::friction_law::gpu {
FrictionSolverDetails::FrictionSolverDetails(
    seissol::initializer::parameters::DRParameters* drParameters)
    : FrictionSolverInterface(drParameters) {}

FrictionSolverDetails::~FrictionSolverDetails() = default;

void FrictionSolverDetails::allocateAuxiliaryMemory() {

  {
    devTimeWeights =
        seissol::memory::allocTyped<double>(ConvergenceOrder, 1, memory::DeviceGlobalMemory);
  }

  {
    data = seissol::memory::allocTyped<FrictionLawData>(1, 1, memory::DeviceGlobalMemory);
  }

  {
    constexpr auto Dim0 = misc::dimSize<init::resample, 0>();
    constexpr auto Dim1 = misc::dimSize<init::resample, 1>();
    const size_t requiredNumBytes = Dim0 * Dim1 * sizeof(real);

    resampleMatrix = seissol::memory::allocTyped<real>(Dim0 * Dim1, 1, memory::DeviceGlobalMemory);
    seissol::memory::memcopyTyped<real>(resampleMatrix,
                                        init::resample::Values,
                                        Dim0 * Dim1,
                                        memory::DeviceGlobalMemory,
                                        memory::Standard);
  }

  {
    devSpaceWeights =
        seissol::memory::allocTyped<real>(misc::NumPaddedPoints, 1, memory::DeviceGlobalMemory);
    seissol::memory::memcopyTyped<real>(devSpaceWeights,
                                        spaceWeights,
                                        misc::NumPaddedPoints,
                                        memory::DeviceGlobalMemory,
                                        memory::Standard);
  }

  {
    const auto data =
        seissol::dr::friction_law::tp::InverseFourierCoefficients<misc::NumTpGridPoints>();
    devTpInverseFourierCoefficients =
        seissol::memory::allocTyped<real>(misc::NumTpGridPoints, 1, memory::DeviceGlobalMemory);
    seissol::memory::memcopyTyped<real>(devTpInverseFourierCoefficients,
                                        data.data().data(),
                                        misc::NumTpGridPoints,
                                        memory::DeviceGlobalMemory,
                                        memory::Standard);
  }

  {
    const auto data = seissol::dr::friction_law::tp::GridPoints<misc::NumTpGridPoints>();
    devTpGridPoints =
        seissol::memory::allocTyped<real>(misc::NumTpGridPoints, 1, memory::DeviceGlobalMemory);
    seissol::memory::memcopyTyped<real>(devTpGridPoints,
                                        data.data().data(),
                                        misc::NumTpGridPoints,
                                        memory::DeviceGlobalMemory,
                                        memory::Standard);
  }

  {
    const auto data = seissol::dr::friction_law::tp::GaussianHeatSource<misc::NumTpGridPoints>();
    devHeatSource =
        seissol::memory::allocTyped<real>(misc::NumTpGridPoints, 1, memory::DeviceGlobalMemory);
    seissol::memory::memcopyTyped<real>(devHeatSource,
                                        data.data().data(),
                                        misc::NumTpGridPoints,
                                        memory::DeviceGlobalMemory,
                                        memory::Standard);
  }
}
} // namespace seissol::dr::friction_law::gpu
