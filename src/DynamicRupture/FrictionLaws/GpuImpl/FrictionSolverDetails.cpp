// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverDetails.h"

#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "GeneratedCode/init.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Typedefs.h"
#include "Memory/MemoryAllocator.h"

namespace seissol::dr::friction_law::gpu {
FrictionSolverDetails::FrictionSolverDetails(
    seissol::initializer::parameters::DRParameters* drParameters)
    : FrictionSolverInterface(drParameters) {}

FrictionSolverDetails::~FrictionSolverDetails() = default;

void FrictionSolverDetails::allocateAuxiliaryMemory(GlobalData* globalData) {
  {
    data_ = seissol::memory::allocTyped<FrictionLawData>(1, 1, memory::Memkind::DeviceGlobalMemory);
  }

  resampleMatrix_ = globalData->resampleMatrix;
  devSpaceWeights_ = globalData->spaceWeights;
  devTpInverseFourierCoefficients_ = globalData->tpInverseFourierCoefficients;
  devHeatSource_ = globalData->heatSource;
  devTpGridPoints_ = globalData->tpGridPoints;
}
} // namespace seissol::dr::friction_law::gpu
