// SPDX-FileCopyrightText: 2022 SeisSol Group
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
#include <GeneratedCode/init.h>
#include <cstddef>

#include "Memory/MemoryAllocator.h"

#include "DynamicRupture/FrictionLaws/TPCommon.h"

namespace seissol::dr::friction_law::gpu {
FrictionSolverDetails::FrictionSolverDetails(
    seissol::initializer::parameters::DRParameters* drParameters)
    : FrictionSolverInterface(drParameters) {}

FrictionSolverDetails::~FrictionSolverDetails() = default;

void FrictionSolverDetails::allocateAuxiliaryMemory(GlobalData* globalData) {
  {
    data = seissol::memory::allocTyped<FrictionLawData<Cfg>>(1, 1, memory::DeviceGlobalMemory);
  }

  resampleMatrix = globalData->resampleMatrix;
  devSpaceWeights = globalData->spaceWeights;
  devTpInverseFourierCoefficients = globalData->tpInverseFourierCoefficients;
  devHeatSource = globalData->heatSource;
  devTpGridPoints = globalData->tpGridPoints;
}
} // namespace seissol::dr::friction_law::gpu
