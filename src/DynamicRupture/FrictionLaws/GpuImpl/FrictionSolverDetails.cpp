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
#include "GeneratedCode/init.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include <Initializer/Typedefs.h>
#include <Memory/GlobalData.h>
#include <cstddef>

#include "Memory/MemoryAllocator.h"

#include "DynamicRupture/FrictionLaws/TPCommon.h"

namespace seissol::dr::friction_law::gpu {
template <typename Cfg>
FrictionSolverDetails<Cfg>::FrictionSolverDetails(
    seissol::initializer::parameters::DRParameters* drParameters)
    : FrictionSolverInterface<Cfg>(drParameters) {}

template <typename Cfg>
FrictionSolverDetails<Cfg>::~FrictionSolverDetails() = default;

template <typename Cfg>
void FrictionSolverDetails<Cfg>::allocateAuxiliaryMemory(const GlobalData& globalData) {
  {
    data = seissol::memory::allocTyped<FrictionLawData<Cfg>>(1, 1, memory::DeviceGlobalMemory);
  }

  const auto& global = globalData.get<Cfg, Executor::Device>();

  resampleMatrix = global.resampleMatrix;
  devSpaceWeights = global.spaceWeights;
  devTpInverseFourierCoefficients = global.tpInverseFourierCoefficients;
  devHeatSource = global.heatSource;
  devTpGridPoints = global.tpGridPoints;
}

#define SEISSOL_CONFIGITER(cfg) template class FrictionSolverDetails<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::dr::friction_law::gpu
