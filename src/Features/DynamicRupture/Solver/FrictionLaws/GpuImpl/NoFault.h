// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_NOFAULT_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/NoFault.h"
#include "Features/DynamicRupture/Solver/Storage.h"

namespace seissol::dr::friction_law::gpu {

/// No friction (see friction_law::NoFault) on the device.
template <typename Cfg>
class NoFault final : public BaseFrictionSolver<Cfg, NoFault<Cfg>>,
                      public friction_law::NoFault<Cfg> {
  public:
  explicit NoFault(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : BaseFrictionSolver<Cfg, NoFault<Cfg>>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {
    friction_law::NoFault<Cfg>::copySpecificStorageDataToLocal(
        data, layerData, seissol::initializer::AllocationPlace::Device);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_NOFAULT_H_
