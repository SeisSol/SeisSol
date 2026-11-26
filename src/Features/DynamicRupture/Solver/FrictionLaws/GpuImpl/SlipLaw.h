// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/SlipLaw.h"
#include "Features/DynamicRupture/Solver/Storage.h"

namespace seissol::dr::friction_law::gpu {

/// The slip law (see friction_law::SlipLaw) on the device.
template <typename Cfg, typename TPMethod>
class SlipLaw final : public BaseFrictionSolver<Cfg, SlipLaw<Cfg, TPMethod>>,
                      public friction_law::SlipLaw<Cfg, TPMethod> {
  public:
  explicit SlipLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : BaseFrictionSolver<Cfg, SlipLaw<Cfg, TPMethod>>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {
    friction_law::SlipLaw<Cfg, TPMethod>::copySpecificStorageDataToLocal(
        data, layerData, seissol::initializer::AllocationPlace::Device);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_SLIPLAW_H_
