// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/AgingLaw.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"

namespace seissol::dr::friction_law::gpu {

/// The aging law (see friction_law::AgingLaw) on the device.
template <typename Cfg, typename TPMethod>
class AgingLaw final : public BaseFrictionSolver<Cfg, AgingLaw<Cfg, TPMethod>>,
                       public friction_law::AgingLaw<Cfg, TPMethod> {
  public:
  explicit AgingLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : BaseFrictionSolver<Cfg, AgingLaw<Cfg, TPMethod>>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {
    friction_law::AgingLaw<Cfg, TPMethod>::copySpecificStorageDataToLocal(
        data, layerData, seissol::initializer::AllocationPlace::Device);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_AGINGLAW_H_
