// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FastVelocityWeakeningLaw.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"

namespace seissol::dr::friction_law::gpu {

/// The rate-and-state law with fast velocity weakening (see friction_law::FastVelocityWeakeningLaw)
/// on the device.
template <typename Cfg, typename TPMethod>
class FastVelocityWeakeningLaw final
    : public BaseFrictionSolver<Cfg, FastVelocityWeakeningLaw<Cfg, TPMethod>>,
      public friction_law::FastVelocityWeakeningLaw<Cfg, TPMethod> {
  public:
  explicit FastVelocityWeakeningLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : BaseFrictionSolver<Cfg, FastVelocityWeakeningLaw<Cfg, TPMethod>>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {
    friction_law::FastVelocityWeakeningLaw<Cfg, TPMethod>::copySpecificStorageDataToLocal(
        data, layerData, seissol::initializer::AllocationPlace::Device);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
