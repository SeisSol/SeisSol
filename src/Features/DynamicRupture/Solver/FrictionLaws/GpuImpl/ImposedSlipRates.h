// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_IMPOSEDSLIPRATES_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/ImposedSlipRates.h"
#include "Features/DynamicRupture/Solver/Storage.h"

namespace seissol::dr::friction_law::gpu {

/// The imposed slip rates (see friction_law::ImposedSlipRates) on the device.
template <typename Cfg, typename STF>
class ImposedSlipRates final : public BaseFrictionSolver<Cfg, ImposedSlipRates<Cfg, STF>>,
                               public friction_law::ImposedSlipRates<Cfg, STF> {
  public:
  explicit ImposedSlipRates(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : BaseFrictionSolver<Cfg, ImposedSlipRates<Cfg, STF>>(drParameters) {}

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData) {
    friction_law::ImposedSlipRates<Cfg, STF>::copySpecificStorageDataToLocal(
        data, layerData, seissol::initializer::AllocationPlace::Device);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_IMPOSEDSLIPRATES_H_
