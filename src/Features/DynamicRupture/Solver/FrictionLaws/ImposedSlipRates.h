// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_IMPOSEDSLIPRATES_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cstdint>

namespace seissol::dr::friction_law {

/**
 * Slip rates set to the values of a source time function `STF` along a prescribed direction, once
 * for the host and the device (see FrictionLawDriver; the solvers are cpu::ImposedSlipRates and
 * gpu::ImposedSlipRates).
 */
template <typename Cfg, typename STF>
class ImposedSlipRates {
  public:
  using real = Real<Cfg>;

  // the slip rate of a sub-step is the one prescribed at the end of its interval
  static constexpr bool SlipRateAtIntervalEnds = true;

  static void copySpecificStorageDataToLocal(FrictionLawData<Cfg>* data,
                                             DynamicRupture::Layer& layerData,
                                             seissol::initializer::AllocationPlace place) {
    data->imposedSlipDirection1 =
        layerData.var<LTSImposedSlipRates::ImposedSlipDirection1>(Cfg(), place);
    data->imposedSlipDirection2 =
        layerData.var<LTSImposedSlipRates::ImposedSlipDirection2>(Cfg(), place);
    STF::copyStorageToLocal(data, layerData, place);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void preHook(FrictionLawContext<Cfg, Type>& /*ctx*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void postHook(FrictionLawContext<Cfg, Type>& /*ctx*/) {}

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void updateFrictionAndSlip(FrictionLawContext<Cfg, Type>& ctx,
                                                       uint32_t timeIndex) {
    const real timeIncrement = ctx.args->deltaT[timeIndex];
    real currentTime = ctx.fullUpdateTime;
    for (uint32_t i = 0; i <= timeIndex; i++) {
      currentTime += ctx.args->deltaT[i];
    }
    const auto ltsFace = ctx.ltsFace;
    const auto& impAndEta = ctx.data->impAndEta[ltsFace];
    const auto& impedanceMatrices = ctx.data->impedanceMatrices[ltsFace];

    // not vectorized: the source time functions branch by the time and call transcendental
    // functions; vectorized, every point would evaluate all of their branches
    forEachPointScalar(ctx, [&](std::uint32_t pointIndex) {
      const real stfEvaluated = STF::evaluate(ctx, currentTime, timeIncrement, pointIndex);

      const auto evalCardinal1 =
          ctx.data->imposedSlipDirection1[ltsFace][pointIndex] * stfEvaluated;
      const auto evalCardinal2 =
          ctx.data->imposedSlipDirection2[ltsFace][pointIndex] * stfEvaluated;

      const auto [tU1, tU2] =
          common::matmulEta<Cfg>(impAndEta, impedanceMatrices, evalCardinal1, evalCardinal2);

      // the prescribed slip rate also changes the fault-normal traction if the impedance couples
      // the normal and the tangential directions (anisotropy); zero otherwise
      const auto tUN =
          common::matmulEtaNormal<Cfg>(impAndEta, impedanceMatrices, evalCardinal1, evalCardinal2);

      const auto traction1 = ctx.faultStresses.traction1[pointIndex] - tU1;
      const auto traction2 = ctx.faultStresses.traction2[pointIndex] - tU2;

      ctx.data->traction1[ltsFace][pointIndex] = traction1;
      ctx.data->traction2[ltsFace][pointIndex] = traction2;

      ctx.data->slipRate1[ltsFace][pointIndex] =
          ctx.data->imposedSlipDirection1[ltsFace][pointIndex] * stfEvaluated;
      ctx.data->slipRate2[ltsFace][pointIndex] =
          ctx.data->imposedSlipDirection2[ltsFace][pointIndex] * stfEvaluated;
      ctx.data->slipRateMagnitude[ltsFace][pointIndex] = misc::magnitude(
          ctx.data->slipRate1[ltsFace][pointIndex], ctx.data->slipRate2[ltsFace][pointIndex]);

      // Update slip
      ctx.data->slip1[ltsFace][pointIndex] +=
          ctx.data->slipRate1[ltsFace][pointIndex] * timeIncrement;
      ctx.data->slip2[ltsFace][pointIndex] +=
          ctx.data->slipRate2[ltsFace][pointIndex] * timeIncrement;
      ctx.data->accumulatedSlipMagnitude[ltsFace][pointIndex] +=
          ctx.data->slipRateMagnitude[ltsFace][pointIndex] * timeIncrement;

      ctx.tractionResults.normalStress[pointIndex] =
          ctx.faultStresses.normalStress[pointIndex] - tUN;
      ctx.tractionResults.traction1[pointIndex] = traction1;
      ctx.tractionResults.traction2[pointIndex] = traction2;
    });
  }

  //! no dynamic stress time
  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static PointValues<Cfg, Type>
      dynamicStressQuantity(FrictionLawContext<Cfg, Type>& /*ctx*/) {
    return PointValues<Cfg, Type>{};
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static void
      saveDynamicStressOutput(FrictionLawContext<Cfg, Type>& /*ctx*/,
                              const PointValues<Cfg, Type>& /*previousQuantity*/,
                              const common::SubStepTimes<Cfg>& /*times*/) {}
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_IMPOSEDSLIPRATES_H_
