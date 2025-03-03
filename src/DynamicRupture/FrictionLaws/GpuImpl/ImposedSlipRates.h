// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_IMPOSEDSLIPRATES_H_

#include "BaseFrictionSolver.h"

namespace seissol::dr::friction_law::gpu {
/**
 * Slip rates are set fixed values
 */
template <typename STF>
class ImposedSlipRates : public BaseFrictionSolver<ImposedSlipRates<STF>> {
  public:
  using BaseFrictionSolver<ImposedSlipRates>::BaseFrictionSolver;

  static void
      copySpecificLtsDataTreeToLocal(FrictionLawData* data,
                                     seissol::initializer::Layer& layerData,
                                     const seissol::initializer::DynamicRupture* const dynRup,
                                     real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSImposedSlipRates*>(dynRup);
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->imposedSlipDirection1 = layerData.var(concreteLts->imposedSlipDirection1, place);
    data->imposedSlipDirection2 = layerData.var(concreteLts->imposedSlipDirection2, place);
    STF::copyLtsTreeToLocal(data, layerData, dynRup, fullUpdateTime);
  }

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& ctx, unsigned timeIndex) {
    const real timeIncrement = ctx.data->deltaT[timeIndex];
    real currentTime = ctx.fullUpdateTime;
    for (int i = 0; i <= timeIndex; i++) {
      currentTime += ctx.data->deltaT[i];
    }

    const auto stfEvaluated = STF::evaluateSTF(ctx, currentTime, timeIncrement);

    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] =
        ctx.faultStresses.traction1[timeIndex] -
        ctx.data->impAndEta[ctx.ltsFace].etaS *
            ctx.data->imposedSlipDirection1[ctx.ltsFace][ctx.pointIndex] * stfEvaluated;
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] =
        ctx.faultStresses.traction2[timeIndex] -
        ctx.data->impAndEta[ctx.ltsFace].etaS *
            ctx.data->imposedSlipDirection2[ctx.ltsFace][ctx.pointIndex] * stfEvaluated;

    ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] =
        ctx.data->imposedSlipDirection1[ctx.ltsFace][ctx.pointIndex] * stfEvaluated;
    ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] =
        ctx.data->imposedSlipDirection2[ctx.ltsFace][ctx.pointIndex] * stfEvaluated;
    ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] =
        misc::magnitude(ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex],
                        ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex]);

    // Update slip
    ctx.data->slip1[ctx.ltsFace][ctx.pointIndex] +=
        ctx.data->slipRate1[ctx.ltsFace][ctx.pointIndex] * timeIncrement;
    ctx.data->slip2[ctx.ltsFace][ctx.pointIndex] +=
        ctx.data->slipRate2[ctx.ltsFace][ctx.pointIndex] * timeIncrement;
    ctx.data->accumulatedSlipMagnitude[ctx.ltsFace][ctx.pointIndex] +=
        ctx.data->slipRateMagnitude[ctx.ltsFace][ctx.pointIndex] * timeIncrement;

    ctx.tractionResults.traction1[timeIndex] = ctx.data->traction1[ctx.ltsFace][ctx.pointIndex];
    ctx.tractionResults.traction2[timeIndex] = ctx.data->traction2[ctx.ltsFace][ctx.pointIndex];
  }

  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& ctx) {}
  SEISSOL_DEVICE static void preHook(FrictionLawContext& ctx) {}
  SEISSOL_DEVICE static void postHook(FrictionLawContext& ctx) {}
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_IMPOSEDSLIPRATES_H_
