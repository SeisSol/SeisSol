// SPDX-FileCopyrightText: 2022 SeisSol Group
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

  static void copySpecificStorageDataToLocal(FrictionLawData* data,
                                             DynamicRupture::Layer& layerData) {
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->imposedSlipDirection1 = layerData.var<LTSImposedSlipRates::ImposedSlipDirection1>(place);
    data->imposedSlipDirection2 = layerData.var<LTSImposedSlipRates::ImposedSlipDirection2>(place);
    STF::copyStorageToLocal(data, layerData);
  }

  SEISSOL_DEVICE static void updateFrictionAndSlip(FrictionLawContext& __restrict ctx,
                                                   uint32_t timeIndex) {
    const real timeIncrement = ctx.args->deltaT[timeIndex];
    real currentTime = ctx.args->fullUpdateTime;
    for (uint32_t i = 0; i <= timeIndex; i++) {
      currentTime += ctx.args->deltaT[i];
    }

    const auto stfEvaluated = STF::evaluateSTF(ctx, currentTime, timeIncrement);

    const auto evalCardinal1 =
        ctx.data->imposedSlipDirection1[ctx.ltsFace][ctx.pointIndex] * stfEvaluated;
    const auto evalCardinal2 =
        ctx.data->imposedSlipDirection2[ctx.ltsFace][ctx.pointIndex] * stfEvaluated;

    const auto [tU1, tU2] = common::matmulEta(ctx.data->impAndEta[ctx.ltsFace],
                                              ctx.data->impedanceMatrices[ctx.ltsFace],
                                              evalCardinal1,
                                              evalCardinal2);

    // the prescribed slip rate also changes the fault-normal traction if the impedance couples the
    // normal and the tangential directions (anisotropy); zero otherwise
    const auto tUN = common::matmulEtaNormal(ctx.data->impAndEta[ctx.ltsFace],
                                             ctx.data->impedanceMatrices[ctx.ltsFace],
                                             evalCardinal1,
                                             evalCardinal2);

    ctx.data->traction1[ctx.ltsFace][ctx.pointIndex] = ctx.faultStresses.traction1[timeIndex] - tU1;
    ctx.data->traction2[ctx.ltsFace][ctx.pointIndex] = ctx.faultStresses.traction2[timeIndex] - tU2;

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

    ctx.tractionResults.normalStress[timeIndex] = ctx.faultStresses.normalStress[timeIndex] - tUN;
    ctx.tractionResults.traction1[timeIndex] = ctx.data->traction1[ctx.ltsFace][ctx.pointIndex];
    ctx.tractionResults.traction2[timeIndex] = ctx.data->traction2[ctx.ltsFace][ctx.pointIndex];
  }

  SEISSOL_DEVICE static void saveDynamicStressOutput(FrictionLawContext& __restrict ctx,
                                                     real time) {}
  SEISSOL_DEVICE static void preHook(FrictionLawContext& __restrict ctx) {}
  SEISSOL_DEVICE static void postHook(FrictionLawContext& __restrict ctx) {}
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_IMPOSEDSLIPRATES_H_
