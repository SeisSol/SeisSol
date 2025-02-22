// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SOURCETIMEFUNCTION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SOURCETIMEFUNCTION_H_

#include "BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "FrictionSolverInterface.h"
#include "ImposedSlipRates.h"
#include <DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h>

#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"

namespace seissol::dr::friction_law::gpu {

class YoffeSTF : public ImposedSlipRates<YoffeSTF> {
  public:
  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSImposedSlipRatesYoffe*>(dynRup);
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->onsetTime = layerData.var(concreteLts->onsetTime, place);
    data->tauS = layerData.var(concreteLts->tauS, place);
    data->tauR = layerData.var(concreteLts->tauR, place);
  }

  SEISSOL_DEVICE static real
      evaluateSTF(FrictionLawContext& ctx, real currentTime, [[maybe_unused]] real timeIncrement) {
    return regularizedYoffe::regularizedYoffe(currentTime -
                                                  ctx.data->onsetTime[ctx.ltsFace][ctx.pointIndex],
                                              ctx.data->tauS[ctx.ltsFace][ctx.pointIndex],
                                              ctx.data->tauR[ctx.ltsFace][ctx.pointIndex]);
  }
};

class GaussianSTF : public ImposedSlipRates<GaussianSTF> {
  public:
  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSImposedSlipRatesGaussian*>(dynRup);
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->onsetTime = layerData.var(concreteLts->onsetTime, place);
    data->riseTime = layerData.var(concreteLts->riseTime, place);
  }

  SEISSOL_DEVICE static real
      evaluateSTF(FrictionLawContext& ctx, real currentTime, real timeIncrement) {
    const real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
        currentTime - ctx.data->onsetTime[ctx.ltsFace][ctx.pointIndex],
        timeIncrement,
        ctx.data->riseTime[ctx.ltsFace][ctx.pointIndex]);
    return smoothStepIncrement / timeIncrement;
  }
};

class DeltaSTF : public ImposedSlipRates<DeltaSTF> {
  public:
  static void copyLtsTreeToLocal(FrictionLawData* data,
                                 seissol::initializer::Layer& layerData,
                                 const seissol::initializer::DynamicRupture* const dynRup,
                                 real fullUpdateTime) {
    const auto* concreteLts =
        dynamic_cast<const seissol::initializer::LTSImposedSlipRatesDelta*>(dynRup);
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->onsetTime = layerData.var(concreteLts->onsetTime, place);
  }

  SEISSOL_DEVICE static real
      evaluateSTF(FrictionLawContext& ctx, real currentTime, real timeIncrement) {
    return deltaPulse::deltaPulse(currentTime - ctx.data->onsetTime[ctx.ltsFace][ctx.pointIndex],
                                  timeIncrement);
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SOURCETIMEFUNCTION_H_
