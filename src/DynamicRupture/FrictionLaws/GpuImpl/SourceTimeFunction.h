// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SOURCETIMEFUNCTION_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SOURCETIMEFUNCTION_H_

#include "BaseFrictionSolver.h"
#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/BaseFrictionSolver.h"
#include "DynamicRupture/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "FrictionSolverInterface.h"
#include "ImposedSlipRates.h"
#include "Numerical/DeltaPulse.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Numerical/RegularizedYoffe.h"

namespace seissol::dr::friction_law::gpu {

template <typename Cfg>
class YoffeSTF : public ImposedSlipRates<Cfg, YoffeSTF<Cfg>> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// a time function scaling the slip of the point (cf. ScriptedSTF)
  static constexpr bool Prescribed = false;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data, DynamicRupture::Layer& layerData) {
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->onsetTime = layerData.var<LTSImposedSlipRatesYoffe::OnsetTime>(Cfg(), place);
    data->tauS = layerData.var<LTSImposedSlipRatesYoffe::TauS>(Cfg(), place);
    data->tauR = layerData.var<LTSImposedSlipRatesYoffe::TauR>(Cfg(), place);
  }

  SEISSOL_DEVICE static real evaluateSTF(FrictionLawContext<Cfg>& __restrict ctx,
                                         real currentTime,
                                         [[maybe_unused]] real timeIncrement) {
    return regularizedYoffe::regularizedYoffe(currentTime -
                                                  ctx.data->onsetTime[ctx.ltsFace][ctx.pointIndex],
                                              ctx.data->tauS[ctx.ltsFace][ctx.pointIndex],
                                              ctx.data->tauR[ctx.ltsFace][ctx.pointIndex]);
  }
};

template <typename Cfg>
class GaussianSTF : public ImposedSlipRates<Cfg, GaussianSTF<Cfg>> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// a time function scaling the slip of the point (cf. ScriptedSTF)
  static constexpr bool Prescribed = false;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data, DynamicRupture::Layer& layerData) {
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->onsetTime = layerData.var<LTSImposedSlipRatesGaussian::OnsetTime>(Cfg(), place);
    data->riseTime = layerData.var<LTSImposedSlipRatesGaussian::RiseTime>(Cfg(), place);
  }

  SEISSOL_DEVICE static real
      evaluateSTF(FrictionLawContext<Cfg>& __restrict ctx, real currentTime, real timeIncrement) {
    const real smoothStepIncrement = gaussianNucleationFunction::smoothStepIncrement(
        currentTime - ctx.data->onsetTime[ctx.ltsFace][ctx.pointIndex],
        timeIncrement,
        ctx.data->riseTime[ctx.ltsFace][ctx.pointIndex]);
    return smoothStepIncrement / timeIncrement;
  }
};

template <typename Cfg>
class DeltaSTF : public ImposedSlipRates<Cfg, DeltaSTF<Cfg>> {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// a time function scaling the slip of the point (cf. ScriptedSTF)
  static constexpr bool Prescribed = false;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data, DynamicRupture::Layer& layerData) {
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->onsetTime = layerData.var<LTSImposedSlipRatesDelta::OnsetTime>(Cfg(), place);
  }

  SEISSOL_DEVICE static real
      evaluateSTF(FrictionLawContext<Cfg>& __restrict ctx, real currentTime, real timeIncrement) {
    return deltaPulse::deltaPulse(currentTime - ctx.data->onsetTime[ctx.ltsFace][ctx.pointIndex],
                                  timeIncrement);
  }
};

/**
 * The slip rates of a script (FL 36), along both directions of the face and per sub-step, as
 * dr::friction_law::SlipRateEvaluator wrote them before the step: see GpuImpl/ScriptedSlipRates.h.
 */
template <typename Cfg>
class ScriptedSTF {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  /// gives the slip rate along both directions itself
  static constexpr bool Prescribed = true;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data, DynamicRupture::Layer& layerData) {
    const auto place = seissol::initializer::AllocationPlace::Device;
    data->scriptedSlipRates = reinterpret_cast<const real*>(
        layerData.var<LTSImposedSlipRatesScript::ScriptSlipRates>(Cfg(), place));
    data->scriptedPoints = layerData.size() * misc::NumPaddedPoints<Cfg>;
  }

  SEISSOL_DEVICE static void slipRates(FrictionLawContext<Cfg>& __restrict ctx,
                                       uint32_t timeIndex,
                                       real& rate1,
                                       real& rate2) {
    const std::size_t point = ctx.ltsFace * misc::NumPaddedPoints<Cfg> + ctx.pointIndex;
    const std::size_t points = ctx.data->scriptedPoints;
    const auto first = static_cast<std::size_t>(2 * timeIndex) * points + point;
    rate1 = ctx.data->scriptedSlipRates[first];
    rate2 = ctx.data->scriptedSlipRates[first + points];
  }
};

} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_SOURCETIMEFUNCTION_H_
