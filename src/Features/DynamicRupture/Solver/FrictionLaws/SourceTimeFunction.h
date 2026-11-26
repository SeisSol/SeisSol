// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SOURCETIMEFUNCTION_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SOURCETIMEFUNCTION_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Numerical/DeltaPulse.h"
#include "Core/Numerical/GaussianNucleationFunction.h"
#include "Core/Numerical/RegularizedYoffe.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cstdint>

/**
 * The source time functions of the imposed slip rates (see ImposedSlipRates), once for the host
 * and the device: the slip rate a point imposes at a time, from its parameters.
 */
namespace seissol::dr::friction_law {

template <typename Cfg>
class YoffeSTF {
  public:
  using real = Real<Cfg>;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data,
                                 DynamicRupture::Layer& layerData,
                                 seissol::initializer::AllocationPlace place) {
    data->onsetTime = layerData.var<LTSImposedSlipRatesYoffe::OnsetTime>(Cfg(), place);
    data->tauS = layerData.var<LTSImposedSlipRatesYoffe::TauS>(Cfg(), place);
    data->tauR = layerData.var<LTSImposedSlipRatesYoffe::TauR>(Cfg(), place);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real evaluate(FrictionLawContext<Cfg, Type>& ctx,
                                          real currentTime,
                                          [[maybe_unused]] real timeIncrement,
                                          std::uint32_t pointIndex) {
    return regularized_yoffe::regularizedYoffe(currentTime -
                                                   ctx.data->onsetTime[ctx.ltsFace][pointIndex],
                                               ctx.data->tauS[ctx.ltsFace][pointIndex],
                                               ctx.data->tauR[ctx.ltsFace][pointIndex]);
  }
};

template <typename Cfg>
class GaussianSTF {
  public:
  using real = Real<Cfg>;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data,
                                 DynamicRupture::Layer& layerData,
                                 seissol::initializer::AllocationPlace place) {
    data->onsetTime = layerData.var<LTSImposedSlipRatesGaussian::OnsetTime>(Cfg(), place);
    data->riseTime = layerData.var<LTSImposedSlipRatesGaussian::RiseTime>(Cfg(), place);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real evaluate(FrictionLawContext<Cfg, Type>& ctx,
                                          real currentTime,
                                          real timeIncrement,
                                          std::uint32_t pointIndex) {
    const real smoothStepIncrement = gaussian_nucleation_function::smoothStepIncrement(
        currentTime - ctx.data->onsetTime[ctx.ltsFace][pointIndex],
        timeIncrement,
        ctx.data->riseTime[ctx.ltsFace][pointIndex]);
    return smoothStepIncrement / timeIncrement;
  }
};

template <typename Cfg>
class DeltaSTF {
  public:
  using real = Real<Cfg>;

  static void copyStorageToLocal(FrictionLawData<Cfg>* data,
                                 DynamicRupture::Layer& layerData,
                                 seissol::initializer::AllocationPlace place) {
    data->onsetTime = layerData.var<LTSImposedSlipRatesDelta::OnsetTime>(Cfg(), place);
  }

  template <common::RangeType Type>
  SEISSOL_HOSTDEVICE static real evaluate(FrictionLawContext<Cfg, Type>& ctx,
                                          real currentTime,
                                          real timeIncrement,
                                          std::uint32_t pointIndex) {
    // Currently, the delta pulse is normalized in time equivalent to FL33 and FL34
    return delta_pulse::deltaPulse(currentTime - ctx.data->onsetTime[ctx.ltsFace][pointIndex],
                                   timeIncrement);
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_SOURCETIMEFUNCTION_H_
