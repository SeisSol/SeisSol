// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "RateAndStateInitializer.h"

#include "Common/ConfigDispatch.h"
#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"
#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::initializer {
void RateAndStateInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  BaseDRInitializer::initializeFault(drStorage);

  const auto rsF0ParamName = faultNameAlternatives({"rs_f0", "RS_f0"});
  const auto rsMuWParamName = faultNameAlternatives({"rs_muw", "RS_muw"});
  const auto rsBParamName = faultNameAlternatives({"rs_b", "RS_b"});

  const auto rsF0Param = !faultProvides(rsF0ParamName);
  const auto rsMuWParam = !faultProvides(rsMuWParamName);
  const auto rsBParam = !faultProvides(rsBParamName);

  logInfo() << "RS parameter source (1 == from parameter file, 0 == from easi file): f0"
            << rsF0Param << "- muW" << rsMuWParam << "- b" << rsBParam;

  for (auto& layer : drStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
      auto* dynStressTimePending = layer.var<LTSRateAndState::DynStressTimePending>(Cfg());
      real(*slipRate1)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::SlipRate1>(Cfg());
      real(*slipRate2)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::SlipRate2>(Cfg());
      real(*mu)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::Mu>(Cfg());

      real(*stateVariable)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSRateAndState::StateVariable>(Cfg());
      const real(*rsSl0)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsSl0>(Cfg());
      const real(*rsA)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsA>(Cfg());

      auto* rsF0 = layer.var<LTSRateAndState::RsF0>(Cfg());
      auto* rsMuW = layer.var<LTSRateAndState::RsMuW>(Cfg());
      auto* rsB = layer.var<LTSRateAndState::RsB>(Cfg());

      auto* convergenceInner = layer.var<LTSRateAndState::ConvergenceInner>(Cfg());
      auto* convergenceOuter = layer.var<LTSRateAndState::ConvergenceOuter>(Cfg());

      // the stress the fault starts out under, which is every source that is in effect at the
      // beginning of the simulation and not only the initial state
      const auto sourceCount = stressSourceCount(*drParameters_);
      const auto* stressSources = layer.var<LTSRateAndState::StressSourceInFaultCS>(Cfg());
      const auto* stressSourceOnset = layer.var<LTSRateAndState::StressSourceOnset>(Cfg());
      const auto* stressSourceRiseTime = layer.var<LTSRateAndState::StressSourceRiseTime>(Cfg());

      const auto initialSlipRate =
          misc::magnitude(drParameters_->rsInitialSlipRate1, drParameters_->rsInitialSlipRate2);

      using namespace dr::misc::quantity_indices;
      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
          dynStressTimePending[ltsFace][pointIndex] = true;
          slipRate1[ltsFace][pointIndex] = drParameters_->rsInitialSlipRate1;
          slipRate2[ltsFace][pointIndex] = drParameters_->rsInitialSlipRate2;

          convergenceInner[ltsFace][pointIndex] = true;
          convergenceOuter[ltsFace][pointIndex] = true;

          if (rsF0Param) {
            rsF0[ltsFace][pointIndex] = drParameters_->rsF0;
          }
          if (rsMuWParam) {
            rsMuW[ltsFace][pointIndex] = drParameters_->muW;
          }
          if (rsBParam) {
            rsB[ltsFace][pointIndex] = drParameters_->rsB;
          }

          // compute initial friction and state
          const auto initialStress = stressAtTime<Cfg>(&stressSources[ltsFace * sourceCount],
                                                       &stressSourceRiseTime[ltsFace * sourceCount],
                                                       &stressSourceOnset[ltsFace * sourceCount],
                                                       sourceCount,
                                                       pointIndex,
                                                       static_cast<real>(0.0));
          const auto stateAndFriction = computeInitialStateAndFriction(initialStress[XY],
                                                                       initialStress[XZ],
                                                                       initialStress[XX],
                                                                       rsA[ltsFace][pointIndex],
                                                                       rsB[ltsFace][pointIndex],
                                                                       rsSl0[ltsFace][pointIndex],
                                                                       drParameters_->rsSr0,
                                                                       rsF0[ltsFace][pointIndex],
                                                                       initialSlipRate);
          stateVariable[ltsFace][pointIndex] = stateAndFriction.stateVariable;
          mu[ltsFace][pointIndex] = stateAndFriction.frictionCoefficient;
        }
      }
    });
  }
}

RateAndStateInitializer::StateAndFriction
    RateAndStateInitializer::computeInitialStateAndFriction(double traction1,
                                                            double traction2,
                                                            double pressure,
                                                            double rsA,
                                                            double rsB,
                                                            double rsSl0,
                                                            double rsSr0,
                                                            double rsF0,
                                                            double initialSlipRate) {
  StateAndFriction result{};
  const double absoluteTraction = misc::magnitude(traction1, traction2);
  const double tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable = rsSl0 / rsSr0 *
                         std::exp((rsA * seissol::dr::friction_law::rs::logsinh(2.0, tmp) - rsF0 -
                                   rsA * std::log(initialSlipRate / rsSr0)) /
                                  rsB);
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const double explog = (rsF0 + rsB * std::log(rsSr0 * result.stateVariable / rsSl0)) / rsA;
  const double expval = seissol::dr::friction_law::rs::computeCExp(explog);
  const double linval = initialSlipRate * 0.5 / rsSr0;
  result.frictionCoefficient =
      rsA * seissol::dr::friction_law::rs::arsinhexp(linval, explog, expval);
  return result;
}

void RateAndStateInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*rsSl0)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsSl0>(Cfg());
    real(*rsA)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsA>(Cfg());

    const auto sl0Name = faultNameAlternatives({"rs_sl0", "RS_sl0"});

    parameterToStorageMap.insert({sl0Name, reinterpret_cast<real*>(rsSl0)});
    parameterToStorageMap.insert({"rs_a", reinterpret_cast<real*>(rsA)});

    const auto insertIfPresent = [&](const auto& name, auto* var) {
      if (faultProvides(name)) {
        parameterToStorageMap.insert({name, reinterpret_cast<real*>(var)});
      }
    };
    insertIfPresent("rs_f0", layer.var<LTSRateAndState::RsF0>(Cfg()));
    insertIfPresent("rs_muw", layer.var<LTSRateAndState::RsMuW>(Cfg()));
    insertIfPresent("rs_b", layer.var<LTSRateAndState::RsB>(Cfg()));
  });
}

RateAndStateInitializer::StateAndFriction
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(double traction1,
                                                                        double traction2,
                                                                        double pressure,
                                                                        double rsA,
                                                                        double /*rsB*/,
                                                                        double /*rsSl0*/,
                                                                        double rsSr0,
                                                                        double /*rsF0*/,
                                                                        double initialSlipRate) {
  StateAndFriction result{};
  const double absoluteTraction = misc::magnitude(traction1, traction2);
  const double tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable =
      rsA * seissol::dr::friction_law::rs::logsinh(2.0 * rsSr0 / initialSlipRate, tmp);
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const double explog = result.stateVariable / rsA;
  const double expval = seissol::dr::friction_law::rs::computeCExp(explog);
  const double linval = initialSlipRate * 0.5 / rsSr0;
  result.frictionCoefficient =
      rsA * seissol::dr::friction_law::rs::arsinhexp(linval, explog, expval);
  return result;
}

void RateAndStateFastVelocityInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, layer);
    real(*rsSrW)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSRateAndStateFastVelocityWeakening::RsSrW>(Cfg());
    parameterToStorageMap.insert({"rs_srW", reinterpret_cast<real*>(rsSrW)});
  });
}

ThermalPressurizationInitializer::ThermalPressurizationInitializer(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
    const std::set<std::string>& faultParameterNames)
    : drParameters_(drParameters), faultParameterNames_(faultParameterNames) {}

void ThermalPressurizationInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  for (auto& layer : drStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
      real(*temperature)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSThermalPressurization::Temperature>(Cfg());
      real(*pressure)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSThermalPressurization::Pressure>(Cfg());
      auto* theta = layer.var<LTSThermalPressurization::Theta>(Cfg());
      auto* sigma = layer.var<LTSThermalPressurization::Sigma>(Cfg());

      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
          temperature[ltsFace][pointIndex] = drParameters_->initialTemperature;
          pressure[ltsFace][pointIndex] = drParameters_->initialPressure;
          for (std::size_t tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
               ++tpGridPointIndex) {
            theta[ltsFace][tpGridPointIndex][pointIndex] = 0.0;
            sigma[ltsFace][tpGridPointIndex][pointIndex] = 0.0;
          }
        }
      }
    });
  }
}

void ThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*halfWidthShearZone)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSThermalPressurization::HalfWidthShearZone>(Cfg());
    real(*hydraulicDiffusivity)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSThermalPressurization::HydraulicDiffusivity>(Cfg());

    const auto halfWidthShearZoneName =
        faultNameAlternatives({"tp_halfWidthShearZone", "TP_half_width_shear_zone"});
    const auto hydraulicDiffusivityName =
        faultNameAlternatives({"tp_hydraulicDiffusivity", "alpha_hy"});

    parameterToStorageMap.insert(
        {halfWidthShearZoneName, reinterpret_cast<real*>(halfWidthShearZone)});
    parameterToStorageMap.insert(
        {hydraulicDiffusivityName, reinterpret_cast<real*>(hydraulicDiffusivity)});
  });
}

RateAndStateThermalPressurizationInitializer::RateAndStateThermalPressurizationInitializer(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
    SeisSol& instance)
    : RateAndStateInitializer(drParameters, instance),
      ThermalPressurizationInitializer(drParameters,
                                       RateAndStateInitializer::faultParameterNames_) {}

void RateAndStateThermalPressurizationInitializer::initializeFault(
    DynamicRupture::Storage& drStorage) {
  RateAndStateInitializer::initializeFault(drStorage);
  ThermalPressurizationInitializer::initializeFault(drStorage);
}

void RateAndStateThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, layer);
  ThermalPressurizationInitializer::addAdditionalParameters(parameterToStorageMap, layer);
}

RateAndStateFastVelocityThermalPressurizationInitializer::
    RateAndStateFastVelocityThermalPressurizationInitializer(
        const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
        SeisSol& instance)
    : RateAndStateFastVelocityInitializer(drParameters, instance),
      ThermalPressurizationInitializer(drParameters,
                                       RateAndStateFastVelocityInitializer::faultParameterNames_) {}

void RateAndStateFastVelocityThermalPressurizationInitializer::initializeFault(
    DynamicRupture::Storage& drStorage) {
  RateAndStateFastVelocityInitializer::initializeFault(drStorage);
  ThermalPressurizationInitializer::initializeFault(drStorage);
}

void RateAndStateFastVelocityThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  RateAndStateFastVelocityInitializer::addAdditionalParameters(parameterToStorageMap, layer);
  ThermalPressurizationInitializer::addAdditionalParameters(parameterToStorageMap, layer);
}

std::string ThermalPressurizationInitializer::faultNameAlternatives(
    const std::vector<std::string>& parameter) {
  for (const auto& name : parameter) {
    if (faultParameterNames_.find(name) != faultParameterNames_.end()) {
      if (name != parameter[0]) {
        logWarning() << "You are using the deprecated fault parameter name" << name << "for"
                     << parameter[0];
      }
      return name;
    }
  }
  return parameter[0];
}

} // namespace seissol::dr::initializer
