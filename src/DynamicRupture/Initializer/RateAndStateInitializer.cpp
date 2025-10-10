// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "RateAndStateInitializer.h"

#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include <Initializer/Parameters/DRParameters.h>
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

  const auto rsF0Param = !faultProvides("rs_f0");
  const auto rsMuWParam = !faultProvides("rs_muw");
  const auto rsBParam = !faultProvides("rs_b");

  logInfo() << "RS parameter source (1 == from parameter file, 0 == from easi file): f0"
            << rsF0Param << "- muW" << rsMuWParam << "- b" << rsBParam;

  for (auto& layer : drStorage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>;
      auto* dynStressTimePending = layer.var<LTSRateAndState::DynStressTimePending>(cfg);
      real(*slipRate1)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::SlipRate1>(cfg);
      real(*slipRate2)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::SlipRate2>(cfg);
      real(*mu)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::Mu>(cfg);

      real(*stateVariable)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSRateAndState::StateVariable>(cfg);
      real(*rsSl0)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsSl0>(cfg);
      real(*rsA)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsA>(cfg);
      auto* initialStressInFaultCS = layer.var<LTSRateAndState::InitialStressInFaultCS>(cfg);

      auto* rsF0 = layer.var<LTSRateAndState::RsF0>(cfg);
      auto* rsMuW = layer.var<LTSRateAndState::RsMuW>(cfg);
      auto* rsB = layer.var<LTSRateAndState::RsB>(cfg);

      const double initialSlipRate =
          misc::magnitude(drParameters->rsInitialSlipRate1, drParameters->rsInitialSlipRate2);

      using namespace dr::misc::quantity_indices;
      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
          dynStressTimePending[ltsFace][pointIndex] = true;
          slipRate1[ltsFace][pointIndex] = drParameters->rsInitialSlipRate1;
          slipRate2[ltsFace][pointIndex] = drParameters->rsInitialSlipRate2;

          if (rsF0Param) {
            rsF0[ltsFace][pointIndex] = drParameters->rsF0;
          }
          if (rsMuWParam) {
            rsMuW[ltsFace][pointIndex] = drParameters->muW;
          }
          if (rsBParam) {
            rsB[ltsFace][pointIndex] = drParameters->rsB;
          }

          // compute initial friction and state
          auto stateAndFriction =
              computeInitialStateAndFriction(initialStressInFaultCS[ltsFace][XY][pointIndex],
                                             initialStressInFaultCS[ltsFace][XZ][pointIndex],
                                             initialStressInFaultCS[ltsFace][XX][pointIndex],
                                             rsA[ltsFace][pointIndex],
                                             drParameters->rsB,
                                             rsSl0[ltsFace][pointIndex],
                                             drParameters->rsSr0,
                                             drParameters->rsF0,
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
                         std::exp((rsA * std::log(std::exp(tmp) - std::exp(-tmp)) - rsF0 -
                                   rsA * std::log(initialSlipRate / rsSr0)) /
                                  rsB);
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const double tmp2 = initialSlipRate * 0.5 / rsSr0 *
                      std::exp((rsF0 + rsB * std::log(rsSr0 * result.stateVariable / rsSl0)) / rsA);
  result.frictionCoefficient = rsA * std::asinh(tmp2);
  return result;
}

void RateAndStateInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;
    real(*rsSl0)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsSl0>(cfg);
    real(*rsA)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSRateAndState::RsA>(cfg);
    const auto sl0Name = faultNameAlternatives({"rs_sl0", "RS_sl0"});

    parameterToStorageMap.insert({sl0Name, reinterpret_cast<real*>(rsSl0)});
    parameterToStorageMap.insert({"rs_a", reinterpret_cast<real*>(rsA)});

    const auto insertIfPresent = [&](const auto& name, auto* var) {
      if (faultProvides(name)) {
        parameterToStorageMap.insert({name, reinterpret_cast<real*>(var)});
      }
    };
    insertIfPresent("rs_f0", layer.var<LTSRateAndState::RsF0>());
    insertIfPresent("rs_muw", layer.var<LTSRateAndState::RsMuW>());
    insertIfPresent("rs_b", layer.var<LTSRateAndState::RsB>());
  });
}

RateAndStateInitializer::StateAndFriction
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(double traction1,
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
  result.stateVariable =
      rsA * std::log(2.0 * rsSr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const double tmp2 = initialSlipRate * 0.5 / rsSr0 * std::exp(result.stateVariable / rsA);
  result.frictionCoefficient = rsA * std::asinh(tmp2);
  return result;
}

void RateAndStateFastVelocityInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;
    RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, layer);
    real(*rsSrW)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSRateAndStateFastVelocityWeakening::RsSrW>(cfg);
    parameterToStorageMap.insert({"rs_srW", reinterpret_cast<real*>(rsSrW)});
  });
}

ThermalPressurizationInitializer::ThermalPressurizationInitializer(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
    const std::set<std::string>& faultParameterNames)
    : drParameters(drParameters), faultParameterNames(faultParameterNames) {}

void ThermalPressurizationInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  for (auto& layer : drStorage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>;
      real(*temperature)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSThermalPressurization::Temperature>(cfg);
      real(*pressure)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSThermalPressurization::Pressure>(cfg);
      auto* theta = layer.var<LTSThermalPressurization::Theta>(cfg);
      auto* sigma = layer.var<LTSThermalPressurization::Sigma>(cfg);

      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
          temperature[ltsFace][pointIndex] = drParameters->initialTemperature;
          pressure[ltsFace][pointIndex] = drParameters->initialPressure;
          for (unsigned tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
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
  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;
    real(*halfWidthShearZone)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSThermalPressurization::HalfWidthShearZone>(cfg);
    real(*hydraulicDiffusivity)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSThermalPressurization::HydraulicDiffusivity>(cfg);
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
      ThermalPressurizationInitializer(drParameters, RateAndStateInitializer::faultParameterNames) {
}

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
                                       RateAndStateFastVelocityInitializer::faultParameterNames) {}

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
    if (faultParameterNames.find(name) != faultParameterNames.end()) {
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
