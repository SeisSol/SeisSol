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
#include <string>
#include <unordered_map>
#include <utils/logger.h>

namespace seissol::dr::initializer {
void RateAndStateInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  BaseDRInitializer::initializeFault(drStorage);

  for (auto& layer : drStorage.leaves(Ghost)) {

    auto* dynStressTimePending = layer.var<LTSRateAndState::DynStressTimePending>();
    real(*slipRate1)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::SlipRate1>();
    real(*slipRate2)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::SlipRate2>();
    real(*mu)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::Mu>();

    real(*stateVariable)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::StateVariable>();
    real(*rsSl0)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::RsSl0>();
    real(*rsA)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::RsA>();
    auto* initialStressInFaultCS = layer.var<LTSRateAndState::InitialStressInFaultCS>();

    const real initialSlipRate =
        misc::magnitude(drParameters->rsInitialSlipRate1, drParameters->rsInitialSlipRate2);

    using namespace dr::misc::quantity_indices;
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = true;
        slipRate1[ltsFace][pointIndex] = drParameters->rsInitialSlipRate1;
        slipRate2[ltsFace][pointIndex] = drParameters->rsInitialSlipRate2;
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
  }
}

RateAndStateInitializer::StateAndFriction
    RateAndStateInitializer::computeInitialStateAndFriction(real traction1,
                                                            real traction2,
                                                            real pressure,
                                                            real rsA,
                                                            real rsB,
                                                            real rsSl0,
                                                            real rsSr0,
                                                            real rsF0,
                                                            real initialSlipRate) {
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
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  real(*rsSl0)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::RsSl0>();
  real(*rsA)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::RsA>();
  parameterToStorageMap.insert({"rs_sl0", reinterpret_cast<real*>(rsSl0)});
  parameterToStorageMap.insert({"rs_a", reinterpret_cast<real*>(rsA)});
}

RateAndStateInitializer::StateAndFriction
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(real traction1,
                                                                        real traction2,
                                                                        real pressure,
                                                                        real rsA,
                                                                        real rsB,
                                                                        real rsSl0,
                                                                        real rsSr0,
                                                                        real rsF0,
                                                                        real initialSlipRate) {
  StateAndFriction result{};
  const real absoluteTraction = misc::magnitude(traction1, traction2);
  const real tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable =
      rsA * std::log(2.0 * rsSr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  if (result.stateVariable < 0) {
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const real tmp2 = initialSlipRate * 0.5 / rsSr0 * std::exp(result.stateVariable / rsA);
  result.frictionCoefficient = rsA * std::asinh(tmp2);
  return result;
}

void RateAndStateFastVelocityInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, layer);
  real(*rsSrW)[misc::NumPaddedPoints] = layer.var<LTSRateAndStateFastVelocityWeakening::RsSrW>();
  parameterToStorageMap.insert({"rs_srW", reinterpret_cast<real*>(rsSrW)});
}

ThermalPressurizationInitializer::ThermalPressurizationInitializer(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters)
    : drParameters(drParameters) {}

void ThermalPressurizationInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  for (auto& layer : drStorage.leaves(Ghost)) {
    real(*temperature)[misc::NumPaddedPoints] = layer.var<LTSThermalPressurization::Temperature>();
    real(*pressure)[misc::NumPaddedPoints] = layer.var<LTSThermalPressurization::Pressure>();
    auto* theta = layer.var<LTSThermalPressurization::Theta>();
    auto* sigma = layer.var<LTSThermalPressurization::Sigma>();

    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters->initialTemperature;
        pressure[ltsFace][pointIndex] = drParameters->initialPressure;
        for (unsigned tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
             ++tpGridPointIndex) {
          theta[ltsFace][tpGridPointIndex][pointIndex] = 0.0;
          sigma[ltsFace][tpGridPointIndex][pointIndex] = 0.0;
        }
      }
    }
  }
}

void ThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  real(*halfWidthShearZone)[misc::NumPaddedPoints] =
      layer.var<LTSThermalPressurization::HalfWidthShearZone>();
  real(*hydraulicDiffusivity)[misc::NumPaddedPoints] =
      layer.var<LTSThermalPressurization::HydraulicDiffusivity>();
  parameterToStorageMap.insert(
      {"tp_halfWidthShearZone", reinterpret_cast<real*>(halfWidthShearZone)});
  parameterToStorageMap.insert(
      {"tp_hydraulicDiffusivity", reinterpret_cast<real*>(hydraulicDiffusivity)});
}

RateAndStateThermalPressurizationInitializer::RateAndStateThermalPressurizationInitializer(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
    SeisSol& instance)
    : RateAndStateInitializer(drParameters, instance),
      ThermalPressurizationInitializer(drParameters) {}

void RateAndStateThermalPressurizationInitializer::initializeFault(
    DynamicRupture::Storage& drStorage) {
  RateAndStateInitializer::initializeFault(drStorage);
  ThermalPressurizationInitializer::initializeFault(drStorage);
}

void RateAndStateThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, layer);
  ThermalPressurizationInitializer::addAdditionalParameters(parameterToStorageMap, layer);
}

RateAndStateFastVelocityThermalPressurizationInitializer::
    RateAndStateFastVelocityThermalPressurizationInitializer(
        const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
        SeisSol& instance)
    : RateAndStateFastVelocityInitializer(drParameters, instance),
      ThermalPressurizationInitializer(drParameters) {}

void RateAndStateFastVelocityThermalPressurizationInitializer::initializeFault(
    DynamicRupture::Storage& drStorage) {
  RateAndStateFastVelocityInitializer::initializeFault(drStorage);
  ThermalPressurizationInitializer::initializeFault(drStorage);
}

void RateAndStateFastVelocityThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  RateAndStateFastVelocityInitializer::addAdditionalParameters(parameterToStorageMap, layer);
  ThermalPressurizationInitializer::addAdditionalParameters(parameterToStorageMap, layer);
}

} // namespace seissol::dr::initializer
