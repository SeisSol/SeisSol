// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "RateAndStateInitializer.h"

#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
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

  const auto rsF0Param = !faultProvides("rs_f0");
  const auto rsMuWParam = !faultProvides("rs_muw");
  const auto rsBParam = !faultProvides("rs_b");

  logInfo() << "RS parameter source (1 == from parameter file, 0 == from easi file): f0"
            << rsF0Param << "- muW" << rsMuWParam << "- b" << rsBParam;

  for (auto& layer : drStorage.leaves(Ghost)) {

    auto* dynStressTimePending = layer.var<LTSRateAndState::DynStressTimePending>();
    real(*slipRate1)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::SlipRate1>();
    real(*slipRate2)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::SlipRate2>();
    real(*mu)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::Mu>();

    real(*stateVariable)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::StateVariable>();
    const real(*rsSl0)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::RsSl0>();
    const real(*rsA)[misc::NumPaddedPoints] = layer.var<LTSRateAndState::RsA>();

    auto* rsF0 = layer.var<LTSRateAndState::RsF0>();
    auto* rsMuW = layer.var<LTSRateAndState::RsMuW>();
    auto* rsB = layer.var<LTSRateAndState::RsB>();

    auto* initialStressInFaultCS = layer.var<LTSRateAndState::InitialStressInFaultCS>();

    const real initialSlipRate =
        misc::magnitude(drParameters_->rsInitialSlipRate1, drParameters_->rsInitialSlipRate2);

    using namespace dr::misc::quantity_indices;
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = true;
        slipRate1[ltsFace][pointIndex] = drParameters_->rsInitialSlipRate1;
        slipRate2[ltsFace][pointIndex] = drParameters_->rsInitialSlipRate2;

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
        const auto stateAndFriction =
            computeInitialStateAndFriction(initialStressInFaultCS[ltsFace][XY][pointIndex],
                                           initialStressInFaultCS[ltsFace][XZ][pointIndex],
                                           initialStressInFaultCS[ltsFace][XX][pointIndex],
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
}

RateAndStateInitializer::StateAndFriction
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(real traction1,
                                                                        real traction2,
                                                                        real pressure,
                                                                        real rsA,
                                                                        real /*rsB*/,
                                                                        real /*rsSl0*/,
                                                                        real rsSr0,
                                                                        real /*rsF0*/,
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
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
    const std::set<std::string>& faultParameterNames)
    : drParameters_(drParameters), faultParameterNames_(faultParameterNames) {}

void ThermalPressurizationInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  for (auto& layer : drStorage.leaves(Ghost)) {
    real(*temperature)[misc::NumPaddedPoints] = layer.var<LTSThermalPressurization::Temperature>();
    real(*pressure)[misc::NumPaddedPoints] = layer.var<LTSThermalPressurization::Pressure>();
    auto* theta = layer.var<LTSThermalPressurization::Theta>();
    auto* sigma = layer.var<LTSThermalPressurization::Sigma>();

    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters_->initialTemperature;
        pressure[ltsFace][pointIndex] = drParameters_->initialPressure;
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

  const auto halfWidthShearZoneName =
      faultNameAlternatives({"tp_halfWidthShearZone", "TP_half_width_shear_zone"});
  const auto hydraulicDiffusivityName =
      faultNameAlternatives({"tp_hydraulicDiffusivity", "alpha_hy"});

  parameterToStorageMap.insert(
      {halfWidthShearZoneName, reinterpret_cast<real*>(halfWidthShearZone)});
  parameterToStorageMap.insert(
      {hydraulicDiffusivityName, reinterpret_cast<real*>(hydraulicDiffusivity)});
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
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
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
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
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
