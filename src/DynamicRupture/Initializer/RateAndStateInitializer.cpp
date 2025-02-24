// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
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
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include <Initializer/Parameters/DRParameters.h>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <string>
#include <unordered_map>
#include <utils/logger.h>

namespace seissol::dr::initializer {
void RateAndStateInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree);
  const auto* concreteLts = dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup);

  for (auto& layer : dynRupTree->leaves(Ghost)) {

    bool(*dynStressTimePending)[misc::NumPaddedPoints] =
        layer.var(concreteLts->dynStressTimePending);
    real(*slipRate1)[misc::NumPaddedPoints] = layer.var(concreteLts->slipRate1);
    real(*slipRate2)[misc::NumPaddedPoints] = layer.var(concreteLts->slipRate2);
    real(*mu)[misc::NumPaddedPoints] = layer.var(concreteLts->mu);

    real(*stateVariable)[misc::NumPaddedPoints] = layer.var(concreteLts->stateVariable);
    real(*rsSl0)[misc::NumPaddedPoints] = layer.var(concreteLts->rsSl0);
    real(*rsA)[misc::NumPaddedPoints] = layer.var(concreteLts->rsA);
    real(*initialStressInFaultCS)[misc::NumPaddedPoints][6] =
        layer.var(concreteLts->initialStressInFaultCS);

    const real initialSlipRate =
        misc::magnitude(drParameters->rsInitialSlipRate1, drParameters->rsInitialSlipRate2);

    using namespace dr::misc::quantity_indices;
    for (unsigned ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = true;
        slipRate1[ltsFace][pointIndex] = drParameters->rsInitialSlipRate1;
        slipRate2[ltsFace][pointIndex] = drParameters->rsInitialSlipRate2;
        // compute initial friction and state
        auto stateAndFriction =
            computeInitialStateAndFriction(initialStressInFaultCS[ltsFace][pointIndex][XY],
                                           initialStressInFaultCS[ltsFace][pointIndex][XZ],
                                           initialStressInFaultCS[ltsFace][pointIndex][XX],
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
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  const auto* concreteLts = dynamic_cast<const seissol::initializer::LTSRateAndState*>(dynRup);
  real(*rsSl0)[misc::NumPaddedPoints] = layer.var(concreteLts->rsSl0);
  real(*rsA)[misc::NumPaddedPoints] = layer.var(concreteLts->rsA);
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
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, layer);
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSRateAndStateFastVelocityWeakening*>(dynRup);
  real(*rsSrW)[misc::NumPaddedPoints] = layer.var(concreteLts->rsSrW);
  parameterToStorageMap.insert({"rs_srW", reinterpret_cast<real*>(rsSrW)});
}

ThermalPressurizationInitializer::ThermalPressurizationInitializer(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters)
    : drParameters(drParameters) {}

void ThermalPressurizationInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::ThermalPressurization*>(dynRup);

  for (auto& layer : dynRupTree->leaves(Ghost)) {
    real(*temperature)[misc::NumPaddedPoints] = layer.var(concreteLts->temperature);
    real(*pressure)[misc::NumPaddedPoints] = layer.var(concreteLts->pressure);
    real(*theta)[misc::NumPaddedPoints][misc::NumTpGridPoints] = layer.var(concreteLts->theta);
    real(*sigma)[misc::NumPaddedPoints][misc::NumTpGridPoints] = layer.var(concreteLts->sigma);
    real(*thetaTmpBuffer)[misc::NumPaddedPoints][misc::NumTpGridPoints] =
        layer.var(concreteLts->thetaTmpBuffer);
    real(*sigmaTmpBuffer)[misc::NumPaddedPoints][misc::NumTpGridPoints] =
        layer.var(concreteLts->sigmaTmpBuffer);

    for (unsigned ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters->initialTemperature;
        pressure[ltsFace][pointIndex] = drParameters->initialPressure;
        for (unsigned tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
             ++tpGridPointIndex) {
          theta[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
          sigma[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
          thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
          sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = 0.0;
        }
      }
    }
  }
}

void ThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::ThermalPressurization*>(dynRup);

  real(*halfWidthShearZone)[misc::NumPaddedPoints] = layer.var(concreteLts->halfWidthShearZone);
  real(*hydraulicDiffusivity)[misc::NumPaddedPoints] = layer.var(concreteLts->hydraulicDiffusivity);
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
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  RateAndStateInitializer::initializeFault(dynRup, dynRupTree);
  ThermalPressurizationInitializer::initializeFault(dynRup, dynRupTree);
}

void RateAndStateThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, layer);
  ThermalPressurizationInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, layer);
}

RateAndStateFastVelocityThermalPressurizationInitializer::
    RateAndStateFastVelocityThermalPressurizationInitializer(
        const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
        SeisSol& instance)
    : RateAndStateFastVelocityInitializer(drParameters, instance),
      ThermalPressurizationInitializer(drParameters) {}

void RateAndStateFastVelocityThermalPressurizationInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  RateAndStateFastVelocityInitializer::initializeFault(dynRup, dynRupTree);
  ThermalPressurizationInitializer::initializeFault(dynRup, dynRupTree);
}

void RateAndStateFastVelocityThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::Layer& layer) {
  RateAndStateFastVelocityInitializer::addAdditionalParameters(
      parameterToStorageMap, dynRup, layer);
  ThermalPressurizationInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, layer);
}

} // namespace seissol::dr::initializer
