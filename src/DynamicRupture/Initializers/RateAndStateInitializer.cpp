#include "RateAndStateInitializer.h"

#include "DynamicRupture/Misc.h"

namespace seissol::dr::initializers {
void RateAndStateInitializer::initializeFault(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree);
  auto* concreteLts = dynamic_cast<seissol::initializers::LTSRateAndState const* const>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    bool(*dynStressTimePending)[misc::numPaddedPoints] = it->var(concreteLts->dynStressTimePending);
    real(*slipRate1)[misc::numPaddedPoints] = it->var(concreteLts->slipRate1);
    real(*slipRate2)[misc::numPaddedPoints] = it->var(concreteLts->slipRate2);
    real(*mu)[misc::numPaddedPoints] = it->var(concreteLts->mu);

    real(*stateVariable)[misc::numPaddedPoints] = it->var(concreteLts->stateVariable);
    real(*rsSl0)[misc::numPaddedPoints] = it->var(concreteLts->rsSl0);
    real(*rsA)[misc::numPaddedPoints] = it->var(concreteLts->rsA);
    real(*initialStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    const real initialSlipRate =
        misc::magnitude(drParameters->rsInitialSlipRate1, drParameters->rsInitialSlipRate2);

    using namespace dr::misc::quantity_indices;
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = drParameters->isDsOutputOn;
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
  StateAndFriction result;
  const double absoluteTraction = misc::magnitude(traction1, traction2);
  const double tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable = rsSl0 / rsSr0 *
                         std::exp((rsA * std::log(std::exp(tmp) - std::exp(-tmp)) - rsF0 -
                                   rsA * std::log(initialSlipRate / rsSr0)) /
                                  rsB);
  if (result.stateVariable < 0) {
    logError() << "Found a negative state variable while initializing the fault. Are you sure your "
                  "setup is correct?";
  }
  const double tmp2 = initialSlipRate * 0.5 / rsSr0 *
                      std::exp((rsF0 + rsB * std::log(rsSr0 * result.stateVariable / rsSl0)) / rsA);
  result.frictionCoefficient = rsA * std::asinh(tmp2);
  return result;
}

void RateAndStateInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTSRateAndState const* const>(dynRup);
  real(*rsSl0)[misc::numPaddedPoints] = it->var(concreteLts->rsSl0);
  real(*rsA)[misc::numPaddedPoints] = it->var(concreteLts->rsA);
  parameterToStorageMap.insert({"rs_sl0", (real*)rsSl0});
  parameterToStorageMap.insert({"rs_a", (real*)rsA});
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
  StateAndFriction result;
  const real absoluteTraction = misc::magnitude(traction1, traction2);
  const real tmp = std::abs(absoluteTraction / (rsA * pressure));
  result.stateVariable =
      rsA * std::log(2.0 * rsSr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  if (result.stateVariable < 0) {
    logError() << "Found a negative state variable while initializing the fault. Are you sure your "
                  "setup is correct?";
  }
  const real tmp2 = initialSlipRate * 0.5 / rsSr0 * std::exp(result.stateVariable / rsA);
  result.frictionCoefficient = rsA * misc::asinh(tmp2);
  return result;
}

void RateAndStateFastVelocityInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSRateAndStateFastVelocityWeakening const* const>(
          dynRup);
  real(*rsSrW)[misc::numPaddedPoints] = it->var(concreteLts->rsSrW);
  parameterToStorageMap.insert({"rs_srW", (real*)rsSrW});
}

void RateAndStateThermalPressurizationInitializer::initializeFault(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  RateAndStateInitializer::initializeFault(dynRup, dynRupTree);

  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSRateAndStateThermalPressurization const* const>(
          dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*temperature)[misc::numPaddedPoints] = it->var(concreteLts->temperature);
    real(*pressure)[misc::numPaddedPoints] = it->var(concreteLts->pressure);
    real(*theta)[misc::numPaddedPoints][misc::numberOfTPGridPoints] = it->var(concreteLts->theta);
    real(*sigma)[misc::numPaddedPoints][misc::numberOfTPGridPoints] = it->var(concreteLts->sigma);
    real(*thetaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints] =
        it->var(concreteLts->thetaTmpBuffer);
    real(*sigmaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints] =
        it->var(concreteLts->sigmaTmpBuffer);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters->initialTemperature;
        pressure[ltsFace][pointIndex] = drParameters->initialPressure;
        for (unsigned tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
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

void RateAndStateThermalPressurizationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateFastVelocityInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);

  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSRateAndStateThermalPressurization const* const>(
          dynRup);

  real(*halfWidthShearZone)[misc::numPaddedPoints] = it->var(concreteLts->halfWidthShearZone);
  real(*hydraulicDiffusivity)[misc::numPaddedPoints] = it->var(concreteLts->hydraulicDiffusivity);
  parameterToStorageMap.insert({"tp_halfWidthShearZone", (real*)halfWidthShearZone});
  parameterToStorageMap.insert({"tp_hydraulicDiffusivity", (real*)hydraulicDiffusivity});
}
} // namespace seissol::dr::initializers
