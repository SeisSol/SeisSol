#include "RateAndStateInitializer.h"

#include "DynamicRupture/Misc.h"

namespace seissol::dr::initializers {
void RateAndStateInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                              seissol::initializers::LTSTree* dynRupTree,
                                              seissol::Interoperability* eInteroperability) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree, eInteroperability);
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    bool(*dynStressTimePending)[misc::numPaddedPoints] = it->var(concreteLts->dynStressTimePending);
    real* averagedSlip = it->var(concreteLts->averagedSlip);
    real(*slipRate1)[misc::numPaddedPoints] = it->var(concreteLts->slipRate1);
    real(*slipRate2)[misc::numPaddedPoints] = it->var(concreteLts->slipRate2);
    real(*mu)[misc::numPaddedPoints] = it->var(concreteLts->mu);

    real(*stateVariable)[misc::numPaddedPoints] = it->var(concreteLts->stateVariable);
    real(*rsSl0)[misc::numPaddedPoints] = it->var(concreteLts->rsSl0);
    real(*rsA)[misc::numPaddedPoints] = it->var(concreteLts->rsA);
    real(*initialStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    real initialSlipRate =
        misc::magnitude(drParameters.rsInitialSlipRate1, drParameters.rsInitialSlipRate2);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = drParameters.isDsOutputOn;
        slipRate1[ltsFace][pointIndex] = drParameters.rsInitialSlipRate1;
        slipRate2[ltsFace][pointIndex] = drParameters.rsInitialSlipRate2;
        // compute initial friction and state
        std::tie(stateVariable[ltsFace][pointIndex], mu[ltsFace][pointIndex]) =
            computeInitialStateAndFriction(initialStressInFaultCS[ltsFace][pointIndex][3],
                                           initialStressInFaultCS[ltsFace][pointIndex][5],
                                           initialStressInFaultCS[ltsFace][pointIndex][0],
                                           rsA[ltsFace][pointIndex],
                                           drParameters.rsB,
                                           rsSl0[ltsFace][pointIndex],
                                           drParameters.rsSr0,
                                           drParameters.rsF0,
                                           initialSlipRate);
      }
      averagedSlip[ltsFace] = 0.0;
    }
    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      eInteroperability->copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRate1, slipRate2, mu);
      eInteroperability->copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVariable);
    }
  }
}

std::pair<real, real>
    RateAndStateInitializer::computeInitialStateAndFriction(real tractionXY,
                                                            real tractionXZ,
                                                            real pressure,
                                                            real rsA,
                                                            real rsB,
                                                            real rsSl0,
                                                            real rsSr0,
                                                            real rsF0,
                                                            real initialSlipRate) {
  real absoluteTraction = misc::magnitude(tractionXY, tractionXZ);
  real tmp = std::abs(absoluteTraction / (rsA * pressure));
  real stateVariable = rsSl0 / rsSr0 *
                       std::exp((rsA * std::log(std::exp(tmp) - std::exp(-tmp)) - rsF0 -
                                 rsA * std::log(initialSlipRate / rsSr0)) /
                                rsB);
  real tmp2 = initialSlipRate * 0.5 / rsSr0 *
              std::exp((rsF0 + rsB * std::log(rsSr0 * stateVariable / rsSl0)) / rsA);
  real mu = rsA * std::asinh(tmp2);
  return {stateVariable, mu};
}

void RateAndStateInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);
  real(*rsSl0)[misc::numPaddedPoints] = it->var(concreteLts->rsSl0);
  real(*rsA)[misc::numPaddedPoints] = it->var(concreteLts->rsA);
  parameterToStorageMap.insert({"rs_sl0", (real*)rsSl0});
  parameterToStorageMap.insert({"rs_a", (real*)rsA});
}

std::pair<real, real>
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(real tractionXY,
                                                                        real tractionXZ,
                                                                        real pressure,
                                                                        real rsA,
                                                                        real rsB,
                                                                        real rsSl0,
                                                                        real rsSr0,
                                                                        real rsF0,
                                                                        real initialSlipRate) {
  real absoluteTraction = misc::magnitude(tractionXY, tractionXZ);
  real tmp = std::abs(absoluteTraction / (rsA * pressure));
  real stateVariable =
      rsA * std::log(2.0 * rsSr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  real tmp2 = initialSlipRate * 0.5 / rsSr0 * std::exp(stateVariable / rsA);
  real mu = rsA * misc::asinh(tmp2);
  return {stateVariable, mu};
}

void RateAndStateFastVelocityInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);
  real(*rsSrW)[misc::numPaddedPoints] = it->var(concreteLts->rsSrW);
  parameterToStorageMap.insert({"rs_srW", (real*)rsSrW});
}

void RateAndStateThermalPressurizationInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* eInteroperability) {
  RateAndStateInitializer::initializeFault(dynRup, dynRupTree, eInteroperability);

  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurization*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*temperature)[misc::numPaddedPoints] = it->var(concreteLts->temperature);
    real(*pressure)[misc::numPaddedPoints] = it->var(concreteLts->pressure);
    real(*theta)[misc::numPaddedPoints][misc::numberOfTPGridPoints] =
        it->var(concreteLts->theta);
    real(*sigma)[misc::numPaddedPoints][misc::numberOfTPGridPoints] =
        it->var(concreteLts->sigma);
    real(*thetaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints] =
        it->var(concreteLts->thetaTmpBuffer);
    real(*sigmaTmpBuffer)[misc::numPaddedPoints][misc::numberOfTPGridPoints] =
        it->var(concreteLts->sigmaTmpBuffer);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters.initialTemperature;
        pressure[ltsFace][pointIndex] = drParameters.initialPressure;
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
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateFastVelocityInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);

  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurization*>(dynRup);

  real(*halfWidthShearZone)[misc::numPaddedPoints] = it->var(concreteLts->halfWidthShearZone);
  real(*hydraulicDiffusivity)[misc::numPaddedPoints] = it->var(concreteLts->hydraulicDiffusivity);
  parameterToStorageMap.insert({"halfWidthShearZone", (real*)halfWidthShearZone});
  parameterToStorageMap.insert({"hydraulicDiffusivity", (real*)hydraulicDiffusivity});
}
} // namespace seissol::dr::initializers