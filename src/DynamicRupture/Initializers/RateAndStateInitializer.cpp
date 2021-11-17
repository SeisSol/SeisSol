#include "RateAndStateInitializer.h"

namespace seissol::dr::initializers {
void RateAndStateInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                              seissol::initializers::LTSTree* dynRupTree,
                                              seissol::Interoperability* e_interoperability) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree, e_interoperability);
  auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    bool(*ds)[numPaddedPoints] = it->var(concreteLts->ds);
    real* averagedSlip = it->var(concreteLts->averagedSlip);
    real(*slipRateStrike)[numPaddedPoints] = it->var(concreteLts->slipRateStrike);
    real(*slipRateDip)[numPaddedPoints] = it->var(concreteLts->slipRateDip);
    real(*mu)[numPaddedPoints] = it->var(concreteLts->mu);

    real(*stateVariable)[numPaddedPoints] = it->var(concreteLts->stateVariable);
    real(*rs_sl0)[numPaddedPoints] = it->var(concreteLts->rs_sl0);
    real(*rs_a)[numPaddedPoints] = it->var(concreteLts->rs_a);
    real(*initialStressInFaultCS)[numPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    real initialSlipRate = std::sqrt(std::pow(drParameters.rs_initialSlipRate1, 2) +
                                     std::pow(drParameters.rs_initialSlipRate2, 2));

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        ds[ltsFace][pointIndex] = drParameters.isDsOutputOn;
        slipRateStrike[ltsFace][pointIndex] = drParameters.rs_initialSlipRate1;
        slipRateDip[ltsFace][pointIndex] = drParameters.rs_initialSlipRate2;
        // compute initial friction and state
        std::tie(stateVariable[ltsFace][pointIndex], mu[ltsFace][pointIndex]) =
            computeInitialStateAndFriction(initialStressInFaultCS[ltsFace][pointIndex][3],
                                           initialStressInFaultCS[ltsFace][pointIndex][5],
                                           initialStressInFaultCS[ltsFace][pointIndex][0],
                                           rs_a[ltsFace][pointIndex],
                                           drParameters.rs_b,
                                           rs_sl0[ltsFace][pointIndex],
                                           drParameters.rs_sr0,
                                           drParameters.rs_f0,
                                           initialSlipRate);
      }
      averagedSlip[ltsFace] = 0.0;
    }
    // can be removed once output is in c++
    for (unsigned int ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      e_interoperability->copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRateStrike, slipRateDip, mu);
      e_interoperability->copyFrictionOutputToFortranStateVar(ltsFace, meshFace, stateVariable);
    }
  }
}

std::pair<real, real>
    RateAndStateInitializer::computeInitialStateAndFriction(real tractionXY,
                                                            real tractionXZ,
                                                            real pressure,
                                                            real rs_a,
                                                            real rs_b,
                                                            real rs_sl0,
                                                            real rs_sr0,
                                                            real rs_f0,
                                                            real initialSlipRate) {
  real absoluteTraction = std::sqrt(std::pow(tractionXY, 2) + std::pow(tractionXZ, 2));
  real tmp = std::abs(absoluteTraction / (rs_a * pressure));
  real stateVariable = rs_sl0 / rs_sr0 *
                       std::exp((rs_a * std::log(std::exp(tmp) - std::exp(-tmp)) - rs_f0 -
                                 rs_a * std::log(initialSlipRate / rs_sr0)) /
                                rs_b);
  real tmp2 = initialSlipRate * 0.5 / rs_sr0 *
              std::exp((rs_f0 + rs_b * std::log(rs_sr0 * stateVariable / rs_sl0)) / rs_a);
  real mu = rs_a * std::asinh(tmp2);
  return {stateVariable, mu};
}

void RateAndStateInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);
  real(*rs_sl0)[numPaddedPoints] = it->var(concreteLts->rs_sl0);
  real(*rs_a)[numPaddedPoints] = it->var(concreteLts->rs_a);
  parameterToStorageMap.insert({"sl0", (real*)rs_sl0});
  parameterToStorageMap.insert({"rs_a", (real*)rs_a});
}

std::pair<real, real>
    RateAndStateFastVelocityInitializer::computeInitialStateAndFriction(real tractionXY,
                                                                        real tractionXZ,
                                                                        real pressure,
                                                                        real rs_a,
                                                                        real rs_b,
                                                                        real rs_sl0,
                                                                        real rs_sr0,
                                                                        real rs_f0,
                                                                        real initialSlipRate) {
  real absoluteTraction = std::sqrt(std::pow(tractionXY, 2) + std::pow(tractionXZ, 2));
  real tmp = std::abs(absoluteTraction / (rs_a * pressure));
  real stateVariable =
      rs_a * std::log(2.0 * rs_sr0 / initialSlipRate * (std::exp(tmp) - std::exp(-tmp)) / 2.0);
  real tmp2 = initialSlipRate * 0.5 / rs_sr0 * std::exp(stateVariable / rs_a);
  // asinh(x)=log(x+sqrt(x^2+1))
  real mu = rs_a * std::log(tmp2 + std::sqrt(std::pow(tmp2, 2) + 1.0));
  return {stateVariable, mu};
}

void RateAndStateFastVelocityInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);
  real(*rs_srW)[numPaddedPoints] = it->var(concreteLts->rs_srW);
  parameterToStorageMap.insert({"rs_srW", (real*)rs_srW});
}

void RateAndStateThermalPressurisationInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* e_interoperability) {
  RateAndStateInitializer::initializeFault(dynRup, dynRupTree, e_interoperability);

  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*temperature)[numPaddedPoints] = it->var(concreteLts->temperature);
    real(*pressure)[numPaddedPoints] = it->var(concreteLts->pressure);
    real(*TP_Theta)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_theta);
    real(*TP_sigma)[numPaddedPoints][TP_grid_nz] = it->var(concreteLts->TP_sigma);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        temperature[ltsFace][pointIndex] = drParameters.initialTemperature;
        pressure[ltsFace][pointIndex] = drParameters.iniPressure;
        for (unsigned iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; ++iTP_grid_nz) {
          TP_Theta[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
          TP_sigma[ltsFace][pointIndex][iTP_grid_nz] = 0.0;
        }
      }
    }
  }
}

void RateAndStateThermalPressurisationInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  RateAndStateFastVelocityInitializer::addAdditionalParameters(parameterToStorageMap, dynRup, it);

  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);

  real(*TP_halfWidthShearZone)[numPaddedPoints] = it->var(concreteLts->TP_halfWidthShearZone);
  real(*alphaHy)[numPaddedPoints] = it->var(concreteLts->alphaHy);
  parameterToStorageMap.insert({"TP_halfWidthShearZone", (real*)TP_halfWidthShearZone});
  parameterToStorageMap.insert({"alphaHy", (real*)alphaHy});
}
} // namespace seissol::dr::initializers