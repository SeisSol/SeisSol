#include "LinearSlipWeakeningInitializer.h"

#include "utils/logger.h"

namespace seissol::dr::initializers {

void LinearSlipWeakeningInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* e_interoperability) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree, e_interoperability);

  auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    bool(*DS)[numPaddedPoints] = it->var(concreteLts->ds);
    real* averagedSlip = it->var(concreteLts->averagedSlip);
    real(*dynStressTime)[numPaddedPoints] = it->var(concreteLts->dynStressTime);
    real(*slipRateStrike)[numPaddedPoints] = it->var(concreteLts->slipRateStrike);
    real(*slipRateDip)[numPaddedPoints] = it->var(concreteLts->slipRateDip);
    real(*mu)[numPaddedPoints] = it->var(concreteLts->mu);
    real(*mu_s)[numPaddedPoints] = it->var(concreteLts->mu_s);
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);

      // initialize padded elements for vectorization
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        dynStressTime[ltsFace][pointIndex] = 0.0;
        DS[ltsFace][pointIndex] = drParameters.isDsOutputOn;
        slipRateStrike[ltsFace][pointIndex] = 0.0;
        slipRateDip[ltsFace][pointIndex] = 0.0;
        // initial friction coefficient is static friction (no slip has yet occurred)
        mu[ltsFace][pointIndex] = mu_s[ltsFace][pointIndex];
      }
      averagedSlip[ltsFace] = 0.0;
      // can be removed once output is in c++
      e_interoperability->copyFrictionOutputToFortranFL2(
          ltsFace, meshFace, averagedSlip, dynStressTime, slipRateStrike, slipRateDip, mu);
    }
  }
}

void LinearSlipWeakeningInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
  real(*d_c)[numPaddedPoints] = it->var(concreteLts->d_c);
  real(*mu_s)[numPaddedPoints] = it->var(concreteLts->mu_s);
  real(*mu_d)[numPaddedPoints] = it->var(concreteLts->mu_d);
  real(*cohesion)[numPaddedPoints] = it->var(concreteLts->cohesion);
  parameterToStorageMap.insert({"d_c", (real*)d_c});
  parameterToStorageMap.insert({"mu_s", (real*)mu_s});
  parameterToStorageMap.insert({"mu_d", (real*)mu_d});
  parameterToStorageMap.insert({"cohesion", (real*)cohesion});
}

void LinearSlipWeakeningForcedRuptureTimeInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* e_interoperability) {
  LinearSlipWeakeningInitializer::initializeFault(dynRup, dynRupTree, e_interoperability);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);
  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real* tn = it->var(concreteLts->tn);
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      tn[ltsFace] = 0.0;
    }
  }
}

void LinearSlipWeakeningForcedRuptureTimeInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);
  real(*forcedRuptureTime)[numPaddedPoints] = it->var(concreteLts->forcedRuptureTime);
  parameterToStorageMap.insert({"forced_rupture_time", (real*)forcedRuptureTime});
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* e_interoperability) {
  LinearSlipWeakeningInitializer::initializeFault(dynRup, dynRupTree, e_interoperability);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*regularisedStrength)[numPaddedPoints] = it->var(concreteLts->regularisedStrength);
    real(*mu)[numPaddedPoints] = it->var(concreteLts->mu);
    real(*initialStressInFaultCS)[numPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      // unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        regularisedStrength[ltsFace][pointIndex] =
            mu[ltsFace][pointIndex] * initialStressInFaultCS[ltsFace][pointIndex][0];
      }
      // can be removed once output is in c++
      e_interoperability->copyFrictionOutputToFortranStrength(
          ltsFace, meshFace, regularisedStrength);
    }
  }
}
} // namespace seissol::dr::initializers
