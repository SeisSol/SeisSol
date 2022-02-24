#include "LinearSlipWeakeningInitializer.h"

#include "utils/logger.h"

namespace seissol::dr::initializers {

void LinearSlipWeakeningInitializer::initializeFault(seissol::initializers::DynamicRupture* dynRup,
                                                     seissol::initializers::LTSTree* dynRupTree,
                                                     seissol::Interoperability* eInteroperability) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree, eInteroperability);

  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    bool(*dynStressTimePending)[misc::numPaddedPoints] = it->var(concreteLts->dynStressTimePending);
    real* averagedSlip = it->var(concreteLts->averagedSlip);
    real(*slipRate1)[misc::numPaddedPoints] = it->var(concreteLts->slipRate1);
    real(*slipRate2)[misc::numPaddedPoints] = it->var(concreteLts->slipRate2);
    real(*mu)[misc::numPaddedPoints] = it->var(concreteLts->mu);
    real(*muS)[misc::numPaddedPoints] = it->var(concreteLts->muS);
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);

      // initialize padded elements for vectorization
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = drParameters.isDsOutputOn;
        slipRate1[ltsFace][pointIndex] = 0.0;
        slipRate2[ltsFace][pointIndex] = 0.0;
        // initial friction coefficient is static friction (no slip has yet occurred)
        mu[ltsFace][pointIndex] = muS[ltsFace][pointIndex];
      }
      averagedSlip[ltsFace] = 0.0;
      // can be removed once output is in c++
      eInteroperability->copyFrictionOutputToFortranSpecific(
          ltsFace, meshFace, averagedSlip, slipRate1, slipRate2, mu);
    }
  }
}

void LinearSlipWeakeningInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
  real(*dC)[misc::numPaddedPoints] = it->var(concreteLts->dC);
  real(*muS)[misc::numPaddedPoints] = it->var(concreteLts->muS);
  real(*muD)[misc::numPaddedPoints] = it->var(concreteLts->muD);
  real(*cohesion)[misc::numPaddedPoints] = it->var(concreteLts->cohesion);
  parameterToStorageMap.insert({"d_c", (real*)dC});
  parameterToStorageMap.insert({"mu_s", (real*)muS});
  parameterToStorageMap.insert({"mu_d", (real*)muD});
  parameterToStorageMap.insert({"cohesion", (real*)cohesion});
}

void LinearSlipWeakeningForcedRuptureTimeInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* eInteroperability) {
  LinearSlipWeakeningInitializer::initializeFault(dynRup, dynRupTree, eInteroperability);
  auto* concreteLts =
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
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);
  real(*forcedRuptureTime)[misc::numPaddedPoints] = it->var(concreteLts->forcedRuptureTime);
  parameterToStorageMap.insert({"forced_rupture_time", (real*)forcedRuptureTime});
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::Interoperability* eInteroperability) {
  LinearSlipWeakeningInitializer::initializeFault(dynRup, dynRupTree, eInteroperability);
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*regularisedStrength)[misc::numPaddedPoints] = it->var(concreteLts->regularisedStrength);
    real(*mu)[misc::numPaddedPoints] = it->var(concreteLts->mu);
    real(*initialStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      const auto& drFaceInformation = it->var(dynRup->faceInformation);
      unsigned meshFace = static_cast<int>(drFaceInformation[ltsFace].meshFace);
      // unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        regularisedStrength[ltsFace][pointIndex] =
            mu[ltsFace][pointIndex] * initialStressInFaultCS[ltsFace][pointIndex][0];
      }
      // can be removed once output is in c++
      eInteroperability->copyFrictionOutputToFortranStrength(
          ltsFace, meshFace, regularisedStrength);
    }
  }
}
} // namespace seissol::dr::initializers
