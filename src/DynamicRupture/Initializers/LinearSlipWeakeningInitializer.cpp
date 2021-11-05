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
    bool(*DS)[numPaddedPoints] = it->var(concreteLts->DS);
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

void seissol::dr::initializers::LinearSlipWeakeningInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {}

void LinearSlipWeakeningInitializer::addAdditionalParameters(
    std::map<std::string, double*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto concreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakening*>(dynRup);
  real(*d_c)[numPaddedPoints] = it->var(concreteLts->d_c);
  real(*mu_s)[numPaddedPoints] = it->var(concreteLts->mu_s);
  real(*mu_d)[numPaddedPoints] = it->var(concreteLts->mu_d);
  parameterToStorageMap.insert({"d_c", (double*)d_c});
  parameterToStorageMap.insert({"mu_s", (double*)mu_s});
  parameterToStorageMap.insert({"mu_d", (double*)mu_d});
}

void LinearSlipWeakeningForcedRuptureTimeInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  LinearSlipWeakeningInitializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {

    real(*forcedRuptureTime)[numPaddedPoints] =
        it->var(concreteLts->forcedRuptureTime); // from faultParameters
    real* tn = it->var(concreteLts->tn);         // = 0

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

      // initialize padded elements for vectorization
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        forcedRuptureTime[ltsFace][pointIndex] = 0.0;
      }
      for (unsigned pointIndex = 0; pointIndex < numberOfPoints; ++pointIndex) {
        if (faultParameters["forced_rupture_time"] != NULL) {
          forcedRuptureTime[ltsFace][pointIndex] = static_cast<real>(
              faultParameters["forced_rupture_time"][(meshFace)*numberOfPoints + pointIndex]);
        }
      }
      tn[ltsFace] = 0.0;
    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}

void LinearSlipWeakeningBimaterialInitializer::initializeFrictionMatrices(
    seissol::initializers::DynamicRupture* dynRup,
    seissol::initializers::LTSTree* dynRupTree,
    seissol::dr::friction_law::BaseFrictionLaw* FrictionLaw,
    std::unordered_map<std::string, double*> faultParameters,
    unsigned* ltsFaceToMeshFace,
    seissol::Interoperability& e_interoperability) {
  LinearSlipWeakeningInitializer::initializeFrictionMatrices(
      dynRup, dynRupTree, FrictionLaw, faultParameters, ltsFaceToMeshFace, e_interoperability);
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);

  unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*regularisedStrength)[numPaddedPoints] = it->var(concreteLts->regularisedStrength);
    real(*mu)[numPaddedPoints] = it->var(concreteLts->mu);
    real(*initialStressInFaultCS)[numPaddedPoints][6] = it->var(dynRup->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      // unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
      for (unsigned pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
        regularisedStrength[ltsFace][pointIndex] =
            mu[ltsFace][pointIndex] * initialStressInFaultCS[ltsFace][pointIndex][0];
      }
    } // lts-face loop
    layerLtsFaceToMeshFace += it->getNumberOfCells();
  } // leaf_iterator loop
}
} // namespace seissol::dr::initializers
