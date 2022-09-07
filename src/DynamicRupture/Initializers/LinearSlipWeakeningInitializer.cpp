#include "LinearSlipWeakeningInitializer.h"

#include "utils/logger.h"

namespace seissol::dr::initializers {

void LinearSlipWeakeningInitializer::initializeFault(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree);

  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSLinearSlipWeakening const* const>(dynRup);
  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    bool(*dynStressTimePending)[misc::numPaddedPoints] = it->var(concreteLts->dynStressTimePending);
    real(*slipRate1)[misc::numPaddedPoints] = it->var(concreteLts->slipRate1);
    real(*slipRate2)[misc::numPaddedPoints] = it->var(concreteLts->slipRate2);
    real(*mu)[misc::numPaddedPoints] = it->var(concreteLts->mu);
    real(*muS)[misc::numPaddedPoints] = it->var(concreteLts->muS);
    real(*forcedRuptureTime)[misc::numPaddedPoints] = it->var(concreteLts->forcedRuptureTime);
    const bool providesForcedRuptureTime = this->faultProvides("forced_rupture_time");
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      // initialize padded elements for vectorization
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = drParameters->isDsOutputOn;
        slipRate1[ltsFace][pointIndex] = 0.0;
        slipRate2[ltsFace][pointIndex] = 0.0;
        // initial friction coefficient is static friction (no slip has yet occurred)
        mu[ltsFace][pointIndex] = muS[ltsFace][pointIndex];
        if (!providesForcedRuptureTime) {
          forcedRuptureTime[ltsFace][pointIndex] = std::numeric_limits<real>::max();
        }
      }
    }
  }
}

void LinearSlipWeakeningInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSLinearSlipWeakening const* const>(dynRup);
  real(*dC)[misc::numPaddedPoints] = it->var(concreteLts->dC);
  real(*muS)[misc::numPaddedPoints] = it->var(concreteLts->muS);
  real(*muD)[misc::numPaddedPoints] = it->var(concreteLts->muD);
  real(*cohesion)[misc::numPaddedPoints] = it->var(concreteLts->cohesion);
  parameterToStorageMap.insert({"d_c", (real*)dC});
  parameterToStorageMap.insert({"mu_s", (real*)muS});
  parameterToStorageMap.insert({"mu_d", (real*)muD});
  parameterToStorageMap.insert({"cohesion", (real*)cohesion});
  if (this->faultProvides("forced_rupture_time")) {
    real(*forcedRuptureTime)[misc::numPaddedPoints] = it->var(concreteLts->forcedRuptureTime);
    parameterToStorageMap.insert({"forced_rupture_time", (real*)forcedRuptureTime});
  }
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(
    seissol::initializers::DynamicRupture const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  LinearSlipWeakeningInitializer::initializeFault(dynRup, dynRupTree);
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSLinearSlipWeakeningBimaterial const* const>(dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*regularisedStrength)[misc::numPaddedPoints] = it->var(concreteLts->regularisedStrength);
    real(*mu)[misc::numPaddedPoints] = it->var(concreteLts->mu);
    real(*cohesion)[misc::numPaddedPoints] = it->var(concreteLts->cohesion);
    real(*initialStressInFaultCS)[misc::numPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
        regularisedStrength[ltsFace][pointIndex] =
            -cohesion[ltsFace][pointIndex] -
            mu[ltsFace][pointIndex] *
                std::min(static_cast<real>(0.0), initialStressInFaultCS[ltsFace][pointIndex][0]);
      }
    }
  }
}
} // namespace seissol::dr::initializers
