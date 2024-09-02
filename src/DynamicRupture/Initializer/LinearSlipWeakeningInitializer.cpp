#include "LinearSlipWeakeningInitializer.h"

#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/tree/LTSInternalNode.hpp"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Kernels/precision.hpp"
#include <algorithm>
#include <limits>
#include <unordered_map>

namespace seissol::dr::initializer {

void LinearSlipWeakeningInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree);

  auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening* const>(dynRup);
  for (auto it = dynRupTree->beginLeaf(seissol::initializer::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    bool(*dynStressTimePending)[misc::NumPaddedPoints] = it->var(concreteLts->dynStressTimePending);
    real(*slipRate1)[misc::NumPaddedPoints] = it->var(concreteLts->slipRate1);
    real(*slipRate2)[misc::NumPaddedPoints] = it->var(concreteLts->slipRate2);
    real(*mu)[misc::NumPaddedPoints] = it->var(concreteLts->mu);
    real(*muS)[misc::NumPaddedPoints] = it->var(concreteLts->muS);
    real(*forcedRuptureTime)[misc::NumPaddedPoints] = it->var(concreteLts->forcedRuptureTime);
    const bool providesForcedRuptureTime = this->faultProvides("forced_rupture_time");
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      // initialize padded elements for vectorization
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = true;
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
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSInternalNode::LeafIterator& it) {
  auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening* const>(dynRup);
  real(*dC)[misc::NumPaddedPoints] = it->var(concreteLts->dC);
  real(*muS)[misc::NumPaddedPoints] = it->var(concreteLts->muS);
  real(*muD)[misc::NumPaddedPoints] = it->var(concreteLts->muD);
  real(*cohesion)[misc::NumPaddedPoints] = it->var(concreteLts->cohesion);
  parameterToStorageMap.insert({"d_c", (real*)dC});
  parameterToStorageMap.insert({"mu_s", (real*)muS});
  parameterToStorageMap.insert({"mu_d", (real*)muD});
  parameterToStorageMap.insert({"cohesion", (real*)cohesion});
  if (this->faultProvides("forced_rupture_time")) {
    real(*forcedRuptureTime)[misc::NumPaddedPoints] = it->var(concreteLts->forcedRuptureTime);
    parameterToStorageMap.insert({"forced_rupture_time", (real*)forcedRuptureTime});
  }
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  LinearSlipWeakeningInitializer::initializeFault(dynRup, dynRupTree);
  auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSLinearSlipWeakeningBimaterial* const>(dynRup);

  for (auto it = dynRupTree->beginLeaf(seissol::initializer::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    real(*regularisedStrength)[misc::NumPaddedPoints] = it->var(concreteLts->regularisedStrength);
    real(*mu)[misc::NumPaddedPoints] = it->var(concreteLts->mu);
    real(*cohesion)[misc::NumPaddedPoints] = it->var(concreteLts->cohesion);
    real(*initialStressInFaultCS)[misc::NumPaddedPoints][6] =
        it->var(concreteLts->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        regularisedStrength[ltsFace][pointIndex] =
            -cohesion[ltsFace][pointIndex] -
            mu[ltsFace][pointIndex] *
                std::min(static_cast<real>(0.0), initialStressInFaultCS[ltsFace][pointIndex][0]);
      }
    }
  }
}
} // namespace seissol::dr::initializer
