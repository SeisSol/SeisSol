#include "LinearSlipWeakeningInitializer.h"

#include "utils/logger.h"

namespace seissol::dr::initializers {

template <typename Config>
void LinearSlipWeakeningInitializer<Config>::initializeFault(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  BaseDRInitializer<Config>::initializeFault(dynRup, dynRupTree);

  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSLinearSlipWeakening<Config> const* const>(dynRup);
  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    bool(*dynStressTimePending)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->dynStressTimePending);
    typename Config::RealT(*slipRate1)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->slipRate1);
    typename Config::RealT(*slipRate2)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->slipRate2);
    typename Config::RealT(*mu)[misc::numPaddedPoints<Config>] = it->var(concreteLts->mu);
    typename Config::RealT(*muS)[misc::numPaddedPoints<Config>] = it->var(concreteLts->muS);
    typename Config::RealT(*forcedRuptureTime)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->forcedRuptureTime);
    const bool providesForcedRuptureTime = this->faultProvides("forced_rupture_time");
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      // initialize padded elements for vectorization
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; ++pointIndex) {
        dynStressTimePending[ltsFace][pointIndex] = drParameters->isDsOutputOn;
        slipRate1[ltsFace][pointIndex] = 0.0;
        slipRate2[ltsFace][pointIndex] = 0.0;
        // initial friction coefficient is static friction (no slip has yet occurred)
        mu[ltsFace][pointIndex] = muS[ltsFace][pointIndex];
        if (!providesForcedRuptureTime) {
          forcedRuptureTime[ltsFace][pointIndex] =
              std::numeric_limits<typename Config::RealT>::max();
        }
      }
    }
  }
}

template <typename Config>
void LinearSlipWeakeningInitializer<Config>::addAdditionalParameters(
    std::unordered_map<std::string, typename Config::RealT*>& parameterToStorageMap,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSLinearSlipWeakening<Config> const* const>(dynRup);
  typename Config::RealT(*dC)[misc::numPaddedPoints<Config>] = it->var(concreteLts->dC);
  typename Config::RealT(*muS)[misc::numPaddedPoints<Config>] = it->var(concreteLts->muS);
  typename Config::RealT(*muD)[misc::numPaddedPoints<Config>] = it->var(concreteLts->muD);
  typename Config::RealT(*cohesion)[misc::numPaddedPoints<Config>] = it->var(concreteLts->cohesion);
  parameterToStorageMap.insert({"d_c", (typename Config::RealT*)dC});
  parameterToStorageMap.insert({"mu_s", (typename Config::RealT*)muS});
  parameterToStorageMap.insert({"mu_d", (typename Config::RealT*)muD});
  parameterToStorageMap.insert({"cohesion", (typename Config::RealT*)cohesion});
  if (this->faultProvides("forced_rupture_time")) {
    typename Config::RealT(*forcedRuptureTime)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->forcedRuptureTime);
    parameterToStorageMap.insert(
        {"forced_rupture_time", (typename Config::RealT*)forcedRuptureTime});
  }
}

template <typename Config>
void LinearSlipWeakeningBimaterialInitializer<Config>::initializeFault(
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    seissol::initializers::LTSTree* const dynRupTree) {
  LinearSlipWeakeningInitializer<Config>::initializeFault(dynRup, dynRupTree);
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSLinearSlipWeakeningBimaterial<Config> const* const>(
          dynRup);

  for (seissol::initializers::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
       it != dynRupTree->endLeaf();
       ++it) {
    typename Config::RealT(*regularisedStrength)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->regularisedStrength);
    typename Config::RealT(*mu)[misc::numPaddedPoints<Config>] = it->var(concreteLts->mu);
    typename Config::RealT(*cohesion)[misc::numPaddedPoints<Config>] =
        it->var(concreteLts->cohesion);
    typename Config::RealT(*initialStressInFaultCS)[misc::numPaddedPoints<Config>][6] =
        it->var(concreteLts->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; ++pointIndex) {
        regularisedStrength[ltsFace][pointIndex] =
            -cohesion[ltsFace][pointIndex] -
            mu[ltsFace][pointIndex] * std::min(static_cast<typename Config::RealT>(0.0),
                                               initialStressInFaultCS[ltsFace][pointIndex][0]);
      }
    }
  }
}
} // namespace seissol::dr::initializers
