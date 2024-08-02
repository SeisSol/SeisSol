#include "AdjointRSFInitializer.h"

#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/tree/LTSInternalNode.hpp"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Kernels/precision.hpp"
#include <cmath>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <utils/logger.h>

namespace seissol::dr::initializer {
void AdjointRSFInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree);
  auto* concreteLts = dynamic_cast<const seissol::initializer::LTSAdjointRSF* const>(dynRup);

  for (seissol::initializer::LTSTree::leaf_iterator it =
           dynRupTree->beginLeaf(seissol::initializer::LayerMask(Ghost));
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

AdjointRSFInitializer::StateAndFriction
    AdjointRSFInitializer::computeInitialStateAndFriction(real traction1,
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
    logWarning()
        << "Found a negative state variable while initializing the fault. Are you sure your "
           "setup is correct?";
  }
  const double tmp2 = initialSlipRate * 0.5 / rsSr0 *
                      std::exp((rsF0 + rsB * std::log(rsSr0 * result.stateVariable / rsSl0)) / rsA);
  result.frictionCoefficient = rsA * std::asinh(tmp2);
  return result;
}

void AdjointRSFInitializer::addAdditionalParameters(
    std::unordered_map<std::string, real*>& parameterToStorageMap,
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSInternalNode::leaf_iterator& it) {
  auto* concreteLts = dynamic_cast<const seissol::initializer::LTSAdjointRSF* const>(dynRup);
  real(*rsSl0)[misc::numPaddedPoints] = it->var(concreteLts->rsSl0);
  real(*rsA)[misc::numPaddedPoints] = it->var(concreteLts->rsA);
  parameterToStorageMap.insert({"rs_sl0", (real*)rsSl0});
  parameterToStorageMap.insert({"rs_a", (real*)rsA});
}
} // namespace seissol::dr::initializer
