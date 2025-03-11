// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LinearSlipWeakeningInitializer.h"

#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include <algorithm>
#include <limits>
#include <unordered_map>

namespace seissol::dr::initializer {

void LinearSlipWeakeningInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  BaseDRInitializer::initializeFault(dynRup, dynRupTree);

  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening*>(dynRup);
  for (auto& layer : dynRupTree->leaves(Ghost)) {
    bool(*dynStressTimePending)[misc::NumPaddedPoints] =
        layer.var(concreteLts->dynStressTimePending);
    real(*slipRate1)[misc::NumPaddedPoints] = layer.var(concreteLts->slipRate1);
    real(*slipRate2)[misc::NumPaddedPoints] = layer.var(concreteLts->slipRate2);
    real(*mu)[misc::NumPaddedPoints] = layer.var(concreteLts->mu);
    real(*muS)[misc::NumPaddedPoints] = layer.var(concreteLts->muS);
    real(*forcedRuptureTime)[misc::NumPaddedPoints] = layer.var(concreteLts->forcedRuptureTime);
    const bool providesForcedRuptureTime = this->faultProvides("forced_rupture_time");
    for (unsigned ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
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
    seissol::initializer::Layer& layer) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSLinearSlipWeakening*>(dynRup);
  real(*dC)[misc::NumPaddedPoints] = layer.var(concreteLts->dC);
  real(*muS)[misc::NumPaddedPoints] = layer.var(concreteLts->muS);
  real(*muD)[misc::NumPaddedPoints] = layer.var(concreteLts->muD);
  real(*cohesion)[misc::NumPaddedPoints] = layer.var(concreteLts->cohesion);
  parameterToStorageMap.insert({"d_c", reinterpret_cast<real*>(dC)});
  parameterToStorageMap.insert({"mu_s", reinterpret_cast<real*>(muS)});
  parameterToStorageMap.insert({"mu_d", reinterpret_cast<real*>(muD)});
  parameterToStorageMap.insert({"cohesion", reinterpret_cast<real*>(cohesion)});
  if (this->faultProvides("forced_rupture_time")) {
    real(*forcedRuptureTime)[misc::NumPaddedPoints] = layer.var(concreteLts->forcedRuptureTime);
    parameterToStorageMap.insert(
        {"forced_rupture_time", reinterpret_cast<real*>(forcedRuptureTime)});
  }
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(
    const seissol::initializer::DynamicRupture* const dynRup,
    seissol::initializer::LTSTree* const dynRupTree) {
  LinearSlipWeakeningInitializer::initializeFault(dynRup, dynRupTree);
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::LTSLinearSlipWeakeningBimaterial*>(dynRup);

  for (auto& layer : dynRupTree->leaves(Ghost)) {
    real(*regularizedStrength)[misc::NumPaddedPoints] = layer.var(concreteLts->regularizedStrength);
    real(*mu)[misc::NumPaddedPoints] = layer.var(concreteLts->mu);
    real(*cohesion)[misc::NumPaddedPoints] = layer.var(concreteLts->cohesion);
    real(*initialStressInFaultCS)[misc::NumPaddedPoints][6] =
        layer.var(concreteLts->initialStressInFaultCS);

    for (unsigned ltsFace = 0; ltsFace < layer.getNumberOfCells(); ++ltsFace) {
      for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        regularizedStrength[ltsFace][pointIndex] =
            -cohesion[ltsFace][pointIndex] -
            mu[ltsFace][pointIndex] *
                std::min(static_cast<real>(0.0), initialStressInFaultCS[ltsFace][pointIndex][0]);
      }
    }
  }
}
} // namespace seissol::dr::initializer
