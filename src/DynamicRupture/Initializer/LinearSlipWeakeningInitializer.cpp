// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LinearSlipWeakeningInitializer.h"

#include "Config.h"
#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <unordered_map>

namespace seissol::dr::initializer {

void LinearSlipWeakeningInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  BaseDRInitializer::initializeFault(drStorage);
  for (auto& layer : drStorage.leaves(Ghost)) {
    bool (*dynStressTimePending)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::DynStressTimePending>(Config());
    real(*slipRate1)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::SlipRate1>(Config());
    real(*slipRate2)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::SlipRate2>(Config());
    real(*mu)[misc::NumPaddedPoints<Config>] = layer.var<LTSLinearSlipWeakening::Mu>(Config());
    const real(*muS)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::MuS>(Config());
    real(*forcedRuptureTime)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::ForcedRuptureTime>(Config());
    const bool providesForcedRuptureTime = this->faultProvides("forced_rupture_time");
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      // initialuint32_t pointIndexts for vectorization
      for (std::size_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Config>; ++pointIndex) {
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
    std::unordered_map<std::string, real*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  real(*dC)[misc::NumPaddedPoints<Config>] = layer.var<LTSLinearSlipWeakening::DC>(Config());
  real(*muS)[misc::NumPaddedPoints<Config>] = layer.var<LTSLinearSlipWeakening::MuS>(Config());
  real(*muD)[misc::NumPaddedPoints<Config>] = layer.var<LTSLinearSlipWeakening::MuD>(Config());
  real(*cohesion)[misc::NumPaddedPoints<Config>] =
      layer.var<LTSLinearSlipWeakening::Cohesion>(Config());
  parameterToStorageMap.insert({"d_c", reinterpret_cast<real*>(dC)});
  parameterToStorageMap.insert({"mu_s", reinterpret_cast<real*>(muS)});
  parameterToStorageMap.insert({"mu_d", reinterpret_cast<real*>(muD)});
  parameterToStorageMap.insert({"cohesion", reinterpret_cast<real*>(cohesion)});
  if (this->faultProvides("forced_rupture_time")) {
    real(*forcedRuptureTime)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::ForcedRuptureTime>(Config());
    parameterToStorageMap.insert(
        {"forced_rupture_time", reinterpret_cast<real*>(forcedRuptureTime)});
  }
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  LinearSlipWeakeningInitializer::initializeFault(drStorage);
  for (auto& layer : drStorage.leaves(Ghost)) {
    real(*regularizedStrength)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(Config());
    const real(*mu)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::Mu>(Config());
    const real(*cohesion)[misc::NumPaddedPoints<Config>] =
        layer.var<LTSLinearSlipWeakening::Cohesion>(Config());
    // the stress the fault starts out under, which is every source that is in effect at the
    // beginning of the simulation and not only the initial state
    const auto sourceCount = stressSourceCount(*drParameters_);
    const auto* stressSources = layer.var<LTSLinearSlipWeakening::StressSourceInFaultCS>(Config());
    const auto* stressSourceOnset = layer.var<LTSLinearSlipWeakening::StressSourceOnset>(Config());
    const auto* stressSourceRiseTime =
        layer.var<LTSLinearSlipWeakening::StressSourceRiseTime>(Config());

    using namespace dr::misc::quantity_indices;
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Config>; ++pointIndex) {
        const auto initialStress =
            stressAtTime<Config>(&stressSources[ltsFace * sourceCount],
                                 &stressSourceRiseTime[ltsFace * sourceCount],
                                 &stressSourceOnset[ltsFace * sourceCount],
                                 sourceCount,
                                 pointIndex,
                                 static_cast<real>(0.0));
        regularizedStrength[ltsFace][pointIndex] =
            -cohesion[ltsFace][pointIndex] -
            mu[ltsFace][pointIndex] * std::min(static_cast<real>(0.0), initialStress[XX]);
      }
    }
  }
}
} // namespace seissol::dr::initializer
