// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LinearSlipWeakeningInitializer.h"

#include "Common/ConfigDispatch.h"
#include "Common/Real.h"
#include "DynamicRupture/Initializer/BaseDRInitializer.h"
#include "DynamicRupture/Misc.h"
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
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
      bool (*dynStressTimePending)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSLinearSlipWeakening::DynStressTimePending>(Cfg());
      real(*slipRate1)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSLinearSlipWeakening::SlipRate1>(Cfg());
      real(*slipRate2)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSLinearSlipWeakening::SlipRate2>(Cfg());
      real(*mu)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSLinearSlipWeakening::Mu>(Cfg());
      const real(*muS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSLinearSlipWeakening::MuS>(Cfg());
      real(*forcedRuptureTime)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSLinearSlipWeakening::ForcedRuptureTime>(Cfg());
      const bool providesForcedRuptureTime = this->faultProvides("forced_rupture_time");
      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        // initialuint32_t pointIndexts for vectorization
        for (std::size_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
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
    });
  }
}

void LinearSlipWeakeningInitializer::addAdditionalParameters(
    std::unordered_map<std::string, void*>& parameterToStorageMap, DynamicRupture::Layer& layer) {
  dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    real(*dC)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSLinearSlipWeakening::DC>(Cfg());
    real(*muS)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSLinearSlipWeakening::MuS>(Cfg());
    real(*muD)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSLinearSlipWeakening::MuD>(Cfg());
    real(*cohesion)[misc::NumPaddedPoints<Cfg>] =
        layer.var<LTSLinearSlipWeakening::Cohesion>(Cfg());
    parameterToStorageMap.insert({"d_c", reinterpret_cast<real*>(dC)});
    parameterToStorageMap.insert({"mu_s", reinterpret_cast<real*>(muS)});
    parameterToStorageMap.insert({"mu_d", reinterpret_cast<real*>(muD)});
    parameterToStorageMap.insert({"cohesion", reinterpret_cast<real*>(cohesion)});
    if (this->faultProvides("forced_rupture_time")) {
      real(*forcedRuptureTime)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSLinearSlipWeakening::ForcedRuptureTime>(Cfg());
      parameterToStorageMap.insert(
          {"forced_rupture_time", reinterpret_cast<real*>(forcedRuptureTime)});
    }
  });
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  LinearSlipWeakeningInitializer::initializeFault(drStorage);
  for (auto& layer : drStorage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
      real(*regularizedStrength)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(Cfg());
      const real(*mu)[misc::NumPaddedPoints<Cfg>] = layer.var<LTSLinearSlipWeakening::Mu>(Cfg());
      const real(*cohesion)[misc::NumPaddedPoints<Cfg>] =
          layer.var<LTSLinearSlipWeakening::Cohesion>(Cfg());
      // the stress the fault starts out under, which is every source that is in effect at the
      // beginning of the simulation and not only the initial state
      const auto sourceCount = stressSourceCount(*drParameters_);
      const auto* stressSources = layer.var<LTSLinearSlipWeakening::StressSourceInFaultCS>(Cfg());
      const auto* stressSourceOnset = layer.var<LTSLinearSlipWeakening::StressSourceOnset>(Cfg());
      const auto* stressSourceRiseTime =
          layer.var<LTSLinearSlipWeakening::StressSourceRiseTime>(Cfg());

      using namespace dr::misc::quantity_indices;
      for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
        for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
          const auto initialStress = stressAtTime<Cfg>(&stressSources[ltsFace * sourceCount],
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
    });
  }
}
} // namespace seissol::dr::initializer
