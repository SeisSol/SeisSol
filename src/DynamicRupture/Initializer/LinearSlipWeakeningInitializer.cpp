// SPDX-FileCopyrightText: 2022 SeisSol Group
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
    bool (*dynStressTimePending)[misc::NumPaddedPoints] =
        layer.var<LTSLinearSlipWeakening::DynStressTimePending>();
    real(*slipRate1)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::SlipRate1>();
    real(*slipRate2)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::SlipRate2>();
    real(*mu)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::Mu>();
    const real(*muS)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::MuS>();
    real(*forcedRuptureTime)[misc::NumPaddedPoints] =
        layer.var<LTSLinearSlipWeakening::ForcedRuptureTime>();
    const bool providesForcedRuptureTime = this->faultProvides("forced_rupture_time");
    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      // initialuint32_t pointIndexts for vectorization
      for (std::size_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
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
  real(*dC)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::DC>();
  real(*muS)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::MuS>();
  real(*muD)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::MuD>();
  real(*cohesion)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::Cohesion>();
  parameterToStorageMap.insert({"d_c", reinterpret_cast<real*>(dC)});
  parameterToStorageMap.insert({"mu_s", reinterpret_cast<real*>(muS)});
  parameterToStorageMap.insert({"mu_d", reinterpret_cast<real*>(muD)});
  parameterToStorageMap.insert({"cohesion", reinterpret_cast<real*>(cohesion)});
  if (this->faultProvides("forced_rupture_time")) {
    real(*forcedRuptureTime)[misc::NumPaddedPoints] =
        layer.var<LTSLinearSlipWeakening::ForcedRuptureTime>();
    parameterToStorageMap.insert(
        {"forced_rupture_time", reinterpret_cast<real*>(forcedRuptureTime)});
  }
}

void LinearSlipWeakeningBimaterialInitializer::initializeFault(DynamicRupture::Storage& drStorage) {
  LinearSlipWeakeningInitializer::initializeFault(drStorage);
  for (auto& layer : drStorage.leaves(Ghost)) {
    real(*regularizedStrength)[misc::NumPaddedPoints] =
        layer.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>();
    const real(*mu)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::Mu>();
    const real(*cohesion)[misc::NumPaddedPoints] = layer.var<LTSLinearSlipWeakening::Cohesion>();
    const auto* initialStressInFaultCS =
        layer.var<LTSLinearSlipWeakening::InitialStressInFaultCS>();

    for (std::size_t ltsFace = 0; ltsFace < layer.size(); ++ltsFace) {
      for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
        regularizedStrength[ltsFace][pointIndex] =
            -cohesion[ltsFace][pointIndex] -
            mu[ltsFace][pointIndex] *
                std::min(static_cast<real>(0.0), initialStressInFaultCS[ltsFace][0][pointIndex]);
      }
    }
  }
}
} // namespace seissol::dr::initializer
