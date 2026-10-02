// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENING_H_

#include "Common/Real.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"

namespace seissol::dr::output {
class LinearSlipWeakening : public ReceiverOutputImpl<LinearSlipWeakening> {
  public:
  template <typename Cfg>
  Real<Cfg> computeLocalStrength(LocalInfo<Cfg>& local) {
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    const auto* const cohesions =
        local.layer->template var<LTSLinearSlipWeakening::Cohesion>(Cfg());
    const auto cohesion = cohesions[local.ltsId][local.gpIndex];

    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return -1.0 * local.frictionCoefficient *
               std::min(effectiveNormalStress, static_cast<real>(0.0)) -
           cohesion;
  }

  template <typename Cfg>
  Real<Cfg> computeLocalStrengthSlope(LocalInfo<Cfg>& local) {
    using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    return effectiveNormalStress < 0 ? local.frictionCoefficient : static_cast<real>(0.0);
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENING_H_
