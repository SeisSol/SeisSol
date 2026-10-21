// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_

#include "Common/Real.h"
#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Memory/Descriptor/DynamicRupture.h"

#include <algorithm>
#include <cmath>

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public ReceiverOutputImpl<LinearSlipWeakeningBimaterial> {
  public:
  template <typename Cfg>
  Real<Cfg> computeLocalStrength(LocalInfo<Cfg>& local) {
    const auto* const regularizedStrengths =
        getCellData<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(local);
    return regularizedStrengths[local.gpIndex];
  }

  template <typename Cfg>
  Real<Cfg> computeLocalStrengthSlope(LocalInfo<Cfg>& local) {
    using real = Real<Cfg>;
    // The Prakash-Clifton regularization low-passes the strength, so only the fraction
    // -expm1(-(V + vStar) dt / prakashLength) of a normal stress change arrives instantaneously --
    // evaluated with the slip rate and the sub time step of the friction solve that produced the
    // stored regularized strength. It scales the slope of the unregularized law.
    const auto* const slipRateMagnitude = getCellData<DynamicRupture::SlipRateMagnitude>(local);
    const auto effectiveNormalStress =
        local.transientNormalTraction + local.iniNormalTraction - local.fluidPressure;
    if (effectiveNormalStress >= 0) {
      return 0.0;
    }

    const auto expval = -(std::max(static_cast<real>(0.0), slipRateMagnitude[local.gpIndex]) +
                          static_cast<real>(drParameters_->vStar)) *
                        static_cast<real>(local.deltaT) /
                        static_cast<real>(drParameters_->prakashLength);
    return local.frictionCoefficient * -std::expm1(expval);
  }

  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = ReceiverOutputImpl::getOutputVariables();
    baseVector.push_back(
        drStorage_->info<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
