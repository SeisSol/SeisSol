// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_

#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public LinearSlipWeakening {
  protected:
  real computeLocalStrength(LocalInfo& local) override {
    const auto* const regularizedStrengths =
        getCellData<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(local);
    return regularizedStrengths[local.gpIndex];
  }

  real computeLocalStrengthSlope(LocalInfo& /*local*/) override {
    // The Prakash-Clifton regularisation low-passes the strength, so only the fraction
    // -expm1(-(V + vStar) dt / prakashLength) of a normal stress change arrives instantaneously.
    // That factor needs the time step of the friction solve, which the receiver output does not
    // see, so the reconstruction leaves the coupling out rather than overstating it.
    return 0.0;
  }

  public:
  [[nodiscard]] std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = LinearSlipWeakening::getOutputVariables();
    baseVector.push_back(
        drStorage_->info<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
