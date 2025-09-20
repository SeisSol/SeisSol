// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_

#include "DynamicRupture/Output/ReceiverBasedOutput.h"
#include <Memory/Descriptor/DynamicRupture.h>

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public LinearSlipWeakening {
  real computeLocalStrength(LocalInfo& local) override {
    const auto* const regularizedStrengths =
        getCellData<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(local);
    return regularizedStrengths[local.gpIndex];
  }

  std::vector<std::size_t> getOutputVariables() const override {
    auto baseVector = LinearSlipWeakening::getOutputVariables();
    baseVector.push_back(
        drStorage->info<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>().index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
