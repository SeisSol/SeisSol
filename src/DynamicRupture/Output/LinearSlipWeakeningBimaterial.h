// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_

#include "DynamicRupture/Output/ReceiverBasedOutput.h"

namespace seissol::dr::output {
class LinearSlipWeakeningBimaterial : public LinearSlipWeakening {
  std::array<real, seissol::multisim::NumSimulations> computeLocalStrength(LocalInfo& local) override {
    using DrLtsDescrType = seissol::initializer::LTSLinearSlipWeakeningBimaterial;
    const auto* const regularizedStrengths =
        getCellData(local, static_cast<DrLtsDescrType*>(drDescr)->regularizedStrength);
    std::array<real, seissol::multisim::NumSimulations> regularizedStrengthsArray{};
    for (unsigned int i = 0; i < seissol::multisim::NumSimulations; ++i) {
      regularizedStrengthsArray[i] = regularizedStrengths[local.nearestGpIndex*seissol::multisim::NumSimulations + i];
    }
    return regularizedStrengthsArray;
  }

  std::vector<std::size_t> getOutputVariables() const override {
    using DrLtsDescrType = seissol::initializer::LTSLinearSlipWeakeningBimaterial;
    auto baseVector = LinearSlipWeakening::getOutputVariables();
    baseVector.push_back(static_cast<DrLtsDescrType*>(drDescr)->regularizedStrength.index);
    return baseVector;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENINGBIMATERIAL_H_
