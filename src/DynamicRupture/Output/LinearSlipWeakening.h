// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENING_H_

#include "DynamicRupture/Output/ReceiverBasedOutput.h"

namespace seissol::dr::output {
class LinearSlipWeakening : public ReceiverOutput {
  protected:
  std::array<real, seissol::multisim::NumSimulations> computeLocalStrength(LocalInfo& local) override {
    using DrLtsDescrType = seissol::initializer::LTSLinearSlipWeakening;
    const auto* const cohesions = local.layer->var(static_cast<DrLtsDescrType*>(drDescr)->cohesion);
    std::array<real, seissol::multisim::NumSimulations> strengthsArray{};

    for (size_t sim = 0; sim < seissol::multisim::NumSimulations; sim++)
    {
      const auto cohesion = cohesions[local.ltsId][local.nearestGpIndex*seissol::multisim::NumSimulations + sim];
      const auto effectiveNormalStress =
        local.transientNormalTraction[sim] + local.iniNormalTraction[sim] - local.fluidPressure[sim];
      strengthsArray[sim] =
          -1.0 * local.frictionCoefficient[sim] *
          std::min(effectiveNormalStress, static_cast<real>(0.0)) - cohesion;
      }
    return strengthsArray;
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_LINEARSLIPWEAKENING_H_
