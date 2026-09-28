// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_NOFAULT_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_NOFAULT_H_

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law::cpu {
/**
 * No friction computation input stress equals output
 */
class NoFault : public BaseFrictionLaw<NoFault> {
  public:
  // the slip rate of a sub-step is the one the friction solve finds at its time point
  static constexpr bool SlipRateAtIntervalEnds = false;
  using BaseFrictionLaw::BaseFrictionLaw;

  static void updateFrictionAndSlip(const FaultStresses<Executor::Host>& faultStresses,
                                    const FaultStresses<Executor::Host>& initialStress,
                                    TractionResults<Executor::Host>& tractionResults,
                                    std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                                    std::array<real, misc::NumPaddedPoints>& strengthBuffer,
                                    std::size_t ltsFace,
                                    uint32_t timeIndex);

  void preHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {
  };
  void postHook(std::array<real, misc::NumPaddedPoints>& stateVariableBuffer, std::size_t ltsFace) {
  };
  //! no dynamic stress time
  [[nodiscard]] const real* dynamicStressQuantity(std::size_t /*ltsFace*/) const { return nullptr; }
  void saveDynamicStressOutput(std::size_t /*ltsFace*/,
                               const real* /*previousQuantity*/,
                               const common::SubStepTimes& /*times*/) {}
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_NOFAULT_H_
