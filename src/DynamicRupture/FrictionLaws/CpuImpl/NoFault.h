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
template <typename Cfg>
class NoFault : public BaseFrictionLaw<Cfg, NoFault<Cfg>> {
  public:
  using real = Real<Cfg>;
  using BaseFrictionLaw<Cfg, NoFault<Cfg>>::BaseFrictionLaw;

  void copyStorageToLocal(DynamicRupture::Layer& layerData) override {}

  static void
      updateFrictionAndSlip(const FaultStresses<Cfg, Executor::Host>& faultStresses,
                            TractionResults<Cfg, Executor::Host>& tractionResults,
                            std::array<real, misc::NumPaddedPoints<Cfg>>& stateVariableBuffer,
                            std::array<real, misc::NumPaddedPoints<Cfg>>& strengthBuffer,
                            std::size_t ltsFace,
                            uint32_t timeIndex);

  void preHook(std::array<real, misc::NumPaddedPoints<Cfg>>& stateVariableBuffer,
               std::size_t ltsFace) {};
  void postHook(std::array<real, misc::NumPaddedPoints<Cfg>>& stateVariableBuffer,
                std::size_t ltsFace) {};
  void saveDynamicStressOutput(std::size_t ltsFace, real time) {};
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_NOFAULT_H_
