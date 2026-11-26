// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/CpuImpl/FrictionSolverHost.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/SevereVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::cpu {

/// The rate-and-state law with severe velocity weakening (see
/// friction_law::SevereVelocityWeakeningLaw) on the host.
template <typename Cfg, typename TPMethod>
class SevereVelocityWeakeningLaw final
    : public FrictionSolverHost<Cfg, SevereVelocityWeakeningLaw<Cfg, TPMethod>>,
      public friction_law::SevereVelocityWeakeningLaw<Cfg, TPMethod> {
  public:
  explicit SevereVelocityWeakeningLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverHost<Cfg, SevereVelocityWeakeningLaw<Cfg, TPMethod>>(drParameters) {}
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
