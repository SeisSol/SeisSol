// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_AGINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_AGINGLAW_H_

#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/AgingLaw.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/CpuImpl/FrictionSolverHost.h"

namespace seissol::dr::friction_law::cpu {

/// The aging law (see friction_law::AgingLaw) on the host.
template <typename Cfg, typename TPMethod>
class AgingLaw final : public FrictionSolverHost<Cfg, AgingLaw<Cfg, TPMethod>>,
                       public friction_law::AgingLaw<Cfg, TPMethod> {
  public:
  explicit AgingLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverHost<Cfg, AgingLaw<Cfg, TPMethod>>(drParameters) {}
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_AGINGLAW_H_
