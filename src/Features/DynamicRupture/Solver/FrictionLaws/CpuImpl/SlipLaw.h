// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SLIPLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SLIPLAW_H_

#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/CpuImpl/FrictionSolverHost.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/SlipLaw.h"

namespace seissol::dr::friction_law::cpu {

/// The slip law (see friction_law::SlipLaw) on the host.
template <typename Cfg, typename TPMethod>
class SlipLaw final : public FrictionSolverHost<Cfg, SlipLaw<Cfg, TPMethod>>,
                      public friction_law::SlipLaw<Cfg, TPMethod> {
  public:
  explicit SlipLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverHost<Cfg, SlipLaw<Cfg, TPMethod>>(drParameters) {}
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_SLIPLAW_H_
