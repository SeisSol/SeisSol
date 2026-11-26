// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_NOFAULT_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_NOFAULT_H_

#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/CpuImpl/FrictionSolverHost.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/NoFault.h"

namespace seissol::dr::friction_law::cpu {

/// No friction (see friction_law::NoFault) on the host.
template <typename Cfg>
class NoFault final : public FrictionSolverHost<Cfg, NoFault<Cfg>>,
                      public friction_law::NoFault<Cfg> {
  public:
  explicit NoFault(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverHost<Cfg, NoFault<Cfg>>(drParameters) {}
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_NOFAULT_H_
