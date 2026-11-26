// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_IMPOSEDSLIPRATES_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_IMPOSEDSLIPRATES_H_

#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/CpuImpl/FrictionSolverHost.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/ImposedSlipRates.h"

namespace seissol::dr::friction_law::cpu {

/// The imposed slip rates (see friction_law::ImposedSlipRates) on the host.
template <typename Cfg, typename STF>
class ImposedSlipRates final : public FrictionSolverHost<Cfg, ImposedSlipRates<Cfg, STF>>,
                               public friction_law::ImposedSlipRates<Cfg, STF> {
  public:
  explicit ImposedSlipRates(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverHost<Cfg, ImposedSlipRates<Cfg, STF>>(drParameters) {}
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_IMPOSEDSLIPRATES_H_
