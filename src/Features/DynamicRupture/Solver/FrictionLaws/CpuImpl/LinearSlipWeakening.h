// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_LINEARSLIPWEAKENING_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_LINEARSLIPWEAKENING_H_

#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/CpuImpl/FrictionSolverHost.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/LinearSlipWeakening.h"

namespace seissol::dr::friction_law::cpu {

template <typename Cfg>
using NoSpecialization = friction_law::NoSpecialization<Cfg>;
template <typename Cfg>
using BiMaterialFault = friction_law::BiMaterialFault<Cfg>;
template <typename Cfg>
using TPApprox = friction_law::TPApprox<Cfg>;

/// The linear slip weakening laws (see friction_law::LinearSlipWeakening) on the host.
template <typename Cfg, typename SpecializationT>
class LinearSlipWeakeningLaw final
    : public FrictionSolverHost<Cfg, LinearSlipWeakeningLaw<Cfg, SpecializationT>>,
      public LinearSlipWeakening<Cfg, SpecializationT> {
  public:
  explicit LinearSlipWeakeningLaw(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverHost<Cfg, LinearSlipWeakeningLaw<Cfg, SpecializationT>>(drParameters) {}
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_CPUIMPL_LINEARSLIPWEAKENING_H_
