// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_AGINGLAW_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_AGINGLAW_H_

#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Features/DynamicRupture/Physics/Dual.h"
#include "Features/DynamicRupture/Physics/RateAndStateCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/SlowVelocityWeakeningLaw.h"

#include <cmath>
#include <cstdint>

namespace seissol::dr::friction_law {

/**
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
template <typename Cfg, typename TPMethod>
class AgingLaw : public SlowVelocityWeakeningLaw<Cfg, AgingLaw<Cfg, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;

  /**
   * Integrates the state variable ODE in time
   * \f[ \frac{\partial \Psi}{\partial t} = 1 - \frac{V}{L} \Psi \f]
   * Analytic solution:
   * \f[\Psi(t) = \frac{L}{V} + \left(\Psi_0 - \frac{L}{V}\right) \exp\left( -\frac{V}{L} \cdot t
   * \right). \f]
   * @param stateVarReference \f$ \Psi_0 \f$
   * @param timeIncrement \f$ t \f$
   * @param localSlipRate \f$ V \f$
   * @return \f$ \Psi(t) \f$
   */
  /// generic over the scalar the slip rate arrives in; see SlowVelocityWeakeningLaw::StateScalar
  template <typename S, common::RangeType Type>
  SEISSOL_HOSTDEVICE static S stateVariableAt(const FrictionLawContext<Cfg, Type>& ctx,
                                              std::uint32_t pointIndex,
                                              real stateVarReference,
                                              real timeIncrement,
                                              S localSlipRate) {
    using std::exp;
    const real localSl0 = ctx.data->sl0[ctx.ltsFace][pointIndex];
    const S preexp1 = -localSlipRate * S(timeIncrement / localSl0);
    // (L / V) (1 - exp(-V t / L)) is t times the mean of the relaxation over the step. Stated that
    // way, neither L / V nor its derivative -L / V^2 appears: at the slip-rate floor the quotient
    // is 2e33 and multiplies a relaxation of 5e-38, and the derivative leaves single precision
    // outright, though the state itself is simply the time step there.
    const S weight = rs::relaxationWeight(-preexp1);
    // The relaxation towards L / V with a single exponential, psi0 + (L / V - psi0) (1 - exp(p)),
    // keeps L / V as its fixed point (see FastVelocityWeakeningLaw::stateVariableAt); with
    // 1 - exp(p) = -p weight, it reads psi0 + weight (t + psi0 p) and leaves L / V out as well.
    // With exp once the step relaxes more than half the way, since mu takes the logarithm of the
    // state and the form with expm1 alone may then round it to zero.
    if (valueOf(-preexp1 * weight) < static_cast<real>(0.5)) {
      return S(stateVarReference) + weight * (S(timeIncrement) + S(stateVarReference) * preexp1);
    }
    return S(stateVarReference) * exp(preexp1) + S(timeIncrement) * weight;
  }
};

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_AGINGLAW_H_
