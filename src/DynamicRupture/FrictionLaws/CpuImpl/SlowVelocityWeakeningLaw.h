// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_

#include "RateAndState.h"

namespace seissol::dr::friction_law::cpu {
template <class Derived, class TPMethod>
class SlowVelocityWeakeningLaw
    : public RateAndStateBase<SlowVelocityWeakeningLaw<Derived, TPMethod>, TPMethod> {
  public:
  SlowVelocityWeakeningLaw() = default;
  using RateAndStateBase<SlowVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyStorageToLocal(DynamicRupture::Layer& layerData) {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

// Note that we need double precision here, since single precision led to NaNs.
#pragma omp declare simd
  template <typename S>
  S updateStateVariable(std::uint32_t pointIndex,
                        std::size_t faceIndex,
                        double stateVarReference,
                        double timeIncrement,
                        S localSlipRate) {
    return static_cast<Derived*>(this)->updateStateVariable(
        pointIndex, faceIndex, stateVarReference, timeIncrement, localSlipRate);
  }

  struct MuDetails {
    std::array<real, misc::NumPaddedPoints> a{};
    std::array<real, misc::NumPaddedPoints> cLin{};
    std::array<real, misc::NumPaddedPoints> cExpLog{};
    std::array<real, misc::NumPaddedPoints> cExp{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const real localA = this->a_[ltsFace][pointIndex];
      const real localSl0 = this->sl0_[ltsFace][pointIndex];
      const real log1 =
          std::log(this->drParameters_.rsSr0 * localStateVariable[pointIndex] / localSl0);
      const real localF0 = this->f0_[ltsFace][pointIndex];
      const real localB = this->b_[ltsFace][pointIndex];

      const real cLin = static_cast<real>(0.5) / this->drParameters_.rsSr0;
      const real cExpLog = (localF0 + localB * log1) / localA;
      const real cExp = rs::computeCExp(cExpLog);

      details.a[pointIndex] = localA;
      details.cLin[pointIndex] = cLin;
      details.cExpLog[pointIndex] = cExpLog;
      details.cExp[pointIndex] = cExp;
    }
    return details;
  }

  /// the precision the state variable of this law is stated in
  using StateScalar = double;

  /// the state variable is a closed-form function of the slip rate, so the inversion can carry it
  /// inside its own iteration instead of relaying it through a fixed point
  static constexpr bool FoldsStateVariable = true;

  /// The friction coefficient at a slip rate, with the state variable evaluated at that very slip
  /// rate. Both dependencies travel through the scalar, so a dual number comes back carrying
  /// d(mu)/dV of the composition.
  template <typename S>
  S updateMuFolded(std::size_t ltsFace,
                   std::uint32_t pointIndex,
                   S slipRate,
                   real stateVarReference,
                   real timeIncrement) {
    using std::log;
    const auto stateVariable =
        static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                         ltsFace,
                                                         static_cast<double>(stateVarReference),
                                                         static_cast<double>(timeIncrement),
                                                         dualCast<StateScalar>(slipRate));

    const S localStateVariable = dualCast<real>(stateVariable);
    const S localA = S(this->a_[ltsFace][pointIndex]);
    const S localSl0 = S(this->sl0_[ltsFace][pointIndex]);
    const S log1 = log(S(this->drParameters_.rsSr0) * localStateVariable / localSl0);
    const S cExpLog =
        (S(this->f0_[ltsFace][pointIndex]) + S(this->b_[ltsFace][pointIndex]) * log1) / localA;
    const S cLin = S(static_cast<real>(0.5) / this->drParameters_.rsSr0);
    return localA * rs::arsinhexp(cLin * slipRate, cExpLog, rs::computeCExp(cExpLog));
  }

  /**
   * Computes the friction coefficient from the state variable and slip rate
   * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Psi
   * / L)}{a} \right)\right).\f]
   * Note that we need double precision here, since single precision led to NaNs.
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Psi \f$
   * @return \f$ \mu \f$
   */
#pragma omp declare simd
  /// generic over the scalar: a dual slip rate carries the derivative out with the value
  template <typename S>
  S updateMu(std::uint32_t pointIndex, S localSlipRateMagnitude, const MuDetails& details) {
    const S lx = S(details.cLin[pointIndex]) * localSlipRateMagnitude;
    return S(details.a[pointIndex]) *
           rs::arsinhexp(lx, S(details.cExpLog[pointIndex]), S(details.cExp[pointIndex]));
  }

  /**
   * Computes the derivative of the friction coefficient with respect to the slip rate.
   * \f[\frac{\partial}{\partial V}\mu = \frac{aC}{\sqrt{(VC)^2 +1}} \text{ with } C =
   * \frac{1}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Psi / L)}{a} \right). \f]
   * Note that we need double precision here, since single precision led to NaNs.
   * @param localSlipRateMagnitude \f$ V \f$
   * @param localStateVariable \f$ \Psi \f$
   * @return \f$ \mu \f$
   */
#pragma omp declare simd

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws, we just copy the buffer into the
   * member variable.
   */
  void resampleStateVar(const std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                        std::size_t ltsFace) const {
#pragma omp simd
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->stateVariable_[ltsFace][pointIndex] = stateVariableBuffer[pointIndex];
    }
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
