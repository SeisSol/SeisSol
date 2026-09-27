// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/Misc.h"
#include "RateAndState.h"

namespace seissol::dr::friction_law::cpu {
template <class TPMethod>
class SevereVelocityWeakeningLaw
    : public RateAndStateBase<SevereVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<SevereVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyStorageToLocal(DynamicRupture::Layer& layerData) {}

// Note that we need double precision here, since single precision led to NaNs.
#pragma omp declare simd
  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S>
  S updateStateVariable(std::uint32_t pointIndex,
                        std::size_t faceIndex,
                        double stateVarReference,
                        double timeIncrement,
                        S localSlipRate) {
    const double localSl0 = this->sl0_[faceIndex][pointIndex];

    const S steadyStateStateVariable = localSlipRate * S(localSl0 / this->drParameters_.rsSr0);

    const double preexp1 = -this->drParameters_.rsSr0 * (timeIncrement / localSl0);
    const double exp1v = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    const S localStateVariable = steadyStateStateVariable * S(exp1m) + S(exp1v * stateVarReference);

    return localStateVariable;
  }

  struct MuDetails {
    std::array<real, misc::NumPaddedPoints> a{};
    std::array<real, misc::NumPaddedPoints> c{};
    std::array<real, misc::NumPaddedPoints> f0{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const real localA = this->a_[ltsFace][pointIndex];
      const real localSl0 = this->sl0_[ltsFace][pointIndex];
      const real c = this->b_[ltsFace][pointIndex] * localStateVariable[pointIndex] /
                     (localStateVariable[pointIndex] + localSl0);

      details.a[pointIndex] = localA;
      details.c[pointIndex] = c;
      details.f0[pointIndex] = this->f0_[ltsFace][pointIndex];
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
    const auto stateVariable = this->updateStateVariable(pointIndex,
                                                         ltsFace,
                                                         static_cast<double>(stateVarReference),
                                                         static_cast<double>(timeIncrement),
                                                         dualCast<StateScalar>(slipRate));

    const S localStateVariable = dualCast<real>(stateVariable);
    const S localSl0 = S(this->sl0_[ltsFace][pointIndex]);
    const S c =
        S(this->b_[ltsFace][pointIndex]) * localStateVariable / (localStateVariable + localSl0);
    return S(this->f0_[ltsFace][pointIndex]) +
           S(this->a_[ltsFace][pointIndex]) * slipRate / (slipRate + S(this->drParameters_.rsSr0)) -
           c;
  }

#pragma omp declare simd
  /// generic over the scalar: a dual slip rate carries the derivative out with the value
  template <typename S>
  S updateMu(std::uint32_t pointIndex, S localSlipRateMagnitude, const MuDetails& details) {
    return S(details.f0[pointIndex]) +
           S(details.a[pointIndex]) * localSlipRateMagnitude /
               (localSlipRateMagnitude + S(this->drParameters_.rsSr0)) -
           S(details.c[pointIndex]);
  }

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

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
