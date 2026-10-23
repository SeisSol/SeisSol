// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_

#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "RateAndState.h"
#include "Solver/MultipleSimulations.h"

#include <cmath>

namespace seissol::dr::friction_law::cpu {

template <typename Cfg, typename TPMethod>
class FastVelocityWeakeningLaw
    : public RateAndStateBase<Cfg, FastVelocityWeakeningLaw<Cfg, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;

  using RateAndStateBase<Cfg, FastVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  void allocateAuxiliaryMemory(GlobalData<Cfg>* globalData) override {
    RateAndStateBase<Cfg, FastVelocityWeakeningLaw, TPMethod>::allocateAuxiliaryMemory(globalData);
    resampleKrnlPrototype_.bindGlobals(*globalData);
  }

  /**
   * Copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyStorageToLocal(DynamicRupture::Layer& layerData) {
    this->srW_ = layerData.var<LTSRateAndStateFastVelocityWeakening::RsSrW>(Cfg());
  }

/**
 * Integrates the state variable ODE in time
 * \f[\frac{\partial \Psi}{\partial t} = - \frac{V}{L}\left(\Psi - \Psi_{ss}(V) \right)\f]
 * with steady state variable \f$\Psi_{ss}\f$.
 * Assume \f$V\f$ is constant through the time interval, then the analytic solution is:
 * \f[ \Psi(t) = \Psi_0 \exp\left( -\frac{V}{L} t \right) + \Psi_{ss} \left( 1 - \exp\left(
 * - \frac{V}{L} t\right) \right).\f]
 * @param stateVarReference \f$ \Psi_0 \f$
 * @param timeIncrement \f$ t \f$
 * @param localSlipRate \f$ V \f$
 * @return \f$ \Psi(t) \f$
 */
#pragma omp declare simd
  /// generic over the scalar the slip rate arrives in, so that the inversion can differentiate
  /// the state variable by the very slip rate it is solving for
  template <typename S>
  [[nodiscard]] S updateStateVariable(std::uint32_t pointIndex,
                                      std::size_t faceIndex,
                                      real stateVarReference,
                                      real timeIncrement,
                                      S localSlipRate) const {
    using std::expm1;
    using std::fmax;
    using std::log;
    using std::pow;
    const real localMuW = this->muW_[faceIndex][pointIndex];
    const real localSrW = this->srW_[faceIndex][pointIndex];
    const real localA = this->a_[faceIndex][pointIndex];
    const real localSl0 = this->sl0_[faceIndex][pointIndex];

    // low-velocity steady state friction coefficient
    const S lowVelocityFriction = fmax(S(static_cast<real>(0)),
                                       S(this->f0_[faceIndex][pointIndex]) -
                                           S(this->b_[faceIndex][pointIndex] - localA) *
                                               log(localSlipRate / S(this->drParameters_.rsSr0)));
    const S steadyStateFrictionCoefficient =
        S(localMuW) +
        (lowVelocityFriction - S(localMuW)) /
            pow(S(static_cast<real>(1.0)) + misc::power<8>(localSlipRate / S(localSrW)),
                static_cast<real>(1.0 / 8.0));
    // TODO: check again, if double precision is necessary here (earlier, there were cancellation
    // issues)
    // stated over V / (2 V_0) rather than under 2 V_0 / V: the reciprocal's derivative is
    // -2 V_0 / V^2, which at the slip-rate floor is 2e64 and leaves single precision, while the
    // state variable it belongs to is perfectly well behaved there
    const S steadyStateStateVariable =
        S(localA) *
        rs::logsinhOver(localSlipRate * S(static_cast<real>(0.5) / this->drParameters_.rsSr0),
                        steadyStateFrictionCoefficient / S(localA));

    // Exact integration over the step with the slip rate held fixed. Written as psi0 + (psi_ss -
    // psi0) * (1 - exp(p)), it keeps psi_ss as its exact fixed point in any precision. The form
    // psi_ss * (1 - exp(p)) + psi0 * exp(p) rounds exp and expm1 independently, and their mismatch
    // moves the fixed point by L / (V dt) times as much, a few thousand here.
    const S preexp1 = -localSlipRate * S(timeIncrement / localSl0);
    const S exp1m = -expm1(preexp1);
    const S localStateVariable =
        S(stateVarReference) + (steadyStateStateVariable - S(stateVarReference)) * exp1m;
    assert((std::isfinite(valueOf(localStateVariable)) ||
            pointIndex >= misc::NumBoundaryGaussPoints<Cfg> * Cfg::NumSimulations) &&
           "Inf/NaN detected");
    return localStateVariable;
  }

  struct MuDetails {
    std::array<real, misc::NumPaddedPoints<Cfg>> a{};
    std::array<real, misc::NumPaddedPoints<Cfg>> cLin{};
    std::array<real, misc::NumPaddedPoints<Cfg>> cExpLog{};
    std::array<real, misc::NumPaddedPoints<Cfg>> cExp{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints<Cfg>>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      const real localA = this->a_[ltsFace][pointIndex];

      const real cLin = static_cast<real>(0.5) / this->drParameters_.rsSr0;
      const real cExpLog = localStateVariable[pointIndex] / localA;
      const real cExp = rs::computeCExp(cExpLog);

      details.a[pointIndex] = localA;
      details.cLin[pointIndex] = cLin;
      details.cExpLog[pointIndex] = cExpLog;
      details.cExp[pointIndex] = cExp;
    }
    return details;
  }

  /// the precision the state variable of this law is stated in
  using StateScalar = real;

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
    const auto stateVariable = this->updateStateVariable(
        pointIndex, ltsFace, stateVarReference, timeIncrement, dualCast<StateScalar>(slipRate));

    const S localStateVariable = dualCast<real>(stateVariable);
    const S localA = S(this->a_[ltsFace][pointIndex]);
    const S cExpLog = localStateVariable / localA;
    const S cLin = S(static_cast<real>(0.5) / this->drParameters_.rsSr0);
    return localA * rs::arsinhexp(cLin * slipRate, cExpLog, rs::computeCExp(cExpLog));
  }

/**
 * Computes the friction coefficient from the state variable and slip rate
 * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp
 * \left(\frac{\Psi}{a}\right)\right). \f]
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
 * \f[\frac{\partial}{\partial V}\mu = \frac{aC}{\sqrt{ (VC)^2 + 1} \text{ with } C =
 * \frac{1}{2V_0} \cdot \exp \left(\frac{\Psi}{a}\right)\right).\f]
 * @param localSlipRateMagnitude \f$ V \f$
 * @param localStateVariable \f$ \Psi \f$
 * @return \f$ \mu \f$
 */
#pragma omp declare simd

  /**
   * Resample the state variable.
   */
  void resampleStateVar(const std::array<real, misc::NumPaddedPoints<Cfg>>& stateVariableBuffer,
                        std::size_t ltsFace) const {
    alignas(Alignment) std::array<real, misc::NumPaddedPoints<Cfg>> deltaStateVar = {0};
    alignas(Alignment) std::array<real, misc::NumPaddedPoints<Cfg>> resampledDeltaStateVar = {0};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      deltaStateVar[pointIndex] =
          stateVariableBuffer[pointIndex] - this->stateVariable_[ltsFace][pointIndex];
    }
    auto resampleKrnl = resampleKrnlPrototype_;
    resampleKrnl.originalQ = deltaStateVar.data();
    resampleKrnl.resampledQ = resampledDeltaStateVar.data();
    resampleKrnl.execute();

#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; pointIndex++) {
      this->stateVariable_[ltsFace][pointIndex] =
          this->stateVariable_[ltsFace][pointIndex] + resampledDeltaStateVar[pointIndex];
    }
  }

  protected:
  real (*__restrict srW_)[misc::NumPaddedPoints<Cfg>]{nullptr};
  dynamicRupture::kernel::resampleParameter<Cfg> resampleKrnlPrototype_;
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_FASTVELOCITYWEAKENINGLAW_H_
