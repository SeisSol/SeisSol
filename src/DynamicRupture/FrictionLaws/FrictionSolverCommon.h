#ifndef SEISSOL_FRICTIONSOLVER_COMMON_H
#define SEISSOL_FRICTIONSOLVER_COMMON_H

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/DynamicRupture.h"
#include "Kernels/DynamicRupture.h"
#include "Numerical_aux/GaussianNucleationFunction.h"
#include <type_traits>

/**
 * Contains common functions required both for CPU and GPU impl.
 * of Dynamic Rupture solvers. The functions placed in
 * this class definition (of the header file) result
 * in the function inlining required for GPU impl.
 */
namespace seissol::dr::friction_law::common {

template <size_t Start, size_t End, size_t Step>
struct ForLoopRange {
  static constexpr size_t start{Start};
  static constexpr size_t end{End};
  static constexpr size_t step{Step};
};

enum class RangeType { CPU, GPU };

template <RangeType Type>
struct NumPoints {
  private:
  using CpuRange = ForLoopRange<0, dr::misc::numPaddedPoints, 1>;
  using GpuRange = ForLoopRange<0, 1, 1>;

  public:
  using Range = typename std::conditional<Type == RangeType::CPU, CpuRange, GpuRange>::type;
};

template <RangeType Type>
struct QInterpolated {
  private:
  using CpuRange = ForLoopRange<0, tensor::QInterpolated::size(), 1>;
  using GpuRange = ForLoopRange<0, tensor::QInterpolated::size(), misc::numPaddedPoints>;

  public:
  using Range = typename std::conditional<Type == RangeType::CPU, CpuRange, GpuRange>::type;
};

/**
 * Asserts whether all relevant arrays are properly aligned
 */
inline void checkAlignmentPreCompute(
    const real qIPlus[CONVERGENCE_ORDER][dr::misc::numQuantities][dr::misc::numPaddedPoints],
    const real qIMinus[CONVERGENCE_ORDER][dr::misc::numQuantities][dr::misc::numPaddedPoints],
    const FaultStresses& faultStresses) {
  using namespace dr::misc::quantity_indices;
  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][U]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][V]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][W]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][N]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T1]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T2]) % ALIGNMENT == 0);

    assert(reinterpret_cast<uintptr_t>(qIMinus[o][U]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][V]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][W]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][N]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T1]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T2]) % ALIGNMENT == 0);

    assert(reinterpret_cast<uintptr_t>(faultStresses.normalStress[o]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(faultStresses.traction1[o]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(faultStresses.traction2[o]) % ALIGNMENT == 0);
  }
}

/**
 * Calculate traction and normal stress at the interface of a face.
 * Using equations (A2) from Pelties et al. 2014
 * Definiton of eta and impedance Z are found in Carsten Uphoff's dissertation on page 47 and in
 * equation (4.51) respectively.
 *
 * @param[out] faultStresses contains normalStress, traction1, traction2
 *             at the 2d face quadrature nodes evaluated at the time
 *             quadrature points
 * @param[in] impAndEta contains eta and impedance values
 * @param[in] qInterpolatedPlus a plus side dofs interpolated at time sub-intervals
 * @param[in] qInterpolatedMinus a minus side dofs interpolated at time sub-intervals
 */
template <RangeType Type = RangeType::CPU>
inline void precomputeStressFromQInterpolated(
    FaultStresses& faultStresses,
    ImpedancesAndEta& impAndEta,
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qStrainInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qStrainInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    unsigned startLoopIndex = 0) {

  static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                "Different number of quadrature points?");

  auto etaP = impAndEta.etaP;
  auto etaS = impAndEta.etaS;
  auto invZp = impAndEta.invZp;
  auto invZs = impAndEta.invZs;
  auto invZpNeig = impAndEta.invZpNeig;
  auto invZsNeig = impAndEta.invZsNeig;

  auto la0P = impAndEta.lambda0P;
  auto mu0P = impAndEta.mu0P;
  auto gaRP = impAndEta.gammaRP;
  auto xi0P = impAndEta.xi0P;
  auto rhoP = impAndEta.rho0P;
  auto la0M = impAndEta.lambda0M;
  auto mu0M = impAndEta.mu0M;
  auto gaRM = impAndEta.gammaRM;
  auto xi0M = impAndEta.xi0M;
  auto rhoM = impAndEta.rho0M;

  using QInterpolatedShapeT = const real(*)[misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
  auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));
  auto* qStrainIPlus = (reinterpret_cast<QInterpolatedShapeT>(qStrainInterpolatedPlus));
  auto* qStrainIMinus = (reinterpret_cast<QInterpolatedShapeT>(qStrainInterpolatedMinus));

  using namespace dr::misc::quantity_indices;
  unsigned DAM = 9;
  unsigned BRE = 10;

#ifndef ACL_DEVICE
  checkAlignmentPreCompute(qIPlus, qIMinus, faultStresses);
#endif

  // Derive averaged field and material properties at t0 of each time step
  real exxP, eyyP, ezzP, exyP, eyzP, ezxP, damP, breP;
  real exxM, eyyM, ezzM, exyM, eyzM, ezxM, damM, breM;

  exxP = 0; eyyP = 0; ezzP = 0; exyP = 0; eyzP = 0; ezxP = 0; damP = 0;
  exxM = 0; eyyM = 0; ezzM = 0; exyM = 0; eyzM = 0; ezxM = 0; damM = 0;

  using Range = typename NumPoints<Type>::Range;
#ifndef ACL_DEVICE
  #pragma omp simd
  /// TODO: test if reduction would work here
  // reduction(+:exxP,eyyP,ezzP,exyP,eyzP,ezxP,damP,exxM,eyyM,ezzM,exyM,eyzM,ezxM,damM)
#endif
  for (auto index = Range::start; index < Range::end; index += Range::step) {
    auto i{startLoopIndex + index};

    exxP += qStrainIPlus[0][XX][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    eyyP += qStrainIPlus[0][YY][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    ezzP += qStrainIPlus[0][ZZ][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    exyP += qStrainIPlus[0][XY][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    eyzP += qStrainIPlus[0][YZ][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    ezxP += qStrainIPlus[0][XZ][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    damP += qStrainIPlus[0][DAM][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    breP += qStrainIPlus[0][BRE][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;

    exxM += qStrainIMinus[0][XX][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    eyyM += qStrainIMinus[0][YY][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    ezzM += qStrainIMinus[0][ZZ][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    exyM += qStrainIMinus[0][XY][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    eyzM += qStrainIMinus[0][YZ][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    ezxM += qStrainIMinus[0][XZ][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    damM += qStrainIMinus[0][DAM][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
    breM += qStrainIMinus[0][BRE][i] * 1.0/seissol::dr::misc::numberOfBoundaryGaussPoints;
  }

  real epsInitxx = 0.0e-4; // eps_xx0
  real epsInityy = 0.0e-3; // eps_yy0
  real epsInitzz = 0.0e-4; // eps_zz0
  real epsInitxy = 0.0e-3; // eps_xx0
  real epsInityz = -0e-1; // eps_yy0
  real epsInitzx = -0e-1; // eps_zz0
  
  //TODO(NONLINEAR) what are these numbers?
  real aB0 = 7.43e9;
  real aB1 = -12.14e9;
  real aB2 = 18.93e9;
  real aB3 = -5.067e9;

  real EspIp = (exxP+epsInitxx) + (eyyP+epsInityy) + (ezzP+epsInitzz);
  real EspIIp = (exxP+epsInitxx)*(exxP+epsInitxx)
    + (eyyP+epsInityy)*(eyyP+epsInityy)
    + (ezzP+epsInitzz)*(ezzP+epsInitzz)
    + 2*(exyP+epsInitxy)*(exyP+epsInitxy)
    + 2*(eyzP+epsInityz)*(eyzP+epsInityz)
    + 2*(ezxP+epsInitzx)*(ezxP+epsInitzx);
  real alphap = damP;
  real xip;
  if (EspIIp > 1e-30){
    xip = EspIp / std::sqrt(EspIIp);
  } else{
    xip = 0.0;
  }

  real EspIm = (exxM+epsInitxx) + (eyyM+epsInityy) + (ezzM+epsInitzz);
  real EspIIm = (exxM+epsInitxx)*(exxM+epsInitxx)
    + (eyyM+epsInityy)*(eyyM+epsInityy)
    + (ezzM+epsInitzz)*(ezzM+epsInitzz)
    + 2*(exyM+epsInitxy)*(exyM+epsInitxy)
    + 2*(eyzM+epsInityz)*(eyzM+epsInityz)
    + 2*(ezxM+epsInitzx)*(ezxM+epsInitzx);
  real alpham = damM;
  real xim;
  if (EspIIm > 1e-30){
    xim = EspIm / std::sqrt(EspIIm);
  } else{
    xim = 0.0;
  }

  // compute laP, muP, laM and muM
  auto muP = (1-breP) * (mu0P
    - alphap*xi0P*gaRP
    - 0.5*alphap*gaRP*xip)
    + breP * (
      (aB0 + 0.5*aB1*xip - 0.5*aB3*xip*xip*xip)
    );
  auto laP = (1-breP) * (la0P
  - alphap*gaRP*(eyyP+epsInityy)/std::sqrt(EspIIp))
  + breP * (
    (2.0*aB2 + 3.0*aB3*xip) + aB1*(eyyP+epsInityy)/std::sqrt(EspIIp)
  );

  auto muM = (1-breM) * (mu0M
    - alpham*xi0M*gaRM
    - 0.5*alpham*gaRM*xim)
    + breM * (
      (aB0 + 0.5*aB1*xim - 0.5*aB3*xim*xim*xim)
    );
  auto laM = (1-breM) * (la0M
  - alpham*gaRM*(eyyM+epsInityy)/std::sqrt(EspIIm))
  + breM * (
    (2.0*aB2 + 3.0*aB3*xim) + aB1*(eyyM+epsInityy)/std::sqrt(EspIIm)
  );

  invZp = 1.0/std::sqrt( rhoP*(laP+2*muP) );
  invZpNeig = 1.0/std::sqrt( rhoM*(laM+2*muM) );
  invZs = 1.0/std::sqrt( rhoP*(muP) );
  invZsNeig = 1.0/std::sqrt( rhoM*(muM) );

  etaP = 1.0/(invZp + invZpNeig);
  etaS = 1.0/(invZs + invZsNeig);

  // Change the values in impAndEta so that they can be used in later calculation
  // of "hat" variables.
  impAndEta.etaP = etaP;
  impAndEta.etaS = etaS;
  impAndEta.invZp = invZp;
  impAndEta.invZs = invZs;
  impAndEta.invZpNeig = invZpNeig;
  impAndEta.invZsNeig = invZsNeig;

  impAndEta.invEtaS = 1.0 / etaS;
  impAndEta.zp = 1.0 / invZp;
  impAndEta.zs = 1.0 / invZs;
  impAndEta.zpNeig = 1.0 / invZpNeig;
  impAndEta.zsNeig = 1.0 / invZsNeig;

  impAndEta.csOcpTZsOZp = std::sqrt( (muP)/rhoP ) / std::sqrt( (laP+2*muP)/rhoP )
                          * impAndEta.zs / impAndEta.zp;
  impAndEta.csOcpTZsOZpNeig = std::sqrt( (muM)/rhoM ) / std::sqrt( (laM+2*muM)/rhoM )
                          * impAndEta.zsNeig / impAndEta.zpNeig;

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
    #pragma omp simd
#endif
    for (auto index = Range::start; index < Range::end; index += Range::step) {
      auto i{startLoopIndex + index};
      faultStresses.normalStress[o][i] =
          etaP * (qIMinus[o][U][i] - qIPlus[o][U][i] + qIPlus[o][N][i] * invZp +
                  qIMinus[o][N][i] * invZpNeig);

      faultStresses.traction1[o][i] =
          etaS * (qIMinus[o][V][i] - qIPlus[o][V][i] + qIPlus[o][T1][i] * invZs +
                  qIMinus[o][T1][i] * invZsNeig);

      faultStresses.traction2[o][i] =
          etaS * (qIMinus[o][W][i] - qIPlus[o][W][i] + qIPlus[o][T2][i] * invZs +
                  qIMinus[o][T2][i] * invZsNeig);
    }
  }
}

/**
 * Asserts whether all relevant arrays are properly aligned
 */
inline void checkAlignmentPostCompute(
    const real qIPlus[CONVERGENCE_ORDER][dr::misc::numQuantities][dr::misc::numPaddedPoints],
    const real qIMinus[CONVERGENCE_ORDER][dr::misc::numQuantities][dr::misc::numPaddedPoints],
    const real imposedStateP[CONVERGENCE_ORDER][dr::misc::numPaddedPoints],
    const real imposedStateM[CONVERGENCE_ORDER][dr::misc::numPaddedPoints],
    const FaultStresses& faultStresses,
    const TractionResults& tractionResults) {
  using namespace dr::misc::quantity_indices;

  assert(reinterpret_cast<uintptr_t>(imposedStateP[U]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[V]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[W]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[N]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[T1]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[T2]) % ALIGNMENT == 0);

  assert(reinterpret_cast<uintptr_t>(imposedStateM[U]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[V]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[W]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[N]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[T1]) % ALIGNMENT == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[T2]) % ALIGNMENT == 0);

  for (size_t o = 0; o < CONVERGENCE_ORDER; ++o) {
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][U]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][V]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][W]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][N]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T1]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T2]) % ALIGNMENT == 0);

    assert(reinterpret_cast<uintptr_t>(qIMinus[o][U]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][V]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][W]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][N]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T1]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T2]) % ALIGNMENT == 0);

    assert(reinterpret_cast<uintptr_t>(faultStresses.normalStress[o]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(tractionResults.traction1[o]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(tractionResults.traction2[o]) % ALIGNMENT == 0);
  }
}

/**
 * Integrate over all Time points with the time weights and calculate the traction for each side
 * according to Carsten Uphoff Thesis: EQ.: 4.60
 *
 * @param[in] faultStresses
 * @param[in] tractionResults
 * @param[in] impAndEta
 * @param[in] qInterpolatedPlus
 * @param[in] qInterpolatedMinus
 * @param[in] timeWeights
 * @param[out] imposedStatePlus
 * @param[out] imposedStateMinus
 */
template <RangeType Type = RangeType::CPU>
inline void postcomputeImposedStateFromNewStress(
    const FaultStresses& faultStresses,
    const TractionResults& tractionResults,
    const ImpedancesAndEta& impAndEta,
    real imposedStatePlus[tensor::QInterpolated::size()],
    real imposedStateMinus[tensor::QInterpolated::size()],
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const double timeWeights[CONVERGENCE_ORDER],
    unsigned startIndex = 0) {

  // set imposed state to zero
  real interPlus[tensor::QInterpolated::size()];
  real interMinus[tensor::QInterpolated::size()];

  using QInterpolatedRange = typename QInterpolated<Type>::Range;
  for (auto index = QInterpolatedRange::start; index < QInterpolatedRange::end;
       index += QInterpolatedRange::step) {
    auto i{startIndex + index};
    imposedStatePlus[i] = static_cast<real>(0.0);
    imposedStateMinus[i] = static_cast<real>(0.0);
    interPlus[i] = static_cast<real>(0.0);
    interMinus[i] = static_cast<real>(0.0);
  }

  const auto invZs = impAndEta.invZs;
  const auto invZp = impAndEta.invZp;
  const auto invZsNeig = impAndEta.invZsNeig;
  const auto invZpNeig = impAndEta.invZpNeig;

  const auto csOcpTZsOZp = impAndEta.csOcpTZsOZp;
  const auto csOcpTZsOZpNeig = impAndEta.csOcpTZsOZpNeig;

  using ImposedStateShapeT = real(*)[misc::numPaddedPoints];
  auto* interP = reinterpret_cast<ImposedStateShapeT>(interPlus);
  auto* interM = reinterpret_cast<ImposedStateShapeT>(interMinus);

  auto* imposedStateP = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
  auto* imposedStateM = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

  using QInterpolatedShapeT = const real(*)[misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  using namespace dr::misc::quantity_indices;
  unsigned DAM = 9;

#ifndef ACL_DEVICE
  checkAlignmentPostCompute(
      qIPlus, qIMinus, imposedStateP, imposedStateM, faultStresses, tractionResults);
#endif

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    auto weight = timeWeights[o];

    using NumPointsRange = typename NumPoints<Type>::Range;
#ifndef ACL_DEVICE
    #pragma omp simd
#endif
    for (auto index = NumPointsRange::start; index < NumPointsRange::end;
         index += NumPointsRange::step) {
      auto i{startIndex + index};

      const auto normalStress = faultStresses.normalStress[o][i];
      const auto traction1 = tractionResults.traction1[o][i];
      const auto traction2 = tractionResults.traction2[o][i];

      interM[N][i] += weight * normalStress;
      interM[T1][i] += weight * traction1;
      interM[T2][i] += weight * traction2;
      interM[U][i] +=
          weight * (qIMinus[o][U][i] - invZpNeig * (normalStress - qIMinus[o][N][i]));
      interM[V][i] +=
          weight * (qIMinus[o][V][i] - invZsNeig * (traction1 - qIMinus[o][T1][i]));
      interM[W][i] +=
          weight * (qIMinus[o][W][i] - invZsNeig * (traction2 - qIMinus[o][T2][i]));

      interP[N][i] += weight * normalStress;
      interP[T1][i] += weight * traction1;
      interP[T2][i] += weight * traction2;
      interP[U][i] += weight * (qIPlus[o][U][i] + invZp * (normalStress - qIPlus[o][N][i]));
      interP[V][i] += weight * (qIPlus[o][V][i] + invZs * (traction1 - qIPlus[o][T1][i]));
      interP[W][i] += weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));
    }
  }
  // Fill in nonlinear Flux term at the end time integratoon point.
  for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints;
        i ++) {
    real vx, vy ,vz; // In global coord.
    vx = impAndEta.faultN[0]*interP[U][i] +
    impAndEta.faultT1[0]*interP[V][i] + impAndEta.faultT2[0]*interP[W][i];
    vy = impAndEta.faultN[1]*interP[U][i] +
    impAndEta.faultT1[1]*interP[V][i] + impAndEta.faultT2[1]*interP[W][i];
    vz = impAndEta.faultN[2]*interP[U][i] +
    impAndEta.faultT1[2]*interP[V][i] + impAndEta.faultT2[2]*interP[W][i];
    imposedStateP[0][i] = - vx * impAndEta.faultN[0]; // minus sign from conservation laws
    imposedStateP[1][i] = - vy * impAndEta.faultN[1];
    imposedStateP[2][i] = - vz * impAndEta.faultN[2];

    imposedStateP[3][i] = - 0.5*(vx*impAndEta.faultN[1] + vy*impAndEta.faultN[0]);
    imposedStateP[4][i] = - 0.5*(vy*impAndEta.faultN[2] + vz*impAndEta.faultN[1]);
    imposedStateP[5][i] = - 0.5*(vx*impAndEta.faultN[2] + vz*impAndEta.faultN[0]);

    // Also need to project traction (N, T1, T2) in rotated coord, back to
    // global coords.
    real trac_x, trac_y ,trac_z; // In global coord.
    trac_x = impAndEta.faultN[0]*interP[N][i] +
    impAndEta.faultT1[0]*interP[T1][i] + impAndEta.faultT2[0]*interP[T2][i];
    trac_y = impAndEta.faultN[1]*interP[N][i] +
    impAndEta.faultT1[1]*interP[T1][i] + impAndEta.faultT2[1]*interP[T2][i];
    trac_z = impAndEta.faultN[2]*interP[N][i] +
    impAndEta.faultT1[2]*interP[T1][i] + impAndEta.faultT2[2]*interP[T2][i];

    // minus sign from conservation laws.
    imposedStateP[6][i] = - trac_x / impAndEta.rho0P;
    imposedStateP[7][i] = - trac_y / impAndEta.rho0P;
    imposedStateP[8][i] = - trac_z / impAndEta.rho0P;
    imposedStateP[9][i] = 0.0;
    imposedStateP[10][i] = 0.0;

    vx = (impAndEta.faultN[0]*interM[U][i] +
    impAndEta.faultT1[0]*interM[V][i] + impAndEta.faultT2[0]*interM[W][i]);
    vy = (impAndEta.faultN[1]*interM[U][i] +
    impAndEta.faultT1[1]*interM[V][i] + impAndEta.faultT2[1]*interM[W][i]);
    vz = (impAndEta.faultN[2]*interM[U][i] +
    impAndEta.faultT1[2]*interM[V][i] + impAndEta.faultT2[2]*interM[W][i]);
    // additional minus sign due to projected face-normal coords
    // are opposite to face normal of the cell "M".
    imposedStateM[0][i] = -1.0*(- vx * impAndEta.faultN[0]);
    imposedStateM[1][i] = -1.0*(- vy * impAndEta.faultN[1]);
    imposedStateM[2][i] = -1.0*(- vz * impAndEta.faultN[2]);

    imposedStateM[3][i] = -1.0*(-0.5*(vx*impAndEta.faultN[1] + vy*impAndEta.faultN[0]));
    imposedStateM[4][i] = -1.0*(-0.5*(vy*impAndEta.faultN[2] + vz*impAndEta.faultN[1]));
    imposedStateM[5][i] = -1.0*(-0.5*(vx*impAndEta.faultN[2] + vz*impAndEta.faultN[0]));

    // minus sign due to traction needs to time (-1,0,0) in roated coords for cell "M"
    trac_x = - (impAndEta.faultN[0]*interM[N][i] +
    impAndEta.faultT1[0]*interM[T1][i] + impAndEta.faultT2[0]*interM[T2][i]);
    trac_y = - (impAndEta.faultN[1]*interM[N][i] +
    impAndEta.faultT1[1]*interM[T1][i] + impAndEta.faultT2[1]*interM[T2][i]);
    trac_z = - (impAndEta.faultN[2]*interM[N][i] +
    impAndEta.faultT1[2]*interM[T1][i] + impAndEta.faultT2[2]*interM[T2][i]);

    // minus sign from conservation laws.
    imposedStateM[6][i] = - trac_x / impAndEta.rho0M;
    imposedStateM[7][i] = - trac_y / impAndEta.rho0M;
    imposedStateM[8][i] = - trac_z / impAndEta.rho0M;
    imposedStateM[9][i] = 0.0;
    imposedStateM[10][i] = 0.0;
  }
}

/**
 * adjusts initial stresses based on the given nucleation ones
 *
 * @param[out] initialStressInFaultCS
 * @param[in] nucleationStressInFaultCS
 * @param[in] t0
 * @param[in] dt
 * @param[in] index - device iteration index
 */
template <RangeType Type = RangeType::CPU,
          typename MathFunctions = seissol::functions::HostStdFunctions>
inline void adjustInitialStress(real initialStressInFaultCS[misc::numPaddedPoints][6],
                                const real nucleationStressInFaultCS[misc::numPaddedPoints][6],
                                real fullUpdateTime,
                                real t0,
                                real dt,
                                unsigned startIndex = 0) {
  if (fullUpdateTime <= t0) {
    const real gNuc =
        gaussianNucleationFunction::smoothStepIncrement<MathFunctions>(fullUpdateTime, dt, t0);

    using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
    #pragma omp simd
#endif
    for (auto index = Range::start; index < Range::end; index += Range::step) {
      auto pointIndex{startIndex + index};
      for (unsigned i = 0; i < 6; i++) {
        initialStressInFaultCS[pointIndex][i] += nucleationStressInFaultCS[pointIndex][i] * gNuc;
      }
    }
  }
}

/**
 * output rupture front, saves update time of the rupture front
 * rupture front is the first registered change in slip rates that exceeds 0.001
 *
 * param[in,out] ruptureTimePending
 * param[out] ruptureTime
 * param[in] slipRateMagnitude
 * param[in] fullUpdateTime
 */
template <RangeType Type = RangeType::CPU>
inline void saveRuptureFrontOutput(bool ruptureTimePending[misc::numPaddedPoints],
                                   real ruptureTime[misc::numPaddedPoints],
                                   const real slipRateMagnitude[misc::numPaddedPoints],
                                   real fullUpdateTime,
                                   unsigned startIndex = 0) {

  using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
  #pragma omp simd
#endif
  for (auto index = Range::start; index < Range::end; index += Range::step) {
    auto pointIndex{startIndex + index};
    constexpr real ruptureFrontThreshold = 0.001;
    if (ruptureTimePending[pointIndex] && slipRateMagnitude[pointIndex] > ruptureFrontThreshold) {
      ruptureTime[pointIndex] = fullUpdateTime;
      ruptureTimePending[pointIndex] = false;
    }
  }
}

/**
 * Save the maximal computed slip rate magnitude in peakSlipRate
 *
 * param[in] slipRateMagnitude
 * param[in, out] peakSlipRate
 */
template <RangeType Type = RangeType::CPU>
inline void savePeakSlipRateOutput(real slipRateMagnitude[misc::numPaddedPoints],
                                   real peakSlipRate[misc::numPaddedPoints],
                                   unsigned startIndex = 0) {

  using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
  #pragma omp simd
#endif
  for (auto index = Range::start; index < Range::end; index += Range::step) {
    auto pointIndex{startIndex + index};
    peakSlipRate[pointIndex] = std::max(peakSlipRate[pointIndex], slipRateMagnitude[pointIndex]);
  }
}

template <RangeType Type = RangeType::CPU>
inline void computeFrictionEnergy(
    DREnergyOutput& energyData,
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const ImpedancesAndEta& impAndEta,
    const double timeWeights[CONVERGENCE_ORDER],
    const real spaceWeights[NUMBER_OF_SPACE_QUADRATURE_POINTS],
    const DRGodunovData& godunovData,
    size_t startIndex = 0) {

  auto* slip = reinterpret_cast<real(*)[misc::numPaddedPoints]>(energyData.slip);
  auto* accumulatedSlip = energyData.accumulatedSlip;
  auto* frictionalEnergy = energyData.frictionalEnergy;
  const double doubledSurfaceArea = godunovData.doubledSurfaceArea;

  using QInterpolatedShapeT = const real(*)[misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  const auto aPlus = impAndEta.etaP * impAndEta.invZp;
  const auto bPlus = impAndEta.etaS * impAndEta.invZs;

  const auto aMinus = impAndEta.etaP * impAndEta.invZpNeig;
  const auto bMinus = impAndEta.etaS * impAndEta.invZsNeig;

  using Range = typename NumPoints<Type>::Range;

  using namespace dr::misc::quantity_indices;
  for (size_t o = 0; o < CONVERGENCE_ORDER; ++o) {
    const auto timeWeight = timeWeights[o];

#ifndef ACL_DEVICE
    #pragma omp simd
#endif
    for (size_t index = Range::start; index < Range::end; index += Range::step) {
      const size_t i{startIndex + index};

      const real interpolatedSlipRate1 = qIMinus[o][U][i] - qIPlus[o][U][i];
      const real interpolatedSlipRate2 = qIMinus[o][V][i] - qIPlus[o][V][i];
      const real interpolatedSlipRate3 = qIMinus[o][W][i] - qIPlus[o][W][i];

      const real interpolatedSlipRateMagnitude =
          misc::magnitude(interpolatedSlipRate1, interpolatedSlipRate2, interpolatedSlipRate3);

      accumulatedSlip[i] += timeWeight * interpolatedSlipRateMagnitude;

      slip[0][i] += timeWeight * interpolatedSlipRate1;
      slip[1][i] += timeWeight * interpolatedSlipRate2;
      slip[2][i] += timeWeight * interpolatedSlipRate3;

      const real interpolatedTraction11 = aPlus * qIMinus[o][XX][i] + aMinus * qIPlus[o][XX][i];
      const real interpolatedTraction12 = bPlus * qIMinus[o][XY][i] + bMinus * qIPlus[o][XY][i];
      const real interpolatedTraction13 = bPlus * qIMinus[o][XZ][i] + bMinus * qIPlus[o][XZ][i];

      const auto spaceWeight = spaceWeights[i];
      const auto weight = -1.0 * timeWeight * spaceWeight * doubledSurfaceArea;
      frictionalEnergy[i] += weight * (interpolatedTraction11 * interpolatedSlipRate1 +
                                       interpolatedTraction12 * interpolatedSlipRate2 +
                                       interpolatedTraction13 * interpolatedSlipRate3);
    }
  }
}

} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_FRICTIONSOLVER_COMMON_H
