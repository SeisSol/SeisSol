#ifndef SEISSOL_FRICTIONSOLVER_COMMON_H
#define SEISSOL_FRICTIONSOLVER_COMMON_H

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/DynamicRupture.h"
#include "Kernels/DynamicRupture.h"

/**
 * Contains common functions required both for CPU and GPU impl.
 * of Dynamic Rupture solvers. The functions placed in
 * this class definition (of the header file) result
 * in the function inlining required for GPU impl.
 */
namespace seissol::dr::friction_law::common {

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
inline void precomputeStressFromQInterpolated(
    FaultStresses& faultStresses,
    const ImpedancesAndEta& impAndEta,
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()]) {

  static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                "Different number of quadrature points?");

  const auto etaP = impAndEta.etaP;
  const auto etaS = impAndEta.etaS;
  const auto invZp = impAndEta.invZp;
  const auto invZs = impAndEta.invZs;
  const auto invZpNeig = impAndEta.invZpNeig;
  const auto invZsNeig = impAndEta.invZsNeig;

  using QInterpolatedShapeT = const real(*)[misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
  auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));

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

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#else
#pragma omp simd
#endif // ACL_DEVICE_OFFLOAD
    for (unsigned i = 0; i < misc::numPaddedPoints; ++i) {
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
inline void postcomputeImposedStateFromNewStress(
    const FaultStresses& faultStresses,
    const TractionResults& tractionResults,
    const ImpedancesAndEta& impAndEta,
    real imposedStatePlus[tensor::QInterpolated::size()],
    real imposedStateMinus[tensor::QInterpolated::size()],
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const double timeWeights[CONVERGENCE_ORDER]) {

  // set imposed state to zero
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#endif // ACL_DEVICE_OFFLOAD
  for (unsigned int i = 0; i < tensor::QInterpolated::size(); i++) {
    imposedStatePlus[i] = static_cast<real>(0.0);
    imposedStateMinus[i] = static_cast<real>(0.0);
  }

  const auto invZs = impAndEta.invZs;
  const auto invZp = impAndEta.invZp;
  const auto invZsNeig = impAndEta.invZsNeig;
  const auto invZpNeig = impAndEta.invZpNeig;

  using ImposedStateShapeT = real(*)[misc::numPaddedPoints];
  auto* imposedStateP = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
  auto* imposedStateM = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

  using QInterpolatedShapeT = const real(*)[misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

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

    assert(reinterpret_cast<uintptr_t>(faultStresses.normalStress[o]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(tractionResults.traction1[o]) % ALIGNMENT == 0);
    assert(reinterpret_cast<uintptr_t>(tractionResults.traction2[o]) % ALIGNMENT == 0);
    auto weight = timeWeights[o];

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#else
#pragma omp simd
#endif // ACL_DEVICE_OFFLOAD
    for (unsigned i = 0; i < misc::numPaddedPoints; ++i) {
      const auto normalStress = faultStresses.normalStress[o][i];
      const auto traction1 = tractionResults.traction1[o][i];
      const auto traction2 = tractionResults.traction2[o][i];

      imposedStateM[N][i] += weight * normalStress;
      imposedStateM[T1][i] += weight * traction1;
      imposedStateM[T2][i] += weight * traction2;
      imposedStateM[U][i] +=
          weight * (qIMinus[o][U][i] - invZpNeig * (normalStress - qIMinus[o][N][i]));
      imposedStateM[V][i] +=
          weight * (qIMinus[o][V][i] - invZsNeig * (traction1 - qIMinus[o][T1][i]));
      imposedStateM[W][i] +=
          weight * (qIMinus[o][W][i] - invZsNeig * (traction2 - qIMinus[o][T2][i]));

      imposedStateP[N][i] += weight * normalStress;
      imposedStateP[T1][i] += weight * traction1;
      imposedStateP[T2][i] += weight * traction2;
      imposedStateP[U][i] += weight * (qIPlus[o][U][i] + invZp * (normalStress - qIPlus[o][N][i]));
      imposedStateP[V][i] += weight * (qIPlus[o][V][i] + invZs * (traction1 - qIPlus[o][T1][i]));
      imposedStateP[W][i] += weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));
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
inline void saveRuptureFrontOutput(bool ruptureTimePending[misc::numPaddedPoints],
                                   real ruptureTime[misc::numPaddedPoints],
                                   const real slipRateMagnitude[misc::numPaddedPoints],
                                   real fullUpdateTime) {
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#else
#pragma omp simd
#endif // ACL_DEVICE_OFFLOAD
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
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
inline void savePeakSlipRateOutput(real slipRateMagnitude[misc::numPaddedPoints],
                                   real peakSlipRate[misc::numPaddedPoints]) {

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#else
#pragma omp simd
#endif // ACL_DEVICE_OFFLOAD
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    peakSlipRate[pointIndex] = std::max(peakSlipRate[pointIndex], slipRateMagnitude[pointIndex]);
  }
}
} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_FRICTIONSOLVER_COMMON_H
