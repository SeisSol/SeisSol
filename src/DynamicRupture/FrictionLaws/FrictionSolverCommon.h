#ifndef SEISSOL_FRICTIONSOLVER_COMMON_H
#define SEISSOL_FRICTIONSOLVER_COMMON_H

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "Initializer/DynamicRupture.h"
#include "Kernels/DynamicRupture.h"

namespace seissol::dr::friction_law {
struct Common {
  /**
   * Contains common functions required both for CPU and GPU impl.
   * of Dynamic Rupture solvers. The functions placed in
   * this class definition (of the header file) result
   * in the function inlining required for GPU impl.
   */

  public:
  /**
   * Calculate traction and normal stress at the interface of a face.
   * Using equations (A2) from Pelties et al. 2014
   * Definiton of eta and impedance Z are found in dissertation of Carsten Uphoff
   *
   * @param[out] faultStresses contains normalStress, traction1, traction2
   *             at the 2d face quadrature nodes evaluated at the time
   *             quadrature points
   * @param[in] impAndEta contains eta and impedance values
   * @param[in] qInterpolatedPlus a plus side dofs interpolated at time sub-intervals
   * @param[in] qInterpolatedMinus a minus side dofs interpolated at time sub-intervals
   */
  static void precomputeStressFromQInterpolated(
      FaultStresses& faultStresses,
      const ImpedancesAndEta& impAndEta,
      const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()]) {

    static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                  "Different number of quadrature points?");

    // this initialization of the kernel could be moved to the initializer,
    // since all inputs outside the j-loop are time independent
    // set inputParam could be extendent for this
    // the kernel then could be a class attribute (but be careful of race conditions since this is
    // computed in parallel!!)

    const auto etaP = impAndEta.etaP;
    const auto etaS = impAndEta.etaS;
    const auto invZp = impAndEta.invZp;
    const auto invZs = impAndEta.invZs;
    const auto invZpNeig = impAndEta.invZpNeig;
    const auto invZsNeig = impAndEta.invZsNeig;

    using QInterpolatedShapeT = const real(*)[misc::numQuantities][misc::numPaddedPoints];
    auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
    auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));

    for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#endif // ACL_DEVICE_OFFLOAD
      for (unsigned i = 0; i < misc::numPaddedPoints; ++i) {
        faultStresses.normalStress[o][i] =
            etaP * (qIMinus[o][6][i] - qIPlus[o][6][i] + qIPlus[o][0][i] * invZp +
                    qIMinus[o][0][i] * invZpNeig);

        faultStresses.traction1[o][i] =
            etaS * (qIMinus[o][7][i] - qIPlus[o][7][i] + qIPlus[o][3][i] * invZs +
                    qIMinus[o][3][i] * invZsNeig);

        faultStresses.traction2[o][i] =
            etaS * (qIMinus[o][8][i] - qIPlus[o][8][i] + qIPlus[o][5][i] * invZs +
                    qIMinus[o][5][i] * invZsNeig);
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
  static void postcomputeImposedStateFromNewStress(
      const FaultStresses& faultStresses,
      const TractionResults& tractionResults,
      const ImpedancesAndEta& impAndEta,
      real imposedStatePlus[tensor::QInterpolated::size()],
      real imposedStateMinus[tensor::QInterpolated::size()],
      const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      const double timeWeights[CONVERGENCE_ORDER]) {

    // this initialization of the kernel could be moved to the initializer
    // set inputParam could be extendent for this (or create own function)
    // the kernel then could be a class attribute and following values are only set once
    //(but be careful of race conditions since this is computed in parallel for each face!!)

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

    for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
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

        imposedStateM[0][i] += weight * normalStress;
        imposedStateM[3][i] += weight * traction1;
        imposedStateM[5][i] += weight * traction2;
        imposedStateM[6][i] +=
            weight * (qIMinus[o][6][i] - invZpNeig * (normalStress - qIMinus[o][0][i]));
        imposedStateM[7][i] +=
            weight * (qIMinus[o][7][i] - invZsNeig * (traction1 - qIMinus[o][3][i]));
        imposedStateM[8][i] +=
            weight * (qIMinus[o][8][i] - invZsNeig * (traction2 - qIMinus[o][5][i]));

        imposedStateP[0][i] += weight * normalStress;
        imposedStateP[3][i] += weight * traction1;
        imposedStateP[5][i] += weight * traction2;
        imposedStateP[6][i] +=
            weight * (qIPlus[o][6][i] + invZp * (normalStress - qIPlus[o][0][i]));
        imposedStateP[7][i] += weight * (qIPlus[o][7][i] + invZs * (traction1 - qIPlus[o][3][i]));
        imposedStateP[8][i] += weight * (qIPlus[o][8][i] + invZs * (traction2 - qIPlus[o][5][i]));
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
  static void saveRuptureFrontOutput(bool ruptureTimePending[misc::numPaddedPoints],
                                     real ruptureTime[misc::numPaddedPoints],
                                     const real slipRateMagnitude[misc::numPaddedPoints],
                                     real fullUpdateTime) {
#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
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
  static void savePeakSlipRateOutput(real slipRateMagnitude[misc::numPaddedPoints],
                                     real peakSlipRate[misc::numPaddedPoints]) {

#ifdef ACL_DEVICE_OFFLOAD
#pragma omp loop bind(parallel)
#endif // ACL_DEVICE_OFFLOAD
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      peakSlipRate[pointIndex] = std::max(peakSlipRate[pointIndex], slipRateMagnitude[pointIndex]);
    }
  }

  /**
   * Compute and store element-averaged slip to determine the magnitude of an earthquake.
   * In calc_seissol.f90 this value will be multiplied by the element surface
   * and the seismic moment is outputted once at the end of the simulation.
   *
   * @param[in] tmpSlip
   * @param[out] averagedSlip
   */
  static void saveAverageSlipOutput(real tmpSlip[misc::numPaddedPoints], real& averagedSlip) {
    real sumOfTmpSlip = 0;
    for (unsigned pointIndex = 0; pointIndex < misc::numberOfBoundaryGaussPoints; pointIndex++) {
      sumOfTmpSlip += tmpSlip[pointIndex];
    }
    averagedSlip += sumOfTmpSlip / misc::numberOfBoundaryGaussPoints;
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_FRICTIONSOLVER_COMMON_H
