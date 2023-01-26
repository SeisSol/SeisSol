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
 * Copes an eigen3 matrix to a 2D yateto tensor
 */
template <typename T, int dim1, int dim2>
void copyEigenToYateto(Eigen::Matrix<T, dim1, dim2> const& matrix,
                       yateto::DenseTensorView<2, T>& tensorView) {
  assert(tensorView.shape(0) == dim1);
  assert(tensorView.shape(1) == dim2);

  tensorView.setZero();
  for (size_t row = 0; row < dim1; ++row) {
    for (size_t col = 0; col < dim2; ++col) {
      tensorView(row, col) = matrix(row, col);
    }
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
    const ImpedanceMatrices& impedanceMatrices,
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    unsigned startLoopIndex = 0) {

  static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                "Different number of quadrature points?");

  seissol::dynamicRupture::kernel::computeTheta krnl;
  krnl.extractVelocities = init::extractVelocities::Values;
  krnl.extractTractions = init::extractTractions::Values;

  // Compute Theta from eq (4.53) in Carsten's thesis
  krnl.Zplus = impedanceMatrices.impedance;
  krnl.Zminus = impedanceMatrices.impedanceNeig;
  krnl.eta = impedanceMatrices.eta;

  alignas(ALIGNMENT) real thetaBuffer[tensor::theta::size()] = {};
  krnl.theta = thetaBuffer;
  auto thetaView = init::theta::view::create(thetaBuffer);

  // TODO: Integrate loop over o into the kernel
  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    krnl.Qplus = qInterpolatedPlus[o];
    krnl.Qminus = qInterpolatedMinus[o];
    krnl.execute();

    // TODO: what about GPUs/Ranges
    for (unsigned i = 0; i < misc::numPaddedPoints; ++i) {
      faultStresses.normalStress[o][i] = thetaView(i, 0);
      faultStresses.traction1[o][i] = thetaView(i, 1);
      faultStresses.traction2[o][i] = thetaView(i, 2);
#ifdef USE_POROELASTIC
      faultStresses.fluidPressure[o][i] = thetaView(i, 3);
#else
      faultStresses.fluidPressure[o][i] = 0.0;
#endif
      faultStresses.effectiveNormalStress[o][i] =
          faultStresses.normalStress[o][i] - faultStresses.fluidPressure[o][i];
    }
  }
}

/**
 * Integrate over all Time points with the time weights and calculate the traction for each side
 * according to Carsten Uphoff Thesis: EQ.: 4.60
 *
 * @param[in] faultStresses
 * @param[in] tractionResults
 * @param[in] impedancenceMatrices
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
    const ImpedanceMatrices& impedanceMatrices,
    real imposedStatePlus[tensor::QInterpolated::size()],
    real imposedStateMinus[tensor::QInterpolated::size()],
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const double timeWeights[CONVERGENCE_ORDER],
    unsigned startIndex = 0) {

  // set imposed state to zero
  using QInterpolatedRange = typename QInterpolated<Type>::Range;
  for (auto index = QInterpolatedRange::start; index < QInterpolatedRange::end;
       index += QInterpolatedRange::step) {
    auto i{startIndex + index};
    imposedStatePlus[i] = static_cast<real>(0.0);
    imposedStateMinus[i] = static_cast<real>(0.0);
  }

  // setup kernel
  seissol::dynamicRupture::kernel::computeImposedStateM krnlM;
  krnlM.extractVelocities = init::extractVelocities::Values;
  krnlM.extractTractions = init::extractTractions::Values;
  krnlM.mapToVelocities = init::mapToVelocities::Values;
  krnlM.mapToTractions = init::mapToTractions::Values;
  krnlM.Zminus = impedanceMatrices.impedanceNeig;
  krnlM.imposedState = imposedStateMinus;

  seissol::dynamicRupture::kernel::computeImposedStateP krnlP;
  krnlP.extractVelocities = init::extractVelocities::Values;
  krnlP.extractTractions = init::extractTractions::Values;
  krnlP.mapToVelocities = init::mapToVelocities::Values;
  krnlP.mapToTractions = init::mapToTractions::Values;
  krnlP.Zplus = impedanceMatrices.impedance;
  krnlP.imposedState = imposedStatePlus;

  alignas(ALIGNMENT) real thetaBuffer[tensor::theta::size()] = {};
  auto thetaView = init::theta::view::create(thetaBuffer);
  krnlM.theta = thetaBuffer;
  krnlP.theta = thetaBuffer;

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    auto weight = timeWeights[o];
    // copy values to yateto dataformat
    for (unsigned i = 0; i < misc::numPaddedPoints; ++i) {
      thetaView(i, 0) = faultStresses.normalStress[o][i];
      thetaView(i, 1) = tractionResults.traction1[o][i];
      thetaView(i, 2) = tractionResults.traction2[o][i];
#ifdef USE_POROELASTIC
      thetaView(i, 3) = faultStresses.fluidPressure[o][i];
#endif
    }
    // execute kernel (and hence update imposedStatePlus/Minus)
    krnlM.Qminus = qInterpolatedMinus[o];
    krnlM.weight = weight;
    krnlM.execute();

    krnlP.Qplus = qInterpolatedPlus[o];
    krnlP.weight = weight;
    krnlP.execute();
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
// See https://github.com/llvm/llvm-project/issues/60163
// NOLINTNEXTLINE
inline void saveRuptureFrontOutput(bool ruptureTimePending[misc::numPaddedPoints],
                                   // See https://github.com/llvm/llvm-project/issues/60163
                                   // NOLINTNEXTLINE
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
inline void savePeakSlipRateOutput(const real slipRateMagnitude[misc::numPaddedPoints],
                                   // See https://github.com/llvm/llvm-project/issues/60163
                                   // NOLINTNEXTLINE
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
