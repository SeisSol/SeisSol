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
  static constexpr size_t size{End - Start};
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
    const ImpedancesAndEta& impAndEta,
    const real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    unsigned startLoopIndex = 0) {

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

#ifndef ACL_DEVICE
  checkAlignmentPreCompute(qIPlus, qIMinus, faultStresses);
#endif

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
  using QInterpolatedRange = typename QInterpolated<Type>::Range;
  for (auto index = QInterpolatedRange::start; index < QInterpolatedRange::end;
       index += QInterpolatedRange::step) {
    auto i{startIndex + index};
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

#ifndef ACL_DEVICE
  checkAlignmentPostCompute(
      qIPlus, qIMinus, imposedStateP, imposedStateM, faultStresses, tractionResults);
#endif

  using NumPointsRange = typename NumPoints<Type>::Range;
  real cachedImposedStateM[dr::misc::numQuantities][NumPointsRange::size];
  real cachedImposedStateP[dr::misc::numQuantities][NumPointsRange::size];

  for (auto index = NumPointsRange::start; index < NumPointsRange::end;
       index += NumPointsRange::step) {
    auto i{startIndex + index};
#pragma unroll
    for (unsigned q = 0; q < dr::misc::numQuantities; ++q) {
      cachedImposedStateM[q][index] = imposedStateM[q][i];
      cachedImposedStateP[q][index] = imposedStateP[q][i];
    }
  }

  for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
    auto weight = static_cast<real>(timeWeights[o]);

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = NumPointsRange::start; index < NumPointsRange::end;
         index += NumPointsRange::step) {
      auto i{startIndex + index};

      const auto normalStress = faultStresses.normalStress[o][i];
      const auto traction1 = tractionResults.traction1[o][i];
      const auto traction2 = tractionResults.traction2[o][i];

      cachedImposedStateM[N][index] += weight * normalStress;
      cachedImposedStateM[T1][index] += weight * traction1;
      cachedImposedStateM[T2][index] += weight * traction2;
      cachedImposedStateM[U][index] +=
          weight * (qIMinus[o][U][i] - invZpNeig * (normalStress - qIMinus[o][N][i]));
      cachedImposedStateM[V][index] +=
          weight * (qIMinus[o][V][i] - invZsNeig * (traction1 - qIMinus[o][T1][i]));
      cachedImposedStateM[W][index] +=
          weight * (qIMinus[o][W][i] - invZsNeig * (traction2 - qIMinus[o][T2][i]));

      cachedImposedStateP[N][index] += weight * normalStress;
      cachedImposedStateP[T1][index] += weight * traction1;
      cachedImposedStateP[T2][index] += weight * traction2;
      cachedImposedStateP[U][index] +=
          weight * (qIPlus[o][U][i] + invZp * (normalStress - qIPlus[o][N][i]));
      cachedImposedStateP[V][index] +=
          weight * (qIPlus[o][V][i] + invZs * (traction1 - qIPlus[o][T1][i]));
      cachedImposedStateP[W][index] +=
          weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));
    }
  }

  for (auto index = NumPointsRange::start; index < NumPointsRange::end;
       index += NumPointsRange::step) {
    auto i{startIndex + index};
#pragma unroll
    for (unsigned q = 0; q < dr::misc::numQuantities; ++q) {
      imposedStateM[q][i] = cachedImposedStateM[q][index];
      imposedStateP[q][i] = cachedImposedStateP[q][index];
    }
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
  real cachedAccumulatedSlip[Range::size];
  real cachedFrictionalEnergy[Range::size];
  real cachedSlip[3][Range::size];

  for (auto index = Range::start; index < Range::end; index += Range::step) {
    auto i{startIndex + index};
    cachedAccumulatedSlip[index] = accumulatedSlip[i];
    cachedFrictionalEnergy[index] = frictionalEnergy[i];
#pragma unroll
    for (unsigned d = 0; d < 3; ++d) {
      cachedSlip[d][index] = slip[d][i];
    }
  }

  using namespace dr::misc::quantity_indices;
  for (size_t o = 0; o < CONVERGENCE_ORDER; ++o) {
    const auto timeWeight = static_cast<real>(timeWeights[o]);

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

      cachedAccumulatedSlip[i] += timeWeight * interpolatedSlipRateMagnitude;

      cachedSlip[0][i] += timeWeight * interpolatedSlipRate1;
      cachedSlip[1][i] += timeWeight * interpolatedSlipRate2;
      cachedSlip[2][i] += timeWeight * interpolatedSlipRate3;

      const real interpolatedTraction11 = aPlus * qIMinus[o][XX][i] + aMinus * qIPlus[o][XX][i];
      const real interpolatedTraction12 = bPlus * qIMinus[o][XY][i] + bMinus * qIPlus[o][XY][i];
      const real interpolatedTraction13 = bPlus * qIMinus[o][XZ][i] + bMinus * qIPlus[o][XZ][i];

      const auto spaceWeight = spaceWeights[i];
      const auto weight = -1.0 * timeWeight * spaceWeight * doubledSurfaceArea;
      cachedFrictionalEnergy[i] += weight * (interpolatedTraction11 * interpolatedSlipRate1 +
                                             interpolatedTraction12 * interpolatedSlipRate2 +
                                             interpolatedTraction13 * interpolatedSlipRate3);
    }
  }

  for (auto index = Range::start; index < Range::end; index += Range::step) {
    auto i{startIndex + index};
    accumulatedSlip[i] = cachedAccumulatedSlip[index];
    frictionalEnergy[i] = cachedFrictionalEnergy[index];
#pragma unroll
    for (unsigned d = 0; d < 3; ++d) {
      slip[d][i] = cachedSlip[d][index];
    }
  }
}

} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_FRICTIONSOLVER_COMMON_H
