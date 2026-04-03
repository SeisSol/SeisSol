// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_

#include "Common/Constants.h"
#include "Common/Executor.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "Initializer/Typedefs.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Solver/MultipleSimulations.h"

#include <cmath>
#include <limits>
#include <type_traits>

/**
 * Contains common functions required both for CPU and GPU impl.
 * of Dynamic Rupture solvers. The functions placed in
 * this class definition (of the header file) result
 * in the function inlining required for GPU impl.
 */
namespace seissol::dr::friction_law::common {

template <uint32_t StartT, uint32_t EndT, uint32_t StepT>
struct ForLoopRange {
  static constexpr uint32_t Start{StartT};
  static constexpr uint32_t End{EndT};
  static constexpr uint32_t Step{StepT};
  static constexpr uint32_t Size{EndT - StartT};
};

enum class RangeType { CPU, GPU };

template <RangeType Type>
struct NumPoints {
  private:
  using CpuRange = ForLoopRange<0, dr::misc::NumPaddedPoints, 1>;
  using GpuRange = ForLoopRange<0, 1, 1>;

  public:
  // Range::Start is 0, and Range::End is seissol::misc::NumPaddedPoints for CPU
  using Range = std::conditional_t<Type == RangeType::CPU, CpuRange, GpuRange>;
};

template <RangeType Type>
struct QInterpolated {
  private:
  using CpuRange = ForLoopRange<0, tensor::QInterpolated::size(), 1>;
  using GpuRange = ForLoopRange<0, tensor::QInterpolated::size(), misc::NumPaddedPoints>;

  public:
  using Range = std::conditional_t<Type == RangeType::CPU, CpuRange, GpuRange>;
};

template <RangeType Type>
struct RangeExecutor;

template <>
struct RangeExecutor<RangeType::CPU> {
  static constexpr Executor Exec = Executor::Host;
};

template <>
struct RangeExecutor<RangeType::GPU> {
  static constexpr Executor Exec = Executor::Device;
};

template <Executor Executor>
struct VariableIndexing;

template <>
struct VariableIndexing<Executor::Host> {
  static constexpr real& index(real (&data)[misc::TimeSteps][misc::NumPaddedPoints], int o, int i) {
    return data[o][i];
  }

  static constexpr real
      index(const real (&data)[misc::TimeSteps][misc::NumPaddedPoints], int o, int i) {
    return data[o][i];
  }
};

template <>
struct VariableIndexing<Executor::Device> {
  static constexpr real& index(real (&data)[misc::TimeSteps], int o, int /*i*/) { return data[o]; }

  static constexpr real index(const real (&data)[misc::TimeSteps], int o, int /*i*/) {
    return data[o];
  }
};

/**
 * Asserts whether all relevant arrays are properly aligned
 */
inline void checkAlignmentPreCompute(
    [[maybe_unused]] const real qIPlus[misc::TimeSteps][dr::misc::NumQuantities]
                                      [dr::misc::NumPaddedPoints],
    [[maybe_unused]] const real qIMinus[misc::TimeSteps][dr::misc::NumQuantities]
                                       [dr::misc::NumPaddedPoints],
    [[maybe_unused]] const FaultStresses<Executor::Host>& faultStresses) {
  using namespace dr::misc::quantity_indices;
  for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][U]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][V]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][W]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][N]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T1]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T2]) % Alignment == 0);

    assert(reinterpret_cast<uintptr_t>(qIMinus[o][U]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][V]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][W]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][N]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T1]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T2]) % Alignment == 0);

    assert(reinterpret_cast<uintptr_t>(faultStresses.normalStress[o]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(faultStresses.traction1[o]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(faultStresses.traction2[o]) % Alignment == 0);
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
 * @param[in] impedanceMatrices contains impedance and eta values, in the poroelastic case, these
 * are non-diagonal matrices
 * @param[in] qInterpolatedPlus a plus side dofs interpolated at time sub-intervals
 * @param[in] qInterpolatedMinus a minus side dofs interpolated at time sub-intervals
 */
template <RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void precomputeStressFromQInterpolated(
    FaultStresses<RangeExecutor<Type>::Exec>& faultStresses,
    const ImpedancesAndEta& impAndEta,
    [[maybe_unused]] const ImpedanceMatrices& impedanceMatrices,
    const real qInterpolatedPlus[misc::TimeSteps][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[misc::TimeSteps][tensor::QInterpolated::size()],
    real etaPDamp,
    uint32_t startLoopIndex = 0) {
  static_assert(tensor::QInterpolated::Shape[seissol::multisim::BasisFunctionDimension] ==
                    tensor::resample::Shape[0],
                "Different number of quadrature points?");

#ifndef USE_POROELASTIC
  const auto etaP = impAndEta.etaP * etaPDamp;
  const auto etaS = impAndEta.etaS;
  const auto invZp = impAndEta.invZp;
  const auto invZs = impAndEta.invZs;
  const auto invZpNeig = impAndEta.invZpNeig;
  const auto invZsNeig = impAndEta.invZsNeig;

  using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
  const auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
  const auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));

  using namespace dr::misc::quantity_indices;

#ifndef ACL_DEVICE
  checkAlignmentPreCompute(qIPlus, qIMinus, faultStresses);
#endif

  for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
    using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto i{startLoopIndex + index};
      VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i) =
          etaP * (qIMinus[o][U][i] - qIPlus[o][U][i] + qIPlus[o][N][i] * invZp +
                  qIMinus[o][N][i] * invZpNeig);

      VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction1, o, i) =
          etaS * (qIMinus[o][V][i] - qIPlus[o][V][i] + qIPlus[o][T1][i] * invZs +
                  qIMinus[o][T1][i] * invZsNeig);

      VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction2, o, i) =
          etaS * (qIMinus[o][W][i] - qIPlus[o][W][i] + qIPlus[o][T2][i] * invZs +
                  qIMinus[o][T2][i] * invZsNeig);
    }
  }
#else
  seissol::dynamicRupture::kernel::computeTheta krnl;
  krnl.extractVelocities = init::extractVelocities::Values;
  krnl.extractTractions = init::extractTractions::Values;

  // Compute Theta from eq (4.53) in Carsten's thesis
  krnl.Zplus = impedanceMatrices.impedance;
  krnl.Zminus = impedanceMatrices.impedanceNeig;
  krnl.eta = impedanceMatrices.eta;

  alignas(Alignment) real thetaBuffer[tensor::theta::size()] = {};
  krnl.theta = thetaBuffer;
  auto thetaView = init::theta::view::create(thetaBuffer);

  for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
    krnl.Qplus = qInterpolatedPlus[o];
    krnl.Qminus = qInterpolatedMinus[o];
    krnl.execute();

    for (uint32_t i = 0; i < misc::NumPaddedPoints; ++i) {
      faultStresses.normalStress[o][i] = thetaView(i, 0);
      faultStresses.traction1[o][i] = thetaView(i, 1);
      faultStresses.traction2[o][i] = thetaView(i, 2);
      faultStresses.fluidPressure[o][i] = thetaView(i, 3);
    }
  }
#endif
}

/**
 * Asserts whether all relevant arrays are properly aligned
 */
inline void checkAlignmentPostCompute(
    [[maybe_unused]] const real qIPlus[misc::TimeSteps][dr::misc::NumQuantities]
                                      [dr::misc::NumPaddedPoints],
    [[maybe_unused]] const real qIMinus[misc::TimeSteps][dr::misc::NumQuantities]
                                       [dr::misc::NumPaddedPoints],
    [[maybe_unused]] const real imposedStateP[misc::TimeSteps][dr::misc::NumPaddedPoints],
    [[maybe_unused]] const real imposedStateM[misc::TimeSteps][dr::misc::NumPaddedPoints],
    [[maybe_unused]] const FaultStresses<Executor::Host>& faultStresses,
    [[maybe_unused]] const TractionResults<Executor::Host>& tractionResults) {
  using namespace dr::misc::quantity_indices;

  assert(reinterpret_cast<uintptr_t>(imposedStateP[U]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[V]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[W]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[N]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[T1]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateP[T2]) % Alignment == 0);

  assert(reinterpret_cast<uintptr_t>(imposedStateM[U]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[V]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[W]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[N]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[T1]) % Alignment == 0);
  assert(reinterpret_cast<uintptr_t>(imposedStateM[T2]) % Alignment == 0);

  for (size_t o = 0; o < misc::TimeSteps; ++o) {
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][U]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][V]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][W]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][N]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T1]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIPlus[o][T2]) % Alignment == 0);

    assert(reinterpret_cast<uintptr_t>(qIMinus[o][U]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][V]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][W]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][N]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T1]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(qIMinus[o][T2]) % Alignment == 0);

    assert(reinterpret_cast<uintptr_t>(faultStresses.normalStress[o]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(tractionResults.traction1[o]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(tractionResults.traction2[o]) % Alignment == 0);
  }
}

/**
 * Integrate over all Time points with the time weights and calculate the traction for each side
 * according to Carsten Uphoff Thesis: EQ.: 4.60
 *
 * @param[in] faultStresses
 * @param[in] tractionResults
 * @param[in] impAndEta
 * @param[in] impedancenceMatrices
 * @param[in] qInterpolatedPlus
 * @param[in] qInterpolatedMinus
 * @param[in] timeWeights
 * @param[out] imposedStatePlus
 * @param[out] imposedStateMinus
 */
template <RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void postcomputeImposedStateFromNewStress(
    const FaultStresses<RangeExecutor<Type>::Exec>& faultStresses,
    const TractionResults<RangeExecutor<Type>::Exec>& tractionResults,
    const ImpedancesAndEta& impAndEta,
    [[maybe_unused]] const ImpedanceMatrices& impedanceMatrices,
    real imposedStatePlus[tensor::QInterpolated::size()],
    real imposedStateMinus[tensor::QInterpolated::size()],
    const real qInterpolatedPlus[misc::TimeSteps][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[misc::TimeSteps][tensor::QInterpolated::size()],
    const real timeWeights[misc::TimeSteps],
    uint32_t startIndex = 0) {

#ifndef USE_POROELASTIC
  const auto invZs = impAndEta.invZs;
  const auto invZp = impAndEta.invZp;
  const auto invZsNeig = impAndEta.invZsNeig;
  const auto invZpNeig = impAndEta.invZpNeig;

  using ImposedStateShapeT = real(*)[misc::NumPaddedPoints];
  auto* imposedStateP = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
  auto* imposedStateM = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

  using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
  const auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  const auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  using namespace dr::misc::quantity_indices;

#ifndef ACL_DEVICE
  checkAlignmentPostCompute(
      qIPlus, qIMinus, imposedStateP, imposedStateM, faultStresses, tractionResults);
#endif

  using NumPointsRange = typename NumPoints<Type>::Range;

  // zero initialize
  real localImposedStateM[dr::misc::NumQuantities][NumPointsRange::Size]{};
  real localImposedStateP[dr::misc::NumQuantities][NumPointsRange::Size]{};

  for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
    const auto weight = timeWeights[o];

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
         index += NumPointsRange::Step) {
      auto i{startIndex + index};

      const auto normalStress =
          VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i);
      const auto traction1 =
          VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.traction1, o, i);
      const auto traction2 =
          VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.traction2, o, i);

      localImposedStateM[N][index] += weight * normalStress;
      localImposedStateM[T1][index] += weight * traction1;
      localImposedStateM[T2][index] += weight * traction2;
      localImposedStateM[U][index] +=
          weight * (qIMinus[o][U][i] - invZpNeig * (normalStress - qIMinus[o][N][i]));
      localImposedStateM[V][index] +=
          weight * (qIMinus[o][V][i] - invZsNeig * (traction1 - qIMinus[o][T1][i]));
      localImposedStateM[W][index] +=
          weight * (qIMinus[o][W][i] - invZsNeig * (traction2 - qIMinus[o][T2][i]));

      localImposedStateP[N][index] += weight * normalStress;
      localImposedStateP[T1][index] += weight * traction1;
      localImposedStateP[T2][index] += weight * traction2;
      localImposedStateP[U][index] +=
          weight * (qIPlus[o][U][i] + invZp * (normalStress - qIPlus[o][N][i]));
      localImposedStateP[V][index] +=
          weight * (qIPlus[o][V][i] + invZs * (traction1 - qIPlus[o][T1][i]));
      localImposedStateP[W][index] +=
          weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));
    }
  }

  for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
       index += NumPointsRange::Step) {
    auto i{startIndex + index};
#pragma unroll
    for (std::uint32_t q = 0; q < dr::misc::NumQuantities; ++q) {
      imposedStateM[q][i] = localImposedStateM[q][index];
      imposedStateP[q][i] = localImposedStateP[q][index];
    }
  }
#else
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

  alignas(Alignment) real thetaBuffer[tensor::theta::size()] = {};
  auto thetaView = init::theta::view::create(thetaBuffer);
  krnlM.theta = thetaBuffer;
  krnlP.theta = thetaBuffer;

  for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
    const auto weight = timeWeights[o];
    // copy values to yateto dataformat
    for (uint32_t i = 0; i < misc::NumPaddedPoints; ++i) {
      thetaView(i, 0) = faultStresses.normalStress[o][i];
      thetaView(i, 1) = tractionResults.traction1[o][i];
      thetaView(i, 2) = tractionResults.traction2[o][i];
      thetaView(i, 3) = faultStresses.fluidPressure[o][i];
    }
    // execute kernel (and hence update imposedStatePlus/Minus)
    krnlM.Qminus = qInterpolatedMinus[o];
    krnlM.weight = weight;
    krnlM.execute();

    krnlP.Qplus = qInterpolatedPlus[o];
    krnlP.weight = weight;
    krnlP.execute();
  }
#endif
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
template <RangeType Type = RangeType::CPU>
// See https://github.com/llvm/llvm-project/issues/60163
// NOLINTNEXTLINE
SEISSOL_HOSTDEVICE inline void
    adjustInitialStress(real initialStressInFaultCS[6][misc::NumPaddedPoints],
                        const real nucleationStressInFaultCS[6][misc::NumPaddedPoints],
                        // See https://github.com/llvm/llvm-project/issues/60163
                        // NOLINTNEXTLINE
                        real initialPressure[misc::NumPaddedPoints],
                        const real nucleationPressure[misc::NumPaddedPoints],
                        real fullUpdateTime,
                        real t0,
                        real s0,
                        real dt,
                        uint32_t startIndex = 0) {
  if (fullUpdateTime <= t0 + s0 && fullUpdateTime >= s0) {
    const real gNuc =
        gaussianNucleationFunction::smoothStepIncrement<real>(fullUpdateTime - s0, dt, t0);

    using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto pointIndex{startIndex + index};
      for (uint32_t i = 0; i < 6; i++) {
        initialStressInFaultCS[i][pointIndex] += nucleationStressInFaultCS[i][pointIndex] * gNuc;
      }
      initialPressure[pointIndex] += nucleationPressure[pointIndex] * gNuc;
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
SEISSOL_HOSTDEVICE inline void
    // See https://github.com/llvm/llvm-project/issues/60163
    // NOLINTNEXTLINE
    saveRuptureFrontOutput(bool ruptureTimePending[misc::NumPaddedPoints],
                           // See https://github.com/llvm/llvm-project/issues/60163
                           // NOLINTNEXTLINE
                           real ruptureTime[misc::NumPaddedPoints],
                           const real slipRateMagnitude[misc::NumPaddedPoints],
                           real fullUpdateTime,
                           uint32_t startIndex = 0) {

  using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto pointIndex{startIndex + index};
    constexpr real RuptureFrontThreshold = 0.001;
    if (ruptureTimePending[pointIndex] && slipRateMagnitude[pointIndex] > RuptureFrontThreshold) {
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
SEISSOL_HOSTDEVICE inline void
    savePeakSlipRateOutput(const real slipRateMagnitude[misc::NumPaddedPoints],
                           // See https://github.com/llvm/llvm-project/issues/60163
                           // NOLINTNEXTLINE
                           real peakSlipRate[misc::NumPaddedPoints],
                           uint32_t startIndex = 0) {

  using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto pointIndex{startIndex + index};
    peakSlipRate[pointIndex] = std::max(peakSlipRate[pointIndex], slipRateMagnitude[pointIndex]);
  }
}
/**
 * update timeSinceSlipRateBelowThreshold (used in Abort Criteria)
 *
 * param[in] slipRateMagnitude
 * param[in] ruptureTimePending
 * param[in, out] timeSinceSlipRateBelowThreshold
 * param[in] dt
 */
template <RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void
    updateTimeSinceSlipRateBelowThreshold(const real slipRateMagnitude[misc::NumPaddedPoints],
                                          const bool ruptureTimePending[misc::NumPaddedPoints],
                                          // See https://github.com/llvm/llvm-project/issues/60163
                                          // NOLINTNEXTLINE
                                          DREnergyOutput& energyData,
                                          const real dt,
                                          const real slipRateThreshold,
                                          uint32_t startIndex = 0) {

  using Range = typename NumPoints<Type>::Range;
  auto* timeSinceSlipRateBelowThreshold = energyData.timeSinceSlipRateBelowThreshold;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto pointIndex{startIndex + index};
    if (not ruptureTimePending[pointIndex]) {
      if (slipRateMagnitude[pointIndex] < slipRateThreshold) {
        timeSinceSlipRateBelowThreshold[pointIndex] += dt;
      } else {
        timeSinceSlipRateBelowThreshold[pointIndex] = 0;
      }
    } else {
      timeSinceSlipRateBelowThreshold[pointIndex] = std::numeric_limits<real>::infinity();
    }
  }
}
template <RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void computeFrictionEnergy(
    DREnergyOutput& energyData,
    const real qInterpolatedPlus[misc::TimeSteps][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[misc::TimeSteps][tensor::QInterpolated::size()],
    const ImpedancesAndEta& impAndEta,
    const real timeWeights[misc::TimeSteps],
    const real spaceWeights[seissol::kernels::NumSpaceQuadraturePoints],
    const DRGodunovData& godunovData,
    const real slipRateMagnitude[misc::NumPaddedPoints],
    const bool energiesFromAcrossFaultVelocities,
    size_t startIndex = 0) {

  auto* slip = reinterpret_cast<real(*)[misc::NumPaddedPoints]>(energyData.slip);
  auto* accumulatedSlip = energyData.accumulatedSlip;
  auto* frictionalEnergy = energyData.frictionalEnergy;
  const real doubledSurfaceAreaN = -static_cast<real>(godunovData.doubledSurfaceArea);

  using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
  const auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  const auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  const auto bPlus = impAndEta.etaS * impAndEta.invZs;
  const auto bMinus = impAndEta.etaS * impAndEta.invZsNeig;

  using Range = typename NumPoints<Type>::Range;
  real localAccumulatedSlip[Range::Size]{};
  real localFrictionalEnergy[Range::Size]{};
  real localSlip[3][Range::Size]{};

  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto i{startIndex + index};
    localAccumulatedSlip[index] = accumulatedSlip[i];
    localFrictionalEnergy[index] = frictionalEnergy[i];
#pragma unroll
    for (uint32_t d = 0; d < 3; ++d) {
      localSlip[d][index] = slip[d][i];
    }
  }

  using namespace dr::misc::quantity_indices;
  for (size_t o = 0; o < misc::TimeSteps; ++o) {
    const auto timeWeight = timeWeights[o];

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (size_t index = Range::Start; index < Range::End; index += Range::Step) {

      const size_t i{startIndex + index}; // startIndex is always 0 for CPU

      const real interpolatedSlipRate1 = qIMinus[o][U][i] - qIPlus[o][U][i];
      const real interpolatedSlipRate2 = qIMinus[o][V][i] - qIPlus[o][V][i];
      const real interpolatedSlipRate3 = qIMinus[o][W][i] - qIPlus[o][W][i];

      if (energiesFromAcrossFaultVelocities) {
        const real interpolatedSlipRateMagnitude =
            misc::magnitude(interpolatedSlipRate1, interpolatedSlipRate2, interpolatedSlipRate3);

        localAccumulatedSlip[index] += timeWeight * interpolatedSlipRateMagnitude;
      } else {
        // we use slipRateMagnitude (computed from slipRate1 and slipRate2 in the friction law)
        // instead of computing the slip rate magnitude from the differences in velocities
        // calculated above (magnitude of the vector (slipRateMagnitudei)). The moment magnitude
        // based on (slipRateMagnitudei) is typically non zero at the end of the earthquake
        // (probably because it incorporates the velocity discontinuities inherent of DG methods,
        // including the contributions of fault normal velocity discontinuity)
        localAccumulatedSlip[index] += timeWeight * slipRateMagnitude[i];
      }

      localSlip[0][index] += timeWeight * interpolatedSlipRate1;
      localSlip[1][index] += timeWeight * interpolatedSlipRate2;
      localSlip[2][index] += timeWeight * interpolatedSlipRate3;

      const real interpolatedTraction12 = bPlus * qIMinus[o][T1][i] + bMinus * qIPlus[o][T1][i];
      const real interpolatedTraction13 = bPlus * qIMinus[o][T2][i] + bMinus * qIPlus[o][T2][i];

      const auto spaceWeight = spaceWeights[i / multisim::NumSimulations];
      const auto weight = timeWeight * spaceWeight * doubledSurfaceAreaN;
      localFrictionalEnergy[index] += weight * (interpolatedTraction12 * interpolatedSlipRate2 +
                                                interpolatedTraction13 * interpolatedSlipRate3);
    }
  }

  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto i{startIndex + index};
    accumulatedSlip[i] = localAccumulatedSlip[index];
    frictionalEnergy[i] = localFrictionalEnergy[index];
#pragma unroll
    for (uint32_t d = 0; d < 3; ++d) {
      slip[d][i] = localSlip[d][index];
    }
  }
}

} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
