// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_

#include <Common/Executor.h>
#include <Equations/Datastructures.h>
#include <Model/CommonDatastructures.h>
#include <cmath>
#include <limits>
#include <type_traits>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "Initializer/Typedefs.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Solver/MultipleSimulations.h"

/**
 * Contains common functions required both for CPU and GPU impl.
 * of Dynamic Rupture solvers. The functions placed in
 * this class definition (of the header file) result
 * in the function inlining required for GPU impl.
 */
namespace seissol::dr::friction_law::common {

template <size_t StartT, size_t EndT, size_t StepT>
struct ForLoopRange {
  static constexpr size_t Start{StartT};
  static constexpr size_t End{EndT};
  static constexpr size_t Step{StepT};
};

enum class RangeType { CPU, GPU };

template <RangeType Type, typename Cfg>
struct NumPoints {
  private:
  using CpuRange = ForLoopRange<0, dr::misc::NumPaddedPoints<Cfg>, 1>;
  using GpuRange = ForLoopRange<0, 1, 1>;

  public:
  // Range::Start is 0, and Range::End is seissol::misc::NumPaddedPoints<Cfg> for CPU
  using Range = std::conditional_t<Type == RangeType::CPU, CpuRange, GpuRange>;
};

template <RangeType Type, typename Cfg>
struct QInterpolated {
  private:
  using CpuRange = ForLoopRange<0, tensor::QInterpolated<Cfg>::size(), 1>;
  using GpuRange = ForLoopRange<0, tensor::QInterpolated<Cfg>::size(), misc::NumPaddedPoints<Cfg>>;

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

template <typename Cfg, Executor Executor>
struct VariableIndexing;

template <typename Cfg>
struct VariableIndexing<Cfg, Executor::Host> {
  template <typename T>
  static constexpr T&
      index(T (&data)[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>], int o, int i) {
    return data[o][i];
  }

  template <typename T>
  static constexpr T
      index(const T (&data)[misc::TimeSteps<Cfg>][misc::NumPaddedPoints<Cfg>], int o, int i) {
    return data[o][i];
  }
};

template <typename Cfg>
struct VariableIndexing<Cfg, Executor::Device> {
  template <typename T>
  static constexpr T& index(T (&data)[misc::TimeSteps<Cfg>], int o, int i) {
    return data[o];
  }

  template <typename T>
  static constexpr T index(const T (&data)[misc::TimeSteps<Cfg>], int o, int i) {
    return data[o];
  }
};

/**
 * Asserts whether all relevant arrays are properly aligned
 */
template <typename Cfg>
inline void checkAlignmentPreCompute(
    const Real<Cfg> qIPlus[misc::TimeSteps<Cfg>][dr::misc::NumQuantities<Cfg>]
                          [dr::misc::NumPaddedPoints<Cfg>],
    const Real<Cfg> qIMinus[misc::TimeSteps<Cfg>][dr::misc::NumQuantities<Cfg>]
                           [dr::misc::NumPaddedPoints<Cfg>],
    const FaultStresses<Cfg, Executor::Host>& faultStresses) {
  using namespace dr::misc::quantity_indices;
  for (unsigned o = 0; o < misc::TimeSteps<Cfg>; ++o) {
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
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void precomputeStressFromQInterpolated(
    FaultStresses<Cfg, RangeExecutor<Type>::Exec>& faultStresses,
    const ImpedancesAndEta<Real<Cfg>>& impAndEta,
    const ImpedanceMatrices<Cfg>& impedanceMatrices,
    const Real<Cfg> qInterpolatedPlus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const Real<Cfg> qInterpolatedMinus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    Real<Cfg> etaPDamp,
    unsigned startLoopIndex = 0) {
  static_assert(tensor::QInterpolated<Cfg>::Shape[seissol::multisim::BasisDim<Cfg>] ==
                    tensor::resample<Cfg>::Shape[0],
                "Different number of quadrature points?");

  if constexpr (model::MaterialTT<Cfg>::Type != model::MaterialType::Poroelastic) {
    const auto etaP = impAndEta.etaP * etaPDamp;
    const auto etaS = impAndEta.etaS;
    const auto invZp = impAndEta.invZp;
    const auto invZs = impAndEta.invZs;
    const auto invZpNeig = impAndEta.invZpNeig;
    const auto invZsNeig = impAndEta.invZsNeig;

    using QInterpolatedShapeT =
        const Real<Cfg>(*)[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>];
    const auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
    const auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));

    using namespace dr::misc::quantity_indices;

#ifndef ACL_DEVICE
    checkAlignmentPreCompute(qIPlus, qIMinus, faultStresses);
#endif

    for (unsigned o = 0; o < misc::TimeSteps<Cfg>; ++o) {
      using Range = typename NumPoints<Type, Cfg>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
      for (auto index = Range::Start; index < Range::End; index += Range::Step) {
        auto i{startLoopIndex + index};
        VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i) =
            etaP * (qIMinus[o][U][i] - qIPlus[o][U][i] + qIPlus[o][N][i] * invZp +
                    qIMinus[o][N][i] * invZpNeig);

        VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.traction1, o, i) =
            etaS * (qIMinus[o][V][i] - qIPlus[o][V][i] + qIPlus[o][T1][i] * invZs +
                    qIMinus[o][T1][i] * invZsNeig);

        VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.traction2, o, i) =
            etaS * (qIMinus[o][W][i] - qIPlus[o][W][i] + qIPlus[o][T2][i] * invZs +
                    qIMinus[o][T2][i] * invZsNeig);
      }
    }
  } else {
    seissol::dynamicRupture::kernel::computeTheta<Cfg> krnl;
    krnl.extractVelocities = init::extractVelocities<Cfg>::Values;
    krnl.extractTractions = init::extractTractions<Cfg>::Values;

    // Compute Theta from eq (4.53) in Carsten's thesis
    krnl.Zplus = impedanceMatrices.impedance;
    krnl.Zminus = impedanceMatrices.impedanceNeig;
    krnl.eta = impedanceMatrices.eta;

    alignas(Alignment) Real<Cfg> thetaBuffer[tensor::theta<Cfg>::size()]{};
    krnl.theta = thetaBuffer;
    auto thetaView = init::theta<Cfg>::view::create(thetaBuffer);

    for (unsigned o = 0; o < Cfg::ConvergenceOrder; ++o) {
      krnl.Qplus = qInterpolatedPlus[o];
      krnl.Qminus = qInterpolatedMinus[o];
      krnl.execute();

      for (unsigned i = 0; i < misc::NumPaddedPoints<Cfg>; ++i) {
        faultStresses.normalStress[o][i] = thetaView(i, 0);
        faultStresses.traction1[o][i] = thetaView(i, 1);
        faultStresses.traction2[o][i] = thetaView(i, 2);
        faultStresses.fluidPressure[o][i] = thetaView(i, 3);
      }
    }
  }
}

/**
 * Asserts whether all relevant arrays are properly aligned
 */
template <typename Cfg>
inline void checkAlignmentPostCompute(
    const Real<Cfg> qIPlus[misc::TimeSteps<Cfg>][dr::misc::NumQuantities<Cfg>]
                          [dr::misc::NumPaddedPoints<Cfg>],
    const Real<Cfg> qIMinus[misc::TimeSteps<Cfg>][dr::misc::NumQuantities<Cfg>]
                           [dr::misc::NumPaddedPoints<Cfg>],
    const Real<Cfg> imposedStateP[misc::TimeSteps<Cfg>][dr::misc::NumPaddedPoints<Cfg>],
    const Real<Cfg> imposedStateM[misc::TimeSteps<Cfg>][dr::misc::NumPaddedPoints<Cfg>],
    const FaultStresses<Cfg, Executor::Host>& faultStresses,
    const TractionResults<Cfg, Executor::Host>& tractionResults) {
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

  for (size_t o = 0; o < misc::TimeSteps<Cfg>; ++o) {
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
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void postcomputeImposedStateFromNewStress(
    const FaultStresses<Cfg, RangeExecutor<Type>::Exec>& faultStresses,
    const TractionResults<Cfg, RangeExecutor<Type>::Exec>& tractionResults,
    const ImpedancesAndEta<Real<Cfg>>& impAndEta,
    const ImpedanceMatrices<Cfg>& impedanceMatrices,
    Real<Cfg> imposedStatePlus[tensor::QInterpolated<Cfg>::size()],
    Real<Cfg> imposedStateMinus[tensor::QInterpolated<Cfg>::size()],
    const Real<Cfg> qInterpolatedPlus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const Real<Cfg> qInterpolatedMinus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const double* timeWeights,
    unsigned startIndex = 0) {

  // set imposed state to zero
  using QInterpolatedRange = typename QInterpolated<Type, Cfg>::Range;
  for (auto index = QInterpolatedRange::Start; index < QInterpolatedRange::End;
       index += QInterpolatedRange::Step) {
    auto i{startIndex + index};
    imposedStatePlus[i] = static_cast<Real<Cfg>>(0.0);
    imposedStateMinus[i] = static_cast<Real<Cfg>>(0.0);
  }
  if constexpr (model::MaterialTT<Cfg>::Type != model::MaterialType::Poroelastic) {
    const auto invZs = impAndEta.invZs;
    const auto invZp = impAndEta.invZp;
    const auto invZsNeig = impAndEta.invZsNeig;
    const auto invZpNeig = impAndEta.invZpNeig;

    using ImposedStateShapeT = Real<Cfg>(*)[misc::NumPaddedPoints<Cfg>];
    auto* imposedStateP = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
    auto* imposedStateM = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

    using QInterpolatedShapeT =
        const Real<Cfg>(*)[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>];
    const auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
    const auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

    using namespace dr::misc::quantity_indices;

#ifndef ACL_DEVICE
    checkAlignmentPostCompute(
        qIPlus, qIMinus, imposedStateP, imposedStateM, faultStresses, tractionResults);
#endif

    for (unsigned o = 0; o < misc::TimeSteps<Cfg>; ++o) {
      auto weight = timeWeights[o];

      using NumPointsRange = typename NumPoints<Type, Cfg>::Range;
#ifndef ACL_DEVICE
#pragma omp simd
#endif
      for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
           index += NumPointsRange::Step) {
        auto i{startIndex + index};

        const auto normalStress = VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(
            faultStresses.normalStress, o, i);
        const auto traction1 = VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(
            tractionResults.traction1, o, i);
        const auto traction2 = VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(
            tractionResults.traction2, o, i);

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
        imposedStateP[U][i] +=
            weight * (qIPlus[o][U][i] + invZp * (normalStress - qIPlus[o][N][i]));
        imposedStateP[V][i] += weight * (qIPlus[o][V][i] + invZs * (traction1 - qIPlus[o][T1][i]));
        imposedStateP[W][i] += weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));
      }
    }
  } else {
    // setup kernel
    seissol::dynamicRupture::kernel::computeImposedStateM<Cfg> krnlM;
    krnlM.extractVelocities = init::extractVelocities<Cfg>::Values;
    krnlM.extractTractions = init::extractTractions<Cfg>::Values;
    krnlM.mapToVelocities = init::mapToVelocities<Cfg>::Values;
    krnlM.mapToTractions = init::mapToTractions<Cfg>::Values;
    krnlM.Zminus = impedanceMatrices.impedanceNeig;
    krnlM.imposedState = imposedStateMinus;

    seissol::dynamicRupture::kernel::computeImposedStateP<Cfg> krnlP;
    krnlP.extractVelocities = init::extractVelocities<Cfg>::Values;
    krnlP.extractTractions = init::extractTractions<Cfg>::Values;
    krnlP.mapToVelocities = init::mapToVelocities<Cfg>::Values;
    krnlP.mapToTractions = init::mapToTractions<Cfg>::Values;
    krnlP.Zplus = impedanceMatrices.impedance;
    krnlP.imposedState = imposedStatePlus;

    alignas(Alignment) Real<Cfg> thetaBuffer[tensor::theta<Cfg>::size()] = {};
    auto thetaView = init::theta<Cfg>::view::create(thetaBuffer);
    krnlM.theta = thetaBuffer;
    krnlP.theta = thetaBuffer;

    for (unsigned o = 0; o < Cfg::ConvergenceOrder; ++o) {
      auto weight = timeWeights[o];
      // copy values to yateto dataformat
      for (unsigned i = 0; i < misc::NumPaddedPoints<Cfg>; ++i) {
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
template <typename Cfg, RangeType Type = RangeType::CPU>
// See https://github.com/llvm/llvm-project/issues/60163
// NOLINTNEXTLINE
SEISSOL_HOSTDEVICE inline void
    adjustInitialStress(Real<Cfg> initialStressInFaultCS[6][misc::NumPaddedPoints<Cfg>],
                        const Real<Cfg> nucleationStressInFaultCS[6][misc::NumPaddedPoints<Cfg>],
                        // See https://github.com/llvm/llvm-project/issues/60163
                        // NOLINTNEXTLINE
                        Real<Cfg> initialPressure[misc::NumPaddedPoints<Cfg>],
                        const Real<Cfg> nucleationPressure[misc::NumPaddedPoints<Cfg>],
                        Real<Cfg> fullUpdateTime,
                        Real<Cfg> t0,
                        Real<Cfg> s0,
                        Real<Cfg> dt,
                        unsigned startIndex = 0) {
  if (fullUpdateTime <= t0 + s0 && fullUpdateTime >= s0) {
    const auto gNuc =
        gaussianNucleationFunction::smoothStepIncrement<Real<Cfg>>(fullUpdateTime - s0, dt, t0);

    using Range = typename NumPoints<Type, Cfg>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto pointIndex{startIndex + index};
      for (unsigned i = 0; i < 6; i++) {
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
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void
    // See https://github.com/llvm/llvm-project/issues/60163
    // NOLINTNEXTLINE
    saveRuptureFrontOutput(bool ruptureTimePending[misc::NumPaddedPoints<Cfg>],
                           // See https://github.com/llvm/llvm-project/issues/60163
                           // NOLINTNEXTLINE
                           Real<Cfg> ruptureTime[misc::NumPaddedPoints<Cfg>],
                           const Real<Cfg> slipRateMagnitude[misc::NumPaddedPoints<Cfg>],
                           Real<Cfg> fullUpdateTime,
                           unsigned startIndex = 0) {

  using Range = typename NumPoints<Type, Cfg>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto pointIndex{startIndex + index};
    constexpr Real<Cfg> RuptureFrontThreshold = 0.001;
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
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void
    savePeakSlipRateOutput(const Real<Cfg> slipRateMagnitude[misc::NumPaddedPoints<Cfg>],
                           // See https://github.com/llvm/llvm-project/issues/60163
                           // NOLINTNEXTLINE
                           Real<Cfg> peakSlipRate[misc::NumPaddedPoints<Cfg>],
                           unsigned startIndex = 0) {

  using Range = typename NumPoints<Type, Cfg>::Range;

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
 * param[in] sumDt
 */
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void updateTimeSinceSlipRateBelowThreshold(
    const Real<Cfg> slipRateMagnitude[misc::NumPaddedPoints<Cfg>],
    const bool ruptureTimePending[misc::NumPaddedPoints<Cfg>],
    // See https://github.com/llvm/llvm-project/issues/60163
    // NOLINTNEXTLINE
    DREnergyOutput<Cfg>& energyData,
    const Real<Cfg> sumDt,
    const Real<Cfg> slipRateThreshold,
    unsigned startIndex = 0) {

  using Range = typename NumPoints<Type, Cfg>::Range;
  auto* timeSinceSlipRateBelowThreshold = energyData.timeSinceSlipRateBelowThreshold;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto pointIndex{startIndex + index};
    if (not ruptureTimePending[pointIndex]) {
      if (slipRateMagnitude[pointIndex] < slipRateThreshold) {
        timeSinceSlipRateBelowThreshold[pointIndex] += sumDt;
      } else {
        timeSinceSlipRateBelowThreshold[pointIndex] = 0;
      }
    } else {
      timeSinceSlipRateBelowThreshold[pointIndex] = std::numeric_limits<Real<Cfg>>::infinity();
    }
  }
}
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void computeFrictionEnergy(
    DREnergyOutput<Cfg>& energyData,
    const Real<Cfg> qInterpolatedPlus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const Real<Cfg> qInterpolatedMinus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const ImpedancesAndEta<Real<Cfg>>& impAndEta,
    const double* timeWeights,
    const Real<Cfg> spaceWeights[seissol::kernels::NumSpaceQuadraturePoints<Cfg>],
    const DRGodunovData<Cfg>& godunovData,
    const Real<Cfg> slipRateMagnitude[misc::NumPaddedPoints<Cfg>],
    const bool energiesFromAcrossFaultVelocities,
    size_t startIndex = 0) {

  auto* slip = reinterpret_cast<Real<Cfg>(*)[misc::NumPaddedPoints<Cfg>]>(energyData.slip);
  auto* accumulatedSlip = energyData.accumulatedSlip;
  auto* frictionalEnergy = energyData.frictionalEnergy;
  const double doubledSurfaceArea = godunovData.doubledSurfaceArea;

  using QInterpolatedShapeT =
      const Real<Cfg>(*)[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>];
  const auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  const auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  const auto bPlus = impAndEta.etaS * impAndEta.invZs;
  const auto bMinus = impAndEta.etaS * impAndEta.invZsNeig;

  using Range = typename NumPoints<Type, Cfg>::Range;

  using namespace dr::misc::quantity_indices;
  for (size_t o = 0; o < misc::TimeSteps<Cfg>; ++o) {
    const auto timeWeight = timeWeights[o];
#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (size_t index = Range::Start; index < Range::End; index += Range::Step) {

      const size_t i{startIndex + index}; // startIndex is always 0 for CPU

      const auto interpolatedSlipRate1 = qIMinus[o][U][i] - qIPlus[o][U][i];
      const auto interpolatedSlipRate2 = qIMinus[o][V][i] - qIPlus[o][V][i];
      const auto interpolatedSlipRate3 = qIMinus[o][W][i] - qIPlus[o][W][i];

      if (energiesFromAcrossFaultVelocities) {
        const auto interpolatedSlipRateMagnitude =
            misc::magnitude(interpolatedSlipRate1, interpolatedSlipRate2, interpolatedSlipRate3);

        accumulatedSlip[i] += timeWeight * interpolatedSlipRateMagnitude;
      } else {
        // we use slipRateMagnitude (computed from slipRate1 and slipRate2 in the friction law)
        // instead of computing the slip rate magnitude from the differences in velocities
        // calculated above (magnitude of the vector (slipRateMagnitudei)). The moment magnitude
        // based on (slipRateMagnitudei) is typically non zero at the end of the earthquake
        // (probably because it incorporates the velocity discontinuities inherent of DG methods,
        // including the contributions of fault normal velocity discontinuity)
        accumulatedSlip[i] += timeWeight * slipRateMagnitude[i];
      }

      slip[0][i] += timeWeight * interpolatedSlipRate1;
      slip[1][i] += timeWeight * interpolatedSlipRate2;
      slip[2][i] += timeWeight * interpolatedSlipRate3;

      const auto interpolatedTraction12 = bPlus * qIMinus[o][T1][i] + bMinus * qIPlus[o][T1][i];
      const auto interpolatedTraction13 = bPlus * qIMinus[o][T2][i] + bMinus * qIPlus[o][T2][i];

      const auto spaceWeight = spaceWeights[i / multisim::NumSimulations<Cfg>];

      const auto weight = -timeWeight * spaceWeight * doubledSurfaceArea;
      frictionalEnergy[i] += weight * (interpolatedTraction12 * interpolatedSlipRate2 +
                                       interpolatedTraction13 * interpolatedSlipRate3);
    }
  }
}

} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
