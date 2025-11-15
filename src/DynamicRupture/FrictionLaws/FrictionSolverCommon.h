// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_

#include <Common/Executor.h>
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
  static constexpr real& index(real (&data)[misc::TimeSteps], int o, int i) { return data[o]; }

  static constexpr real index(const real (&data)[misc::TimeSteps], int o, int i) { return data[o]; }
};

/**
 * Asserts whether all relevant arrays are properly aligned
 */
inline void checkAlignmentPreCompute(
    const real qIPlus[misc::TimeSteps][dr::misc::NumQuantities][dr::misc::NumPaddedPoints],
    const real qIMinus[misc::TimeSteps][dr::misc::NumQuantities][dr::misc::NumPaddedPoints],
    const FaultStresses<Executor::Host>& faultStresses) {
  using namespace dr::misc::quantity_indices;
  for (unsigned o = 0; o < misc::TimeSteps; ++o) {
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
    const ImpedanceMatrices& impedanceMatrices,
    const real qInterpolatedPlus[misc::TimeSteps][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[misc::TimeSteps][tensor::QInterpolated::size()],
    real etaPDamp,
    unsigned startLoopIndex = 0) {
  static_assert(tensor::QInterpolated::Shape[seissol::multisim::BasisFunctionDimension] ==
                    tensor::resample::Shape[0],
                "Different number of quadrature points?");

  if constexpr (model::MaterialT::Type != model::MaterialType::Poroelastic) {
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

    for (unsigned o = 0; o < misc::TimeSteps; ++o) {
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
  } else {
    // poroelastic kernel (for CPU+GPU)
    // TODO: generalize and unify with the above (probably using either templates or Yateto)
    // (the v1.1.0-1.3.1 Yateto+selector matrix based kernel was removed since GPU support was
    // missing)

    using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
    const auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
    const auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));

    using namespace dr::misc::quantity_indices;

    for (unsigned o = 0; o < misc::TimeSteps; ++o) {
      using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
      for (auto index = Range::Start; index < Range::End; index += Range::Step) {
        auto i{startLoopIndex + index};

        // Compute Theta from eq (4.53) in Carsten's thesis

        real velDiff[4]{};
        velDiff[0] = qIMinus[o][U][i] - qIPlus[o][U][i];
        velDiff[1] = qIMinus[o][V][i] - qIPlus[o][V][i];
        velDiff[2] = qIMinus[o][W][i] - qIPlus[o][W][i];
        velDiff[3] = qIMinus[o][FU][i] - qIPlus[o][FU][i];

        real strP[4]{};
        real strM[4]{};
        const auto rowCompute = [&](auto linear, auto index) {
#pragma unroll
          for (std::uint32_t j = 0; j < 4; ++j) {
            strP[j] += impedanceMatrices.impedance[linear * 4 + j] * qIPlus[o][index][i];
            strM[j] += impedanceMatrices.impedanceNeig[linear * 4 + j] * qIMinus[o][index][i];
          }
        };
        rowCompute(0, N);
        rowCompute(1, T1);
        rowCompute(2, T2);
        rowCompute(3, FP);

        real res[4]{};
#pragma unroll
        for (std::uint32_t k = 0; k < 4; ++k) {
#pragma unroll
          for (std::uint32_t j = 0; j < 4; ++j) {
            res[j] += impedanceMatrices.eta[k * 4 + j] * (velDiff[k] + strP[k] + strM[k]);
          }
        }

        VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i) =
            res[0];
        VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction1, o, i) = res[1];
        VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction2, o, i) = res[2];
        VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.fluidPressure, o, i) =
            res[3];
      }
    }
  }
}

/**
 * Asserts whether all relevant arrays are properly aligned
 */
inline void checkAlignmentPostCompute(
    const real qIPlus[misc::TimeSteps][dr::misc::NumQuantities][dr::misc::NumPaddedPoints],
    const real qIMinus[misc::TimeSteps][dr::misc::NumQuantities][dr::misc::NumPaddedPoints],
    const real imposedStateP[misc::TimeSteps][dr::misc::NumPaddedPoints],
    const real imposedStateM[misc::TimeSteps][dr::misc::NumPaddedPoints],
    const FaultStresses<Executor::Host>& faultStresses,
    const TractionResults<Executor::Host>& tractionResults) {
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
    const ImpedanceMatrices& impedanceMatrices,
    real imposedStatePlus[tensor::QInterpolated::size()],
    real imposedStateMinus[tensor::QInterpolated::size()],
    const real qInterpolatedPlus[misc::TimeSteps][tensor::QInterpolated::size()],
    const real qInterpolatedMinus[misc::TimeSteps][tensor::QInterpolated::size()],
    const double timeWeights[misc::TimeSteps],
    unsigned startIndex = 0) {

  // set imposed state to zero
  using QInterpolatedRange = typename QInterpolated<Type>::Range;
  for (auto index = QInterpolatedRange::Start; index < QInterpolatedRange::End;
       index += QInterpolatedRange::Step) {
    auto i{startIndex + index};
    imposedStatePlus[i] = static_cast<real>(0.0);
    imposedStateMinus[i] = static_cast<real>(0.0);
  }
  if constexpr (model::MaterialT::Type != model::MaterialType::Poroelastic) {
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

    for (std::uint32_t o = 0; o < misc::TimeSteps; ++o) {
      auto weight = timeWeights[o];

      using NumPointsRange = typename NumPoints<Type>::Range;
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
    // poroelastic kernel (for CPU+GPU)
    // TODO: generalize and unify with the above (probably using either templates or Yateto)
    // (the v1.1.0-1.3.1 Yateto+selector matrix based kernel was removed since GPU support was
    // missing)

    using ImposedStateShapeT = real(*)[misc::NumPaddedPoints];
    auto* imposedStateP = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
    auto* imposedStateM = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

    using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
    const auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
    const auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

    using namespace dr::misc::quantity_indices;

    for (std::uint32_t o = 0; o < misc::TimeSteps; ++o) {
      auto weight = timeWeights[o];

      using NumPointsRange = typename NumPoints<Type>::Range;
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
        const auto fluidPressure =
            VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.fluidPressure, o, i);

        const auto handleSide = [&](auto& imposedState, const auto& qI, const auto& mZ, real sign) {
          constexpr std::uint32_t Count = 4;

          imposedState[N][i] += weight * normalStress;
          imposedState[T1][i] += weight * traction1;
          imposedState[T2][i] += weight * traction2;
          imposedState[FP][i] += weight * fluidPressure;

          real diff[Count]{};
          diff[0] = (normalStress - qI[o][N][i]) * sign;
          diff[1] = (traction1 - qI[o][T1][i]) * sign;
          diff[2] = (traction2 - qI[o][T2][i]) * sign;
          diff[3] = (fluidPressure - qI[o][FP][i]) * sign;

          const auto handleEntry = [&](auto linear, auto index) {
            real acc = 0;
#pragma unroll
            for (int k = 0; k < Count; ++k) {
              acc += mZ[Count * k + linear] * diff[k];
            }
            imposedState[index][i] += weight * (qI[o][index][i] + acc);
          };

          handleEntry(0, U);
          handleEntry(1, V);
          handleEntry(2, W);
          handleEntry(3, FU);
        };

        handleSide(imposedStateM, qIMinus, impedanceMatrices.impedanceNeig, -1);
        handleSide(imposedStateP, qIPlus, impedanceMatrices.impedance, 1);
      }
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
                        unsigned startIndex = 0) {
  if (fullUpdateTime <= t0 + s0 && fullUpdateTime >= s0) {
    const real gNuc =
        gaussianNucleationFunction::smoothStepIncrement<real>(fullUpdateTime - s0, dt, t0);

    using Range = typename NumPoints<Type>::Range;

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
                           unsigned startIndex = 0) {

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
                           unsigned startIndex = 0) {

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
 * param[in] sumDt
 */
template <RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void
    updateTimeSinceSlipRateBelowThreshold(const real slipRateMagnitude[misc::NumPaddedPoints],
                                          const bool ruptureTimePending[misc::NumPaddedPoints],
                                          // See https://github.com/llvm/llvm-project/issues/60163
                                          // NOLINTNEXTLINE
                                          DREnergyOutput& energyData,
                                          const real sumDt,
                                          const real slipRateThreshold,
                                          unsigned startIndex = 0) {

  using Range = typename NumPoints<Type>::Range;
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
    const double timeWeights[misc::TimeSteps],
    const real spaceWeights[seissol::kernels::NumSpaceQuadraturePoints],
    const DRGodunovData& godunovData,
    const real slipRateMagnitude[misc::NumPaddedPoints],
    const bool energiesFromAcrossFaultVelocities,
    size_t startIndex = 0) {

  auto* slip = reinterpret_cast<real(*)[misc::NumPaddedPoints]>(energyData.slip);
  auto* accumulatedSlip = energyData.accumulatedSlip;
  auto* frictionalEnergy = energyData.frictionalEnergy;
  const double doubledSurfaceArea = godunovData.doubledSurfaceArea;

  using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
  const auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  const auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  const auto bPlus = impAndEta.etaS * impAndEta.invZs;
  const auto bMinus = impAndEta.etaS * impAndEta.invZsNeig;

  using Range = typename NumPoints<Type>::Range;

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

      const real interpolatedTraction12 = bPlus * qIMinus[o][T1][i] + bMinus * qIPlus[o][T1][i];
      const real interpolatedTraction13 = bPlus * qIMinus[o][T2][i] + bMinus * qIPlus[o][T2][i];

      const auto spaceWeight = spaceWeights[i / multisim::NumSimulations];

      const auto weight = -timeWeight * spaceWeight * doubledSurfaceArea;
      frictionalEnergy[i] += weight * (interpolatedTraction12 * interpolatedSlipRate2 +
                                       interpolatedTraction13 * interpolatedSlipRate3);
    }
  }
}

} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
