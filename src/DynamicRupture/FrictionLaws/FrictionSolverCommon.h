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

  if constexpr (model::MaterialT::Type == model::MaterialType::Elastic ||
                model::MaterialT::Type == model::MaterialType::Viscoelastic) {
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
  } else {
    // poroelastic kernel (for CPU+GPU)
    // TODO: generalize and unify with the above (probably using either templates or Yateto)
    // (the v1.1.0-1.3.1 Yateto+selector matrix based kernel was removed since GPU support was
    // missing)
    //
    // etaPDamp scales the fault-normal response only. The elastic branch does that by damping
    // etaP; here eta is not diagonal, but since res[0] = sum_k eta(0,k) * x_k, damping row 0 of
    // eta and damping res[0] are the same thing -- so it is applied at the very end. The
    // poroelastic fluid pressure res[3] is deliberately left undamped, matching the elastic
    // behaviour of only touching the normal traction.

    using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
    const auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
    const auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));

    using namespace dr::misc::quantity_indices;

    for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
      using Range = typename NumPoints<Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
      for (auto index = Range::Start; index < Range::End; index += Range::Step) {
        auto i{startLoopIndex + index};

        constexpr uint32_t Count =
            model::MaterialT::Type == model::MaterialType::Poroelastic ? 4 : 3;

        // Compute Theta from eq (4.53) in Carsten's thesis

        real velDiff[Count]{};
        velDiff[0] = qIMinus[o][U][i] - qIPlus[o][U][i];
        velDiff[1] = qIMinus[o][V][i] - qIPlus[o][V][i];
        velDiff[2] = qIMinus[o][W][i] - qIPlus[o][W][i];

        real strP[Count]{};
        real strM[Count]{};
        const auto rowCompute = [&](auto linear, auto index) {
#pragma unroll
          for (std::uint32_t j = 0; j < Count; ++j) {
            strP[j] += impedanceMatrices.impedance[linear * Count + j] * qIPlus[o][index][i];
            strM[j] += impedanceMatrices.impedanceNeig[linear * Count + j] * qIMinus[o][index][i];
          }
        };
        rowCompute(0, N);
        rowCompute(1, T1);
        rowCompute(2, T2);

        if constexpr (model::MaterialT::Type == model::MaterialType::Poroelastic) {
          velDiff[3] = qIMinus[o][FU][i] - qIPlus[o][FU][i];
          rowCompute(3, FP);
        }

        real res[Count]{};
#pragma unroll
        for (std::uint32_t k = 0; k < Count; ++k) {
#pragma unroll
          for (std::uint32_t j = 0; j < Count; ++j) {
            res[j] += impedanceMatrices.eta[k * Count + j] * (velDiff[k] + strP[k] + strM[k]);
          }
        }

        VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i) =
            res[0] * etaPDamp;
        VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction1, o, i) = res[1];
        VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.traction2, o, i) = res[2];
        if constexpr (model::MaterialT::Type == model::MaterialType::Poroelastic) {
          VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.fluidPressure, o, i) =
              res[3];
        }
      }
    }
  }
}

/**
 * Seeds the post-solve tractions with the trial ("stick") state.
 *
 * At the moment this only concerns the fault-normal stress: for every impedance that is block
 * diagonal in the fault-normal direction versus the two tangential ones -- i.e. everything except
 * anisotropy -- a frictional interface with no opening leaves the normal stress at its trial value,
 * so the friction laws never touch it. Seeding it here keeps those friction laws free of
 * boilerplate and makes postcomputeImposedStateFromNewStress read a single, well-defined source.
 *
 * @param[in] faultStresses trial stresses from precomputeStressFromQInterpolated
 * @param[out] tractionResults
 */
template <RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void
    initializeTractionResults(const FaultStresses<RangeExecutor<Type>::Exec>& faultStresses,
                              TractionResults<RangeExecutor<Type>::Exec>& tractionResults,
                              uint32_t startIndex = 0) {
  using Range = typename NumPoints<Type>::Range;

  for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      const auto i{startIndex + index};
      VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.normalStress, o, i) =
          VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, o, i);
    }
  }
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
    assert(reinterpret_cast<uintptr_t>(tractionResults.normalStress[o]) % Alignment == 0);
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

  using NumPointsRange = typename NumPoints<Type>::Range;

  // zero initialize
  real localImposedStateM[dr::misc::NumQuantities][NumPointsRange::Size]{};
  real localImposedStateP[dr::misc::NumQuantities][NumPointsRange::Size]{};

  using ImposedStateShapeT = real(*)[misc::NumPaddedPoints];
  auto* imposedStateP = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
  auto* imposedStateM = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

  using QInterpolatedShapeT = const real(*)[misc::NumQuantities][misc::NumPaddedPoints];
  const auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  const auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  // set imposed state to zero
  if constexpr (model::MaterialT::Type == model::MaterialType::Elastic ||
                model::MaterialT::Type == model::MaterialType::Viscoelastic) {
    const auto invZs = impAndEta.invZs;
    const auto invZp = impAndEta.invZp;
    const auto invZsNeig = impAndEta.invZsNeig;
    const auto invZpNeig = impAndEta.invZpNeig;

    using namespace dr::misc::quantity_indices;

#ifndef ACL_DEVICE
    checkAlignmentPostCompute(
        qIPlus, qIMinus, imposedStateP, imposedStateM, faultStresses, tractionResults);
#endif

    for (uint32_t o = 0; o < misc::TimeSteps; ++o) {
      const auto weight = timeWeights[o];

#ifndef ACL_DEVICE
#pragma omp simd
#endif
      for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
           index += NumPointsRange::Step) {
        auto i{startIndex + index};

        const auto normalStress =
            VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.normalStress, o, i);
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
  } else {
    // poroelastic kernel (for CPU+GPU)
    // TODO: generalize and unify with the above (probably using either templates or Yateto)
    // (the v1.1.0-1.3.1 Yateto+selector matrix based kernel was removed since GPU support was
    // missing)

    using namespace dr::misc::quantity_indices;

    for (std::uint32_t o = 0; o < misc::TimeSteps; ++o) {
      const auto weight = timeWeights[o];

#ifndef ACL_DEVICE
#pragma omp simd
#endif
      for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
           index += NumPointsRange::Step) {
        auto i{startIndex + index};

        const auto normalStress =
            VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.normalStress, o, i);
        const auto traction1 =
            VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.traction1, o, i);
        const auto traction2 =
            VariableIndexing<RangeExecutor<Type>::Exec>::index(tractionResults.traction2, o, i);
        const auto fluidPressure =
            VariableIndexing<RangeExecutor<Type>::Exec>::index(faultStresses.fluidPressure, o, i);

        const auto handleSide = [&](auto& imposedState, const auto& qI, const auto& mZ, real sign) {
          constexpr std::uint32_t Count =
              model::MaterialT::Type == model::MaterialType::Poroelastic ? 4 : 3;

          imposedState[N][i] += weight * normalStress;
          imposedState[T1][i] += weight * traction1;
          imposedState[T2][i] += weight * traction2;

          real diff[Count]{};
          diff[0] = (normalStress - qI[o][N][i]) * sign;
          diff[1] = (traction1 - qI[o][T1][i]) * sign;
          diff[2] = (traction2 - qI[o][T2][i]) * sign;

          if constexpr (model::MaterialT::Type == model::MaterialType::Poroelastic) {
            imposedState[FP][i] += weight * fluidPressure;
            diff[3] = (fluidPressure - qI[o][FP][i]) * sign;
          }

          const auto handleEntry = [&](auto linear, auto index) {
            real acc = 0;
#pragma unroll
            for (std::uint32_t k = 0; k < Count; ++k) {
              acc += mZ[Count * k + linear] * diff[k];
            }
            imposedState[index][i] += weight * (qI[o][index][i] + acc);
          };

          handleEntry(0, U);
          handleEntry(1, V);
          handleEntry(2, W);
          if constexpr (model::MaterialT::Type == model::MaterialType::Poroelastic) {
            handleEntry(3, FU);
          }
        };

        handleSide(localImposedStateM, qIMinus, impedanceMatrices.impedanceNeig, -1);
        handleSide(localImposedStateP, qIPlus, impedanceMatrices.impedance, 1);
      }
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

  using namespace dr::misc::quantity_indices;

  real bPlus11{};
  real bPlus12{};
  real bPlus21{};
  real bPlus22{};
  real bMinus11{};
  real bMinus12{};
  real bMinus21{};
  real bMinus22{};

  if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
    constexpr auto Rows = 3;
    bPlus11 = godunovData.tractionPlusMatrix[Rows * 1 + 1];
    bPlus12 = godunovData.tractionPlusMatrix[Rows * 1 + 2];
    bPlus21 = godunovData.tractionPlusMatrix[Rows * 2 + 1];
    bPlus22 = godunovData.tractionPlusMatrix[Rows * 2 + 2];
    bMinus11 = godunovData.tractionMinusMatrix[Rows * 1 + 1];
    bMinus12 = godunovData.tractionMinusMatrix[Rows * 1 + 2];
    bMinus21 = godunovData.tractionMinusMatrix[Rows * 2 + 1];
    bMinus22 = godunovData.tractionMinusMatrix[Rows * 2 + 2];
  } else {
    bPlus11 = impAndEta.etaS * impAndEta.invZs;
    bPlus12 = 0;
    bPlus21 = 0;
    bPlus22 = impAndEta.etaS * impAndEta.invZs;
    bMinus11 = impAndEta.etaS * impAndEta.invZsNeig;
    bMinus12 = 0;
    bMinus21 = 0;
    bMinus22 = impAndEta.etaS * impAndEta.invZsNeig;
  }

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

      const auto qIPlusT1 = qIPlus[o][T1][i];
      const auto qIPlusT2 = qIPlus[o][T2][i];
      const auto qIMinusT1 = qIMinus[o][T1][i];
      const auto qIMinusT2 = qIMinus[o][T2][i];

      const real interpolatedTraction12 =
          bPlus11 * qIMinusT1 + bPlus12 * qIMinusT2 + bMinus11 * qIPlusT1 + bMinus12 * qIPlusT2;
      const real interpolatedTraction13 =
          bPlus21 * qIMinusT1 + bPlus22 * qIMinusT2 + bMinus21 * qIPlusT1 + bMinus22 * qIPlusT2;

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

/**
  Anisotropy projection handling.
  Has no effect for isotropy.

  Returns {etaProj, invEtaProj}
 */
SEISSOL_HOSTDEVICE inline std::pair<real, real>
    projectEta(const ImpedancesAndEta& impAndEta,
               [[maybe_unused]] const ImpedanceMatrices& impedanceMatrices,
               [[maybe_unused]] real t1,
               [[maybe_unused]] real t2,
               [[maybe_unused]] real tmag) {
  if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
    constexpr std::uint32_t Count =
        model::MaterialT::Type == model::MaterialType::Poroelastic ? 4 : 3;

    const real n1 = (tmag > 0) ? (t1 / tmag) : static_cast<real>(1.0);
    const real n2 = (tmag > 0) ? (t2 / tmag) : static_cast<real>(0.0);

    const real etaProj = impedanceMatrices.eta[Count * 1 + 1] * n1 * n1 +
                         impedanceMatrices.eta[Count * 1 + 2] * n1 * n2 +
                         impedanceMatrices.eta[Count * 2 + 1] * n2 * n1 +
                         impedanceMatrices.eta[Count * 2 + 2] * n2 * n2;

    return {etaProj, static_cast<real>(1.0) / etaProj};
  } else {
    return {impAndEta.etaS, impAndEta.invEtaS};
  }
}

/**
  Anisotropy normal/shear coupling handling.
  Has no effect for isotropy.

  Returns the coefficient

    c = (eta * n)_n = eta_ns * n1 + eta_nd * n2

  for the unit shear direction n = (t1, t2) / tmag. With it, the fault-normal stress after the
  friction solve is

    sigma(V) = sigma_stick - V * c,

  i.e. shear slip changes the normal stress whenever eta is not block diagonal in the fault-normal
  direction. c is identically zero for every isotropic material, so the whole correction disappears
  there.
 */
SEISSOL_HOSTDEVICE inline real
    projectEtaNormal([[maybe_unused]] const ImpedancesAndEta& impAndEta,
                     [[maybe_unused]] const ImpedanceMatrices& impedanceMatrices,
                     [[maybe_unused]] real t1,
                     [[maybe_unused]] real t2,
                     [[maybe_unused]] real tmag) {
  if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
    // the anisotropic block is always 3x3 (no fluid pressure component)
    constexpr std::uint32_t Count = 3;

    const real n1 = (tmag > 0) ? (t1 / tmag) : static_cast<real>(1.0);
    const real n2 = (tmag > 0) ? (t2 / tmag) : static_cast<real>(0.0);

    // eta is a dense, column-major tensor: eta[col * Count + row]
    return impedanceMatrices.eta[Count * 1 + 0] * n1 + impedanceMatrices.eta[Count * 2 + 0] * n2;
  } else {
    return static_cast<real>(0.0);
  }
}

/**
  Anisotropy direction handling.
  Has no effect for isotropy.

  Returns a 2-element vector
 */
SEISSOL_HOSTDEVICE inline std::pair<real, real>
    matmulEta(const ImpedancesAndEta& impAndEta,
              [[maybe_unused]] const ImpedanceMatrices& impedanceMatrices,
              real v1,
              real v2) {
  if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
    constexpr std::uint32_t Count =
        model::MaterialT::Type == model::MaterialType::Poroelastic ? 4 : 3;

    const real w1 =
        impedanceMatrices.eta[Count * 1 + 1] * v1 + impedanceMatrices.eta[Count * 1 + 2] * v2;

    const real w2 =
        impedanceMatrices.eta[Count * 2 + 1] * v1 + impedanceMatrices.eta[Count * 2 + 2] * v2;

    return {w1, w2};
  } else {
    return {impAndEta.etaS * v1, impAndEta.etaS * v2};
  }
}

/**
  Anisotropy normal/shear coupling, for a known slip rate vector.
  Has no effect for isotropy.

  Returns the fault-normal component of eta * (0, v1, v2), i.e.

    (eta * v)_n = eta_ns * v1 + eta_nd * v2,

  so that the fault-normal traction after the solve is sigma = sigma_stick - (eta * v)_n. This is
  the counterpart of projectEtaNormal for the case where the slip rate vector -- and not just its
  direction -- is known.
 */
SEISSOL_HOSTDEVICE inline real
    matmulEtaNormal([[maybe_unused]] const ImpedancesAndEta& impAndEta,
                    [[maybe_unused]] const ImpedanceMatrices& impedanceMatrices,
                    [[maybe_unused]] real v1,
                    [[maybe_unused]] real v2) {
  if constexpr (model::MaterialT::Type == model::MaterialType::Anisotropic) {
    // the anisotropic block is always 3x3 (no fluid pressure component)
    constexpr std::uint32_t Count = 3;

    // eta is a dense, column-major tensor: eta[col * Count + row]
    return impedanceMatrices.eta[Count * 1 + 0] * v1 + impedanceMatrices.eta[Count * 2 + 0] * v2;
  } else {
    return static_cast<real>(0.0);
  }
}

} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
