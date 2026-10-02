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
#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "GeneratedCode/tensor.h"
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

template <typename Cfg, RangeType Type>
struct NumPoints {
  private:
  using CpuRange = ForLoopRange<0, dr::misc::NumPaddedPoints<Cfg>, 1>;
  using GpuRange = ForLoopRange<0, 1, 1>;

  public:
  // Range::Start is 0, and Range::End is seissol::misc::NumPaddedPoints for CPU
  using Range = std::conditional_t<Type == RangeType::CPU, CpuRange, GpuRange>;
};

template <typename Cfg, RangeType Type>
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
  static constexpr Real<Cfg>& index(Real<Cfg> (&data)[misc::NumPaddedPoints<Cfg>], int i) {
    return data[i];
  }

  static constexpr Real<Cfg> index(const Real<Cfg> (&data)[misc::NumPaddedPoints<Cfg>], int i) {
    return data[i];
  }
};

template <typename Cfg>
struct VariableIndexing<Cfg, Executor::Device> {
  static constexpr Real<Cfg>& index(Real<Cfg>& data, int /*i*/) { return data; }

  static constexpr Real<Cfg> index(const Real<Cfg>& data, int /*i*/) { return data; }
};

/**
 * Calculate traction and normal stress at the interface of a face.
 * Using equations (A2) from Pelties et al. 2014
 * Definiton of eta and impedance Z are found in Carsten Uphoff's dissertation on page 47 and in
 * equation (4.51) respectively.
 * Only handles a single timestep
 *
 * @param[out] faultStresses contains normalStress, traction1, traction2
 *             at the 2d face quadrature nodes evaluated at the selected time
 *             quadrature point
 * @param[in] impAndEta contains eta and impedance values
 * @param[in] impedanceMatrices contains impedance and eta values, in the poroelastic case, these
 * are non-diagonal matrices
 * @param[in] qInterpolatedPlus a plus side dofs interpolated at time sub-intervals
 * @param[in] qInterpolatedMinus a minus side dofs interpolated at time sub-intervals
 * @param[in] step the timestep to handle
 */
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void precomputeStressFromQInterpolated(
    FaultStresses<Cfg, RangeExecutor<Type>::Exec>& __restrict faultStresses,
    const ImpedancesAndEta<Cfg>& __restrict impAndEta,
    [[maybe_unused]] const ImpedanceMatrices<Cfg>& __restrict impedanceMatrices,
    const Real<Cfg> qInterpolatedPlus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const Real<Cfg> qInterpolatedMinus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    Real<Cfg> etaPDamp,
    uint32_t step,
    uint32_t startLoopIndex = 0) {
  static_assert(tensor::QInterpolated<Cfg>::Shape[multisim::BasisDim<Cfg>] ==
                    tensor::resample<Cfg>::Shape[0],
                "Different number of quadrature points?");

  const auto o = step;

  using QInterpolatedShapeT =
      const Real<Cfg>(*__restrict)[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>];
  const auto* __restrict qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
  const auto* __restrict qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));

  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Elastic ||
                model::MaterialOf<Cfg>::Type == model::MaterialType::Viscoelastic) {
    const auto etaP = impAndEta.etaP * etaPDamp;
    const auto etaS = impAndEta.etaS;
    const auto invZp = impAndEta.invZp;
    const auto invZs = impAndEta.invZs;
    const auto invZpNeig = impAndEta.invZpNeig;
    const auto invZsNeig = impAndEta.invZsNeig;

    using namespace dr::misc::quantity_indices;

    using Range = typename NumPoints<Cfg, Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto i{startLoopIndex + index};
      VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, i) =
          etaP * (qIMinus[o][U][i] - qIPlus[o][U][i] + qIPlus[o][N][i] * invZp +
                  qIMinus[o][N][i] * invZpNeig);

      VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.traction1, i) =
          etaS * (qIMinus[o][V][i] - qIPlus[o][V][i] + qIPlus[o][T1][i] * invZs +
                  qIMinus[o][T1][i] * invZsNeig);

      VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.traction2, i) =
          etaS * (qIMinus[o][W][i] - qIPlus[o][W][i] + qIPlus[o][T2][i] * invZs +
                  qIMinus[o][T2][i] * invZsNeig);
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

    using namespace dr::misc::quantity_indices;

    using Range = typename NumPoints<Cfg, Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      auto i{startLoopIndex + index};

      constexpr uint32_t Count =
          model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic ? 4 : 3;

      // Compute Theta from eq (4.53) in Carsten's thesis

      Real<Cfg> velDiff[Count]{};
      velDiff[0] = qIMinus[o][U][i] - qIPlus[o][U][i];
      velDiff[1] = qIMinus[o][V][i] - qIPlus[o][V][i];
      velDiff[2] = qIMinus[o][W][i] - qIPlus[o][W][i];

      Real<Cfg> strP[Count]{};
      Real<Cfg> strM[Count]{};
      const auto rowCompute = [&](auto linear, auto qindex) {
#pragma unroll
        for (std::uint32_t j = 0; j < Count; ++j) {
          strP[j] += impedanceMatrices.impedance[linear * Count + j] * qIPlus[o][qindex][i];
          strM[j] += impedanceMatrices.impedanceNeig[linear * Count + j] * qIMinus[o][qindex][i];
        }
      };
      rowCompute(0, N);
      rowCompute(1, T1);
      rowCompute(2, T2);

      if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic) {
        velDiff[3] = qIMinus[o][FU][i] - qIPlus[o][FU][i];
        rowCompute(3, FP);
      }

      Real<Cfg> res[Count]{};
#pragma unroll
      for (std::uint32_t k = 0; k < Count; ++k) {
#pragma unroll
        for (std::uint32_t j = 0; j < Count; ++j) {
          res[j] += impedanceMatrices.eta[k * Count + j] * (velDiff[k] + strP[k] + strM[k]);
        }
      }

      VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, i) =
          res[0] * etaPDamp;
      VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.traction1, i) = res[1];
      VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.traction2, i) = res[2];
      if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic) {
        VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.fluidPressure, i) =
            res[3];
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
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void initializeTractionResults(
    const FaultStresses<Cfg, RangeExecutor<Type>::Exec>& __restrict faultStresses,
    TractionResults<Cfg, RangeExecutor<Type>::Exec>& __restrict tractionResults,
    uint32_t startIndex = 0) {
  using Range = typename NumPoints<Cfg, Type>::Range;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    const auto i{startIndex + index};
    VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(tractionResults.normalStress, i) =
        VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(faultStresses.normalStress, i);
  }
}

/**
 * Integrate over all Time points with the time weights and calculate the traction for each side
 * according to Carsten Uphoff Thesis: EQ.: 4.60
 *
 * @param[inout] state
 * @param[in] faultStresses
 * @param[in] tractionResults
 * @param[in] impAndEta
 * @param[in] impedancenceMatrices
 * @param[in] qInterpolatedPlus
 * @param[in] qInterpolatedMinus
 * @param[in] step
 * @param[in] weight
 */
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void postcomputeImposedStateFromNewStress(
    ImposedState<Cfg, RangeExecutor<Type>::Exec>& __restrict state,
    [[maybe_unused]] const FaultStresses<Cfg, RangeExecutor<Type>::Exec>& __restrict faultStresses,
    const TractionResults<Cfg, RangeExecutor<Type>::Exec>& __restrict tractionResults,
    const ImpedancesAndEta<Cfg>& __restrict impAndEta,
    [[maybe_unused]] const ImpedanceMatrices<Cfg>& __restrict impedanceMatrices,
    const Real<Cfg> qInterpolatedPlus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const Real<Cfg> qInterpolatedMinus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    uint32_t step,
    Real<Cfg> weight,
    uint32_t startIndex = 0) {

  using NumPointsRange = typename NumPoints<Cfg, Type>::Range;

  const auto o = step;

  using Acc = VariableIndexing<Cfg, RangeExecutor<Type>::Exec>;

  using QInterpolatedShapeT =
      const Real<Cfg>(*__restrict)[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>];
  const auto* __restrict qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  const auto* __restrict qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Elastic ||
                model::MaterialOf<Cfg>::Type == model::MaterialType::Viscoelastic) {
    const auto invZs = impAndEta.invZs;
    const auto invZp = impAndEta.invZp;
    const auto invZsNeig = impAndEta.invZsNeig;
    const auto invZpNeig = impAndEta.invZpNeig;

    using namespace dr::misc::quantity_indices;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
         index += NumPointsRange::Step) {
      auto i{startIndex + index};

      const auto normalStress = Acc::index(tractionResults.normalStress, i);
      const auto traction1 = Acc::index(tractionResults.traction1, i);
      const auto traction2 = Acc::index(tractionResults.traction2, i);

      Acc::index(state.minus[N], i) += weight * normalStress;
      Acc::index(state.minus[T1], i) += weight * traction1;
      Acc::index(state.minus[T2], i) += weight * traction2;
      Acc::index(state.minus[U], i) +=
          weight * (qIMinus[o][U][i] - invZpNeig * (normalStress - qIMinus[o][N][i]));
      Acc::index(state.minus[V], i) +=
          weight * (qIMinus[o][V][i] - invZsNeig * (traction1 - qIMinus[o][T1][i]));
      Acc::index(state.minus[W], i) +=
          weight * (qIMinus[o][W][i] - invZsNeig * (traction2 - qIMinus[o][T2][i]));

      Acc::index(state.plus[N], i) += weight * normalStress;
      Acc::index(state.plus[T1], i) += weight * traction1;
      Acc::index(state.plus[T2], i) += weight * traction2;
      Acc::index(state.plus[U], i) +=
          weight * (qIPlus[o][U][i] + invZp * (normalStress - qIPlus[o][N][i]));
      Acc::index(state.plus[V], i) +=
          weight * (qIPlus[o][V][i] + invZs * (traction1 - qIPlus[o][T1][i]));
      Acc::index(state.plus[W], i) +=
          weight * (qIPlus[o][W][i] + invZs * (traction2 - qIPlus[o][T2][i]));
    }
  } else {
    // poroelastic kernel (for CPU+GPU)
    // TODO: generalize and unify with the above (probably using either templates or Yateto)
    // (the v1.1.0-1.3.1 Yateto+selector matrix based kernel was removed since GPU support was
    // missing)

    using namespace dr::misc::quantity_indices;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
         index += NumPointsRange::Step) {
      auto i{startIndex + index};

      const auto normalStress = Acc::index(tractionResults.normalStress, i);
      const auto traction1 = Acc::index(tractionResults.traction1, i);
      const auto traction2 = Acc::index(tractionResults.traction2, i);
      const auto fluidPressure = Acc::index(faultStresses.fluidPressure, i);

      const auto handleSide =
          [&](auto& imposedState, const auto& qI, const auto& mZ, Real<Cfg> sign) {
            constexpr std::uint32_t Count =
                model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic ? 4 : 3;

            Acc::index(imposedState[N], i) += weight * normalStress;
            Acc::index(imposedState[T1], i) += weight * traction1;
            Acc::index(imposedState[T2], i) += weight * traction2;

            Real<Cfg> diff[Count]{};
            diff[0] = (normalStress - qI[o][N][i]) * sign;
            diff[1] = (traction1 - qI[o][T1][i]) * sign;
            diff[2] = (traction2 - qI[o][T2][i]) * sign;

            if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic) {
              Acc::index(imposedState[FP], i) += weight * fluidPressure;
              diff[3] = (fluidPressure - qI[o][FP][i]) * sign;
            }

            const auto handleEntry = [&](auto linear, auto qindex) {
              Real<Cfg> acc = 0;
#pragma unroll
              for (std::uint32_t k = 0; k < Count; ++k) {
                acc += mZ[Count * k + linear] * diff[k];
              }
              Acc::index(imposedState[qindex], i) += weight * (qI[o][qindex][i] + acc);
            };

            handleEntry(0, U);
            handleEntry(1, V);
            handleEntry(2, W);
            if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Poroelastic) {
              handleEntry(3, FU);
            }
          };

      handleSide(state.minus, qIMinus, impedanceMatrices.impedanceNeig, -1);
      handleSide(state.plus, qIPlus, impedanceMatrices.impedance, 1);
    }
  }
}

/**
 * Store the imposed state to memory. (and potentially run some last accumulation steps)
 *
 * @param[in] state
 * @param[out] imposedStatePlus
 * @param[out] imposedStateMinus
 */
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void
    finalizeImposedState(const ImposedState<Cfg, RangeExecutor<Type>::Exec>& __restrict state,
                         Real<Cfg> imposedStatePlus[tensor::QInterpolated<Cfg>::size()],
                         Real<Cfg> imposedStateMinus[tensor::QInterpolated<Cfg>::size()],
                         uint32_t startIndex = 0) {

  using NumPointsRange = typename NumPoints<Cfg, Type>::Range;

  using ImposedStateShapeT = Real<Cfg>(*__restrict)[misc::NumPaddedPoints<Cfg>];
  auto* __restrict imposedStateP = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
  auto* __restrict imposedStateM = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

  for (auto index = NumPointsRange::Start; index < NumPointsRange::End;
       index += NumPointsRange::Step) {
    auto i{startIndex + index};
#pragma unroll
    for (std::uint32_t q = 0; q < dr::misc::NumQuantities<Cfg>; ++q) {
      imposedStateM[q][i] =
          VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(state.minus[q], i);
      imposedStateP[q][i] =
          VariableIndexing<Cfg, RangeExecutor<Type>::Exec>::index(state.plus[q], i);
    }
  }
}

/**
 * The initial stress in effect at the given time, in the layout of the fault stresses it is added
 * to: every stress source of the face at the fraction of it that has been applied so far.
 *
 * The initial state is the source without a rise time, so it enters the sum like any other and
 * needs no case of its own.
 *
 * @param[out] initialStress
 * @param[in] stressSourceInFaultCS the stress of every source of this face
 * @param[in] stressSourcePressure
 * @param[in] stressSourceOnset the onset of every source of this face, per point
 * @param[in] stressSourceRiseTime the rise time of every source of this face, per point
 * @param[in] sourceCount
 * @param[in] fullUpdateTime
 */
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void computeInitialStress(
    FaultStresses<Cfg, RangeExecutor<Type>::Exec>& __restrict initialStress,
    const Real<Cfg> (*__restrict stressSourceInFaultCS)[6][misc::NumPaddedPoints<Cfg>],
    const Real<Cfg> (*__restrict stressSourcePressure)[misc::NumPaddedPoints<Cfg>],
    const Real<Cfg> (*__restrict stressSourceOnset)[misc::NumPaddedPoints<Cfg>],
    const Real<Cfg> (*__restrict stressSourceRiseTime)[misc::NumPaddedPoints<Cfg>],
    std::uint32_t sourceCount,
    Real<Cfg> fullUpdateTime,
    uint32_t startIndex = 0) {
  constexpr auto Exec = RangeExecutor<Type>::Exec;
  using Range = typename NumPoints<Cfg, Type>::Range;

  // the components of the stress tensor which take part in the fault-normal Riemann problem
  constexpr std::size_t NormalIndex = 0;
  constexpr std::size_t Traction1Index = 3;
  constexpr std::size_t Traction2Index = 5;

#ifndef ACL_DEVICE
#pragma omp simd
#endif
  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    const auto i{startIndex + index};
    VariableIndexing<Cfg, Exec>::index(initialStress.normalStress, i) = static_cast<Real<Cfg>>(0.0);
    VariableIndexing<Cfg, Exec>::index(initialStress.traction1, i) = static_cast<Real<Cfg>>(0.0);
    VariableIndexing<Cfg, Exec>::index(initialStress.traction2, i) = static_cast<Real<Cfg>>(0.0);
    VariableIndexing<Cfg, Exec>::index(initialStress.fluidPressure, i) =
        static_cast<Real<Cfg>>(0.0);
  }

  for (std::uint32_t source = 0; source < sourceCount; ++source) {
#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (auto index = Range::Start; index < Range::End; index += Range::Step) {
      const auto i{startIndex + index};
      const auto fraction = stressSourceFraction<Real<Cfg>>(
          fullUpdateTime, stressSourceRiseTime[source][i], stressSourceOnset[source][i]);
      VariableIndexing<Cfg, Exec>::index(initialStress.normalStress, i) +=
          stressSourceInFaultCS[source][NormalIndex][i] * fraction;
      VariableIndexing<Cfg, Exec>::index(initialStress.traction1, i) +=
          stressSourceInFaultCS[source][Traction1Index][i] * fraction;
      VariableIndexing<Cfg, Exec>::index(initialStress.traction2, i) +=
          stressSourceInFaultCS[source][Traction2Index][i] * fraction;
      VariableIndexing<Cfg, Exec>::index(initialStress.fluidPressure, i) +=
          stressSourcePressure[source][i] * fraction;
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
                           uint32_t startIndex = 0) {

  using Range = typename NumPoints<Cfg, Type>::Range;

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
                           uint32_t startIndex = 0) {

  using Range = typename NumPoints<Cfg, Type>::Range;

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
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void updateTimeSinceSlipRateBelowThreshold(
    const Real<Cfg> slipRateMagnitude[misc::NumPaddedPoints<Cfg>],
    const bool ruptureTimePending[misc::NumPaddedPoints<Cfg>],
    // See https://github.com/llvm/llvm-project/issues/60163
    // NOLINTNEXTLINE
    DREnergyOutput<Cfg>& __restrict energyData,
    const Real<Cfg> dt,
    const Real<Cfg> slipRateThreshold,
    uint32_t startIndex = 0) {

  using Range = typename NumPoints<Cfg, Type>::Range;
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
      timeSinceSlipRateBelowThreshold[pointIndex] = std::numeric_limits<Real<Cfg>>::infinity();
    }
  }
}
template <typename Cfg, RangeType Type = RangeType::CPU>
SEISSOL_HOSTDEVICE inline void computeFrictionEnergy(
    DREnergyOutput<Cfg>& __restrict energyData,
    const Real<Cfg> qInterpolatedPlus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const Real<Cfg> qInterpolatedMinus[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()],
    const ImpedancesAndEta<Cfg>& __restrict impAndEta,
    const Real<Cfg> timeWeights[misc::TimeSteps<Cfg>],
    const Real<Cfg> spaceWeights[seissol::kernels::NumSpaceQuadraturePoints<Cfg>],
    const DRGodunovData<Cfg>& __restrict godunovData,
    const Real<Cfg> slipRateMagnitude[misc::NumPaddedPoints<Cfg>],
    const bool energiesFromAcrossFaultVelocities,
    size_t startIndex = 0) {

  auto* slip = reinterpret_cast<Real<Cfg>(*)[misc::NumPaddedPoints<Cfg>]>(energyData.slip);
  auto* accumulatedSlip = energyData.accumulatedSlip;
  auto* frictionalEnergy = energyData.frictionalEnergy;
  const Real<Cfg> doubledSurfaceAreaN = -static_cast<Real<Cfg>>(godunovData.doubledSurfaceArea);

  using QInterpolatedShapeT =
      const Real<Cfg>(*)[misc::NumQuantities<Cfg>][misc::NumPaddedPoints<Cfg>];
  const auto* __restrict qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  const auto* __restrict qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);

  using namespace dr::misc::quantity_indices;

  Real<Cfg> bPlus11{};
  Real<Cfg> bPlus12{};
  Real<Cfg> bPlus21{};
  Real<Cfg> bPlus22{};
  Real<Cfg> bMinus11{};
  Real<Cfg> bMinus12{};
  Real<Cfg> bMinus21{};
  Real<Cfg> bMinus22{};
  // the fault-normal column: with an anisotropic impedance the normal traction contributes to the
  // interpolated *shear* traction as well
  Real<Cfg> bPlus10{};
  Real<Cfg> bPlus20{};
  Real<Cfg> bMinus10{};
  Real<Cfg> bMinus20{};

  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
    constexpr auto Rows = 3;
    bPlus10 = godunovData.tractionPlusMatrix[Rows * 1 + 0];
    bPlus11 = godunovData.tractionPlusMatrix[Rows * 1 + 1];
    bPlus12 = godunovData.tractionPlusMatrix[Rows * 1 + 2];
    bPlus20 = godunovData.tractionPlusMatrix[Rows * 2 + 0];
    bPlus21 = godunovData.tractionPlusMatrix[Rows * 2 + 1];
    bPlus22 = godunovData.tractionPlusMatrix[Rows * 2 + 2];
    bMinus10 = godunovData.tractionMinusMatrix[Rows * 1 + 0];
    bMinus11 = godunovData.tractionMinusMatrix[Rows * 1 + 1];
    bMinus12 = godunovData.tractionMinusMatrix[Rows * 1 + 2];
    bMinus20 = godunovData.tractionMinusMatrix[Rows * 2 + 0];
    bMinus21 = godunovData.tractionMinusMatrix[Rows * 2 + 1];
    bMinus22 = godunovData.tractionMinusMatrix[Rows * 2 + 2];
  } else {
    bPlus10 = 0;
    bPlus11 = impAndEta.etaS * impAndEta.invZs;
    bPlus12 = 0;
    bPlus20 = 0;
    bPlus21 = 0;
    bPlus22 = impAndEta.etaS * impAndEta.invZs;
    bMinus10 = 0;
    bMinus11 = impAndEta.etaS * impAndEta.invZsNeig;
    bMinus12 = 0;
    bMinus20 = 0;
    bMinus21 = 0;
    bMinus22 = impAndEta.etaS * impAndEta.invZsNeig;
  }

  using Range = typename NumPoints<Cfg, Type>::Range;
  Real<Cfg> localAccumulatedSlip[Range::Size]{};
  Real<Cfg> localFrictionalEnergy[Range::Size]{};
  Real<Cfg> localSlip[3][Range::Size]{};

  for (auto index = Range::Start; index < Range::End; index += Range::Step) {
    auto i{startIndex + index};
    localAccumulatedSlip[index] = accumulatedSlip[i];
    localFrictionalEnergy[index] = frictionalEnergy[i];
#pragma unroll
    for (uint32_t d = 0; d < 3; ++d) {
      localSlip[d][index] = slip[d][i];
    }
  }

  for (size_t o = 0; o < misc::TimeSteps<Cfg>; ++o) {
    const auto timeWeight = timeWeights[o];

#ifndef ACL_DEVICE
#pragma omp simd
#endif
    for (size_t index = Range::Start; index < Range::End; index += Range::Step) {

      const size_t i{startIndex + index}; // startIndex is always 0 for CPU

      const Real<Cfg> interpolatedSlipRate1 = qIMinus[o][U][i] - qIPlus[o][U][i];
      const Real<Cfg> interpolatedSlipRate2 = qIMinus[o][V][i] - qIPlus[o][V][i];
      const Real<Cfg> interpolatedSlipRate3 = qIMinus[o][W][i] - qIPlus[o][W][i];

      if (energiesFromAcrossFaultVelocities) {
        const Real<Cfg> interpolatedSlipRateMagnitude =
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

      const auto qIPlusN = qIPlus[o][N][i];
      const auto qIPlusT1 = qIPlus[o][T1][i];
      const auto qIPlusT2 = qIPlus[o][T2][i];
      const auto qIMinusN = qIMinus[o][N][i];
      const auto qIMinusT1 = qIMinus[o][T1][i];
      const auto qIMinusT2 = qIMinus[o][T2][i];

      // tau* = b+ tau+ + b- tau-, i.e. b+ pairs with the *plus* side -- matching the
      // computeTractionInterpolated kernel in EnergyOutput, which contracts tractionPlusMatrix
      // with QInterpolatedPlus. Only relevant for a bimaterial interface, where b+ != b-.
      const Real<Cfg> interpolatedTraction12 = bPlus10 * qIPlusN + bPlus11 * qIPlusT1 +
                                               bPlus12 * qIPlusT2 + bMinus10 * qIMinusN +
                                               bMinus11 * qIMinusT1 + bMinus12 * qIMinusT2;
      const Real<Cfg> interpolatedTraction13 = bPlus20 * qIPlusN + bPlus21 * qIPlusT1 +
                                               bPlus22 * qIPlusT2 + bMinus20 * qIMinusN +
                                               bMinus21 * qIMinusT1 + bMinus22 * qIMinusT2;

      const auto spaceWeight = spaceWeights[i / Cfg::NumSimulations];
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
template <typename Cfg>
SEISSOL_HOSTDEVICE inline std::pair<Real<Cfg>, Real<Cfg>>
    projectEta(const ImpedancesAndEta<Cfg>& impAndEta,
               [[maybe_unused]] const ImpedanceMatrices<Cfg>& impedanceMatrices,
               [[maybe_unused]] Real<Cfg> t1,
               [[maybe_unused]] Real<Cfg> t2,
               [[maybe_unused]] Real<Cfg> tmag) {
  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
    // the anisotropic block is always 3x3 (no fluid pressure component)
    constexpr std::uint32_t Count = 3;

    const Real<Cfg> n1 = (tmag > 0) ? (t1 / tmag) : static_cast<Real<Cfg>>(1.0);
    const Real<Cfg> n2 = (tmag > 0) ? (t2 / tmag) : static_cast<Real<Cfg>>(0.0);

    const Real<Cfg> etaProj = impedanceMatrices.eta[Count * 1 + 1] * n1 * n1 +
                              impedanceMatrices.eta[Count * 1 + 2] * n1 * n2 +
                              impedanceMatrices.eta[Count * 2 + 1] * n2 * n1 +
                              impedanceMatrices.eta[Count * 2 + 2] * n2 * n2;

    return {etaProj, static_cast<Real<Cfg>>(1.0) / etaProj};
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
template <typename Cfg>
SEISSOL_HOSTDEVICE inline Real<Cfg>
    projectEtaNormal([[maybe_unused]] const ImpedancesAndEta<Cfg>& impAndEta,
                     [[maybe_unused]] const ImpedanceMatrices<Cfg>& impedanceMatrices,
                     [[maybe_unused]] Real<Cfg> t1,
                     [[maybe_unused]] Real<Cfg> t2,
                     [[maybe_unused]] Real<Cfg> tmag) {
  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
    // the anisotropic block is always 3x3 (no fluid pressure component)
    constexpr std::uint32_t Count = 3;

    const Real<Cfg> n1 = (tmag > 0) ? (t1 / tmag) : static_cast<Real<Cfg>>(1.0);
    const Real<Cfg> n2 = (tmag > 0) ? (t2 / tmag) : static_cast<Real<Cfg>>(0.0);

    // eta is a dense, column-major tensor: eta[col * Count + row]
    return impedanceMatrices.eta[Count * 1 + 0] * n1 + impedanceMatrices.eta[Count * 2 + 0] * n2;
  } else {
    return static_cast<Real<Cfg>>(0.0);
  }
}

/**
  Anisotropy direction handling.
  Has no effect for isotropy.

  Returns a 2-element vector
 */
template <typename Cfg>
SEISSOL_HOSTDEVICE inline std::pair<Real<Cfg>, Real<Cfg>>
    matmulEta(const ImpedancesAndEta<Cfg>& impAndEta,
              [[maybe_unused]] const ImpedanceMatrices<Cfg>& impedanceMatrices,
              Real<Cfg> v1,
              Real<Cfg> v2) {
  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
    // the anisotropic block is always 3x3 (no fluid pressure component)
    constexpr std::uint32_t Count = 3;

    // eta is a dense, column-major tensor: eta[col * Count + row]
    const Real<Cfg> w1 =
        impedanceMatrices.eta[Count * 1 + 1] * v1 + impedanceMatrices.eta[Count * 2 + 1] * v2;

    const Real<Cfg> w2 =
        impedanceMatrices.eta[Count * 1 + 2] * v1 + impedanceMatrices.eta[Count * 2 + 2] * v2;

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
template <typename Cfg>
SEISSOL_HOSTDEVICE inline Real<Cfg>
    matmulEtaNormal([[maybe_unused]] const ImpedancesAndEta<Cfg>& impAndEta,
                    [[maybe_unused]] const ImpedanceMatrices<Cfg>& impedanceMatrices,
                    [[maybe_unused]] Real<Cfg> v1,
                    [[maybe_unused]] Real<Cfg> v2) {
  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
    // the anisotropic block is always 3x3 (no fluid pressure component)
    constexpr std::uint32_t Count = 3;

    // eta is a dense, column-major tensor: eta[col * Count + row]
    return impedanceMatrices.eta[Count * 1 + 0] * v1 + impedanceMatrices.eta[Count * 2 + 0] * v2;
  } else {
    return static_cast<Real<Cfg>>(0.0);
  }
}

/**
  Anisotropy: updates the slip direction.
  Has no effect for isotropy.

  With a non-isotropic shear impedance block, the slip rate and the trial traction are no longer
  collinear. The exact condition is

    tau0 = (S * I + V * eta_ss) n ,   |n| = 1

  hence n = normalize((S * I + V * eta_ss)^-1 tau0). The 2x2 inverse is written out through the
  adjugate; its determinant cancels when normalising, so it is never formed and there is no
  singularity unless S and V both vanish. For eta_ss = eta * I this returns tau0 / |tau0| exactly,
  i.e. the isotropic case is untouched.

  Sweeping n -> V -> n twice reaches machine precision; a single sweep already brings the relative
  error of |V| from ~1e-4 down to ~1e-8.

  @param[in] strength the fault strength belonging to the current slip rate. Note that it can
                      always be recovered as S = n^T tau0 - V * eta_proj, without evaluating the
                      friction law again.
  @param[in] slipRate the current slip rate magnitude V
  @param[in] t1, t2, tmag the trial (stick) shear traction and its magnitude
 */
template <typename Cfg>
SEISSOL_HOSTDEVICE inline std::pair<Real<Cfg>, Real<Cfg>>
    updateSlipDirection([[maybe_unused]] const ImpedancesAndEta<Cfg>& impAndEta,
                        [[maybe_unused]] const ImpedanceMatrices<Cfg>& impedanceMatrices,
                        [[maybe_unused]] Real<Cfg> strength,
                        [[maybe_unused]] Real<Cfg> slipRate,
                        Real<Cfg> t1,
                        Real<Cfg> t2,
                        Real<Cfg> tmag) {
  if constexpr (model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic) {
    // the anisotropic block is always 3x3 (no fluid pressure component)
    constexpr std::uint32_t Count = 3;

    // the very same 2x2 block, in the same convention, that matmulEta and projectEta use
    const Real<Cfg> e11 = impedanceMatrices.eta[Count * 1 + 1];
    const Real<Cfg> e12 = impedanceMatrices.eta[Count * 1 + 2];
    const Real<Cfg> e21 = impedanceMatrices.eta[Count * 2 + 1];
    const Real<Cfg> e22 = impedanceMatrices.eta[Count * 2 + 2];

    // adjugate of (S * I + V * eta_ss), applied to tau0
    const Real<Cfg> u1 = (strength + slipRate * e22) * t1 - slipRate * e12 * t2;
    const Real<Cfg> u2 = -slipRate * e21 * t1 + (strength + slipRate * e11) * t2;

    const Real<Cfg> umag = misc::magnitude(u1, u2);
    if (umag > 0) {
      const Real<Cfg> inv = static_cast<Real<Cfg>>(1.0) / umag;
      return {u1 * inv, u2 * inv};
    }
  }

  const Real<Cfg> n1 = (tmag > 0) ? (t1 / tmag) : static_cast<Real<Cfg>>(1.0);
  const Real<Cfg> n2 = (tmag > 0) ? (t2 / tmag) : static_cast<Real<Cfg>>(0.0);
  return {n1, n2};
}

/**
 * Slip rate magnitude and slip direction of a strength that is affine in the fault-normal traction.
 */
template <typename Cfg>
struct SlipRateSolution {
  Real<Cfg> slipRate{};
  Real<Cfg> direction1{};
  Real<Cfg> direction2{};
  /// the trial traction along the converged slip direction; equals strength + etaEff * slipRate
  Real<Cfg> projectedTraction{};
  /// the divisor the slip rate was obtained with, eta + slope * (eta n)_n
  Real<Cfg> etaEff{};
};

/**
 * Solves
 *
 *   tau0 = (S I + V eta_ss) n ,   S = strength + strengthSlope * V * (eta n)_n ,   |n| = 1
 *
 * for the slip rate V and the slip direction n. With an anisotropic impedance the slip is not
 * parallel to the trial traction, and the strength follows the fault-normal traction, which in turn
 * follows the slip rate. Sweeping n -> V -> n twice resolves both: the strength being affine in the
 * normal traction keeps the closed form for V exact, so only the direction has to be iterated, and
 * the first sweep alone reproduces the isotropic formula.
 *
 * Every projection is a no-op for an isotropic impedance, where n is the direction of tau0 and the
 * result reduces to V = (|tau0| - strength) / eta.
 */
template <typename Cfg>
SEISSOL_HOSTDEVICE inline SlipRateSolution<Cfg>
    solveSlipRate(const ImpedancesAndEta<Cfg>& impAndEta,
                  const ImpedanceMatrices<Cfg>& impedanceMatrices,
                  Real<Cfg> traction1,
                  Real<Cfg> traction2,
                  Real<Cfg> tractionMagnitude,
                  Real<Cfg> strength,
                  Real<Cfg> strengthSlope) {
  const Real<Cfg> invAbsolute = (tractionMagnitude > 0)
                                    ? static_cast<Real<Cfg>>(1.0) / tractionMagnitude
                                    : static_cast<Real<Cfg>>(0.0);
  Real<Cfg> n1 = traction1 * invAbsolute;
  Real<Cfg> n2 = traction2 * invAbsolute;
  Real<Cfg> projectedTraction = tractionMagnitude;
  Real<Cfg> eta =
      projectEta(impAndEta, impedanceMatrices, traction1, traction2, tractionMagnitude).first;
  Real<Cfg> etaNormal =
      projectEtaNormal(impAndEta, impedanceMatrices, traction1, traction2, tractionMagnitude);
  Real<Cfg> slipRate{};
  Real<Cfg> etaEff{};

  // the sweep can only move the direction where the shear block of eta is not a multiple of the
  // identity, so one pass is the exact closed form for every other material
  constexpr std::uint32_t DirectionSweeps =
      model::MaterialOf<Cfg>::Type == model::MaterialType::Anisotropic ? 2 : 1;
  for (std::uint32_t sweep = 0; sweep < DirectionSweeps; ++sweep) {
    // S(V) = S0 + slope * (eta * n)_n * V is exact, so the closed form survives
    etaEff = eta + strengthSlope * etaNormal;
    // a pathologically large coupling must never flip the sign of the divisor
    etaEff = (etaEff > 0) ? etaEff : eta;
    slipRate = std::max(static_cast<Real<Cfg>>(0.0), (projectedTraction - strength) / etaEff);

    if (sweep + 1 == DirectionSweeps) {
      break;
    }

    const Real<Cfg> localStrength = projectedTraction - slipRate * eta;
    const auto [d1, d2] = updateSlipDirection(impAndEta,
                                              impedanceMatrices,
                                              localStrength,
                                              slipRate,
                                              traction1,
                                              traction2,
                                              tractionMagnitude);
    n1 = d1;
    n2 = d2;
    projectedTraction = n1 * traction1 + n2 * traction2;
    eta = projectEta(impAndEta, impedanceMatrices, n1, n2, static_cast<Real<Cfg>>(1.0)).first;
    etaNormal = projectEtaNormal(impAndEta, impedanceMatrices, n1, n2, static_cast<Real<Cfg>>(1.0));
  }

  return {slipRate, n1, n2, projectedTraction, etaEff};
}

} // namespace seissol::dr::friction_law::common

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVERCOMMON_H_
