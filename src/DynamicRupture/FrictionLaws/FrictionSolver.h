// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Parallel/Runtime/Stream.h"

#include <functional>
#include <memory>
#include <vector>

namespace seissol::dr::friction_law {
/**
 * Abstract Base for friction solver class with the public interface
 * Only needed to be able to store a unique_ptr<FrictionSolver> in TimeCluster.
 * BaseFrictionLaw has a template argument for CRTP, hence, we can't store a pointer to any
 * BaseFrictionLaw. It also does not depend on the configuration the friction law computes in.
 *
 * Note: this class, FrictionSolver, must be trivially copyable. It is (or was?) important for GPU
 * offloading.
 */
class FrictionSolver {
  public:
  virtual ~FrictionSolver() = default;

  struct FrictionTime {
    std::vector<double> deltaT;
  };

  virtual void setupLayer(DynamicRupture::Layer& layerData,
                          seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void evaluate(double fullUpdateTime,
                        const FrictionTime& frictionTime,
                        const double* timeWeights,
                        seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  /**
   * compute the DeltaT from the current timePoints of the configuration `Cfg`; call this function
   * before evaluate to set the correct DeltaT
   */
  template <typename Cfg>
  static FrictionTime computeDeltaT(const std::vector<double>& timePoints);

  virtual seissol::initializer::AllocationPlace allocationPlace() {
    return seissol::initializer::AllocationPlace::Host;
  }
};

/**
 * The data of a friction solver of the configuration `Cfg`, taken from the storage of a layer.
 */
template <typename Cfg>
class FrictionSolverImpl : public FrictionSolver {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  explicit FrictionSolverImpl(const FrictionLawParameters<Real<Cfg>>& userDRParameters)
      : drParameters_(userDRParameters) {}

  /**
   * copies all common parameters from the DynamicRupture LTS to the local attributes
   */
  void copyStorageToLocal(DynamicRupture::Layer& layerData);

  virtual void allocateAuxiliaryMemory(GlobalData<Cfg>* globalData) {
    spaceWeights_ = globalData->*init::quadweights<Cfg>::PoolMember;
  }

  protected:
  /**
   * Adjust initial stress by adding nucleation stress * nucleation function
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf.
   */
  real deltaT_[misc::TimeSteps<Cfg>] = {};

  FrictionLawParameters<Real<Cfg>> drParameters_;
  ImpedancesAndEta<Cfg>* __restrict impAndEta_{};
  ImpedanceMatrices<Cfg>* __restrict impedanceMatrices_{};
  real fullUpdateTime_{};
  // CS = coordinate system
  real (*__restrict stressSourceInFaultCS_)[6][misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict cohesion_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict mu_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict accumulatedSlipMagnitude_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slip1_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slip2_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRateMagnitude_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRate1_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict slipRate2_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict ruptureTime_)[misc::NumPaddedPoints<Cfg>]{};
  bool (*__restrict ruptureTimePending_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict peakSlipRate_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict traction1_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict traction2_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict imposedStatePlus_)[tensor::QInterpolated<Cfg>::size()]{};
  real (*__restrict imposedStateMinus_)[tensor::QInterpolated<Cfg>::size()]{};
  const real* __restrict spaceWeights_{};
  DREnergyOutput<Cfg>* __restrict energyData_{};
  DRGodunovData<Cfg>* __restrict godunovData_{};
  real (*__restrict stressSourcePressure_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict stressSourceOnset_)[misc::NumPaddedPoints<Cfg>]{};
  real (*__restrict stressSourceRiseTime_)[misc::NumPaddedPoints<Cfg>]{};

  // be careful only for some FLs initialized:
  real (*__restrict dynStressTime_)[misc::NumPaddedPoints<Cfg>]{};
  bool (*__restrict dynStressTimePending_)[misc::NumPaddedPoints<Cfg>]{};

  real (*__restrict qInterpolatedPlus_)[misc::TimeSteps<Cfg>][tensor::QInterpolated<Cfg>::size()]{};
  real (*__restrict qInterpolatedMinus_)[misc::TimeSteps<Cfg>]
                                        [tensor::QInterpolated<Cfg>::size()]{};
};

/// Creates a friction solver for the configuration with the given id.
using FrictionSolverFactory = std::function<std::unique_ptr<FrictionSolver>(ConfigId)>;

/// The friction solver that `factory` creates for the configuration `Cfg`.
template <typename Cfg>
std::unique_ptr<FrictionSolverImpl<Cfg>> makeFrictionSolver(const FrictionSolverFactory& factory) {
  // the factory creates the solver of the configuration it is asked for
  return std::unique_ptr<FrictionSolverImpl<Cfg>>(
      static_cast<FrictionSolverImpl<Cfg>*>(factory(configIdOf<Cfg>()).release()));
}

/// The factory of the friction solvers `T<Cfg>`, with the parameters `parameters`.
template <template <typename> typename T>
FrictionSolverFactory
    makeFrictionSolverFactory(const seissol::initializer::parameters::DRParameters& parameters) {
  return [parameters](ConfigId config) {
    return dispatchConfig(config, [&](auto cfg) -> std::unique_ptr<FrictionSolver> {
      using Cfg = decltype(cfg);
      return std::make_unique<T<Cfg>>(FrictionLawParameters<Real<Cfg>>(parameters));
    });
  };
}

} // namespace seissol::dr::friction_law

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_FRICTIONSOLVER_H_
