// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONSOLVER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONSOLVER_H_

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Real.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"

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
   * Makes the device read the time at the start of each step from the given clock, instead of
   * the time passed to evaluate().
   */
  virtual void setClock(const double* /*deviceClock*/) {}

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
  using real = Real<Cfg>;

  explicit FrictionSolverImpl(const FrictionLawParameters<Real<Cfg>>& userDRParameters)
      : drParameters_(userDRParameters) {}

  virtual void allocateAuxiliaryMemory(GlobalData<Cfg>* globalData) {
    spaceWeights_ = globalData->*init::quadweights<Cfg>::PoolMember;
  }

  protected:
  FrictionLawParameters<Real<Cfg>> drParameters_;
  const real* __restrict spaceWeights_{};
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

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_FRICTIONSOLVER_H_
