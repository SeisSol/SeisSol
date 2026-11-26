// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_

#include "Core/Common/Constants.h"
#include "Core/Common/Marker.h"
#include "Core/Common/Real.h"
#include "Core/Numerical/Functions.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawContext.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawDriver.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolverCommon.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/FrictionSolverDetails.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/FrictionSolverInterface.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "FrictionSolverInterface.h"
#include "GeneratedCode/init.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <algorithm>

namespace seissol::dr::friction_law::gpu {

// what a friction law works on at the point of a thread (see friction_law::FrictionLawContext)
template <typename Cfg>
using FrictionLawContext = friction_law::FrictionLawContext<Cfg, common::RangeType::GPU>;

template <typename Cfg, typename Derived>
class BaseFrictionSolver : public FrictionSolverDetails<Cfg> {
  private:
  // only the CRTP subclass `Derived` constructs this base
  explicit BaseFrictionSolver(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverDetails<Cfg>(drParameters) {}
  friend Derived;

  public:
  using real = Real<Cfg>;

  ~BaseFrictionSolver() override = default;

  SEISSOL_DEVICE static void evaluatePoint(FrictionLawContext<Cfg>& __restrict ctx) {
    FrictionLawDriver<Cfg, Derived>::evaluate(ctx);
  }

  void setupLayer(DynamicRupture::Layer& layerData,
                  seissol::parallel::runtime::StreamRuntime& runtime) override {
    this->currLayerSize_ = layerData.size();
    FrictionSolverInterface<Cfg>::copyStorageToLocal(&this->dataHost_, layerData);
    Derived::copySpecificStorageDataToLocal(&this->dataHost_, layerData);
    this->dataHost_.drParameters = this->drParameters_;
    device::DeviceInstance::instance().api().copyToAsync(
        this->data_, &this->dataHost_, sizeof(FrictionLawData<Cfg>), runtime.stream());
  }

  void evaluateKernel(seissol::parallel::runtime::StreamRuntime& runtime,
                      double fullUpdateTime,
                      const double* timeWeights,
                      const FrictionSolver::FrictionTime& frictionTime);

  void evaluate(double fullUpdateTime,
                const FrictionSolver::FrictionTime& frictionTime,
                const double* timeWeights,
                seissol::parallel::runtime::StreamRuntime& runtime) override {
    if (this->currLayerSize_ == 0) {
      return;
    }

    if constexpr (model::MaterialOf<Cfg>::SupportsDR) {
      evaluateKernel(runtime, fullUpdateTime, timeWeights, frictionTime);
    } else {
      logError() << "The material" << model::MaterialOf<Cfg>::Text
                 << "does not support DR friction law computations.";
    }
  }
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_BASEFRICTIONSOLVER_H_
