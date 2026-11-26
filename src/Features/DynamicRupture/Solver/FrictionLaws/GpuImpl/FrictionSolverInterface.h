// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_

#include "Core/Common/Real.h"
#include "Core/Memory/Tree/Layer.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLawData.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "GeneratedCode/tensor.h"

// A sycl-independent interface is required for interacting with the wp solver
// which, in its turn, is not supposed to know anything about SYCL
namespace seissol::dr::friction_law::gpu {
template <typename Cfg>
class FrictionSolverInterface : public FrictionSolverImpl<Cfg> {
  public:
  using real = Real<Cfg>;

  explicit FrictionSolverInterface(const FrictionLawParameters<Real<Cfg>>& drParameters)
      : FrictionSolverImpl<Cfg>(drParameters) {}
  ~FrictionSolverInterface() override = default;

  seissol::initializer::AllocationPlace allocationPlace() override {
    return seissol::initializer::AllocationPlace::Device;
  }

  void setClock(const double* deviceClock) override { clock_ = deviceClock; }

  static void copyStorageToLocal(FrictionLawData<Cfg>* data, DynamicRupture::Layer& layerData) {
    copyCommonStorageToLocal(data, layerData, seissol::initializer::AllocationPlace::Device);
  }

  protected:
  FrictionLawData<Cfg> dataHost_;

  // the time at the start of each step on the device; null to take the time passed to evaluate()
  const double* clock_{nullptr};
};
} // namespace seissol::dr::friction_law::gpu

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_FRICTIONLAWS_GPUIMPL_FRICTIONSOLVERINTERFACE_H_
