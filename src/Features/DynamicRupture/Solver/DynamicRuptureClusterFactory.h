// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURECLUSTERFACTORY_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURECLUSTERFACTORY_H_

#include "Core/Common/Executor.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Compute/FaceCluster.h"
#include "Core/TimeStepping/Compute/FaceClusterFactory.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "Features/DynamicRupture/Output/OutputManager.h"
#include "Features/DynamicRupture/Solver/DynamicRuptureCluster.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::dr {

/**
 * Makes the clusters of the dynamic rupture faces: one per layer of the faces that has any, with
 * the friction law of the run, and the on-fault receivers of its output.
 */
class DynamicRuptureClusterFactory : public solver::FaceClusterFactory {
  public:
  /// the faces in `storage`, with the global data `globalData`, the friction laws for the host and
  /// the device, and the output `outputManager`; the run ends at `endTime`
  DynamicRuptureClusterFactory(RunContext& run,
                               DynamicRupture::Storage& storage,
                               GlobalDataStore& globalData,
                               const friction_law::FrictionSolverFactory& frictionLaw,
                               const friction_law::FrictionSolverFactory& frictionLawDevice,
                               output::OutputManager* outputManager,
                               double endTime);

  [[nodiscard]] std::size_t faceCount(std::size_t layerId) const override;

  /// decides when the receivers of each output on the fault take their samples: at the multiples
  /// of its interval, which may follow the time step of the first time cluster with faces on any
  /// rank
  void prepare(const initializer::ClusterLayout& clusterLayout) override;

  std::unique_ptr<solver::FaceCluster>
      makeCluster(std::size_t layerId,
                  const initializer::ClusterLayout& clusterLayout,
                  Executor executor,
                  const solver::ClusterMonitoring& monitoring) override;

  private:
  RunContext& run_;
  DynamicRupture::Storage& storage_;
  GlobalDataStore& globalData_;
  const friction_law::FrictionSolverFactory& frictionLaw_;
  const friction_law::FrictionSolverFactory& frictionLawDevice_;
  output::OutputManager* outputManager_;
  double endTime_;

  //! the outputs of receivers on the fault, and when they take their samples
  std::vector<solver::PickpointPlan> pickpoints_;
};

} // namespace seissol::dr

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURECLUSTERFACTORY_H_
