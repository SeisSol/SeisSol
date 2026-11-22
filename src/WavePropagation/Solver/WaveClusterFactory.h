// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_WAVECLUSTERFACTORY_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_WAVECLUSTERFACTORY_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Executor.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Compute/CellCluster.h"
#include "Core/TimeStepping/Compute/CellClusterFactory.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "WavePropagation/InitialField.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace seissol::solver {

/// Makes the clusters of the cells whose waves the discontinuous Galerkin scheme propagates.
class WaveClusterFactory : public CellClusterFactory {
  public:
  /// With the global data `globalData` of the configurations `configs`, the initial conditions of
  /// the cells of each configuration (by its id), the gravitational acceleration of the free
  /// surface, and where a cluster of the given number of cells computes (`executionPlace`).
  WaveClusterFactory(
      RunContext& run,
      GlobalDataStore& globalData,
      std::vector<ConfigId> configs,
      const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& initialConditions,
      double gravitationalAcceleration,
      std::function<Executor(std::size_t)> executionPlace);

  std::unique_ptr<CellCluster> makeCluster(LTS::Layer& layer,
                                           const initializer::ClusterLayout& clusterLayout,
                                           bool printProgress,
                                           const ClusterMonitoring& monitoring) override;

  private:
  RunContext& run_;
  GlobalDataStore& globalData_;
  std::vector<ConfigId> configs_;
  const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& initialConditions_;
  double gravitationalAcceleration_;
  std::function<Executor(std::size_t)> executionPlace_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_WAVECLUSTERFACTORY_H_
