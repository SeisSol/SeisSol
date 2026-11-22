// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRuptureClusterFactory.h"

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Executor.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/MPI.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Compute/FaceCluster.h"
#include "Core/TimeStepping/Compute/FaceClusterFactory.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "Features/DynamicRupture/Output/OutputManager.h"
#include "Features/DynamicRupture/Solver/DynamicRuptureCluster.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"

#include <cassert>
#include <cstddef>
#include <limits>
#include <memory>
#include <mpi.h>
#include <vector>

namespace seissol::dr {

DynamicRuptureClusterFactory::DynamicRuptureClusterFactory(
    RunContext& run,
    DynamicRupture::Storage& storage,
    GlobalDataStore& globalData,
    const friction_law::FrictionSolverFactory& frictionLaw,
    const friction_law::FrictionSolverFactory& frictionLawDevice,
    output::OutputManager* outputManager,
    double endTime)
    : run_(run), storage_(storage), globalData_(globalData), frictionLaw_(frictionLaw),
      frictionLawDevice_(frictionLawDevice), outputManager_(outputManager), endTime_(endTime) {}

std::size_t DynamicRuptureClusterFactory::faceCount(std::size_t layerId) const {
  return storage_.layer(layerId).size();
}

void DynamicRuptureClusterFactory::prepare(const initializer::ClusterLayout& clusterLayout) {
  std::vector<std::size_t> facesPerCluster(clusterLayout.globalClusterCount);
  for (const auto& layer : storage_.leaves()) {
    facesPerCluster[layer.getIdentifier().lts] += layer.size();
  }

  std::size_t outputCluster = std::numeric_limits<std::size_t>::max();
  for (std::size_t clusterId = 0; clusterId < facesPerCluster.size(); ++clusterId) {
    if (facesPerCluster[clusterId] > 0) {
      outputCluster = clusterId;
      break;
    }
  }
  outputCluster = Mpi::mpi.allreduce(outputCluster, MPI_MIN);

  const auto outputTimestep = outputCluster == std::numeric_limits<std::size_t>::max()
                                  ? std::numeric_limits<double>::infinity()
                                  : clusterLayout.timestepRate(outputCluster);

  // the on-fault receivers sample at the multiples of their interval, counted from time zero
  pickpointSampling_.emplace(outputManager_->pickpointInterval(outputTimestep), endTime_);
}

std::unique_ptr<solver::FaceCluster>
    DynamicRuptureClusterFactory::makeCluster(std::size_t layerId,
                                              const initializer::ClusterLayout& clusterLayout,
                                              Executor executor,
                                              const solver::ClusterMonitoring& monitoring) {
  assert(pickpointSampling_.has_value());
  auto& layer = storage_.layer(layerId);
  const auto clusterId = layer.getIdentifier().lts;
  return dispatchConfig(layer.getIdentifier().config,
                        [&](auto cfg) -> std::unique_ptr<solver::FaceCluster> {
                          using Cfg = decltype(cfg);
                          return std::make_unique<solver::DynamicRuptureCluster<Cfg>>(
                              clusterLayout.timestepRate(clusterId),
                              clusterLayout.clusterRate(clusterId),
                              executor,
                              layerId,
                              pickpointSampling_.value(),
                              globalData_.get<Cfg>(),
                              &layer,
                              frictionLaw_,
                              frictionLawDevice_,
                              outputManager_,
                              run_,
                              monitoring.loopStatistics,
                              monitoring.deviceTimeline,
                              monitoring.actorStateStatistics);
                        });
}

} // namespace seissol::dr
