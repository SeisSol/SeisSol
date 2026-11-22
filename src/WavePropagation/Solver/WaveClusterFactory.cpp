// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WaveClusterFactory.h"

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Executor.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Compute/CellCluster.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "WavePropagation/InitialField.h"
#include "WavePropagation/Solver/WaveCluster.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace seissol::solver {

WaveClusterFactory::WaveClusterFactory(
    RunContext& run,
    GlobalDataStore& globalData,
    std::vector<ConfigId> configs,
    const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& initialConditions,
    double gravitationalAcceleration,
    std::function<Executor(std::size_t)> executionPlace)
    : run_(run), globalData_(globalData), configs_(std::move(configs)),
      initialConditions_(initialConditions), gravitationalAcceleration_(gravitationalAcceleration),
      executionPlace_(std::move(executionPlace)) {}

std::unique_ptr<CellCluster>
    WaveClusterFactory::makeCluster(LTS::Layer& layer,
                                    const initializer::ClusterLayout& clusterLayout,
                                    bool printProgress,
                                    const ClusterMonitoring& monitoring) {
  const auto clusterId = layer.getIdentifier().lts;
  return dispatchConfig(
      layer.getIdentifier().config, [&](auto cfg) -> std::unique_ptr<CellCluster> {
        using Cfg = decltype(cfg);
        const typename WaveCluster<Cfg>::Setup setup{
            .run = run_,
            .globalData = globalData_,
            .configs = configs_,
            .initialConditions = &initialConditions_.at(configIdOf<Cfg>()),
            .gravitationalAcceleration = gravitationalAcceleration_};
        return std::make_unique<WaveCluster<Cfg>>(clusterId,
                                                  clusterId,
                                                  layer.id(),
                                                  clusterLayout.timestepRate(clusterId),
                                                  clusterLayout.clusterRate(clusterId),
                                                  executionPlace_(layer.size()),
                                                  printProgress,
                                                  &layer,
                                                  setup,
                                                  monitoring.loopStatistics,
                                                  monitoring.deviceTimeline,
                                                  monitoring.actorStateStatistics);
      });
}

} // namespace seissol::solver
