// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <Solver/Clustering/Computation/TimeCluster.h>
namespace seissol::solver::clustering::computation {

class CopyCluster : public TimeCluster {
  public:
  CopyCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              bool usePlasticity,
              double maxTimeStepSize,
              long timeStepRate,
              bool printProgress,
              CompoundGlobalData globalData,
              seissol::initializer::Layer* clusterData,
              seissol::initializer::LTS* lts,
              seissol::SeisSol& seissolInstance,
              LoopStatistics* loopStatistics,
              ActorStateStatistics* actorStateStatistics,
              std::shared_ptr<communication::SendNeighborCluster> neighbor,
              const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor)
      : TimeCluster(clusterId,
                    globalClusterId,
                    profilingId,
                    usePlasticity,
                    maxTimeStepSize,
                    timeStepRate,
                    printProgress,
                    globalData,
                    clusterData,
                    lts,
                    seissolInstance,
                    loopStatistics,
                    actorStateStatistics,
                    cpuExecutor),
        neighbor(neighbor) {}
  LayerType getLayerType() const override { return Copy; }

  protected:
  void start() override { dataSent = false; }

  void runCompute(ComputeStep step) override {
    if (step == ComputeStep::Predict && dataSent) {
      neighbor->stopTo(streamRuntime);
    }

    TimeCluster::runCompute(step);

    if (step == ComputeStep::Predict) {
      neighbor->startFrom(streamRuntime);
      dataSent = true;
    }
  }

  void reset() override {
    if (dataSent) {
      neighbor->stopTo(streamRuntime);
      dataSent = false;
    }
  }

  private:
  std::shared_ptr<communication::SendNeighborCluster> neighbor;
  bool dataSent{false};
};

} // namespace seissol::solver::clustering::computation
