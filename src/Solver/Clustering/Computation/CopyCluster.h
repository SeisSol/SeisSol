// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMPUTATION_COPYCLUSTER_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMPUTATION_COPYCLUSTER_H_

#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <Solver/Clustering/Computation/TimeCluster.h>

#include <utility>
namespace seissol::solver::clustering::computation {

class CopyCluster : public TimeCluster {
  public:
  CopyCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              bool usePlasticity,
              double maxTimeStepSize,
              long timeStepRate,
              CompoundGlobalData globalData,
              seissol::initializer::Layer* clusterData,
              seissol::initializer::LTS* lts,
              seissol::SeisSol& seissolInstance,
              LoopStatistics* loopStatistics,
              ActorStateStatistics* actorStateStatistics,
              std::shared_ptr<communication::SendNeighborCluster> neighbor,
              const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
              double priority)
      : TimeCluster(clusterId,
                    globalClusterId,
                    profilingId,
                    usePlasticity,
                    maxTimeStepSize,
                    timeStepRate,
                    false,
                    globalData,
                    clusterData,
                    lts,
                    seissolInstance,
                    loopStatistics,
                    actorStateStatistics,
                    cpuExecutor,
                    priority),
        neighbor(std::move(neighbor)) {}
  LayerType getLayerType() const override { return Copy; }

  protected:
  void start() override {
    TimeCluster::start();
    dataSent = false;
  }

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
    TimeCluster::reset();
    if (dataSent) {
      neighbor->stopTo(streamRuntime);
      dataSent = false;
    }
  }

  std::string description() const override { return "copy-cell"; }

  private:
  std::shared_ptr<communication::SendNeighborCluster> neighbor;
  bool dataSent{false};
};

} // namespace seissol::solver::clustering::computation
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMPUTATION_COPYCLUSTER_H_
