#pragma once

#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <memory>
namespace seissol::solver::clustering::computation {

class GhostCluster : public AbstractTimeCluster {
  public:
  GhostCluster(double maxTimeStepSize,
               long timeStepRate,
               Executor executor,
               std::shared_ptr<communication::RecvNeighborCluster> neighbor)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate, executor), neighbor(neighbor) {}

  LayerType getLayerType() const override { return Ghost; }

  protected:
  bool emptyStep(ComputeStep step) const override { return step != ComputeStep::Predict; }

  void start() override {}

  void runCompute(ComputeStep step) override {
    neighbor->startFrom(streamRuntime);
    neighbor->stopTo(streamRuntime);
  }

  private:
  std::shared_ptr<communication::RecvNeighborCluster> neighbor;
};

} // namespace seissol::solver::clustering::computation
