#pragma once

#include "LtsWeights.h"
#include <functional>

namespace seissol::initializers::time_stepping {
class WeightModel {
  public:
  LtsWeights& ltsWeights;

  WeightModel(LtsWeights& ltsWeights) : ltsWeights(ltsWeights) {}
  virtual ~WeightModel() {}
};

class NodeWeightModel : public WeightModel {
  public:
  NodeWeightModel(LtsWeights& ltsWeights) : WeightModel(ltsWeights) {}
  virtual int evaluateNumberOfConstraints() = 0;
  virtual void setVertexWeights() = 0;
  virtual void setAllowedImbalances() = 0;
  virtual ~NodeWeightModel() {}
};

class EdgeWeightModel : public WeightModel {
  public:
  EdgeWeightModel(LtsWeights& ltsWeights) : WeightModel(ltsWeights) {}
  virtual void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                         const std::vector<idx_t>&>& graph) = 0;

  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph,
                      std::function<int(idx_t, idx_t)>& factor);

  virtual ~EdgeWeightModel() {}
};

class ExponentialWeights : public NodeWeightModel {
  public:
  explicit ExponentialWeights(LtsWeights& ltsWeights) : NodeWeightModel(ltsWeights) {}
  ~ExponentialWeights() override final {}

  int evaluateNumberOfConstraints() override final { return 1; }
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class ExponentialBalancedWeights : public NodeWeightModel {
  public:
  explicit ExponentialBalancedWeights(LtsWeights& ltsWeights) : NodeWeightModel(ltsWeights) {}
  ~ExponentialBalancedWeights() override final {}

  int evaluateNumberOfConstraints() override final { return 2; }
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class EncodedBalancedWeights : public NodeWeightModel {
  public:
  explicit EncodedBalancedWeights(LtsWeights& ltsWeights) : NodeWeightModel(ltsWeights) {}
  ~EncodedBalancedWeights() override final {}

  int evaluateNumberOfConstraints() override final;
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class Naive : public EdgeWeightModel {
  public:
  Naive(LtsWeights& ltsWeights) : EdgeWeightModel(ltsWeights) {}
  ~Naive() override final {}

  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph) override final;
};

class PenalizeBetweenClusters : public EdgeWeightModel {
  public:
  PenalizeBetweenClusters(LtsWeights& ltsWeights) : EdgeWeightModel(ltsWeights) {}
  ~PenalizeBetweenClusters() override final {}

  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph) override final;
};

class ApproximateCommunication : public EdgeWeightModel {
  public:
  ApproximateCommunication(LtsWeights& ltsWeights) : EdgeWeightModel(ltsWeights) {}
  ~ApproximateCommunication() override final {}

  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph) override final;
};
} // namespace seissol::initializers::time_stepping
