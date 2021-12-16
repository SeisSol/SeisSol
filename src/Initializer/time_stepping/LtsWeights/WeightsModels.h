#pragma once

#include "LtsWeights.h"
#include <functional>

namespace seissol::initializers::time_stepping {
class WeightModel {
  public:
  LtsWeights& ltsWeights;

  WeightModel(LtsWeights& ltsWeights) noexcept : ltsWeights(ltsWeights) {}
  virtual ~WeightModel() noexcept {}
};

class NodeWeightModel : public WeightModel {
  public:
  NodeWeightModel(LtsWeights& ltsWeights) noexcept : WeightModel(ltsWeights) {}
  virtual int evaluateNumberOfConstraints() const = 0;
  virtual void setVertexWeights() = 0;
  virtual void setAllowedImbalances() = 0;
  virtual ~NodeWeightModel() noexcept {}
};

class EdgeWeightModel : public WeightModel {
  public:
  EdgeWeightModel(LtsWeights& ltsWeights) noexcept : WeightModel(ltsWeights) {}
  virtual void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                         const std::vector<idx_t>&>& graph) = 0;

  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph,
                      std::function<int(idx_t, idx_t)>& factor);

  virtual ~EdgeWeightModel() noexcept {}
};

class ExponentialWeights : public NodeWeightModel {
  public:
  explicit ExponentialWeights(LtsWeights& ltsWeights) noexcept : NodeWeightModel(ltsWeights) {}
  ~ExponentialWeights() noexcept override final {}

  int evaluateNumberOfConstraints() const override final { return 1; }
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class ExponentialBalancedWeights : public NodeWeightModel {
  public:
  explicit ExponentialBalancedWeights(LtsWeights& ltsWeights) noexcept : NodeWeightModel(ltsWeights) {}
  ~ExponentialBalancedWeights() noexcept override final {}

  int evaluateNumberOfConstraints() const override final { return 2; }
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class EncodedBalancedWeights : public NodeWeightModel {
  public:
  explicit EncodedBalancedWeights(LtsWeights& ltsWeights) noexcept: NodeWeightModel(ltsWeights) {}
  ~EncodedBalancedWeights() noexcept override final {}

  int evaluateNumberOfConstraints() const override final;
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
  PenalizeBetweenClusters(LtsWeights& ltsWeights) noexcept : EdgeWeightModel(ltsWeights) {}
  ~PenalizeBetweenClusters() noexcept override final {}

  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph) override final;
};

class ApproximateCommunication : public EdgeWeightModel {
  public:
  ApproximateCommunication(LtsWeights& ltsWeights) noexcept: EdgeWeightModel(ltsWeights) {}
  ~ApproximateCommunication() noexcept override final {}

  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph) override final;
};
} // namespace seissol::initializers::time_stepping
