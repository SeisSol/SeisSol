#ifndef SEISSOL_LTSWEIGHTSMODELS_H
#define SEISSOL_LTSWEIGHTSMODELS_H

#include "LtsWeights.h"

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
  virtual void setEdgeWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph) = 0;

  void setEdgeWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph,
      std::function<int(idx_t, idx_t)>& factor);

  void setBalancedMessagingWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph);

  void setMessageCountWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph,
      std::function<int(idx_t, idx_t)>& factor);

  virtual ~EdgeWeightModel() noexcept {}
};

class ExponentialWeights final : public NodeWeightModel {
  public:
  explicit ExponentialWeights(LtsWeights& ltsWeights) noexcept : NodeWeightModel(ltsWeights) {}
  ~ExponentialWeights() noexcept override final {}

  int evaluateNumberOfConstraints() const override final { return 1; }
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class ExponentialBalancedWeights final : public NodeWeightModel {
  public:
  explicit ExponentialBalancedWeights(LtsWeights& ltsWeights) noexcept
      : NodeWeightModel(ltsWeights) {}
  ~ExponentialBalancedWeights() noexcept override final {}

  int evaluateNumberOfConstraints() const override final { return 2; }
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class ExponentialBalancedWeightsWithBalancedMessaging final : public NodeWeightModel {
  public:
  explicit ExponentialBalancedWeightsWithBalancedMessaging(LtsWeights& ltsWeights) noexcept
      : NodeWeightModel(ltsWeights) {}
  ~ExponentialBalancedWeightsWithBalancedMessaging() noexcept override final {}

  int evaluateNumberOfConstraints() const override final;
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class ExponentialBalancedWeightsWithMessageCount final : public NodeWeightModel {
  public:
  explicit ExponentialBalancedWeightsWithMessageCount(LtsWeights& ltsWeights) noexcept
      : NodeWeightModel(ltsWeights) {}
  ~ExponentialBalancedWeightsWithMessageCount() noexcept override final {}

  int evaluateNumberOfConstraints() const override final { return 3; }
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class EncodedBalancedWeights final : public NodeWeightModel {
  public:
  explicit EncodedBalancedWeights(LtsWeights& ltsWeights) noexcept : NodeWeightModel(ltsWeights) {}
  ~EncodedBalancedWeights() noexcept override final {}

  int evaluateNumberOfConstraints() const override final;
  void setVertexWeights() override final;
  void setAllowedImbalances() override final;
};

class Naive final : public EdgeWeightModel {
  public:
  Naive(LtsWeights& ltsWeights) : EdgeWeightModel(ltsWeights) {}
  ~Naive() override final {}

  void setEdgeWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph) override final;
};

class ApproximateCommunication final : public EdgeWeightModel {
  public:
  ApproximateCommunication(LtsWeights& ltsWeights) noexcept : EdgeWeightModel(ltsWeights) {}
  ~ApproximateCommunication() noexcept override final {}

  void setEdgeWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph) override final;
};

class ApproximateCommunicationWithBalancedMessaging final : public EdgeWeightModel {
  public:
  ApproximateCommunicationWithBalancedMessaging(LtsWeights& ltsWeights) noexcept
      : EdgeWeightModel(ltsWeights) {}
  ~ApproximateCommunicationWithBalancedMessaging() noexcept override final {}

  void setEdgeWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph) override final;
};

class ApproximateCommunicationWithMessageCount final : public EdgeWeightModel {
  public:
  ApproximateCommunicationWithMessageCount(LtsWeights& ltsWeights) noexcept
      : EdgeWeightModel(ltsWeights) {}
  ~ApproximateCommunicationWithMessageCount() noexcept override final {}

  void setEdgeWeights(
      std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
          graph) override final;
};

} // namespace seissol::initializers::time_stepping

#endif // SEISSOL_LTSWEIGHTSMODELS_H
