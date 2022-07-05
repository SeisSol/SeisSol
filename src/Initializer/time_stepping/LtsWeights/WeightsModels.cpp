#include <PUML/Downward.h>
#include <PUML/PUML.h>
#include <PUML/Upward.h>

#include "WeightsModels.h"

#include <Initializer/ParameterDB.h>
#include <Initializer/typedefs.hpp>
#include <Parallel/MPI.h>
#include <set>
#include <generated_code/init.h>
#include <stdio.h>
#include <stdlib.h>

namespace seissol::initializers::time_stepping {

void ExponentialWeights::setVertexWeights() {
  auto m_ncon = ltsWeights.getNcon();
  auto& m_details = ltsWeights.getDetails();
  auto& m_cellCosts = ltsWeights.getCellCosts();
  auto m_rate = ltsWeights.getRate();
  auto& m_clusterIds = ltsWeights.getClusterIds();
  auto& m_vertexWeights = ltsWeights.getVertexWeights();

  assert(m_ncon == 1 && "single constraint partitioning");
  int maxCluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];
  }
}

void ExponentialWeights::setAllowedImbalances() {
  auto m_ncon = ltsWeights.getNcon();
  auto& m_imbalances = ltsWeights.getImbalances();

  assert(m_ncon == 1 && "single constraint partitioning");
  m_imbalances.resize(m_ncon);

  constexpr double tinyLtsWeightImbalance{1.01};
  m_imbalances[0] = tinyLtsWeightImbalance;
}

void ExponentialBalancedWeights::setVertexWeights() {
  auto m_ncon = ltsWeights.getNcon();
  auto& m_details = ltsWeights.getDetails();
  auto& m_cellCosts = ltsWeights.getCellCosts();
  auto m_rate = ltsWeights.getRate();
  auto& m_clusterIds = ltsWeights.getClusterIds();
  auto& m_vertexWeights = ltsWeights.getVertexWeights();

  assert(m_ncon == 2 && "binary constraints partitioning");
  int maxCluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];

    constexpr int memoryWeight{1};
    m_vertexWeights[m_ncon * cell + 1] = memoryWeight;
  }
}

void ExponentialBalancedWeightsWithMessageCount::setVertexWeights() {
  auto m_ncon = ltsWeights.getNcon();
  auto& m_details = ltsWeights.getDetails();
  auto& m_cellCosts = ltsWeights.getCellCosts();
  auto m_rate = ltsWeights.getRate();
  auto& m_clusterIds = ltsWeights.getClusterIds();
  auto& m_vertexWeights = ltsWeights.getVertexWeights();

  int maxCluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];

    constexpr int memoryWeight{1};
    m_vertexWeights[m_ncon * cell + 1] = memoryWeight;

    // to be set later in edge weights
    m_vertexWeights[m_ncon * cell + 2] = 0;
  }
}

void ExponentialBalancedWeights::setAllowedImbalances() {
  auto m_ncon = ltsWeights.getNcon();
  auto& m_imbalances = ltsWeights.getImbalances();

  assert(m_ncon == 2 && "binary constaints partitioning");
  m_imbalances.resize(m_ncon);

  constexpr double tinyLtsWeightImbalance{1.01};
  m_imbalances[0] = tinyLtsWeightImbalance;

  constexpr double mediumLtsMemoryImbalance{1.05};
  m_imbalances[1] = mediumLtsMemoryImbalance;
}

void ExponentialBalancedWeightsWithMessageCount::setAllowedImbalances() {
  auto m_ncon = ltsWeights.getNcon();
  auto& m_imbalances = ltsWeights.getImbalances();

  assert(m_ncon == 3 && "ternary constaints partitioning");
  m_imbalances.resize(m_ncon);

  constexpr double tinyLtsWeightImbalance{1.01};
  m_imbalances[0] = tinyLtsWeightImbalance;

  constexpr double mediumLtsMemoryImbalance{1.05};
  m_imbalances[1] = mediumLtsMemoryImbalance;

  constexpr double mediumLtsMCImbalance{1.01};
  m_imbalances[2] = mediumLtsMCImbalance;
}

int ExponentialBalancedWeightsWithBalancedMessaging::evaluateNumberOfConstraints() const {
  const auto& m_details = ltsWeights.getDetails();
  const auto m_rate = ltsWeights.getRate();

  int maxCluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, m_rate);

  // 1 for memory, 1 for computation time, and 1 every time cluster (for the neighbors of the node
  // we can send messages to)
  //  + 2 + 1 -> extra 1 needed due to the smallest time cluster being cluster 0
  return maxCluster + 3;
}

void ExponentialBalancedWeightsWithBalancedMessaging::setAllowedImbalances() {
  const auto m_ncon = ltsWeights.getNcon();
  auto& m_imbalances = ltsWeights.getImbalances();

  m_imbalances.resize(m_ncon);

  constexpr double mediumLtsWeightImbalance{1.05};
  for (int i = 0; i < m_ncon; ++i) {
    m_imbalances[i] = mediumLtsWeightImbalance;
  }
}

void ExponentialBalancedWeightsWithBalancedMessaging::setVertexWeights() {
  auto m_ncon = ltsWeights.getNcon();
  auto& m_details = ltsWeights.getDetails();
  auto& m_cellCosts = ltsWeights.getCellCosts();
  auto m_rate = ltsWeights.getRate();
  auto& m_clusterIds = ltsWeights.getClusterIds();
  auto& m_vertexWeights = ltsWeights.getVertexWeights();

  int maxCluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];

    constexpr int memoryWeight{1};
    m_vertexWeights[m_ncon * cell + 1] = memoryWeight;

    // leave the rest 0 as they will be set in edgeweights
    for (int i = 2; i < m_ncon; i++) {
      m_vertexWeights[m_ncon * cell + i] = 0;
    }
  }
}

int EncodedBalancedWeights::evaluateNumberOfConstraints() const {
  const auto& m_details = ltsWeights.getDetails();
  const auto m_rate = ltsWeights.getRate();

  int maxCluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, m_rate);
  return maxCluster + 1;
}

void EncodedBalancedWeights::setVertexWeights() {
  const auto m_ncon = ltsWeights.getNcon();
  const auto& m_cellCosts = ltsWeights.getCellCosts();
  const auto& m_clusterIds = ltsWeights.getClusterIds();
  auto& m_vertexWeights = ltsWeights.getVertexWeights();

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    for (int i = 0; i < m_ncon; ++i) {
      m_vertexWeights[m_ncon * cell + i] = 0;
    }
    m_vertexWeights[m_ncon * cell + m_clusterIds[cell]] = m_cellCosts[cell];
  }
}

void EncodedBalancedWeights::setAllowedImbalances() {
  const auto m_ncon = ltsWeights.getNcon();
  auto& m_imbalances = ltsWeights.getImbalances();

  m_imbalances.resize(m_ncon);

  constexpr double mediumLtsWeightImbalance{1.05};
  for (int i = 0; i < m_ncon; ++i) {
    m_imbalances[i] = mediumLtsWeightImbalance;
  }
}

void Naive::setEdgeWeights(
    std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
        graph) {
  // Naive strategy means we will fill all of them with the same value

  std::vector<int>& m_edgeWeights = ltsWeights.getEdgeWeights();
  std::fill(m_edgeWeights.begin(), m_edgeWeights.end(), 100);
}

inline int calc_offset(int rank, int i) { return (rank > i) ? i : i - 1; }

void EdgeWeightModel::setEdgeWeights(
    std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
        graph,
    std::function<int(idx_t, idx_t)>& factor) {
  std::vector<int>& edgeWeights = ltsWeights.getEdgeWeights();

  ltsWeights.apply_constraints(graph, edgeWeights, factor, OffsetType::edgeWeight);
}

int ipow(int x, int y) {
  assert(y >= 0);

  if (y == 0) {
    return 1;
  }
  int result = x;
  while (--y) {
    result *= x;
  }
  return result;
}

void ApproximateCommunication::setEdgeWeights(
    std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
        graph) {
  unsigned rate = ltsWeights.getRate();
  const auto& m_details = ltsWeights.getDetails();
  int global_max_cluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, rate);

  std::function<int(idx_t, idx_t)> factor = [global_max_cluster, rate](idx_t cluster1,
                                                                       idx_t cluster2) {
    int rt = ipow(rate, global_max_cluster - cluster1);
    return rt;
  };

  EdgeWeightModel::setEdgeWeights(graph, factor);
}

void ApproximateCommunicationWithBalancedMessaging::setEdgeWeights(
    std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
        graph) {
  unsigned rate = ltsWeights.getRate();
  const auto& m_details = ltsWeights.getDetails();
  int global_max_cluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, rate);

  std::function<int(idx_t, idx_t)> factor = [global_max_cluster, rate](idx_t cluster1,
                                                                       idx_t cluster2) {
    int rt = ipow(rate, global_max_cluster - cluster1);
    return rt;
  };

  EdgeWeightModel::setEdgeWeights(graph, factor);
  EdgeWeightModel::setBalancedMessagingWeights(graph);
}

void ApproximateCommunicationWithMessageCount::setEdgeWeights(
    std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
        graph) {
  unsigned rate = ltsWeights.getRate();
  const auto& m_details = ltsWeights.getDetails();
  int global_max_cluster =
      ltsWeights.getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, rate);

  std::function<int(idx_t, idx_t)> factor = [global_max_cluster, rate](idx_t cluster1,
                                                                       idx_t cluster2) {
    int rt = ipow(rate, global_max_cluster - cluster1);
    return rt;
  };

  EdgeWeightModel::setEdgeWeights(graph, factor);
  EdgeWeightModel::setMessageCountWeights(graph, factor);
}

void EdgeWeightModel::setBalancedMessagingWeights(
    std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
        graph) {
  auto& m_vertexWeights = ltsWeights.getVertexWeights();

  std::function<int(idx_t, idx_t)> factor = [](idx_t cluster1, idx_t cluster2) { return 1; };

  ltsWeights.apply_constraints(graph, m_vertexWeights, factor, OffsetType::balancedMsg);
}

void EdgeWeightModel::setMessageCountWeights(
    std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&, const std::vector<idx_t>&>&
        graph,
    std::function<int(idx_t, idx_t)>& factor) {
  auto& m_vertexWeights = ltsWeights.getVertexWeights();

  ltsWeights.apply_constraints(graph, m_vertexWeights, factor, OffsetType::minMsg);
}

} // namespace seissol::initializers::time_stepping