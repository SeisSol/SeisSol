#include <PUML/PUML.h>
#include <PUML/Downward.h>
#include <PUML/Upward.h>

#include "WeightsModels.h"

#include <Initializer/typedefs.hpp>
#include <Initializer/ParameterDB.h>
#include <Parallel/MPI.h>

#include <generated_code/init.h>


namespace seissol::initializers::time_stepping {

void ExponentialWeights::setVertexWeights() {
  assert(m_ncon == 1 && "single constraint partitioning");
  int maxCluster = getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, wiggleFactor, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];
  }
}

void ExponentialWeights::setAllowedImbalances() {
  assert(m_ncon == 1 && "single constraint partitioning");
  m_imbalances.resize(m_ncon);

  constexpr double tinyLtsWeightImbalance{1.01};
  m_imbalances[0] = tinyLtsWeightImbalance;
}


void ExponentialBalancedWeights::setVertexWeights() {
  assert(m_ncon == 2 && "binary constaints partitioning");
  int maxCluster = getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, wiggleFactor, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];

    constexpr int memoryWeight{1};
    m_vertexWeights[m_ncon * cell + 1] = memoryWeight;
  }
}


void ExponentialBalancedWeights::setAllowedImbalances() {
  assert(m_ncon == 2 && "binary constaints partitioning");
  m_imbalances.resize(m_ncon);

  constexpr double tinyLtsWeightImbalance{1.01};
  m_imbalances[0] = tinyLtsWeightImbalance;

  constexpr double mediumLtsMemoryImbalance{1.05};
  m_imbalances[1] = mediumLtsMemoryImbalance;
}


int EncodedBalancedWeights::evaluateNumberOfConstraints() {
  int maxCluster = getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, wiggleFactor, m_rate);
  return maxCluster + 1;
}


void EncodedBalancedWeights::setVertexWeights() {
  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    for (int i = 0; i < m_ncon; ++i) {
      m_vertexWeights[m_ncon * cell + i] = 0;
    }
    m_vertexWeights[m_ncon * cell + m_clusterIds[cell]] = m_cellCosts[cell];
  }
}


void EncodedBalancedWeights::setAllowedImbalances() {
  m_imbalances.resize(m_ncon);

  constexpr double mediumLtsWeightImbalance{1.05};
  for (int i = 0; i < m_ncon; ++i) {
    m_imbalances[i] = mediumLtsWeightImbalance;
  }
}
}
