// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WeightsModels.h"

#include "generated_code/init.h"
#include <Initializer/TimeStepping/LtsWeights/LtsWeights.h>
#include <cassert>

namespace seissol::initializer::time_stepping {

void ExponentialWeights::setVertexWeights() {
  assert(m_ncon == 1 && "single constraint partitioning");
  const int maxCluster =
      getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, wiggleFactor, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    const int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];
  }
}

void ExponentialWeights::setAllowedImbalances() {
  assert(m_ncon == 1 && "single constraint partitioning");
  m_imbalances.resize(m_ncon);

  constexpr double TinyLtsWeightImbalance{1.01};
  m_imbalances[0] = TinyLtsWeightImbalance;
}

void ExponentialBalancedWeights::setVertexWeights() {
  assert(m_ncon == 2 && "binary constaints partitioning");
  const int maxCluster =
      getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, wiggleFactor, m_rate);

  for (unsigned cell = 0; cell < m_cellCosts.size(); ++cell) {
    const int factor = LtsWeights::ipow(m_rate, maxCluster - m_clusterIds[cell]);
    m_vertexWeights[m_ncon * cell] = factor * m_cellCosts[cell];

    constexpr int MemoryWeight{1};
    m_vertexWeights[m_ncon * cell + 1] = MemoryWeight;
  }
}

void ExponentialBalancedWeights::setAllowedImbalances() {
  assert(m_ncon == 2 && "binary constaints partitioning");
  m_imbalances.resize(m_ncon);

  constexpr double TinyLtsWeightImbalance{1.01};
  m_imbalances[0] = TinyLtsWeightImbalance;

  constexpr double MediumLtsMemoryImbalance{1.05};
  m_imbalances[1] = MediumLtsMemoryImbalance;
}

int EncodedBalancedWeights::evaluateNumberOfConstraints() {
  const int maxCluster =
      getCluster(m_details.globalMaxTimeStep, m_details.globalMinTimeStep, wiggleFactor, m_rate);
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

  constexpr double MediumLtsWeightImbalance{1.05};
  for (int i = 0; i < m_ncon; ++i) {
    m_imbalances[i] = MediumLtsWeightImbalance;
  }
}
} // namespace seissol::initializer::time_stepping
