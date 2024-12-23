// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_LTSWEIGHTS_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_LTSWEIGHTS_H_

#include "Geometry/PUMLReader.h"
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "Initializer/TimeStepping/GlobalTimestep.h"

#ifndef PUML_PUML_H
namespace PUML {
class TETPUML;
}
#endif // PUML_PUML_H

namespace seissol {
class SeisSol;
namespace initializer::time_stepping {
struct LtsWeightsConfig {
  seissol::initializer::parameters::BoundaryFormat boundaryFormat;
  std::string velocityModel;
  unsigned rate{};
  int vertexWeightElement{};
  int vertexWeightDynamicRupture{};
  int vertexWeightFreeSurfaceWithGravity{};
};

double computeLocalCostOfClustering(const std::vector<int>& clusterIds,
                                    const std::vector<int>& cellCosts,
                                    unsigned int rate,
                                    double wiggleFactor,
                                    double minimalTimestep);

double computeGlobalCostOfClustering(const std::vector<int>& clusterIds,
                                     const std::vector<int>& cellCosts,
                                     unsigned int rate,
                                     double wiggleFactor,
                                     double minimalTimestep,
                                     MPI_Comm comm);

std::vector<int> enforceMaxClusterId(const std::vector<int>& clusterIds, int maxClusterId);

int computeMaxClusterIdAfterAutoMerge(const std::vector<int>& clusterIds,
                                      const std::vector<int>& cellCosts,
                                      unsigned int rate,
                                      double maximalAdmissibleCost,
                                      double wiggleFactor,
                                      double minimalTimestep);

class LtsWeights {
  public:
  LtsWeights(const LtsWeightsConfig& config, seissol::SeisSol& seissolInstance);

  virtual ~LtsWeights() = default;
  void computeWeights(PUML::TETPUML const& mesh, double maximumAllowedTimeStep);

  [[nodiscard]] const int* vertexWeights() const;
  [[nodiscard]] const double* imbalances() const;
  [[nodiscard]] int nWeightsPerVertex() const;

  private:
  seissol::SeisSol& seissolInstance;

  protected:
  seissol::initializer::GlobalTimestep m_details;

  seissol::initializer::GlobalTimestep collectGlobalTimeStepDetails(double maximumAllowedTimeStep);
  void computeMaxTimesteps(const std::vector<double>& pWaveVel,
                           std::vector<double>& timeSteps,
                           double maximumAllowedTimeStep);
  int getCluster(double timestep, double globalMinTimestep, double wiggleFactor, unsigned rate);
  int getBoundaryCondition(const void* boundaryCond, size_t cell, unsigned face);
  std::vector<int> computeClusterIds(double curWiggleFactor);
  // returns number of reductions for maximum difference
  int computeClusterIdsAndEnforceMaximumDifferenceCached(double curWiggleFactor);
  int enforceMaximumDifference();
  int enforceMaximumDifferenceLocal(int maxDifference = 1);
  std::vector<int> computeCostsPerTimestep();

  static int ipow(int x, int y);

  virtual void setVertexWeights() = 0;
  virtual void setAllowedImbalances() = 0;
  virtual int evaluateNumberOfConstraints() = 0;

  std::string m_velocityModel;
  unsigned m_rate{};
  std::vector<int> m_vertexWeights;
  std::vector<double> m_imbalances;
  std::vector<int> m_cellCosts;
  int m_vertexWeightElement{};
  int m_vertexWeightDynamicRupture{};
  int m_vertexWeightFreeSurfaceWithGravity{};
  int m_ncon{std::numeric_limits<int>::infinity()};
  const PUML::TETPUML* m_mesh{nullptr};
  std::vector<int> m_clusterIds;
  double wiggleFactor = 1.0;
  std::map<double, decltype(m_clusterIds)> clusteringCache; // Maps wiggle factor to clustering
  seissol::initializer::parameters::BoundaryFormat boundaryFormat;
  struct ComputeWiggleFactorResult {
    int maxClusterId;
    double wiggleFactor;
    double cost;
  };
  ComputeWiggleFactorResult computeBestWiggleFactor(std::optional<double> baselineCost,
                                                    bool isAutoMergeUsed);
};
} // namespace initializer::time_stepping
} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_LTSWEIGHTS_H_
