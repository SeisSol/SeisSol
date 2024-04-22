/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT tum.de, https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017 - 2020, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Class for calculating weights for load balancing
 **/

#ifndef INITIALIZER_TIMESTEPPING_LTSWEIGHTS_H_
#define INITIALIZER_TIMESTEPPING_LTSWEIGHTS_H_

#include <Geometry/PUMLReader.h>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "Initializer/time_stepping/GlobalTimestep.hpp"

#ifndef PUML_PUML_H
namespace PUML { class TETPUML; }
#endif // PUML_PUML_H


namespace seissol {
  class SeisSol;
  namespace initializer::time_stepping {
struct LtsWeightsConfig {
  seissol::initializer::parameters::BoundaryFormat boundaryFormat;
  std::string velocityModel{};
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

  const int *vertexWeights() const;
  const double *imbalances() const;
  int nWeightsPerVertex() const;

private:
  seissol::SeisSol& seissolInstance;
protected:
  seissol::initializer::GlobalTimestep m_details;

  seissol::initializer::GlobalTimestep collectGlobalTimeStepDetails(double maximumAllowedTimeStep);
  void computeMaxTimesteps(std::vector<double> const &pWaveVel, std::vector<double> &timeSteps, double maximumAllowedTimeStep);
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

  seissol::initializer::parameters::BoundaryFormat boundaryFormat;
  std::string m_velocityModel{};
  unsigned m_rate{};
  std::vector<int> m_vertexWeights{};
  std::vector<double> m_imbalances{};
  std::vector<int> m_cellCosts{};
  int m_vertexWeightElement{};
  int m_vertexWeightDynamicRupture{};
  int m_vertexWeightFreeSurfaceWithGravity{};
  int m_ncon{std::numeric_limits<int>::infinity()};
  const PUML::TETPUML * m_mesh{nullptr};
  std::vector<int> m_clusterIds{};
  double wiggleFactor = 1.0;
  std::map<double, decltype(m_clusterIds)> clusteringCache; // Maps wiggle factor to clustering
  struct ComputeWiggleFactorResult {
    int maxClusterId;
    double wiggleFactor;
    double cost;
  };
  ComputeWiggleFactorResult computeBestWiggleFactor(std::optional<double> baselineCost,
                                                    bool isAutoMergeUsed);
};
}
}

#endif
