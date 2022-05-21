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

#include <string>
#include <vector>
#include <limits>

// needed for the idx_t type declaration
#include <parmetis.h>
// for std pair
#include <utility>


#ifndef PUML_PUML_H
namespace PUML { class TETPUML; }
#endif // PUML_PUML_H


namespace seissol::initializers::time_stepping {
struct LtsWeightsConfig {
  std::string velocityModel{};
  unsigned rate{};
  int vertexWeightElement{};
  int vertexWeightDynamicRupture{};
  int vertexWeightFreeSurfaceWithGravity{};
};


class NodeWeightModel;
class EdgeWeightModel;

class LtsWeights {
private:
  struct GlobalTimeStepDetails {
    double globalMinTimeStep{};
    double globalMaxTimeStep{};
    std::vector<double> timeSteps{};
  } m_details;

  std::string m_velocityModel{};
  unsigned m_rate{};
  std::vector<int> m_vertexWeights{};
  std::vector<int> m_edgeWeights{};
  std::vector<double> m_imbalances{};
  std::vector<int> m_cellCosts{};
  int m_vertexWeightElement{};
  int m_vertexWeightDynamicRupture{};
  int m_vertexWeightFreeSurfaceWithGravity{};
  int m_ncon{std::numeric_limits<int>::infinity()};
  const PUML::TETPUML* m_mesh{nullptr};
  std::vector<int> m_clusterIds{};

  NodeWeightModel* m_nodeWeightModel;
  EdgeWeightModel* m_edgeWeightModel;

  public:
  LtsWeights(const LtsWeightsConfig& config) noexcept;

  ~LtsWeights() noexcept;
  void computeNodeWeights(PUML::TETPUML const& mesh, double maximumAllowedTimeStep);
  void computeEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                     const std::vector<idx_t>&>& graph);

  const int* vertexWeights() const;
  const double* imbalances() const;
  int nWeightsPerVertex() const;

  const int* edgeWeights() const;
  int edgeCount() const;

  GlobalTimeStepDetails collectGlobalTimeStepDetails(double maximumAllowedTimeStep);
  void computeMaxTimesteps(std::vector<double> const& pWaveVel, std::vector<double>& timeSteps,
                           double maximumAllowedTimeStep);
  int getCluster(double timestep, double globalMinTimestep, unsigned rate);
  int getBoundaryCondition(int const* boundaryCond, unsigned cell, unsigned face);
  std::vector<int> computeClusterIds();
  int enforceMaximumDifference();
  int enforceMaximumDifferenceLocal(int maxDifference = 1);
  std::vector<int> computeCostsPerTimestep();

  static int ipow(int x, int y);

  void setVertexWeights();
  void setAllowedImbalances();
  int evaluateNumberOfConstraints() const;
  void setEdgeWeights(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                 const std::vector<idx_t>&>& graph);

  const GlobalTimeStepDetails& getDetails() const;
  const std::string& getVelocityModel() const;
  unsigned getRate() const;
  std::vector<int>& getVertexWeights();
  std::vector<int>& getEdgeWeights();
  std::vector<double>& getImbalances();
  const std::vector<int>& getCellCosts() const;
  int getVertexWeightElement() const;
  int getVertexWeightDynamicRupture() const;
  int getVertexWeightFreeSurfaceWithGravity() const;
  int getNcon() const;
  const PUML::TETPUML* getMesh() const;
  const std::vector<int>& getClusterIds() const;

  void addWeightModels(NodeWeightModel* nwm, EdgeWeightModel* ewm);

  static int find_rank(const std::vector<idx_t> &vrtxdist, idx_t elemId);
  
  std::vector<std::unordered_map<idx_t, int>>
  exchangeGhostLayer(std::tuple<const std::vector<idx_t>&, const std::vector<idx_t>&,
                                const std::vector<idx_t>&>& graph);
};
}

#endif
