// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_LTSWEIGHTS_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_LTSWEIGHTS_H_

#include "Geometry/PUMLReader.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"

#include <limits>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace seissol {
class SeisSol;
namespace initializer::time_stepping {
struct LtsWeightsConfig {
  seissol::initializer::parameters::BoundaryFormat boundaryFormat;
  std::vector<uint64_t> rate;
  int vertexWeightElement{};
  int vertexWeightDynamicRupture{};
  int vertexWeightFreeSurfaceWithGravity{};
};

double computeLocalCostOfClustering(const std::vector<int>& clusterIds,
                                    const std::vector<int>& cellCosts,
                                    const std::vector<uint64_t>& rate,
                                    double wiggleFactor,
                                    double minimalTimestep);

double computeGlobalCostOfClustering(const std::vector<int>& clusterIds,
                                     const std::vector<int>& cellCosts,
                                     const std::vector<uint64_t>& rate,
                                     double wiggleFactor,
                                     double minimalTimestep,
                                     MPI_Comm comm);

std::vector<int> enforceMaxClusterId(const std::vector<int>& clusterIds, int maxClusterId);

int computeMaxClusterIdAfterAutoMerge(const std::vector<int>& clusterIds,
                                      const std::vector<int>& cellCosts,
                                      const std::vector<uint64_t>& rate,
                                      double maximalAdmissibleCost,
                                      double wiggleFactor,
                                      double minimalTimestep);

std::uint64_t ratepow(const std::vector<std::uint64_t>& rate, std::uint64_t a, std::uint64_t b);

class LtsWeights {
  public:
  LtsWeights(const LtsWeightsConfig& config, seissol::SeisSol& seissolInstance);

  virtual ~LtsWeights() = default;
  void computeWeights(const seissol::geometry::PumlMesh& meshTopology,
                      const seissol::geometry::PumlMesh& meshGeometry);

  [[nodiscard]] const int* vertexWeights() const;
  [[nodiscard]] const double* imbalances() const;
  [[nodiscard]] const std::vector<int>& clusterIds() const;
  [[nodiscard]] const std::vector<double>& timesteps() const;
  [[nodiscard]] int nWeightsPerVertex() const;
  [[nodiscard]] double getWiggleFactor() const;

  private:
  seissol::SeisSol& seissolInstance_;

  protected:
  seissol::initializer::GlobalTimestep details_;

  seissol::initializer::GlobalTimestep collectGlobalTimeStepDetails();
  std::uint64_t getCluster(double timestep,
                           double globalMinTimestep,
                           double wiggleFactor,
                           const std::vector<uint64_t>& rate);
  FaceType getBoundaryCondition(const void* boundaryCond, size_t cell, unsigned face);
  std::vector<int> computeClusterIds(double curWiggleFactor);
  // returns number of reductions for maximum difference
  int computeClusterIdsAndEnforceMaximumDifferenceCached(double curWiggleFactor);
  int enforceMaximumDifference();
  int enforceMaximumDifferenceLocal(int maxDifference = 1);
  std::vector<int> computeCostsPerTimestep();

  virtual void setVertexWeights() = 0;
  virtual void setAllowedImbalances() = 0;
  virtual int evaluateNumberOfConstraints() = 0;

  std::vector<uint64_t> rate_;
  std::vector<int> vertexWeights_;
  std::vector<double> imbalances_;
  std::vector<int> cellCosts_;
  int vertexWeightElement_{};
  int vertexWeightDynamicRupture_{};
  int vertexWeightFreeSurfaceWithGravity_{};
  int ncon_{std::numeric_limits<int>::infinity()};
  const geometry::PumlMesh* meshTopology_{nullptr};
  const geometry::PumlMesh* meshGeometry_{nullptr};
  std::vector<int> clusterIds_;
  double wiggleFactor_ = 1.0;
  std::map<double, decltype(clusterIds_), std::greater<>>
      clusteringCache_; // Maps wiggle factor to clustering
  seissol::initializer::parameters::BoundaryFormat boundaryFormat_;
  struct ComputeWiggleFactorResult {
    int maxClusterId;
    double wiggleFactor;
    double cost;
  };
  ComputeWiggleFactorResult computeBestWiggleFactor(std::optional<double> baselineCost,
                                                    bool isAutoMergeUsed);
  void prepareDifferenceEnforcement();

  std::vector<std::pair<int, std::vector<std::size_t>>> rankToSharedFaces_;
  std::unordered_map<std::size_t, std::size_t> localFaceIdToLocalCellId_;
  std::unordered_map<std::size_t, std::pair<std::size_t, std::size_t>> sharedFaceToExchangeId_;
  std::vector<std::size_t> boundaryCells_;
};
} // namespace initializer::time_stepping
} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_LTSWEIGHTS_H_
