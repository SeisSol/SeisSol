// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_LTSPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_LTSPARAMETERS_H_

#include "Initializer/Clustering/ClusterCostModel.h"
#include "ParameterReader.h"

#include <cstdint>

namespace seissol::initializer::parameters {

enum class LtsWeightsTypes : int {
  ExponentialWeights = 0,
  ExponentialBalancedWeights,
  EncodedBalancedWeights,
  Count
};

struct VertexWeightParameters {
  std::uint64_t weightElement;
  std::uint64_t weightDynamicRupture;
  std::uint64_t weightFreeSurfaceWithGravity;
};

enum class AutoMergeCostBaseline {
  // Use cost without wiggle and cluster merge as baseline
  MaxWiggleFactor,
  // First find best wiggle factor (without merge) and use this as baseline
  BestWiggleFactor,
};

AutoMergeCostBaseline parseAutoMergeCostBaseline(std::string str);

enum class LtsClusteringSearch {
  // Sweep the wiggle factor on a grid, merging clusters from the top. The rate vector is
  // taken as given.
  Legacy,
  // Choose the rate vector as well, by shortest path on the divisor lattice.
  Lattice,
};

LtsClusteringSearch parseLtsClusteringSearch(std::string str);

class LtsParameters {
  private:
  std::vector<uint64_t> rate_;
  double wiggleFactorMinimum_{};
  double wiggleFactorStepsize_{};
  bool wiggleFactorEnforceMaximumDifference_{};
  int maxNumberOfClusters_{std::numeric_limits<int>::max() - 1};
  bool autoMergeClusters_{};
  double allowedPerformanceLossRatioAutoMerge_{};
  AutoMergeCostBaseline autoMergeCostBaseline_{AutoMergeCostBaseline::BestWiggleFactor};
  LtsWeightsTypes ltsWeightsType_{LtsWeightsTypes::ExponentialWeights};
  LtsClusteringSearch clusteringSearch_{LtsClusteringSearch::Legacy};
  double costUpdate_{1.0};
  double costLaunch_{0.0};
  double costFill_{0.0};
  std::uint64_t maxRate_{0};

  public:
  [[nodiscard]] std::vector<uint64_t> getRate() const;
  [[nodiscard]] bool isWiggleFactorUsed() const;
  [[nodiscard]] double getWiggleFactorMinimum() const;
  [[nodiscard]] double getWiggleFactorStepsize() const;
  [[nodiscard]] bool getWiggleFactorEnforceMaximumDifference() const;
  [[nodiscard]] int getMaxNumberOfClusters() const;
  [[nodiscard]] bool isAutoMergeUsed() const;
  [[nodiscard]] double getAllowedPerformanceLossRatioAutoMerge() const;
  [[nodiscard]] AutoMergeCostBaseline getAutoMergeCostBaseline() const;
  [[nodiscard]] LtsWeightsTypes getLtsWeightsType() const;
  [[nodiscard]] LtsClusteringSearch getClusteringSearch() const;
  /// Upper bound on a single ratio between adjacent clusters; 0 means unbounded. Only the
  /// lattice search can honor this, because the legacy search does not choose the rates.
  [[nodiscard]] std::uint64_t getMaxRate() const;
  [[nodiscard]] ClusterCostModel getClusterCostModel() const;

  LtsParameters() = default;

  LtsParameters(const std::vector<uint64_t>& rates,
                double wiggleFactorMinimum,
                double wiggleFactorStepsize,
                bool wigleFactorEnforceMaximumDifference,
                int maxNumberOfClusters,
                bool ltsAutoMergeClusters,
                double allowedPerformanceLossRatioAutoMerge,
                AutoMergeCostBaseline autoMergeCostBaseline,
                LtsWeightsTypes ltsWeightsType,
                LtsClusteringSearch clusteringSearch = LtsClusteringSearch::Legacy,
                double costUpdate = 1.0,
                double costLaunch = 0.0,
                double costFill = 0.0,
                std::uint64_t maxRate = 0);
};

struct TimeSteppingParameters {
  VertexWeightParameters vertexWeight{};
  double cfl{};
  double maxTimestepWidth{};
  double endTime{};
  LtsParameters lts;

  TimeSteppingParameters() = default;

  TimeSteppingParameters(VertexWeightParameters vertexWeight,
                         double cfl,
                         double maxTimestepWidth,
                         double endTime,
                         LtsParameters lts);
};

LtsParameters readLtsParameters(ParameterReader* baseReader);
TimeSteppingParameters readTimeSteppingParameters(ParameterReader* baseReader);

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_LTSPARAMETERS_H_
