// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_LTSPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_LTSPARAMETERS_H_

#include "ParameterReader.h"

namespace seissol::initializer::parameters {

enum class LtsWeightsTypes : int {
  ExponentialWeights = 0,
  ExponentialBalancedWeights,
  EncodedBalancedWeights,
  Count
};

struct VertexWeightParameters {
  int weightElement;
  int weightDynamicRupture;
  int weightFreeSurfaceWithGravity;
};

enum class AutoMergeCostBaseline {
  // Use cost without wiggle and cluster merge as baseline
  MaxWiggleFactor,
  // First find best wiggle factor (without merge) and use this as baseline
  BestWiggleFactor,
};

AutoMergeCostBaseline parseAutoMergeCostBaseline(std::string str);

class LtsParameters {
  private:
  unsigned int rate{};
  double wiggleFactorMinimum{};
  double wiggleFactorStepsize{};
  bool wiggleFactorEnforceMaximumDifference{};
  int maxNumberOfClusters = std::numeric_limits<int>::max() - 1;
  bool autoMergeClusters{};
  double allowedPerformanceLossRatioAutoMerge{};
  AutoMergeCostBaseline autoMergeCostBaseline = AutoMergeCostBaseline::BestWiggleFactor;
  LtsWeightsTypes ltsWeightsType;
  double finalWiggleFactor = 1.0;

  public:
  [[nodiscard]] unsigned int getRate() const;
  [[nodiscard]] bool isWiggleFactorUsed() const;
  [[nodiscard]] double getWiggleFactorMinimum() const;
  [[nodiscard]] double getWiggleFactorStepsize() const;
  [[nodiscard]] bool getWiggleFactorEnforceMaximumDifference() const;
  [[nodiscard]] int getMaxNumberOfClusters() const;
  [[nodiscard]] bool isAutoMergeUsed() const;
  [[nodiscard]] double getAllowedPerformanceLossRatioAutoMerge() const;
  [[nodiscard]] AutoMergeCostBaseline getAutoMergeCostBaseline() const;
  [[nodiscard]] double getWiggleFactor() const;
  [[nodiscard]] LtsWeightsTypes getLtsWeightsType() const;
  void setWiggleFactor(double factor);
  void setMaxNumberOfClusters(int numClusters);

  LtsParameters() = default;

  LtsParameters(unsigned int rate,
                double wiggleFactorMinimum,
                double wiggleFactorStepsize,
                bool wigleFactorEnforceMaximumDifference,
                int maxNumberOfClusters,
                bool ltsAutoMergeClusters,
                double allowedPerformanceLossRatioAutoMerge,
                AutoMergeCostBaseline autoMergeCostBaseline,
                LtsWeightsTypes ltsWeightsType);
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
