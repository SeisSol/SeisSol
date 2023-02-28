#ifndef SEISSOL_LTSCONFIGURATION_H
#define SEISSOL_LTSCONFIGURATION_H

#include <memory>
#include <yaml-cpp/yaml.h>

namespace seissol::initializers::time_stepping {

enum class AutoMergeCostBaseline {
  // Use cost without wiggle and cluster merge as baseline
  MaxWiggleFactor,
  // First find best wiggle factor (without merge) and use this as baseline
  BestWiggleFactor,
};

AutoMergeCostBaseline parseAutoMergeCostBaseline(std::string str);

class LtsParameters {
  private:
  unsigned int rate;
  double wiggleFactorMinimum;
  double wiggleFactorStepsize;
  bool wiggleFactorEnforceMaximumDifference;
  unsigned int maxNumberOfClusters;
  bool autoMergeClusters;
  double allowedPerformanceLossRatioAutoMerge;
  AutoMergeCostBaseline autoMergeCostBaseline = AutoMergeCostBaseline::BestWiggleFactor;

  public:
  [[nodiscard]] unsigned int getRate() const;
  [[nodiscard]] bool isWiggleFactorUsed() const;
  [[nodiscard]] double getWiggleFactorMinimum() const;
  [[nodiscard]] double getWiggleFactorStepsize() const;
  [[nodiscard]] bool getWiggleFactorEnforceMaximumDifference() const;
  [[nodiscard]] int getMaxNumberOfClusters() const;
  [[nodiscard]] bool isAutoMergeUsed() const;
  [[nodiscard]] double getAllowedPerformanceLossRatioAutoMerge() const;
  [[nodiscard]] AutoMergeCostBaseline getClusterMergingBaseline() const;

  LtsParameters(unsigned int rate,
                double wiggleFactorMinimum,
                double wiggleFactorStepsize,
                bool wigleFactorEnforceMaximumDifference,
                int maxNumberOfClusters,
                bool ltsAutoMergeClusters,
                double allowedPerformanceLossRatioAutoMerge,
                AutoMergeCostBaseline autoMergeCostBaseline);
};

LtsParameters readLtsParametersFromYaml(std::shared_ptr<YAML::Node>& params);

} // namespace seissol::initializers::time_stepping

#endif // SEISSOL_LTSCONFIGURATION_H
