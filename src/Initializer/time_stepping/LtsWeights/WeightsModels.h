#ifndef SEISSOL_LTSWEIGHTSMODELS_H
#define SEISSOL_LTSWEIGHTSMODELS_H

#include "LtsWeights.h"


namespace seissol::initializers::time_stepping {

class ExponentialWeights : public LtsWeights {
public:
  explicit ExponentialWeights(const LtsWeightsConfig& config, const LtsParameters* ltsParameters)
      : LtsWeights(config, ltsParameters) {}
  ~ExponentialWeights() override = default;

protected:
  int evaluateNumberOfConstraints() final { return 1; }
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};


class ExponentialBalancedWeights : public LtsWeights {
public:
  explicit ExponentialBalancedWeights(const LtsWeightsConfig& config,
                                      const LtsParameters* ltsParameters)
      : LtsWeights(config, ltsParameters) {
  }
  ~ExponentialBalancedWeights() override = default;

protected:
  int evaluateNumberOfConstraints() final { return 2; }
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};


class EncodedBalancedWeights : public LtsWeights {
public:
  explicit EncodedBalancedWeights(const LtsWeightsConfig& config,
                                  const LtsParameters* ltsParameters)
      : LtsWeights(config, ltsParameters) {}
  ~EncodedBalancedWeights() override = default;

protected:
  int evaluateNumberOfConstraints() final;
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};
}

#endif //SEISSOL_LTSWEIGHTSMODELS_H
