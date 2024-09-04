#ifndef SEISSOL_LTSWEIGHTSMODELS_H
#define SEISSOL_LTSWEIGHTSMODELS_H

#include "LtsWeights.h"

namespace seissol {
  class SeisSol;

namespace initializer::time_stepping {

class ExponentialWeights : public LtsWeights {
public:
  explicit ExponentialWeights(const LtsWeightsConfig& config, seissol::SeisSol& seissolInstance)
      : LtsWeights(config, seissolInstance) {}
  ~ExponentialWeights() override = default;

protected:
  int evaluateNumberOfConstraints() final { return 1; }
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};


class ExponentialBalancedWeights : public LtsWeights {
public:
  explicit ExponentialBalancedWeights(const LtsWeightsConfig& config,
                                      seissol::SeisSol& seissolInstance) 
      : LtsWeights(config, seissolInstance) {
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
                                  seissol::SeisSol& seissolInstance)
      : LtsWeights(config, seissolInstance) {}
  ~EncodedBalancedWeights() override = default;

protected:
  int evaluateNumberOfConstraints() final;
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};
}
}

#endif //SEISSOL_LTSWEIGHTSMODELS_H
