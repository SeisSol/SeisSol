// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_WEIGHTSMODELS_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_WEIGHTSMODELS_H_

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
      : LtsWeights(config, seissolInstance) {}
  ~ExponentialBalancedWeights() override = default;

  protected:
  int evaluateNumberOfConstraints() final { return 2; }
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};

class EncodedBalancedWeights : public LtsWeights {
  public:
  explicit EncodedBalancedWeights(const LtsWeightsConfig& config, seissol::SeisSol& seissolInstance)
      : LtsWeights(config, seissolInstance) {}
  ~EncodedBalancedWeights() override = default;

  protected:
  int evaluateNumberOfConstraints() final;
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};
} // namespace initializer::time_stepping
} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSWEIGHTS_WEIGHTSMODELS_H_
