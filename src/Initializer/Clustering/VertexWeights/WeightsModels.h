// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_WEIGHTSMODELS_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_WEIGHTSMODELS_H_

#include "Initializer/Clustering/VertexWeights/VertexWeightModel.h"

namespace seissol::initializer {

class ExponentialWeights : public VertexWeightModel {
  protected:
  int evaluateNumberOfConstraints() final { return 1; }
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};

class ExponentialBalancedWeights : public VertexWeightModel {
  protected:
  int evaluateNumberOfConstraints() final { return 2; }
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};

class EncodedBalancedWeights : public VertexWeightModel {
  protected:
  int evaluateNumberOfConstraints() final;
  void setVertexWeights() final;
  void setAllowedImbalances() final;
};
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_WEIGHTSMODELS_H_
