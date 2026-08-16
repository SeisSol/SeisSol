// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_INTEGRATIONDATA_H_

#include "Common/Constants.h"
#include "Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"

namespace seissol::tensor {
struct ET;
} // namespace seissol::tensor

namespace seissol::model {

// TODO: remove zeroLengthArrayHandler when only initialized where relevant

struct PoroelasticLocalData {
  real sourceMatrix[zeroLengthArrayHandler(kernels::size<tensor::ET>())]{};
  // NOLINTNEXTLINE
  real G[PoroElasticMaterial::NumQuantities]{};
  // NOLINTNEXTLINE
  real Zinv[PoroElasticMaterial::NumQuantities][ConvergenceOrder * ConvergenceOrder]{};

  // preferrably double; will be compared closely against the "default" timestep width almost all
  // the time
  double typicalTimeStepWidth{};
};
struct PoroelasticNeighborData {};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_INTEGRATIONDATA_H_
