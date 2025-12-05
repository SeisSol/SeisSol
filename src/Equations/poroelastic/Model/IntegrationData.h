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
#include "Kernels/Precision.h"

namespace seissol::model {

template <typename Cfg>
struct PoroelasticLocalData {
  using real = Real<Cfg>;

  real sourceMatrix[seissol::tensor::ET<Cfg>::size()]{};
  real G[PoroElasticMaterial::NumQuantities]{};
  real typicalTimeStepWidth{};
  real Zinv[PoroElasticMaterial::NumQuantities][Cfg::ConvergenceOrder * Cfg::ConvergenceOrder]{};
};
struct PoroelasticNeighborData {};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_POROELASTIC_MODEL_INTEGRATIONDATA_H_
