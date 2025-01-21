// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_

#include "generated_code/tensor.h"

namespace seissol::model {

struct ViscoElasticLocalData {
  real sourceMatrix[seissol::tensor::ET::size()];
};
struct ViscoElasticNeighborData {};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_
