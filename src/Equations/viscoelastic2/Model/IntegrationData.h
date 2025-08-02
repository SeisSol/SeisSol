// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_INTEGRATIONDATA_H_

#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"

namespace seissol::model {

struct ViscoElasticLocalData {
  real E[tensor::E<Cfg>::size()];
  real w[tensor::w<Cfg>::size()];
  real W[tensor::W<Cfg>::size()];
};
struct ViscoElasticNeighborData {
  real w[tensor::w<Cfg>::size()];
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_INTEGRATIONDATA_H_
