// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_

#include "Common/Constants.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"

namespace seissol::tensor {
struct ET;
struct E;
struct w;
struct W;
} // namespace seissol::tensor

namespace seissol::model {

// TODO: remove zeroLengthArrayHandler when only initialized where relevant

struct ViscoElasticLocalDataExtend {
  real sourceMatrix[zeroLengthArrayHandler(kernels::size<tensor::ET>())]{};
};
struct ViscoElasticNeighborDataExtend {};

struct ViscoElasticLocalDataSplit {
  // NOLINTNEXTLINE
  real E[zeroLengthArrayHandler(kernels::size<tensor::E>())]{};
  real w[zeroLengthArrayHandler(kernels::size<tensor::w>())]{};
  // NOLINTNEXTLINE
  real W[zeroLengthArrayHandler(kernels::size<tensor::W>())]{};
};
struct ViscoElasticNeighborDataSplit {
  real w[zeroLengthArrayHandler(kernels::size<tensor::w>())]{};
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_
