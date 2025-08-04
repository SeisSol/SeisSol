// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_

#include "GeneratedCode/tensor.h"
#include <Kernels/Precision.h>

namespace seissol::tensor {
template <typename>
class ET;
} // namespace seissol::tensor

namespace seissol::model {

template <typename Cfg>
struct ViscoElasticQELocalData {
  Real<Cfg> sourceMatrix[seissol::tensor::ET<Cfg>::size()];
};
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC_MODEL_INTEGRATIONDATA_H_
