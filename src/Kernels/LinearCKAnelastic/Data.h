// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_

#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"

#include <cstddef>
#include <yateto/InitTools.h>

namespace seissol::tensor {
struct ET;
struct E;
struct w;
struct W;
} // namespace seissol::tensor

namespace seissol::kernels::solver::linearckanelastic {

// TODO: maybe at some point remove the zeroLengthArrayHandlers

struct AnelasticLocalData {
  // NOLINTNEXTLINE
  real E[zeroLengthArrayHandler(kernels::size<tensor::E>())]{};
  real w[zeroLengthArrayHandler(kernels::size<tensor::w>())]{};
  // NOLINTNEXTLINE
  real W[zeroLengthArrayHandler(kernels::size<tensor::W>())]{};
};

struct AnelasticNeighborData {
  real w[zeroLengthArrayHandler(kernels::size<tensor::w>())]{};
};

} // namespace seissol::kernels::solver::linearckanelastic
#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_
