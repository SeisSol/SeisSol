// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_

#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"

#include <cstddef>
#include <yateto/InitTools.h>

namespace seissol::kernels::solver::linearckanelastic {

// TODO: maybe at some point remove the zeroGuards

struct AnelasticLocalData {
  // NOLINTNEXTLINE
  real E[zeroGuard(kernels::size<tensor::E<Config>>())]{};
  real w[zeroGuard(kernels::size<tensor::w<Config>>())]{};
  // NOLINTNEXTLINE
  real W[zeroGuard(kernels::size<tensor::W<Config>>())]{};
};

struct AnelasticNeighborData {
  real w[zeroGuard(kernels::size<tensor::w<Config>>())]{};
};

} // namespace seissol::kernels::solver::linearckanelastic
#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_
