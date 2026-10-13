// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_

#include "Common/Real.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"

#include <cstddef>
#include <yateto/InitTools.h>

namespace seissol::kernels::solver::linearckanelastic {

// TODO: maybe at some point remove the zeroGuards

template <typename Cfg>
struct AnelasticLocalData {
  // NOLINTNEXTLINE
  Real<Cfg> E[zeroGuard(kernels::size<tensor::E<Cfg>>())]{};
  Real<Cfg> w[zeroGuard(kernels::size<tensor::w<Cfg>>())]{};
  // NOLINTNEXTLINE
  Real<Cfg> W[zeroGuard(kernels::size<tensor::W<Cfg>>())]{};
};

template <typename Cfg>
struct AnelasticNeighborData {
  Real<Cfg> w[zeroGuard(kernels::size<tensor::w<Cfg>>())]{};
};

} // namespace seissol::kernels::solver::linearckanelastic
#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_DATA_H_
