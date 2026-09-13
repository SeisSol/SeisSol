// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_DATA_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_DATA_H_

#include "Common/Constants.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"

#include <array>

namespace seissol::tensor {
struct epsInit;
} // namespace seissol::tensor

namespace seissol::kernels::solver::nonlinearck {

/// What the nonlinear kernels need about a cell beyond its material.
struct NonLinearLocalData {
  /// Strain the cell carries before the first timestep, in the precision the
  /// kernels work in. The material holds it as parameters; the step kernel
  /// reads it once per time node and adds it to every strain it forms.
  // NOLINTNEXTLINE
  alignas(Alignment) real epsInit[zeroGuard(kernels::size<tensor::epsInit>())]{};

  /// Largest wave speed over the cell's nodes, from the damage-dependent
  /// moduli. The Rusanov dissipation is scaled with it.
  double maxWaveSpeed{};
};

/// The same, for the far side of each of the four faces. Only the wave speed
/// crosses: the neighbour's stress arrives through the integrals, so its
/// material does not have to.
struct NonLinearNeighborData {
  std::array<double, 4> maxWaveSpeed{};
};

} // namespace seissol::kernels::solver::nonlinearck

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_DATA_H_
