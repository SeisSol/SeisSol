// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_DATA_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_DATA_H_

#include <array>

namespace seissol::kernels::solver::nonlinearck {

/// What the nonlinear kernels need about a cell beyond its material.
struct NonLinearLocalData {
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
