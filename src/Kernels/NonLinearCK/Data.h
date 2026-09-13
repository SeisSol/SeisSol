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
struct materialParameters;
} // namespace seissol::tensor

namespace seissol::kernels::solver::nonlinearck {

/// What the nonlinear kernels need about a cell beyond its material.
struct NonLinearLocalData {
  /// Strain the cell carries before the first timestep, in the precision the
  /// kernels work in. The material holds it as parameters; the step kernel
  /// reads it once per time node and adds it to every strain it forms.
  // NOLINTNEXTLINE
  alignas(Alignment) real epsInit[zeroGuard(kernels::size<tensor::epsInit>())]{};

  /// Largest wave speed the material can ever reach, from alpha in [0, 1] and
  /// the invariant ratio in [-sqrt(3), sqrt(3)]. It bounds the dissipation a
  /// face needs from the far side without anyone transporting a state.
  real maxWaveSpeedBound{};

  /// The material as the kernels read it: in the precision they work in,
  /// converted once rather than per timestep, and in one array rather than as
  /// arguments -- a kernel argument that is a scalar is uniform over a batch,
  /// and a material is not. Which entry is which parameter is the codegen's
  /// order, and the generated names are what fills it.
  // NOLINTNEXTLINE
  alignas(Alignment) real parameters[zeroGuard(kernels::size<tensor::materialParameters>())]{};
};

/// The same, for the far side of each of the four faces. Only the wave speed
/// crosses: the neighbour's stress arrives through the integrals, so its
/// material does not have to.
struct NonLinearNeighborData {
  std::array<real, 4> maxWaveSpeedBound{};
};

} // namespace seissol::kernels::solver::nonlinearck

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_DATA_H_
