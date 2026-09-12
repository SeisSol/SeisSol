// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_

#include "GeneratedCode/tensor.h"

#include <cstddef>
#include <variant>
#include <yateto/InitTools.h>

namespace seissol::numerical {
template <typename>
class MonomialBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::nonlinearck {

struct NonLinearLocalData;
struct NonLinearNeighborData;

class Spacetime;
class Time;
class Local;
class Neighbor;

struct Solver {
  using SpacetimeKernelT = Spacetime;
  using TimeKernelT = Time;
  using LocalKernelT = Local;
  using NeighborKernelT = Neighbor;

  template <typename RealT>
  using TimeBasis = seissol::numerical::MonomialBasis<RealT>;

  /// A cell hands its neighbours the time-integrated state and the
  /// time-integrated stress. The stress is carried rather than recomputed
  /// because the flux is nonlinear in time: the integral of the stress over a
  /// timestep is not the stress of the integrated strain once the internal
  /// variables move within the step. Carrying it also means neither side of a
  /// face ever needs the other's material.
  static constexpr std::size_t IntegralsSize = tensor::I::size() + tensor::sigmaI::size();
  static constexpr std::size_t DerivativesSize = yateto::computeFamilySize<tensor::dQ>();

  using LocalData = NonLinearLocalData;
  using NeighborData = NonLinearNeighborData;
};

} // namespace seissol::kernels::solver::nonlinearck
#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_
