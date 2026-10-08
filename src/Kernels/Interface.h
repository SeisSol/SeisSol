// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_INTERFACE_H_
#define SEISSOL_SRC_KERNELS_INTERFACE_H_

#include "Common/Constants.h"
#include "Common/Real.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::kernels {
/// The thread-local scratch data of the integration of a cell of the configuration `Cfg`.
template <typename Cfg>
struct LocalTmp {
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  alignas(Alignment) real timeIntegratedAne[zeroGuard(kernels::size<tensor::Iane<Cfg>>())]{};
  GravitationalFreeSurfaceBc<Cfg> gravitationalFreeSurfaceBc;
  alignas(Alignment) std::array<
      real,
      tensor::averageNormalDisplacement<Cfg>::size()> nodalAvgDisplacements[Cell::NumFaces]{};
  // the time derivatives of the cell in its current step, where a boundary reads its state at a
  // time (the nonlinear Dirichlet boundary); null otherwise
  const real* timeDerivatives{nullptr};
  explicit LocalTmp(double graviationalAcceleration)
      : gravitationalFreeSurfaceBc(graviationalAcceleration) {};
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_INTERFACE_H_
