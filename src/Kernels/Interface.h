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
#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::tensor {
struct Iane;
struct sourceI;
} // namespace seissol::tensor

namespace seissol::kernels {
struct LocalTmp {
  alignas(Alignment) real timeIntegratedAne[zeroGuard(kernels::size<tensor::Iane>())]{};
  /// Time-integrated source of the quantities that carry no flux. It is
  /// produced by the predictor and consumed by the corrector of the same cell
  /// in the same timestep, and no neighbour ever asks for it.
  alignas(Alignment) real sourceIntegral[zeroGuard(kernels::size<tensor::sourceI>())]{};
  /// What the derivative recursion transports by, for one cell and one step.
  /// Assembled by the predictor and consumed by the recursion right after, so
  /// it never outlives the call and no cluster can read another's.
  alignas(Alignment) real transport[3][zeroGuard(kernels::size<tensor::transport>(0))]{};
  GravitationalFreeSurfaceBc gravitationalFreeSurfaceBc;
  alignas(Alignment)
      std::array<real,
                 tensor::averageNormalDisplacement::size()> nodalAvgDisplacements[Cell::NumFaces]{};
  explicit LocalTmp(double graviationalAcceleration)
      : gravitationalFreeSurfaceBc(graviationalAcceleration) {};
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_INTERFACE_H_
