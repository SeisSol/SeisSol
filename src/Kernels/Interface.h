// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_INTERFACE_H_
#define SEISSOL_SRC_KERNELS_INTERFACE_H_

#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Memory/Descriptor/LTS.h"
#include <Common/Constants.h>

namespace seissol::tensor {
class Iane;
} // namespace seissol::tensor

namespace seissol::kernels {
struct LocalTmp {
  alignas(Alignment) real
      timeIntegratedAne[zeroLengthArrayHandler(kernels::size<tensor::Iane>())]{};
  alignas(Alignment)
      std::array<real, tensor::averageNormalDisplacement::size()> nodalAvgDisplacements[4];
  GravitationalFreeSurfaceBc gravitationalFreeSurfaceBc;
  LocalTmp(double graviationalAcceleration)
      : gravitationalFreeSurfaceBc(graviationalAcceleration) {};
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_INTERFACE_H_
