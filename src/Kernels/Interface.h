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
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::kernels {
struct LocalTmp {
  alignas(Alignment) real timeIntegratedAne[zeroGuard(kernels::size<tensor::Iane<Config>>())]{};
  GravitationalFreeSurfaceBc gravitationalFreeSurfaceBc;
  alignas(Alignment) std::array<
      real,
      tensor::averageNormalDisplacement<Config>::size()> nodalAvgDisplacements[Cell::NumFaces]{};
  explicit LocalTmp(double graviationalAcceleration)
      : gravitationalFreeSurfaceBc(graviationalAcceleration) {};
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_INTERFACE_H_
