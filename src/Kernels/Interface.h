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

namespace seissol::kernels {

template <typename Cfg>
struct LocalTmp {
  alignas(Alignment)
      Real<Cfg> timeIntegratedAne[zeroLengthArrayHandler(kernels::size<tensor::Iane<Cfg>>())]{};
  GravitationalFreeSurfaceBc<Cfg> gravitationalFreeSurfaceBc;
  alignas(Alignment) std::array<
      Real<Cfg>,
      tensor::averageNormalDisplacement<Cfg>::size()> nodalAvgDisplacements[Cell::NumFaces];
  LocalTmp(double graviationalAcceleration)
      : gravitationalFreeSurfaceBc(graviationalAcceleration) {};
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_INTERFACE_H_
