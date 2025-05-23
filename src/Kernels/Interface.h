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
#include "Memory/Tree/InterfaceHelper.h"
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
#ifndef ACL_DEVICE
LTSTREE_GENERATE_INTERFACE_GETTERED(LocalData,
                                    initializer::LTS,
                                    cellInformation,
                                    localIntegration,
                                    neighboringIntegration,
                                    dofs,
                                    dofsAne,
                                    faceDisplacements,
                                    boundaryMapping,
                                    material)
LTSTREE_GENERATE_INTERFACE_GETTERED(
    NeighborData, initializer::LTS, cellInformation, neighboringIntegration, dofs, dofsAne)
#else
LTSTREE_GENERATE_INTERFACE_GETTERED(LocalData,
                                    initializer::LTS,
                                    cellInformation,
                                    localIntegration,
                                    neighboringIntegration,
                                    dofs,
                                    dofsAne,
                                    faceDisplacements,
                                    faceDisplacementsDevice,
                                    plasticity,
                                    boundaryMapping,
                                    boundaryMappingDevice,
                                    material)
LTSTREE_GENERATE_INTERFACE_GETTERED(
    NeighborData, initializer::LTS, cellInformation, neighboringIntegration, dofs, dofsAne)
#endif
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_INTERFACE_H_
