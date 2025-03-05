// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_INTERFACE_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_INTERFACE_H_

#include "Equations/elastic/Kernels/GravitationalFreeSurfaceBC.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/InterfaceHelper.h"

namespace seissol::kernels {
struct LocalTmp {
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
                                    faceDisplacements,
                                    boundaryMapping,
                                    material)
LTSTREE_GENERATE_INTERFACE_GETTERED(
    NeighborData, initializer::LTS, cellInformation, neighboringIntegration, dofs)
#else
LTSTREE_GENERATE_INTERFACE_GETTERED(LocalData,
                                    initializer::LTS,
                                    cellInformation,
                                    localIntegration,
                                    neighboringIntegration,
                                    dofs,
                                    faceDisplacements,
                                    faceDisplacementsDevice,
                                    plasticity,
                                    boundaryMapping,
                                    boundaryMappingDevice,
                                    material)
LTSTREE_GENERATE_INTERFACE_GETTERED(
    NeighborData, initializer::LTS, cellInformation, neighboringIntegration, dofs)
#endif
} // namespace seissol::kernels

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_INTERFACE_H_
