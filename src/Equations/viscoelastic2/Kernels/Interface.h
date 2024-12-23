// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 **/

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_INTERFACE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_INTERFACE_H_

#include "Equations/elastic/Kernels/GravitationalFreeSurfaceBC.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/InterfaceHelper.h"
#include "Kernels/Precision.h"

namespace seissol::kernels {
struct alignas(Alignment) LocalTmp {
  alignas(Alignment) real timeIntegratedAne[tensor::Iane::size()]{};
  alignas(Alignment)
      std::array<real, tensor::averageNormalDisplacement::size()> nodalAvgDisplacements[4]{};
  GravitationalFreeSurfaceBc gravitationalFreeSurfaceBc;
  LocalTmp(double gravitationalAcceleration)
      : gravitationalFreeSurfaceBc(gravitationalAcceleration) {};
};
LTSTREE_GENERATE_INTERFACE_GETTERED(LocalData,
                                    initializer::LTS,
                                    cellInformation,
                                    localIntegration,
                                    dofs,
                                    dofsAne,
                                    faceDisplacements)
LTSTREE_GENERATE_INTERFACE_GETTERED(
    NeighborData, initializer::LTS, cellInformation, neighboringIntegration, dofs, dofsAne)
} // namespace seissol::kernels

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_INTERFACE_H_
