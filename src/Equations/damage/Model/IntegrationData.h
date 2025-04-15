// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_INTEGRATIONDATA_H_

// #include "Datastructures.h"
#include "Kernels/Precision.h"
#include "generated_code/tensor.h"
#include <Geometry/MeshDefinition.h>

#include<array>

namespace seissol::model {

struct DamageLocalData {
    // std::vector<VrtxCoords>* globalVertPtr;
    Vertex localVertices[4];
    real localVolume;
    real localSurfaces[4];
    std::array<std::array<double, 3>,4> localNormal;
    std::array<std::array<double, 3>,4> localTangent1;
    std::array<std::array<double, 3>,4> localTangent2;
    // unsigned int globalMeshId;
};
struct DamageNeighborData {};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DAMAGE_MODEL_INTEGRATIONDATA_H_
