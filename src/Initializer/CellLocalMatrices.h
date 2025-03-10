// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_CELLLOCALMATRICES_H_
#define SEISSOL_SRC_INITIALIZER_CELLLOCALMATRICES_H_

#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/Boundary.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Lut.h"

#include <array>


namespace seissol::initializer {
class EasiBoundary;
/**
 * Computes the star matrices A*, B*, and C*, and solves the Riemann problems at the interfaces.
 **/
void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTSTree* ltsTree,
                                 LTS* lts,
                                 Lut* ltsLut,
                                 const TimeStepping& timeStepping,
                                 const parameters::ModelParameters& modelParameters);

void initializeBoundaryMappings(const seissol::geometry::MeshReader& meshReader,
                                const EasiBoundary* easiBoundary,
                                LTSTree* ltsTree,
                                LTS* lts,
                                Lut* ltsLut);

void initializeDynamicRuptureMatrices(
    const seissol::geometry::MeshReader& i_meshReader,
    LTSTree* io_ltsTree,
    LTS* i_lts,
    Lut* i_ltsLut,
    std::array<LTSTree*, MULTIPLE_SIMULATIONS> dynRupTree,
    std::array<std::shared_ptr<DynamicRupture>, MULTIPLE_SIMULATIONS> dynRup,
    unsigned* ltsFaceToMeshFace,
    const GlobalData& global,
    double etaHack);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CELLLOCALMATRICES_H_
