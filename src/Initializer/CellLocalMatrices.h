// SPDX-FileCopyrightText: 2015 SeisSol Group
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
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Memory/Tree/Backmap.h>

namespace seissol::initializer {
class EasiBoundary;
/**
 * Computes the star matrices A*, B*, and C*, and solves the Riemann problems at the interfaces.
 **/
void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTS::Tree* ltsTree,
                                 const ClusterLayout& clusterLayout,
                                 const parameters::ModelParameters& modelParameters);

void initializeBoundaryMappings(const seissol::geometry::MeshReader& meshReader,
                                const EasiBoundary* easiBoundary,
                                LTS::Tree* ltsTree);

void initializeDynamicRuptureMatrices(const seissol::geometry::MeshReader& meshReader,
                                      LTS::Tree* ltsTree,
                                      LTS::Backmap* backmap,
                                      DynamicRupture::Tree* dynRupTree,
                                      unsigned* ltsFaceToMeshFace,
                                      const GlobalData& global,
                                      double etaHack);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CELLLOCALMATRICES_H_
