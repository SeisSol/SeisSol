// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_CELLLOCALMATRICES_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_CELLLOCALMATRICES_H_

#include "Geometry/MeshReader.h"
#include "Initializer/Model/FaceSamples.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer {

/**
 * Computes the star matrices A*, B*, and C*, and solves the Riemann problems at the interfaces.
 * Where the flux varies along a face, it is built from the material both sides see at the points
 * of the face, which `faceSamples` holds.
 **/
void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTS::Storage& ltsStorage,
                                 const ClusterLayout& clusterLayout,
                                 const parameters::ModelParameters& modelParameters,
                                 const GlobalData& global,
                                 const FaceSamples& faceSamples);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_CELLLOCALMATRICES_H_
