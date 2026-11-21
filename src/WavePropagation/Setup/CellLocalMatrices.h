// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SETUP_CELLLOCALMATRICES_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SETUP_CELLLOCALMATRICES_H_

#include "Core/Geometry/MeshReader.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "WavePropagation/WavePropagationParameters.h"

namespace seissol::initializer {

/**
 * Computes the star matrices A*, B*, and C*, and solves the Riemann problems at the interfaces.
 **/
void initializeCellLocalMatrices(const seissol::geometry::MeshReader& meshReader,
                                 LTS::Storage& ltsStorage,
                                 const ClusterLayout& clusterLayout,
                                 const parameters::WavePropagationParameters& parameters);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_WAVEPROPAGATION_SETUP_CELLLOCALMATRICES_H_
