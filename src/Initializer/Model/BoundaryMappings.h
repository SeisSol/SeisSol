// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_BOUNDARYMAPPINGS_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_BOUNDARYMAPPINGS_H_

#include "Common/Constants.h"
#include "Geometry/MeshReader.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"

#include <array>
#include <cstddef>
#include <optional>

namespace seissol::initializer {

class EasiBoundary;

/**
 * Maps the face nodes (nodes2D) onto face side of the tetrahedron with the given vertices and
 * writes three global coordinates per node to nodes, as BoundaryFaceInformation::nodes stores them.
 **/
void computeBoundaryNodes(const std::array<const double*, Cell::NumVertices>& vertices,
                          std::size_t side,
                          real* nodes);

void initializeBoundaryMappings(const seissol::geometry::MeshReader& meshReader,
                                const std::optional<EasiBoundary>& easiBoundary,
                                LTS::Storage& ltsStorage);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_BOUNDARYMAPPINGS_H_
