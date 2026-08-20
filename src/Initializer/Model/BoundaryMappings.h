// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_BOUNDARYMAPPINGS_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_BOUNDARYMAPPINGS_H_

#include "Geometry/MeshReader.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer {

class EasiBoundary;
void initializeBoundaryMappings(const seissol::geometry::MeshReader& meshReader,
                                const std::optional<EasiBoundary>& easiBoundary,
                                LTS::Storage& ltsStorage);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_BOUNDARYMAPPINGS_H_
