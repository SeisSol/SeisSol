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

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace seissol::initializer {

/// The node `node` of the face basis of the configuration `Cfg`, in the reference coordinates of a
/// face.
template <typename Cfg>
std::array<double, 2> referenceFaceNode(std::size_t node);

class DirichletCondition;
/// Sets up the faces of the cells of `ltsStorage` that need face data. A Dirichlet face takes the
/// conditions `dirichletConditions`, of which `conditionOfSimulation` names the one of every
/// simulation of the run; the simulations without one have none.
void initializeBoundaryMappings(
    const seissol::geometry::MeshReader& meshReader,
    const std::vector<std::optional<DirichletCondition>>& dirichletConditions,
    const std::vector<std::size_t>& conditionOfSimulation,
    LTS::Storage& ltsStorage);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_BOUNDARYMAPPINGS_H_
