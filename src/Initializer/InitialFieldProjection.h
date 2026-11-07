// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_
#define SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_

#include "Geometry/MeshReader.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Physics/InitialField.h"

#include <memory>
#include <string>
#include <vector>

namespace seissol::initializer {
/// Projects the initial conditions `iniFields`, given for the cells of each configuration by its
/// id, onto the cells of `storage`.
void projectInitialField(
    const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& iniFields,
    const seissol::geometry::MeshReader& meshReader,
    LTS::Storage& storage);

/// The values of the scripted field `path` at the quadrature points of every element, as the
/// configuration `Cfg` projects them, for each of its fused simulations: by element, point,
/// simulation and component of the state (model::stateNames, memory variables included).
///
/// The script reads x, y, z, t, sim, the group of the element, and the parameters of the material
/// of its cell in `storage` (model::scriptParameterNames), which are NaN for an element of another
/// configuration. What else it reads is an error.
template <typename Cfg>
std::vector<double> projectScriptFields(const std::string& path,
                                        double time,
                                        const seissol::geometry::MeshReader& meshReader,
                                        const LTS::Storage& storage,
                                        bool needsTime);

void projectScriptInitialField(const std::string& path,
                               const seissol::geometry::MeshReader& meshReader,
                               LTS::Storage& storage,
                               bool needsTime);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_
