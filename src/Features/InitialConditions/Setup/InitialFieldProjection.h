// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_FEATURES_INITIALCONDITIONS_SETUP_INITIALFIELDPROJECTION_H_
#define SEISSOL_SRC_FEATURES_INITIALCONDITIONS_SETUP_INITIALFIELDPROJECTION_H_

#include "Core/Data/SimulationFiles.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Typedefs.h"
#include "Features/InitialConditions/Physics/InitialField.h"

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

/// The values of the easi fields `iniFields` at the quadrature points of every element, as the
/// configuration `Cfg` projects them: by element, point, quantity and field.
template <typename Cfg>
std::vector<double> projectEasiFields(const std::vector<std::string>& iniFields,
                                      double time,
                                      const seissol::geometry::MeshReader& meshReader,
                                      bool needsTime);

/// Projects the easi fields `iniFields`, given for every simulation of the run, onto the cells of
/// `storage`; each file is sampled once.
void projectEasiInitialField(const parameters::SimulationFiles& iniFields,
                             const seissol::geometry::MeshReader& meshReader,
                             LTS::Storage& storage,
                             bool needsTime);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_SETUP_INITIALFIELDPROJECTION_H_
