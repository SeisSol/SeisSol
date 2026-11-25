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
#include "Core/IO/Easi/EasiContext.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Typedefs.h"
#include "WavePropagation/Physics/InitialField.h"

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

/// Projects the easi fields `iniFields`, given for every simulation of the run, onto the cells of
/// `storage`; each file is sampled once.
void projectEasiInitialField(const parameters::SimulationFiles& iniFields,
                             const seissol::geometry::MeshReader& meshReader,
                             LTS::Storage& storage,
                             bool needsTime,
                             const data::EasiContext& easi);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_FEATURES_INITIALCONDITIONS_SETUP_INITIALFIELDPROJECTION_H_
