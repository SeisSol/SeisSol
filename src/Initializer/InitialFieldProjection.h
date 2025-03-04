// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_
#define SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_

#include <memory>
#include <vector>

#include "Geometry/MeshReader.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Lut.h"
#include "Physics/InitialField.h"

namespace seissol::initializer {
void projectInitialField(const std::vector<std::unique_ptr<physics::InitialField>>& iniFields,
                         const GlobalData& globalData,
                         const seissol::geometry::MeshReader& meshReader,
                         seissol::initializer::MemoryManager& memoryManager,
                         LTS const& lts,
                         const Lut& ltsLut);

std::vector<double> projectEasiFields(const std::vector<std::string>& iniFields,
                                      double time,
                                      const seissol::geometry::MeshReader& meshReader,
                                      bool needsTime);

void projectEasiInitialField(const std::vector<std::string>& iniFields,
                             const GlobalData& globalData,
                             const seissol::geometry::MeshReader& meshReader,
                             seissol::initializer::MemoryManager& memoryManager,
                             LTS const& lts,
                             const Lut& ltsLut,
                             bool needsTime);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_
