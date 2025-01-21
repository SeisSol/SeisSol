// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_SOURCETERM_MANAGER_H_
#define SEISSOL_SRC_SOURCETERM_MANAGER_H_

#include "Geometry/MeshReader.h"
#include "Memory/Tree/Lut.h"
#include "Solver/time_stepping/TimeManager.h"
#include <cstdarg>

namespace seissol::sourceterm {

class Manager {
  public:
  Manager() = default;
  ~Manager() = default;

  static void loadSources(seissol::initializer::parameters::PointSourceType sourceType,
                          const char* fileName,
                          const seissol::geometry::MeshReader& mesh,
                          seissol::initializer::LTSTree* ltsTree,
                          seissol::initializer::LTS* lts,
                          seissol::initializer::Lut* ltsLut,
                          time_stepping::TimeManager& timeManager);
};

} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_MANAGER_H_
