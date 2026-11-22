// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_FEATURES_SOURCETERM_SETUP_MANAGER_H_
#define SEISSOL_SRC_FEATURES_SOURCETERM_SETUP_MANAGER_H_

#include "Core/Geometry/MeshReader.h"
#include "Core/TimeStepping/Settings.h"
#include "Core/TimeStepping/TimeManager.h"
#include "Features/SourceTerm/SourceParameters.h"

#include <cstdarg>

namespace seissol::sourceterm {

class Manager {
  public:
  Manager() = default;
  ~Manager() = default;

  static void loadSources(seissol::initializer::parameters::PointSourceType sourceType,
                          const char* fileName,
                          const seissol::geometry::MeshReader& mesh,
                          LTS::Storage& ltsStorage,
                          LTS::Backmap& backmap,
                          const SimulationSettings& settings,
                          solver::TimeManager& timeManager);
};

} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_FEATURES_SOURCETERM_SETUP_MANAGER_H_
