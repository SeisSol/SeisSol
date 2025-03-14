// SPDX-FileCopyrightText: 2015 SeisSol Group
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
#include <Memory/Descriptor/LTS.h>
#include <Initializer/Parameters/SourceParameters.h>
#include "Solver/Clustering/TimeManager.h"
#include <cstdarg>

namespace seissol::solver::clustering {
class TimeManager;
} // namespace seissol::solver::clustering

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
                          seissol::solver::clustering::TimeManager& timeManager);
};

} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_MANAGER_H_
