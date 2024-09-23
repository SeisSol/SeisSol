// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_SOURCETERM_MANAGER_H_
#define SEISSOL_SRC_SOURCETERM_MANAGER_H_

#include "Geometry/MeshReader.h"
#include "Initializer/Tree/Lut.h"
#include "Solver/time_stepping/TimeManager.h"
#include <cstdarg>

namespace seissol::sourceterm {

class Manager {
  public:
  Manager() = default;
  ~Manager() = default;

  void loadSources(seissol::initializer::parameters::PointSourceType sourceType,
                   const char* fileName,
                   const seissol::geometry::MeshReader& mesh,
                   seissol::initializer::LTSTree* ltsTree,
                   seissol::initializer::LTS* lts,
                   seissol::initializer::Lut* ltsLut,
                   time_stepping::TimeManager& timeManager);
};

} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_MANAGER_H_
