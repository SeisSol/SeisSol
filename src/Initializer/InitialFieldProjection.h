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

#ifndef SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_
#define SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_

#include <memory>

#include "Geometry/MeshReader.h"
#include "Initializer/LTS.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Tree/Lut.h"
#include "Initializer/Typedefs.h"
#include "Physics/InitialField.h"

namespace seissol::initializer {
void projectInitialField(const std::vector<std::unique_ptr<physics::InitialField>>& iniFields,
                         const GlobalData& globalData,
                         const seissol::geometry::MeshReader& meshReader,
                         seissol::initializer::MemoryManager& memoryManager,
                         LTS const& lts,
                         const Lut& ltsLut);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_INITIALFIELDPROJECTION_H_
