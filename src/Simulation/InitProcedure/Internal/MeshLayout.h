// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_SIMULATION_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
#define SEISSOL_SRC_SIMULATION_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_

#include "Core/Common/BasicTypedefs.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Memory/Tree/Colormap.h"
#include "Core/TimeStepping/Layout/Halo.h"

#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace seissol::initializer::internal {

std::vector<ClusterMap> layoutCells(const std::vector<std::size_t>& color,
                                    const std::vector<std::size_t>& ghostColor,
                                    const LTSColorMap& colormap,
                                    const geometry::MeshReader& meshReader);

std::vector<std::vector<std::size_t>> layoutDR(const std::vector<std::size_t>& color,
                                               const std::vector<std::size_t>& ghostColor,
                                               const LTSColorMap& colormap,
                                               const geometry::MeshReader& meshReader);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_SIMULATION_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
