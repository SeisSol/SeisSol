// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_

#include <Geometry/MeshReader.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/TimeStepping/Halo.h>
#include <Memory/Tree/Colormap.h>
#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace seissol::initializer::internal {

struct ClusterLayout {
  std::vector<std::size_t> cellMap;
  std::vector<RemoteCellRegion> regions;
};

std::vector<ClusterLayout> layoutCells(const std::vector<std::size_t>& color,
                                       const std::vector<std::size_t>& ghostColor,
                                       const LTSColorMap& colormap,
                                       const geometry::MeshReader& meshReader);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
