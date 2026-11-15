// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_

#include "Geometry/MeshReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Tree/Colormap.h"

#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace seissol::initializer::internal {

std::vector<ClusterMap> layoutCells(const std::vector<std::size_t>& color,
                                    const std::vector<std::size_t>& ghostColor,
                                    const LTSColorMap& colormap,
                                    const geometry::MeshReader& meshReader);

/**
 * What combineCopyInterior did: per color, whether the copy layer of that color took the interior
 * cells of its time cluster, and for how many clusters that moved cells.
 */
struct CombinedLayers {
  std::vector<bool> copyLayers;
  std::size_t moved{0};
};

/**
 * Moves the interior cells of each time cluster behind its copy cells, into the copy layer, where
 * the two layers hold fewer than `threshold` cells together: then the cluster computes them in one
 * go instead of two smaller ones, and gives up overlapping its interior with the halo exchange.
 */
CombinedLayers combineCopyInterior(std::vector<ClusterMap>& map,
                                   const LTSColorMap& colormap,
                                   std::size_t threshold);

std::vector<std::vector<std::size_t>> layoutDR(const std::vector<std::size_t>& color,
                                               const std::vector<std::size_t>& ghostColor,
                                               const LTSColorMap& colormap,
                                               const geometry::MeshReader& meshReader);

/**
 * Per color, whether it is the copy layer of a time cluster whose copy and interior layers hold
 * fewer than `threshold` fault faces together.
 */
std::vector<bool> faceLayersBelow(const std::vector<std::vector<std::size_t>>& map,
                                  const LTSColorMap& colormap,
                                  std::size_t threshold);

/**
 * Moves the fault faces of the copy layer of each time cluster into its interior layer where
 * `combine` is set for the color of the copy layer; the interior faces then wait for the ghost
 * cells as well. Returns the number of clusters whose copy faces moved.
 */
std::size_t combineCopyInteriorDR(std::vector<std::vector<std::size_t>>& map,
                                  const LTSColorMap& colormap,
                                  const std::vector<bool>& combine);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
