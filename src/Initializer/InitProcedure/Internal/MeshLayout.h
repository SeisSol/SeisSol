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
 * Moves the interior cells of each time cluster behind its copy cells, into the copy layer, where
 * the two layers hold fewer than `threshold` cells together: then the cluster computes them in one
 * go instead of two smaller ones, and gives up overlapping its interior with the halo exchange.
 * Returns the number of clusters whose interior cells moved.
 */
std::size_t combineCopyInterior(std::vector<ClusterMap>& map,
                                const LTSColorMap& colormap,
                                std::size_t threshold);

std::vector<std::vector<std::size_t>> layoutDR(const std::vector<std::size_t>& color,
                                               const std::vector<std::size_t>& ghostColor,
                                               const LTSColorMap& colormap,
                                               const geometry::MeshReader& meshReader);

/**
 * Moves the fault faces of the copy layer of each time cluster into its interior layer where the
 * two hold fewer than `threshold` faces together; the interior faces then wait for the ghost cells
 * as well. Returns the number of clusters whose copy faces moved.
 */
std::size_t combineCopyInteriorDR(std::vector<std::vector<std::size_t>>& map,
                                  const LTSColorMap& colormap,
                                  std::size_t threshold);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
