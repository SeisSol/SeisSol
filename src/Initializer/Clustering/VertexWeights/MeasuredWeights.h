// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_MEASUREDWEIGHTS_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_MEASUREDWEIGHTS_H_

#include "Initializer/Clustering/VertexWeights/MeasuredCellCosts.h"
#include "Initializer/Clustering/VertexWeights/VertexWeightModel.h"

#include <cstddef>
#include <mpi.h>
#include <optional>
#include <string>
#include <vector>

namespace seissol {
class LoopStatistics;
namespace geometry {
class MeshReader;
} // namespace geometry
} // namespace seissol

namespace seissol::initializer {

class MemoryManager;

/**
 * @brief Splits the time the loops took over the cells that they computed.
 *
 * The times are those LoopStatistics recorded per layer (see LoopStatistics::recordLayers). The
 * cells of a time cluster (in a configuration) share the time of its cell loops evenly, those of
 * the interior and the copy layer alike: what the cluster costs on this rank, including the fixed
 * costs of its loops, which matter for the small clusters. The faces of dynamic rupture of a
 * cluster share the time of its face loops the same way; a cell gets the share of each face it
 * lies at (a half, if the cell on the other side is on this rank as well).
 *
 * Nothing if no loop was timed at all on some rank (as within recordings of super-timesteps, or on
 * own queues without DEVICE_OWN_QUEUE_TIMING). Collective.
 */
std::optional<MeasuredCellCosts> measureCellCosts(const LoopStatistics& statistics,
                                                  MemoryManager& memoryManager,
                                                  const geometry::MeshReader& meshReader,
                                                  double simulatedTime,
                                                  MPI_Comm comm);

/**
 * @brief Writes the measured costs of all ranks to the HDF5 file `file`, as the datasets `ids` and
 * `costs` (in no particular order of the cells). Collective.
 */
void writeCellCosts(const std::string& file,
                    const MeasuredCellCosts& costs,
                    double simulatedTime,
                    MPI_Comm comm);

/**
 * @brief Reads the costs of the cells `firstCell` ... `firstCell + count - 1` (by their ids in the
 * mesh file) from a file writeCellCosts() wrote. Collective.
 */
std::vector<double>
    readCellCosts(const std::string& file, std::size_t firstCell, std::size_t count, MPI_Comm comm);

/**
 * @brief Brings the measured costs of the cells to the ranks that read the cells
 * `firstCell` ... `firstCell + count - 1`, as readCellCosts() does from a file. Collective.
 */
std::vector<double> distributeCellCosts(const MeasuredCellCosts& costs,
                                        std::size_t firstCell,
                                        std::size_t count,
                                        MPI_Comm comm);

/**
 * @brief Vertex weights from measured costs: a cell weighs what it cost (see measureCellCosts()),
 * scaled such that the weights of all cells sum to those of ExponentialWeights.
 *
 * The cells a rank partitions are those it read from the mesh file, which the ranks read in runs,
 * one after the other; so the costs are those of the cells `firstCell` ... with `firstCell` the
 * number of cells the ranks before read.
 */
class MeasuredWeights : public VertexWeightModel {
  public:
  //! with the costs from `file`, see writeCellCosts()
  explicit MeasuredWeights(std::string file);
  //! with the costs of the cells of this rank, as measured in the run before
  explicit MeasuredWeights(const MeasuredCellCosts& costs);

  protected:
  std::size_t evaluateNumberOfConstraints() final { return 1; }
  void setVertexWeights() final;
  void setAllowedImbalances() final;

  private:
  std::string file_;
  const MeasuredCellCosts* costs_{nullptr};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_VERTEXWEIGHTS_MEASUREDWEIGHTS_H_
