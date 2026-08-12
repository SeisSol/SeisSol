// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERSMOOTHER_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERSMOOTHER_H_

#include "Geometry/PUMLReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Parameters/MeshParameters.h"

#include <cstddef>
#include <mpi.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::initializer {

/// Decodes the boundary condition of a face into a FaceType.
FaceType decodeFaceType(const void* boundaryCond,
                        std::size_t cell,
                        std::uint8_t face,
                        parameters::BoundaryFormat boundaryFormat,
                        const FaceMap& faceMap);

/// Which cluster id differences are admissible across a face.
///
/// The bound is on the cluster *index*, not on the timestep ratio. That is what the ghost
/// cluster construction requires: TimeManager::addClusters() asserts that a cluster only
/// ever exchanges with its immediate index neighbors, regardless of the rate vector.
struct SmoothingRule {
  std::size_t maxDifference{1};
  std::size_t dynamicRuptureDifference{0};

  [[nodiscard]] std::size_t differenceFor(FaceType faceType) const {
    return faceType == FaceType::DynamicRupture ? dynamicRuptureDifference : maxDifference;
  }
};

/// Enforces a bound on the cluster id difference between neighboring cells.
///
/// This is a min-plus relaxation: a cell is only ever demoted to a finer cluster, never
/// promoted, so the iteration is monotone and reaches a unique fixed point. The sweeps run
/// under OpenMP and read neighbor ids while other threads write them; the resulting race is
/// benign for the fixed point but does make the reported number of demotions vary between
/// runs. That behavior is inherited from the original implementation and kept deliberately.
///
/// The halo exchange pattern and its buffers are built once, because the fixed point is
/// re-entered for every wiggle factor candidate.
class ClusterSmoother {
  public:
  ClusterSmoother(const geometry::PumlMesh& mesh,
                  parameters::BoundaryFormat boundaryFormat,
                  const FaceMap& faceMap);

  /// One local sweep plus one halo exchange. Returns the number of demotions on this rank.
  int relaxOnce(std::vector<std::size_t>& clusterIds, const SmoothingRule& rule);

  /// Iterates to the global fixed point. Returns the total number of demotions.
  int relax(std::vector<std::size_t>& clusterIds, const SmoothingRule& rule, MPI_Comm comm);

  /// Number of ranks this rank exchanges cluster ids with.
  [[nodiscard]] std::size_t exchangeCount() const { return rankToSharedFaces_.size(); }

  /// Cells owning at least one shared face.
  [[nodiscard]] const std::vector<std::size_t>& boundaryCells() const { return boundaryCells_; }

  private:
  const geometry::PumlMesh* mesh_{nullptr};
  parameters::BoundaryFormat boundaryFormat_;
  const FaceMap* faceMap_;

  std::vector<std::pair<int, std::vector<std::size_t>>> rankToSharedFaces_;
  std::unordered_map<std::size_t, std::size_t> localFaceIdToLocalCellId_;
  std::unordered_map<std::size_t, std::pair<std::size_t, std::size_t>> sharedFaceToExchangeId_;
  std::vector<std::size_t> boundaryCells_;

  // reused across sweeps rather than reallocated on every fixed point iteration
  std::vector<MPI_Request> requests_;
  std::vector<std::vector<std::size_t>> ghost_;
  std::vector<std::vector<std::size_t>> copy_;
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_CLUSTERSMOOTHER_H_
