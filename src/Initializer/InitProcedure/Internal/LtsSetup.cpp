// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "LtsSetup.h"

#include "Common/Constants.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <vector>

namespace seissol::initializer::internal {

/**
 * Derives the storage requirements of a single cell from its face types and the time cluster IDs
 * of its face neighbors.
 *
 * The result is encoded in the LtsSetup bitmap. The field positions follow from BufferCountBits,
 * Cell::NumFaces and BufferCount (see LtsSetup.h); with the current values the layout is:
 *
 *   bits  0 - 7:  the BufferType supplied by each face neighbor, two bits per face
 *   bits  8 - 11: one flag per face, set iff that neighbor runs at the same time step
 *   bits 12 - 14: one flag per BufferType, set iff this cell stores data of that type
 *
 *     Example: a cell with GTS neighbors over faces 0 and 1, a neighbor in a coarser cluster
 *     over face 2, and a free-surface boundary over face 3.
 *
 *     [ 15 | 14 13 12 | 11 10  9  8 |  7  6  5  4 |  3  2  1  0 ]
 *     [  - |  1  0  1 |  0  0  1  1 |  0  0  0  1 |  0  0  0  0 ]
 *
 *  Faces 0 and 1 supply StepIntegrals and are marked as same-timestep. Face 2 supplies
 *  Derivatives, since this cell integrates over its own sub-interval. Face 3 is a boundary face
 *  and keeps the default. The cell itself stores StepIntegrals for its GTS neighbors and
 *  AccumulatedIntegrals for the coarser one, but no Derivatives.
 *
 * @return LTS setup of the cell
 * @param ownPrimary primary cell information struct of the cell in consideration
 * @param ownSecondary secondary cell information struct of the cell in consideration
 * @param neighborClusters face-neighbor LTS cluster IDs
 **/
LtsSetup getLtsSetup(const CellLocalInformation& ownPrimary,
                     const SecondaryCellLocalInformation& ownSecondary,
                     const std::array<std::uint64_t, Cell::NumFaces>& neighborClusters) {
  // reset the LTS setup
  LtsSetup ltsSetup{};

  // A solver whose face flux needs both traces at once reads the cell's own
  // integrals in the neighbouring integration, where they are not passed in
  // and cannot be rebuilt from the derivatives. Such a cell keeps them
  // whatever its faces would ask for -- including a cell that has no face
  // with a neighbour at all, which would otherwise get no buffer.
  if constexpr (kernels::Solver::RequiresOwnIntegrals) {
    ltsSetup.setHasBuffer(true, BufferType::StepIntegrals);
  }

  // iterate over the faces
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {

    const auto bcType = getBCType(ownPrimary.faceTypes[face]);

    if (bcType == BCType::External) {
      // continue for external boundary conditions without Neighbor kernel usage
      continue;
    } else if (ownPrimary.faceTypes[face] == FaceType::DynamicRupture) {
      // dynamic rupture faces are always global time stepping but operate on derivatives

      // face-neighbor provides GTS+derivatives
      ltsSetup.setNeighborBuffer(face, BufferType::Derivatives);
      ltsSetup.setNeighborGTSRelation(face, true);

      // cell is required to provide derivatives for dynamic rupture
      ltsSetup.setHasBuffer(true, BufferType::Derivatives);
    }
    // derive the LTS setup based on the cluster ids
    else {
      // neighboring cluster has a larger time step than this cluster
      if (ownSecondary.clusterId < neighborClusters[face]) {
        // neighbor delivers time derivatives
        ltsSetup.setNeighborBuffer(face, BufferType::Derivatives);

        // the cell-local buffer is used in LTS-fashion
        ltsSetup.setHasBuffer(true, BufferType::AccumulatedIntegrals);
      } else if (ownSecondary.clusterId == neighborClusters[face]) {
        // GTS relation
        ltsSetup.setNeighborGTSRelation(face, true);
        ltsSetup.setHasBuffer(true, BufferType::StepIntegrals);
        ltsSetup.setNeighborBuffer(face, BufferType::StepIntegrals);
      } else if (ownSecondary.clusterId > neighborClusters[face]) {
        // cell is required to provide derivatives
        ltsSetup.setHasBuffer(true, BufferType::Derivatives);
        ltsSetup.setNeighborBuffer(face, BufferType::AccumulatedIntegrals);
      }
    }
  }

  return ltsSetup;
}

/**
 * Derives the lts setups of all given cells.
 **/
void deriveLtsSetups(const MeshLayout& layout, LTS::Storage& storage) {
  MPI_Datatype ghostElementType = MPI_DATATYPE_NULL;
  MPI_Datatype ghostElementTypePre = MPI_DATATYPE_NULL;

  // cf. partially https://stackoverflow.com/a/33624425
  const int datatypeCount = 1;
  const std::vector<int> datatypeBlocklen{1};
  const std::vector<MPI_Aint> datatypeDisplacement{offsetof(CellLocalInformation, ltsSetup)};
  const std::vector<MPI_Datatype> datatypeDatatype{MPI_UINT16_T};
  static_assert(sizeof(uint16_t) == sizeof(LtsSetup));

  MPI_Type_create_struct(datatypeCount,
                         datatypeBlocklen.data(),
                         datatypeDisplacement.data(),
                         datatypeDatatype.data(),
                         &ghostElementTypePre);
  const MPI_Aint lb = 0;
  const MPI_Aint extent = sizeof(CellLocalInformation);
  MPI_Type_create_resized(ghostElementTypePre, lb, extent, &ghostElementType);
  MPI_Type_commit(&ghostElementType);

  // iterate over time clusters
  for (auto& layer : storage.leaves(Ghost)) {
    auto* primaryInformationLocal = layer.var<LTS::CellInformation>();
    const auto* secondaryInformationLocal = layer.var<LTS::SecondaryInformation>();
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      std::array<uint64_t, Cell::NumFaces> neighborClusters{};
      for (std::size_t face = 0; face < Cell::NumFaces; face++) {
        // only continue for non-boundary faces
        if (isInternalFaceType(primaryInformationLocal[cell].faceTypes[face])) {
          // get neighboring cell id
          const auto& neighbor = secondaryInformationLocal[cell].faceNeighbors[face];

          // get neighboring setup
          neighborClusters[face] = storage.lookup<LTS::SecondaryInformation>(neighbor).clusterId;
        }
      }

      // set the lts setup for this cell
      primaryInformationLocal[cell].ltsSetup = getLtsSetup(
          primaryInformationLocal[cell], secondaryInformationLocal[cell], neighborClusters);

      // assert that the cell operates at least on one buffer
      assert(primaryInformationLocal[cell].ltsSetup.hasAnyBuffer());
    }
  }

  // get setup in the ghost layer
  haloCommunication<LTS::CellInformation>(layout, storage, ghostElementType);

  // we won't need the ghost element type after this anymore
  MPI_Type_free(&ghostElementType);
}

} // namespace seissol::initializer::internal
