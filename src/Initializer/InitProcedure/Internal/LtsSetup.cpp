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
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <vector>

namespace seissol::initializer::internal {

namespace {

/**
 * Gets the lts setup in relation to the four face neighbors.
 *   Remark: Remember to perform the required normalization step.
 *
 * -------------------------------------------------------------------------------
 *
 *  0 in one of the first four bits: Face neighboring data are buffers.
 *  1 in one of the first four bits: Face neighboring data are derivatives.
 *
 *     Example 1:
 *     [           12 rem. bits               | buf/der bits ]
 *     [  -  -  -  -  -  -  -  -  -  -  -  -  |  0  1  1  0  ]
 *     [ 15 14 13 12 11 10  9  8  7  6  5  4  |  3  2  1  0  ]
 *  In Example 1 the data for face neighbors 0 and 3 are buffers and for 1 and 2 derivatives.
 *
 *  0 in one of bits 4 - 7: No global time stepping
 *  1 in one of bits 4 - 7: The  current cell has a global time stepping relation with the face
 *neighbor.
 *
 *     Example 2:
 *     [       8 rem. bits       |   GTS bits  | buf/der bits ]
 *     [  -  -  -  -  -  -  -  - | 0  0  1  1  |  0  1  1  0  ]
 *     [ 15 14 13 12 11 10  9  8 | 7  6  5  4  |  3  2  1  0  ]
 *  In Example 2 the data of face neighbors 0 and 3 are buffers, 1 and 2 deliver derivatives
 *  Face neighbor 0 has a GTS-relation and this cell works directly on the delivered buffer.
 *  Face neighbor 1 has a GTS-relation, but delivers derivatives -> The derivatives have to
 *translated to time integrated DOFs first. Face neighbor 2 has a LTS-relation and receives
 *derivatives from its neighbor -> The derivates have to be used for a partial time integration.
 *  Face neighbor 3 has a LTS-relation and can operate on the buffers directly.
 *
 * -------------------------------------------------------------------------------
 *
 *  1 in the eigth bit: the cell is required to work on time integration buffers.
 *  1 in the nineth bit: the cell is required to compute time derivatives.
 *
 *     Example 3:
 *     [     remaining     | der. buf. |       first 8 bits       ]
 *     [  -  -  -  -  -  - |  0    1   |  -  -  -  -  -  -  -  -  ]
 *     [ 15 14 13 12 11 10 |  9    8   |  7  6  5  4  3  2  1  0  ]
 *  In Example 3 only a buffer is stored as for example in global time stepping.
 *
 *     Example 4:
 *     [     remaining     | der. buf. |       first 8 bits       ]
 *     [  -  -  -  -  -  - |  1    1   |  -  -  -  -  -  -  -  -  ]
 *     [ 15 14 13 12 11 10 |  9    8   |  7  6  5  4  3  2  1  0  ]
 *  In Example 4 both (buffer+derivative) is stored.
 *
 * -------------------------------------------------------------------------------
 *
 *  1 in the tenth bit: the cell local buffer is a LTS buffer (reset on request only).
 *
 *     Example 5:
 *     [   remaining    | LTS buf. |          first 10 bits         ]
 *     [  -  -  -  -  - |     1    |  -  1  -  -  -  -  -  -  -  -  ]
 *     [ 15 14 13 12 11 |    10    |  9  8  7  6  5  4  3  2  1  0  ]
 *  In Example 5 the buffer is a LTS buffer (reset on request only). GTS buffers are updated in
 *every time step.
 *
 * @return LTS setup (without correction)
 * @param ownPrimary primary cell information struct of the cell in consideration
 * @param ownSecondary secondary cell information struct of the cell in consideration
 * @param neighborClusters face-neighbor LTS cluster IDs
 * @param copy true if the cell is part of the copy layer (only required for correctness in dynamic
 *rupture computations).
 **/
LtsSetup getLtsSetup(const CellLocalInformation& ownPrimary,
                     const SecondaryCellLocalInformation& ownSecondary,
                     const std::array<uint64_t, Cell::NumFaces>& neighborClusters,
                     bool copy = false) {
  // reset the LTS setup
  LtsSetup ltsSetup{};

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

      if (copy) {
        // set the buffer invalid in copy layers
        // TODO: Minor improvements possible: Non-DR MPI-neighbor for example
        ltsSetup.setHasBuffer(true, BufferType::AccumulatedIntegrals);
      }
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

} // namespace

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
    const auto isCopy = layer.getIdentifier().halo == HaloType::Copy;
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
      primaryInformationLocal[cell].ltsSetup = LtsSetup(getLtsSetup(primaryInformationLocal[cell],
                                                                    secondaryInformationLocal[cell],
                                                                    neighborClusters,
                                                                    isCopy));

      // assert that the cell operates at least on buffers or derivatives
      assert(primaryInformationLocal[cell].ltsSetup.hasBuffer(BufferType::StepIntegrals) ||
             primaryInformationLocal[cell].ltsSetup.hasBuffer(BufferType::AccumulatedIntegrals) ||
             primaryInformationLocal[cell].ltsSetup.hasBuffer(BufferType::Derivatives));
    }
  }

  // get setup in the ghost layer
  haloCommunication<LTS::CellInformation>(layout, storage, ghostElementType);

  // we won't need the ghost element type after this anymore
  MPI_Type_free(&ghostElementType);
}

} // namespace seissol::initializer::internal
