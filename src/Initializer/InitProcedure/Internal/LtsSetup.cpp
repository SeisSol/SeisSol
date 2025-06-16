// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "LtsSetup.h"
#include <Initializer/BasicTypedefs.h>
#include <Initializer/CellLocalInformation.h>
#include <Memory/MemoryContainer.h>
#include <Memory/Tree/Layer.h>
#include <array>
#include <cstddef>

namespace {

using namespace seissol;

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
 *  Proposed extension (WIP): bits 11-14 specify if we need to project to a nodal basis on the faces
 **/
unsigned short determineLtsSetup(const CellLocalInformation* cellInformation,
                                 const SecondaryCellLocalInformation* secondaryInformation,
                                 std::size_t self,
                                 std::array<std::size_t, 4> neighbors,
                                 bool copyLayer) {
  // reset the LTS setup
  unsigned short ltsSetup = 0;

  const auto& ownPrimary = cellInformation[self];
  const auto& ownSecondary = secondaryInformation[self];

  // iterate over the faces
  for (unsigned int face = 0; face < 4; face++) {
    // continue for boundary conditions
    if (ownPrimary.faceTypes[face] == FaceType::Outflow) {
      continue;
    }
    // fake neighbors are GTS
    else if (ownPrimary.faceTypes[face] == FaceType::FreeSurface ||
             ownPrimary.faceTypes[face] == FaceType::FreeSurfaceGravity ||
             ownPrimary.faceTypes[face] == FaceType::Dirichlet ||
             ownPrimary.faceTypes[face] == FaceType::Analytical) {
      ltsSetup |= (1 << (face + 4));
    }
    // dynamic rupture faces are always global time stepping but operate on derivatives
    else if (ownPrimary.faceTypes[face] == FaceType::DynamicRupture) {
      // face-neighbor provides GTS derivatives
      // face-neighbor provides derivatives
      ltsSetup |= (1 << face);
      ltsSetup |= (1 << (face + 4));

      // cell is required to provide derivatives for dynamic rupture
      ltsSetup |= (1 << 9);

      if (copyLayer) { // set the buffer invalid in copy layers
        // TODO(unknown): Minor improvements possible: Non-DR MPI-neighbor for example
        ltsSetup |= (1 << 10);
      }
    }
    // derive the LTS setup based on the cluster ids
    else {
      // neighboring cluster has a larger time step than this cluster
      if (ownSecondary.clusterId < secondaryInformation[neighbors[face]].clusterId) {
        // neighbor delivers time derivatives
        ltsSetup |= (1 << face);

        // the cell-local buffer is used in LTS-fashion
        ltsSetup |= (1 << 10);
      }
      // GTS relation
      else if (ownSecondary.clusterId == secondaryInformation[neighbors[face]].clusterId) {
        ltsSetup |= (1 << (face + 4));
      }

      // cell is required to provide derivatives
      if (ownSecondary.clusterId > secondaryInformation[neighbors[face]].clusterId) {
        ltsSetup |= (1 << 9);
      }
      // cell is required to provide a buffer
      else {
        ltsSetup |= (1 << 8);
      }
    }

    // true lts buffer with gts required derivatives
    if ((ltsSetup >> 10) % 2 == 1 && (ltsSetup >> 4) % 16 != 0) {
      ltsSetup |= (1 << 9);
    }
  }

  /*
   * Normalize for special case "free surface/dirichlet on derivatives":
   *   If a cell provides either buffers in a LTS fashion or derivatives only,
   *   the neighboring contribution of the boundary intergral is required to work on the cells
   * derivatives. Free surface/dirichlet boundary conditions work on the cells DOFs in the
   * neighboring contribution: Enable cell local derivatives in this case and mark that the "fake
   * neighbor" provides derivatives.
   */
  for (int face = 0; face < 4; ++face) {
    // check for special case free-surface/dirichlet requirements
    const bool isSpecialCase = ownPrimary.faceTypes[face] == FaceType::FreeSurface ||
                               ownPrimary.faceTypes[face] == FaceType::FreeSurfaceGravity ||
                               ownPrimary.faceTypes[face] == FaceType::Dirichlet ||
                               ownPrimary.faceTypes[face] == FaceType::Analytical;
    if (isSpecialCase &&              // special case face
        ((ltsSetup >> 10) % 2 == 1 || // lts fashion buffer
         (ltsSetup >> 8) % 2 == 0)) { // no buffer at all
      ltsSetup |= (1 << 9);           // enable derivatives computation
      ltsSetup |= (1 << face);        // enable derivatives for the fake face neighbor
    }
  }

  return ltsSetup;
}

/**
 * Normalizes the LTS setup for the special case "GTS on derivatives":
 *   If a face neighbor provides true buffers to cells with larger time steps,
 *   the local cell is required to operate on derivatives of this face neighbor.
 *
 *   Example:
 *        | own |  fn 1 |  fn 2 | fn 3 | fn 4 |
 *   local|  dt |    dt | 0.5dt |   dt |   dt |
 *   fn 4 |  dt | 0.5dt |   2dt |   dt |   dt |
 *         -----------------------------------
 *   In the example the local cell is connected via face 4 with to a GTS neighbor.
 *   Face neighbor 4 is required to deliver true buffers to its second face neighbor.
 *   It follows that the local cell has to operate on the derivatives of face neighbor 4.
 **/
unsigned short correctLtsSetup(const CellLocalInformation* cellInformation,
                               std::size_t self,
                               std::array<std::size_t, 4> neighbors) {

  const auto& ownPrimary = cellInformation[self];
  auto localLtsSetup = ownPrimary.ltsSetup;
  // iterate over the face neighbors
  for (int face = 0; face < 4; ++face) {
    // enforce derivatives if this is a "GTS on derivatives" relation
    if ((((localLtsSetup >> (face + 4)) % 2) != 0) &&
        (cellInformation[neighbors[face]].ltsSetup >> 10) % 2 == 1) {
      localLtsSetup |= (1 << face);
    }
  }
  return localLtsSetup;
}

// removes all non-ghost-layer-relevant cell info
unsigned short correctGhostLayer(const CellLocalInformation* cellInformation,
                                 std::size_t self,
                                 std::array<std::size_t, 4> neighbors) {
  /*auto otherLtsSetup = cellInformation[neighbors[face]].ltsSetup;
  // TODO:
  return otherLtsSetup;*/
  return cellInformation[self].ltsSetup;
}

} // namespace

namespace seissol::initializer::internal {
void handleLtsSetup(memory::MemoryContainer& container) {
  const auto* secondaryCellInformation =
      container.volume.var(container.wpdesc.secondaryInformation);
  auto* cellInformation = container.volume.var(container.wpdesc.cellInformation);
  // TODO: tree to layer mapping
  for (auto& layer : container.volume.leaves(Ghost)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < layer.size(); ++i) {
      const std::size_t self = secondaryCellInformation[i].linearId;
      std::array<std::size_t, 4> neighbors;

      for (int face = 0; face < 4; ++face) {
        neighbors[face] = secondaryCellInformation[i].faceNeighborGlobalIds[face];
      }
      cellInformation[i].ltsSetup = determineLtsSetup(cellInformation,
                                                      secondaryCellInformation,
                                                      self,
                                                      neighbors,
                                                      layer.getIdentifier().halo == HaloType::Copy);
    }
  }

  // do sync

  // pass 4: correct LTS setup
  for (auto& layer : container.volume.leaves(Ghost)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < layer.size(); ++i) {
      const std::size_t self = secondaryCellInformation[i].linearId;
      std::array<std::size_t, 4> neighbors;

      for (int face = 0; face < 4; ++face) {
        neighbors[face] = secondaryCellInformation[i].faceNeighborGlobalIds[face];
      }
      cellInformation[i].ltsSetup = correctLtsSetup(cellInformation, self, neighbors);
    }
  }

  // do sync
  for (auto& layer : container.volume.leaves(Interior | Copy)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t i = 0; i < layer.size(); ++i) {
      const std::size_t self = secondaryCellInformation[i].linearId;
      std::array<std::size_t, 4> neighbors;
      for (int face = 0; face < 4; ++face) {
        neighbors[face] = secondaryCellInformation[i].faceNeighborGlobalIds[face];
      }
      cellInformation[i].ltsSetup = correctGhostLayer(cellInformation, self, neighbors);
    }
  }
}
} // namespace seissol::initializer::internal
