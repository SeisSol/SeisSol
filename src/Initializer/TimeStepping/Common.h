// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOSRC_INITIALIZER_TIMESTEPPING_COMMON_H_
#define SEISSOSRC_INITIALIZER_TIMESTEPPING_COMMON_H_

#include <Initializer/CellLocalInformation.h>
#include <Initializer/LtsSetup.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/LTSTree.h>
#include <cassert>
#include <mpi.h>
#include <set>

#include "Parallel/MPI.h"
#include "Initializer/Typedefs.h"



namespace seissol::initializer::time_stepping {
  
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
 *  1 in one of bits 4 - 7: The  current cell has a global time stepping relation with the face neighbor.
 *
 *     Example 2:
 *     [       8 rem. bits       |   GTS bits  | buf/der bits ]
 *     [  -  -  -  -  -  -  -  - | 0  0  1  1  |  0  1  1  0  ]
 *     [ 15 14 13 12 11 10  9  8 | 7  6  5  4  |  3  2  1  0  ]
 *  In Example 2 the data of face neighbors 0 and 3 are buffers, 1 and 2 deliver derivatives
 *  Face neighbor 0 has a GTS-relation and this cell works directly on the delivered buffer.
 *  Face neighbor 1 has a GTS-relation, but delivers derivatives -> The derivatives have to translated to time integrated DOFs first.
 *  Face neighbor 2 has a LTS-relation and receives derivatives from its neighbor -> The derivates have to be used for a partial time integration.
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
 *  In Example 5 the buffer is a LTS buffer (reset on request only). GTS buffers are updated in every time step.
 *
 * @return lts setup.
 * @param localCluster global id of the cluster to which this cell belongs.
 * @param neighboringClusterIds global ids of the clusters the face neighbors belong to (if present).
 * @param faceTypes types of the four faces.
 * @param faceNeighborIds face neighbor ids.
 * @param copy true if the cell is part of the copy layer (only required for correctness in dynamic rupture computations).
 **/
static LtsSetup getLtsSetup(const CellLocalInformation& ownPrimary,
  const SecondaryCellLocalInformation& ownSecondary,
  const std::array<uint64_t, Cell::NumFaces>& neighborClusters,
                                  bool copy = false ) {
  // reset the LTS setup
  LtsSetup ltsSetup{};

  // iterate over the faces
  for( std::size_t face = 0; face < Cell::NumFaces; face++ ) {
    // continue for boundary conditions
    if (ownPrimary.faceTypes[face] == FaceType::Outflow) {
      continue;
    }
    // fake neighbors are GTS
    else if(isExternalBoundaryFaceType(ownPrimary.faceTypes[face])) {
      ltsSetup.setNeighborGTS(face, true);
    }
    // dynamic rupture faces are always global time stepping but operate on derivatives
    else if( ownPrimary.faceTypes[face] == FaceType::DynamicRupture ) {
      // face-neighbor provides GTS derivatives
      // face-neighbor provides derivatives
      ltsSetup.setNeighborHasDerivatives(face, true);
      ltsSetup.setNeighborGTS(face, true);

      // cell is required to provide derivatives for dynamic rupture
      ltsSetup.setHasDerivatives(true);

      if( copy ) { // set the buffer invalid in copy layers
                     // TODO: Minor improvements possible: Non-DR MPI-neighbor for example
        ltsSetup.setCacheBuffers(true);
      }
    }
    // derive the LTS setup based on the cluster ids
    else {
      // neighboring cluster has a larger time step than this cluster
      if( ownSecondary.clusterId < neighborClusters[face] ) {
        // neighbor delivers time derivatives
        ltsSetup.setNeighborHasDerivatives(face, true);

        // the cell-local buffer is used in LTS-fashion
        ltsSetup.setCacheBuffers(true);
      }
      // GTS relation
      else if( ownSecondary.clusterId == neighborClusters[face] ) {
        ltsSetup.setNeighborGTS(face, true);
      }

      // cell is required to provide derivatives
      if( ownSecondary.clusterId > neighborClusters[face] ) {
        ltsSetup.setHasDerivatives(true);
      }
      // cell is required to provide a buffer
      else {
        ltsSetup.setHasBuffers(true);
      }
    }

    // true lts buffer with gts required derivatives
    bool hasGTS = false;
    for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
      hasGTS |= ltsSetup.neighborGTS(i);
    }
    if(ltsSetup.cacheBuffers() && hasGTS) {
      ltsSetup.setHasDerivatives(true);
    }
  }

  /*
   * Normalize for special case "free surface/dirichlet on derivatives":
   *   If a cell provides either buffers in a LTS fashion or derivatives only,
   *   the neighboring contribution of the boundary intergral is required to work on the cells derivatives.
   *   Free surface/dirichlet boundary conditions work on the cells DOFs in the neighboring contribution:
   *   Enable cell local derivatives in this case and mark that the "fake neighbor" provides derivatives.
   */
  for( std::size_t face = 0; face < Cell::NumFaces; face++ ) {
    // check for special case free-surface/dirichlet requirements
    const bool isSpecialCase = isExternalBoundaryFaceType(ownPrimary.faceTypes[face]);
    if (isSpecialCase &&       // special case face
       ( ltsSetup.cacheBuffers() ||             // lts fashion buffer
         !ltsSetup.hasBuffers() )         ) {  // no buffer at all
      ltsSetup.setHasDerivatives(true);       // enable derivatives computation
      ltsSetup.setNeighborHasDerivatives(face, true);  // enable derivatives for the fake face neighbor
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
 *
 * @param neighboringSetups local time stepping setups for the neighboring cells, set to GTS (240) if not defined (e.g. in case of boundary conditions).
 * @param localLtsSetup local time stepping setup of the local cell.
 **/
static LtsSetup normalizeLtsSetup( const LtsSetup &localLtsSetup,
                                    const std::array<bool, Cell::NumFaces>&  neighborCache) {
  LtsSetup output(localLtsSetup);
  // iterate over the face neighbors
  for( std::size_t face = 0; face < Cell::NumFaces; face++ ) {
    // enforce derivatives if this is a "GTS on derivatives" relation
    if( localLtsSetup.neighborGTS(face) && neighborCache[face] ) {
      output.setNeighborHasDerivatives(face, true);
    }
  }
  return output;
}

/**
 * Synchronizes the LTS setups in the ghost layers.
 *
 * @param numberOfClusters number of clusters
 * @param meshStructure mesh structure.
 * @param cellLocalInformation cell local information which lts setup will be written using present face and cluster.
 domain.
 */
static void synchronizeLtsSetups( unsigned int                 numberOfClusters,
                                  struct MeshStructure        *meshStructure,
                                  struct CellLocalInformation *cellLocalInformation ) {
  // linear lts setups in the ghost regions
  // 0: local cluster
  // 1: ghost region
  std::vector<std::vector<uint16_t>> sendBuffer(numberOfClusters);
  std::vector<std::vector<uint16_t>> receiveBuffer(numberOfClusters);

  std::vector<MPI_Request> sendRequestsRaw;
  std::vector<MPI_Request> recvRequestsRaw;
  std::vector<MPI_Request*> sendRequests;
  std::vector<MPI_Request*> receiveRequests;
  for( unsigned int cluster = 0; cluster < numberOfClusters; cluster++ ) {
    for( unsigned int region = 0; region < meshStructure[cluster].numberOfRegions; region++ ) {
      sendRequestsRaw.push_back(MPI_REQUEST_NULL);
      recvRequestsRaw.push_back(MPI_REQUEST_NULL);
    }
  }
  std::size_t point = 0;
  for( unsigned int cluster = 0; cluster < numberOfClusters; cluster++ ) {
    sendRequests.push_back(sendRequestsRaw.data() + point);
    receiveRequests.push_back(recvRequestsRaw.data() + point);
    point += meshStructure[cluster].numberOfRegions;
  }

  std::size_t cell = 0;

  for( unsigned int cluster = 0; cluster < numberOfClusters; cluster++ ) {
    sendBuffer[cluster].resize(meshStructure[cluster].numberOfCopyCells);
    receiveBuffer[cluster].resize(meshStructure[cluster].numberOfGhostCells);

    // jump over ghost layer
    cell += meshStructure[cluster].numberOfGhostCells;

    // fill the linear buffer
    for( unsigned int copyCell = 0; copyCell < meshStructure[cluster].numberOfCopyCells; copyCell++ ) {
      sendBuffer[cluster][copyCell] = cellLocalInformation[cell].ltsSetup.unwrap();

      cell++;
    }

    unsigned int copyRegionOffset  = 0;
    unsigned int ghostRegionOffset = 0;
    for( unsigned int region = 0; region < meshStructure[cluster].numberOfRegions; region++ ) {
      // post the communication requests
      MPI_Isend( sendBuffer[cluster].data()+copyRegionOffset,                    // buffer
                 meshStructure[cluster].numberOfCopyRegionCells[region], // size
                 MPI_UINT16_T,                                            // data type
                 meshStructure[cluster].neighboringClusters[region][0],  // destination
                 meshStructure[cluster].sendIdentifiers[region],         // message tag
                 seissol::MPI::mpi.comm(),                                      // communicator
                 sendRequests[cluster]+region );           // mpi request

      copyRegionOffset += meshStructure[cluster].numberOfCopyRegionCells[region];

      MPI_Irecv( receiveBuffer[cluster].data()+ghostRegionOffset,                 // buffer
                 meshStructure[cluster].numberOfGhostRegionCells[region], // size
                 MPI_UINT16_T,                                             // data type
                 meshStructure[cluster].neighboringClusters[region][0],   // source
                 meshStructure[cluster].receiveIdentifiers[region],       // message tag
                 seissol::MPI::mpi.comm(),                                       // communicator
                 receiveRequests[cluster]+region );         // mpi request

      ghostRegionOffset += meshStructure[cluster].numberOfGhostRegionCells[region];
    }

    // jump over interior
    cell +=  meshStructure[cluster].numberOfInteriorCells;
  }

  // wait for communication
  for( unsigned int cluster = 0; cluster < numberOfClusters; cluster++ ) {
    MPI_Waitall( meshStructure[cluster].numberOfRegions, // size
                 sendRequests[cluster], // array of requests
                 MPI_STATUS_IGNORE );                         // mpi status
    MPI_Waitall( meshStructure[cluster].numberOfRegions, // size
                 receiveRequests[cluster],    // array of requests
                 MPI_STATUS_IGNORE );                         // mpi status
  }

  // update ghost cells with received information
  cell = 0;
  for( unsigned int cluster = 0; cluster < numberOfClusters; cluster++ ) {
    for( unsigned int ghostCell = 0; ghostCell < meshStructure[cluster].numberOfGhostCells; ghostCell++ ) {
      cellLocalInformation[cell].ltsSetup = LtsSetup(receiveBuffer[cluster][ghostCell]);

      // assert at least derivatives or buffers are offered by the ghost cells
      assert( ( ( cellLocalInformation[cell].ltsSetup.unwrap() >> 8 ) % 2 ||
                ( cellLocalInformation[cell].ltsSetup.unwrap() >> 9 ) % 2 )
                == true );

      cell++;
    }

    // jump over copy layer and interior
    cell +=   meshStructure[cluster].numberOfCopyCells
              + meshStructure[cluster].numberOfInteriorCells;
  }
}

/**
 * Derives the lts setups of all given cells.
 *
 * @param numberOfCluster number of clusters.
 * @param meshStructure mesh structure.
 * @param cellLocalInformation cell local information which lts setup will be written using present face and cluster.
 *
 * @todo Is inline a good idea? Static does not work!
 **/
inline void deriveLtsSetups( unsigned int                 numberOfClusters,
                             struct MeshStructure        *meshStructure,
                            LTS::Storage& storage  ) {

  auto* primaryInformationGlobal = storage.var<LTS::CellInformation>();
  const auto* secondaryInformationGlobal = storage.var<LTS::SecondaryInformation>();

  // iterate over time clusters
  for(auto& layer : storage.leaves(Ghost)) {
    const auto isCopy = layer.getIdentifier().halo == HaloType::Copy;
    auto* primaryInformationLocal = layer.var<LTS::CellInformation>();
    const auto* secondaryInformationLocal = layer.var<LTS::SecondaryInformation>();
    for( std::size_t cell = 0; cell < layer.size(); ++cell ) {
      std::array<uint64_t, Cell::NumFaces> neighborClusters{};
      for( std::size_t face = 0; face < Cell::NumFaces; face++ ) {
        // only continue for non-boundary faces
        if( isInternalFaceType(primaryInformationLocal[cell].faceTypes[face]) ) {
          // get neighboring cell id
          const auto neighbor = secondaryInformationLocal[cell].faceNeighborIds[face];

          // get neighboring setup
          neighborClusters[face] = secondaryInformationGlobal[neighbor].clusterId;
        }
      }

      // set the lts setup for this cell
      primaryInformationLocal[cell].ltsSetup = LtsSetup(getLtsSetup( primaryInformationLocal[cell], secondaryInformationLocal[cell],
                                                              neighborClusters,
                                                              isCopy ));
      // assert that the cell operates at least on buffers or derivatives
      assert (primaryInformationLocal[cell].ltsSetup.hasBuffers() || primaryInformationLocal[cell].ltsSetup.hasDerivatives());
    }
  }

// exchange ltsSetup of the ghost layer for the normalization step
  synchronizeLtsSetups( numberOfClusters,
                        meshStructure,
                        primaryInformationGlobal );

  // iterate over cells and normalize the setups
  for(auto& layer : storage.leaves(Ghost)) {
    const auto isCopy = layer.getIdentifier().halo == HaloType::Copy;
    auto* primaryInformationLocal = layer.var<LTS::CellInformation>();
    const auto* secondaryInformationLocal = layer.var<LTS::SecondaryInformation>();
    for( std::size_t cell = 0; cell < layer.size(); ++cell ) {
      std::array<bool, Cell::NumFaces> neighborCache{};

      // collect lts setups
      for( std::size_t face = 0; face < Cell::NumFaces; face++ ) {
        // only continue for non-boundary faces
        if( isInternalFaceType(primaryInformationLocal[cell].faceTypes[face]) ) {
          const auto neighbor = secondaryInformationLocal[cell].faceNeighborIds[face];
          neighborCache[face] = primaryInformationGlobal[neighbor].ltsSetup.cacheBuffers();
        }
      }

      primaryInformationLocal[cell].ltsSetup = normalizeLtsSetup(primaryInformationLocal[cell].ltsSetup,  neighborCache );

      // assert that the cell operates at least on buffers or derivatives
      assert (primaryInformationLocal[cell].ltsSetup.hasBuffers() || primaryInformationLocal[cell].ltsSetup.hasDerivatives());
    }
  }

// get final setup in the ghost layer (after normalization)
  synchronizeLtsSetups( numberOfClusters,
                        meshStructure,
                        primaryInformationGlobal );
}

} // namespace seissol::initializer::time_stepping


#endif // SEISSOSRC_INITIALIZER_TIMESTEPPING_COMMON_H_
