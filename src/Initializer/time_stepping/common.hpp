/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Common functions for the setups of all lts strategies.
 **/

#ifndef COMMON_HPP
#define COMMON_HPP

#include <cassert>
#include <set>

#include "Parallel/MPI.h"
#include "Initializer/typedefs.hpp"

namespace seissol {
namespace initializers {
namespace time_stepping {
  
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
 * @param i_localCluster global id of the cluster to which this cell belongs.
 * @param i_neighboringClusterIds global ids of the clusters the face neighbors belong to (if present).
 * @param i_faceTypes types of the four faces.
 * @param i_faceNeighborIds face neighbor ids.
 * @param i_copy true if the cell is part of the copy layer (only required for correctness in dynamic rupture computations).
 **/
static unsigned short getLtsSetup(unsigned int i_localClusterId,
                                  unsigned int i_neighboringClusterIds[4],
                                  const FaceType i_faceTypes[4],
                                  const unsigned int i_faceNeighborIds[4], // TODO: Remove, outdated
                                  bool i_copy = false ) {
  // reset the LTS setup
  unsigned short l_ltsSetup = 0;

  // iterate over the faces
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // continue for boundary conditions
    if (i_faceTypes[l_face] == FaceType::outflow) {
      continue;
    }
    // fake neighbors are GTS
    else if(i_faceTypes[l_face] == FaceType::freeSurface ||
	    i_faceTypes[l_face] == FaceType::freeSurfaceGravity ||
	    i_faceTypes[l_face] == FaceType::dirichlet ||
	    i_faceTypes[l_face] == FaceType::analytical) {
      l_ltsSetup |= (1 << (l_face+4) );
    }
    // dynamic rupture faces are always global time stepping but operate on derivatives
    else if( i_faceTypes[l_face] == FaceType::dynamicRupture ) {
      // face-neighbor provides GTS derivatives
      // face-neighbor provides derivatives
      l_ltsSetup |= ( 1 <<  l_face      );
      l_ltsSetup |= ( 1 << (l_face + 4) );

      // cell is required to provide derivatives for dynamic rupture
      l_ltsSetup |= ( 1 << 9 );

      if( i_copy ) { // set the buffer invalid in copy layers
                     // TODO: Minor improvements possible: Non-DR MPI-neighbor for example
        l_ltsSetup |= ( 1 << 10 );
      }
    }
    // derive the LTS setup based on the cluster ids
    else {
      // neighboring cluster has a larger time step than this cluster
      if( i_localClusterId < i_neighboringClusterIds[l_face] ) {
        // neighbor delivers time derivatives
        l_ltsSetup |= ( 1 << l_face );

        // the cell-local buffer is used in LTS-fashion
        l_ltsSetup |= ( 1 << 10     );
      }
      // GTS relation
      else if( i_localClusterId == i_neighboringClusterIds[l_face] ) {
        l_ltsSetup |= ( 1 << (l_face + 4) );
      }

      // cell is required to provide derivatives
      if( i_localClusterId > i_neighboringClusterIds[l_face] ) {
        l_ltsSetup |= ( 1 << 9 );
      }
      // cell is required to provide a buffer
      else {
        l_ltsSetup |= ( 1 << 8 );
      }
    }

    // true lts buffer with gts required derivatives
    if( (l_ltsSetup >> 10)%2 == 1 && (l_ltsSetup >> 4)%16 != 0 ) {
      l_ltsSetup |= ( 1 << 9 );
    }
  }

  /*
   * Normalize for special case "free surface/dirichlet on derivatives":
   *   If a cell provides either buffers in a LTS fashion or derivatives only,
   *   the neighboring contribution of the boundary intergral is required to work on the cells derivatives.
   *   Free surface/dirichlet boundary conditions work on the cells DOFs in the neighboring contribution:
   *   Enable cell local derivatives in this case and mark that the "fake neighbor" provides derivatives.
   */
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // check for special case free-surface/dirichlet requirements
    const bool isSpecialCase = i_faceTypes[l_face] == FaceType::freeSurface ||
      i_faceTypes[l_face] == FaceType::freeSurfaceGravity ||
      i_faceTypes[l_face] == FaceType::dirichlet ||
      i_faceTypes[l_face] == FaceType::analytical;
    if (isSpecialCase &&       // special case face
       ( (l_ltsSetup >> 10) % 2 == 1 ||             // lts fashion buffer
         (l_ltsSetup >> 8 ) % 2 == 0 )         ) {  // no buffer at all
      l_ltsSetup |= ( 1 << 9 );       // enable derivatives computation
      l_ltsSetup |= ( 1 << l_face);  // enable derivatives for the fake face neighbor
    }
  }

  return l_ltsSetup;
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
 * @param i_neighboringSetups local time stepping setups for the neighboring cells, set to GTS (240) if not defined (e.g. in case of boundary conditions).
 * @param io_localLtsSetup local time stepping setup of the local cell.
 **/
static void normalizeLtsSetup( unsigned short  i_neighboringLtsSetups[4],
                               unsigned short &io_localLtsSetup ) {
  // iterate over the face neighbors
  for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
    // enforce derivatives if this is a "GTS on derivatives" relation
    if( (io_localLtsSetup >> (l_face + 4))%2 && (i_neighboringLtsSetups[l_face] >> 10)%2 == 1 ) {
      io_localLtsSetup |= (1 << l_face);
    }
  }
}

/**
 * Synchronizes the LTS setups in the ghost layers.
 *
 * @param i_numberOfClusters number of clusters
 * @param io_meshStructure mesh structure.
 * @param io_cellLocalInformation cell local information which lts setup will be written using present face and cluster.
 domain.
 */
#ifdef USE_MPI
static void synchronizeLtsSetups( unsigned int                 i_numberOfClusters,
                                  struct MeshStructure        *io_meshStructure,
                                  struct CellLocalInformation *io_cellLocalInformation ) {
  // linear lts setups in the ghost regions
  // 0: local cluster
  // 1: ghost region
  unsigned short **l_sendBuffer    = new unsigned short*[ i_numberOfClusters ];
  unsigned short **l_receiveBuffer = new unsigned short*[ i_numberOfClusters ];

  unsigned l_cell = 0;

  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    l_sendBuffer[l_cluster]    = new unsigned short[ io_meshStructure[l_cluster].numberOfCopyCells  ];
    l_receiveBuffer[l_cluster] = new unsigned short[ io_meshStructure[l_cluster].numberOfGhostCells ];

    // jump over ghost layer
    l_cell += io_meshStructure[l_cluster].numberOfGhostCells;

    // fill the linear buffer
    for( unsigned int l_copyCell = 0; l_copyCell < io_meshStructure[l_cluster].numberOfCopyCells; l_copyCell++ ) {
      l_sendBuffer[l_cluster][l_copyCell] = io_cellLocalInformation[l_cell].ltsSetup;

      l_cell++;
    }

    unsigned int l_copyRegionOffset  = 0;
    unsigned int l_ghostRegionOffset = 0;
    for( unsigned int l_region = 0; l_region < io_meshStructure[l_cluster].numberOfRegions; l_region++ ) {
      // post the communication requests
      MPI_Isend( l_sendBuffer[l_cluster]+l_copyRegionOffset,                    // buffer
                 io_meshStructure[l_cluster].numberOfCopyRegionCells[l_region], // size
                 MPI_UNSIGNED_SHORT,                                            // data type
                 io_meshStructure[l_cluster].neighboringClusters[l_region][0],  // destination
                 io_meshStructure[l_cluster].sendIdentifiers[l_region],         // message tag
                 seissol::MPI::mpi.comm(),                                      // communicator
                 io_meshStructure[l_cluster].sendRequests+l_region );           // mpi request

      l_copyRegionOffset += io_meshStructure[l_cluster].numberOfCopyRegionCells[l_region];

      MPI_Irecv( l_receiveBuffer[l_cluster]+l_ghostRegionOffset,                 // buffer
                 io_meshStructure[l_cluster].numberOfGhostRegionCells[l_region], // size
                 MPI_UNSIGNED_SHORT,                                             // data type
                 io_meshStructure[l_cluster].neighboringClusters[l_region][0],   // source
                 io_meshStructure[l_cluster].receiveIdentifiers[l_region],       // message tag
                 seissol::MPI::mpi.comm(),                                       // communicator
                 io_meshStructure[l_cluster].receiveRequests+l_region );         // mpi request

      l_ghostRegionOffset += io_meshStructure[l_cluster].numberOfGhostRegionCells[l_region];
    }

    // jump over interior
    l_cell +=  io_meshStructure[l_cluster].numberOfInteriorCells;
  }

  // wait for communication
  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    MPI_Waitall( io_meshStructure[l_cluster].numberOfRegions, // size
                 io_meshStructure[l_cluster].receiveRequests, // array of requests
                 MPI_STATUS_IGNORE );                         // mpi status
    MPI_Waitall( io_meshStructure[l_cluster].numberOfRegions, // size
                 io_meshStructure[l_cluster].sendRequests,    // array of requests
                 MPI_STATUS_IGNORE );                         // mpi status
  }

  // update ghost cells with received information
  l_cell = 0;
  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    for( unsigned int l_ghostCell = 0; l_ghostCell < io_meshStructure[l_cluster].numberOfGhostCells; l_ghostCell++ ) {
      io_cellLocalInformation[l_cell].ltsSetup = l_receiveBuffer[l_cluster][l_ghostCell];

      // assert at least derivatives or buffers are offered by the ghost cells
      assert( ( ( io_cellLocalInformation[l_cell].ltsSetup >> 8 ) % 2 ||
                ( io_cellLocalInformation[l_cell].ltsSetup >> 9 ) % 2 )
                == true );

      l_cell++;
    }

    // free memory
    delete[] l_receiveBuffer[l_cluster];
    delete[] l_sendBuffer[l_cluster];

    // jump over copy layer and interior
    l_cell +=   io_meshStructure[l_cluster].numberOfCopyCells
              + io_meshStructure[l_cluster].numberOfInteriorCells;
  }

  delete[] l_sendBuffer;
  delete[] l_receiveBuffer;
}
#endif

/**
 * Derives the lts setups of all given cells.
 *
 * @param i_numberOfCluster number of clusters.
 * @param io_meshStructure mesh structure.
 * @param io_cellLocalInformation cell local information which lts setup will be written using present face and cluster.
 *
 * @todo Is inline a good idea? Static does not work!
 **/
inline void deriveLtsSetups( unsigned int                 i_numberOfClusters,
                             struct MeshStructure        *io_meshStructure,
                             struct CellLocalInformation *io_cellLocalInformation ) {
  unsigned int l_cell = 0;

  // iterate over time clusters
  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    // jump over ghost layer which has invalid information for the face neighbors
    l_cell += io_meshStructure[l_cluster].numberOfGhostCells;

    unsigned int l_numberOfClusterCells = io_meshStructure[l_cluster].numberOfCopyCells  +
                                          io_meshStructure[l_cluster].numberOfInteriorCells;

    // iterate over copy and interior
    for( unsigned int l_clusterCell = 0; l_clusterCell < l_numberOfClusterCells; l_clusterCell++ ) {
      // cluster ids of the four face neighbors
      unsigned int l_neighboringClusterIds[4] = {0};
      // collect cluster ids
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // only continue for valid faces
        if (io_cellLocalInformation[l_cell].faceTypes[l_face] == FaceType::regular ||
           io_cellLocalInformation[l_cell].faceTypes[l_face] == FaceType::periodic ||
           io_cellLocalInformation[l_cell].faceTypes[l_face] == FaceType::dynamicRupture) {
	  // get neighboring cell id
	  unsigned int l_neighbor = io_cellLocalInformation[l_cell].faceNeighborIds[l_face];

          // set the cluster id
          l_neighboringClusterIds[l_face] = io_cellLocalInformation[l_neighbor].clusterId;
        }
      }

      // set the lts setup for this cell
      io_cellLocalInformation[l_cell].ltsSetup = getLtsSetup( io_cellLocalInformation[l_cell].clusterId,
                                                              l_neighboringClusterIds,
                                                              io_cellLocalInformation[l_cell].faceTypes,
                                                              io_cellLocalInformation[l_cell].faceNeighborIds,
                                                              (l_clusterCell < io_meshStructure[l_cluster].numberOfCopyCells) );
      // assert that the cell operates at least on buffers or derivatives
      assert( ( ( io_cellLocalInformation[l_cell].ltsSetup >> 8 ) % 2 ||
                ( io_cellLocalInformation[l_cell].ltsSetup >> 9 ) % 2 )
                == true );

      l_cell++;
    }
  }

// exchange ltsSetup of the ghost layer for the normalization step
#ifdef USE_MPI
  synchronizeLtsSetups( i_numberOfClusters,
                        io_meshStructure,
                        io_cellLocalInformation );
#endif

  // iterate over cells and normalize the setups
  l_cell = 0;
  for( unsigned int l_cluster = 0; l_cluster < i_numberOfClusters; l_cluster++ ) {
    // jump over ghost layer
    l_cell += io_meshStructure[l_cluster].numberOfGhostCells;

    unsigned int l_numberOfClusterCells = io_meshStructure[l_cluster].numberOfCopyCells  +
                                          io_meshStructure[l_cluster].numberOfInteriorCells;

    for( unsigned int l_clusterCell = 0; l_clusterCell < l_numberOfClusterCells; l_clusterCell++ ) {
      unsigned short l_neighboringSetups[4];

      // reset to gts
      l_neighboringSetups[0] = l_neighboringSetups[1] = l_neighboringSetups[2] = l_neighboringSetups[3] = 240;

      // collect lts setups
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // only continue for non-boundary faces
        if( io_cellLocalInformation[l_cell].faceTypes[l_face] == FaceType::regular  ||
            io_cellLocalInformation[l_cell].faceTypes[l_face] == FaceType::periodic ||
            io_cellLocalInformation[l_cell].faceTypes[l_face] == FaceType::dynamicRupture ) {
          // get neighboring cell id
          unsigned int l_neighbor = io_cellLocalInformation[l_cell].faceNeighborIds[l_face];

          // get neighboring setup
          l_neighboringSetups[l_face] = io_cellLocalInformation[l_neighbor].ltsSetup;
        }
      }

      // do the normalization
      normalizeLtsSetup( l_neighboringSetups, io_cellLocalInformation[l_cell].ltsSetup );

      // assert that the cell operates at least on buffers or derivatives
      assert( ( ( io_cellLocalInformation[l_cell].ltsSetup >> 8 ) % 2 ||
                ( io_cellLocalInformation[l_cell].ltsSetup >> 9 ) % 2 )
                == true );

      l_cell++;
    }
  }

// get final setup in the ghost layer (after normalization)
#ifdef USE_MPI
  synchronizeLtsSetups( i_numberOfClusters,
                        io_meshStructure,
                        io_cellLocalInformation );
#endif
}

}}}

#endif
