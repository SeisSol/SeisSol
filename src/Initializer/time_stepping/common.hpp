/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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

#include <set>
#include <Kernels/Time.h>

namespace seissol {
namespace initializers {
namespace time_stepping {

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
                 MPI_COMM_WORLD,                                                // communicator
                 io_meshStructure[l_cluster].sendRequests+l_region );           // mpi request

      l_copyRegionOffset += io_meshStructure[l_cluster].numberOfCopyRegionCells[l_region];

      MPI_Irecv( l_receiveBuffer[l_cluster]+l_ghostRegionOffset,                 // buffer
                 io_meshStructure[l_cluster].numberOfGhostRegionCells[l_region], // size
                 MPI_UNSIGNED_SHORT,                                             // data type
                 io_meshStructure[l_cluster].neighboringClusters[l_region][0],   // source
                 io_meshStructure[l_cluster].receiveIdentifiers[l_region],       // message tag
                 MPI_COMM_WORLD,                                                 // communicator
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
      unsigned int l_neighboringClusterIds[4];
      // collect cluster ids
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // only continue for valid faces
        if( io_cellLocalInformation[l_cell].faceTypes[l_face] == regular  ||
            io_cellLocalInformation[l_cell].faceTypes[l_face] == periodic ||
            io_cellLocalInformation[l_cell].faceTypes[l_face] == dynamicRupture ) {
          // get neighboring cell id
          unsigned int l_neighbor = io_cellLocalInformation[l_cell].faceNeighborIds[l_face];

          // set the cluster id
          l_neighboringClusterIds[l_face] = io_cellLocalInformation[l_neighbor].clusterId;
        }
      }

      // set the lts setup for this cell
      io_cellLocalInformation[l_cell].ltsSetup = seissol::kernels::Time::getLtsSetup( io_cellLocalInformation[l_cell].clusterId,
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
        if( io_cellLocalInformation[l_cell].faceTypes[l_face] == regular  ||
            io_cellLocalInformation[l_cell].faceTypes[l_face] == periodic ||
            io_cellLocalInformation[l_cell].faceTypes[l_face] == dynamicRupture ) {
          // get neighboring cell id
          unsigned int l_neighbor = io_cellLocalInformation[l_cell].faceNeighborIds[l_face];

          // get neighboring setup
          l_neighboringSetups[l_face] = io_cellLocalInformation[l_neighbor].ltsSetup;
        }
      }

      // do the normalization
      kernels::Time::normalizeLtsSetup( l_neighboringSetups, io_cellLocalInformation[l_cell].ltsSetup );

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

/**
 * Derives the lts setups of all given cells.
 *
 * @param i_numberOfCells number of cells to derive the lts setup for.
 * @param io_cellLocalInformation cell local information which lts setup will be written using present face and cluster.
 * @param i_ghostLayer true if the cells are part of the ghost layer, information is required for corner cases where every cell-neighbor is part of the computational domain.
 *
 * @todo Is inline a good idea? Static does not work!
 **/
inline void deriveLtsSetups(  unsigned int          i_numberOfCells,
                              CellLocalInformation *io_cellLocalInformation,
                              bool                  i_ghostLayer = false ) {
  // iterate over all cells
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell ++ ) {
    // cluster ids of the four face neighbors
    unsigned int l_neighboringClusterIds[4];

    // collect cluster ids
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      // only continue for valid faces
      if( io_cellLocalInformation[l_cell].faceTypes[l_face] == regular  ||
          io_cellLocalInformation[l_cell].faceTypes[l_face] == periodic ||
          io_cellLocalInformation[l_cell].faceTypes[l_face] == dynamicRupture ) {
        // get neighboring cell id
        unsigned int l_neighbor = io_cellLocalInformation[l_cell].faceNeighborIds[l_face];

        // set the cluster id
        l_neighboringClusterIds[l_face] = io_cellLocalInformation[l_neighbor].clusterId;
      }
    }

    // set the lts setup for this cell
    io_cellLocalInformation[l_cell].ltsSetup = seissol::kernels::Time::getLtsSetup( io_cellLocalInformation[l_cell].clusterId,
                                                                                    l_neighboringClusterIds,
                                                                                    io_cellLocalInformation[l_cell].faceTypes,
                                                                                    io_cellLocalInformation[l_cell].faceNeighborIds,
                                                                                    i_ghostLayer );
  }

  // iterate over all cells and normalize the setups
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell ++ ) {
    unsigned short l_neighboringSetups[4];

    // reset to gts
    l_neighboringSetups[0] = l_neighboringSetups[1] = l_neighboringSetups[2] = l_neighboringSetups[3] = 240;

    // collect lts setups
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      // only continue for non-boundary faces
      if( io_cellLocalInformation[l_cell].faceTypes[l_face] == regular  ||
          io_cellLocalInformation[l_cell].faceTypes[l_face] == periodic ||
          io_cellLocalInformation[l_cell].faceTypes[l_face] == dynamicRupture ) {
        // get neighboring cell id
        unsigned int l_neighbor = io_cellLocalInformation[l_cell].faceNeighborIds[l_face];

        // get neighboring setup
        l_neighboringSetups[l_face] = io_cellLocalInformation[l_neighbor].ltsSetup;
      }
    }

    // do the normalization
    kernels::Time::normalizeLtsSetup( l_neighboringSetups, io_cellLocalInformation[l_cell].ltsSetup );
  }
}

}}}

#endif
