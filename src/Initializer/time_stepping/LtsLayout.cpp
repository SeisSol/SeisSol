/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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
 * Layoyt the LTS schemes compute on.
 **/

#include "Parallel/MPI.h"

#include "utils/logger.h"

#include "LtsLayout.h"
#include "MultiRate.hpp"
#include <iterator>

seissol::initializers::time_stepping::LtsLayout::LtsLayout():
 m_cellTimeStepWidths(       NULL ),
 m_cellClusterIds(           NULL ),
 m_globalTimeStepWidths(     NULL ),
 m_globalTimeStepRates(      NULL ),
 m_plainCopyRegions(         NULL ),
 m_numberOfPlainGhostCells(  NULL ),
 m_plainGhostCellClusterIds( NULL ) {}

seissol::initializers::time_stepping::LtsLayout::~LtsLayout() {
  // free memory of member variables
  if( m_cellTimeStepWidths       != NULL ) delete[] m_cellTimeStepWidths;
  if( m_cellClusterIds           != NULL ) delete[] m_cellClusterIds;
  if( m_globalTimeStepWidths     != NULL ) delete[] m_globalTimeStepWidths;
  if( m_globalTimeStepRates      != NULL ) delete[] m_globalTimeStepRates;
  if( m_numberOfPlainGhostCells  != NULL ) delete[] m_numberOfPlainGhostCells;
  if( m_plainGhostCellClusterIds != NULL ) for( unsigned int l_rank = 0; l_rank < m_plainNeighboringRanks.size(); l_rank++ ) delete[] m_plainGhostCellClusterIds[l_rank];
  if( m_plainGhostCellClusterIds != NULL ) delete[] m_plainGhostCellClusterIds;
  if( m_plainCopyRegions         != NULL ) delete[] m_plainCopyRegions;
}

void seissol::initializers::time_stepping::LtsLayout::setMesh( const MeshReader &i_mesh ) {
  // TODO: remove the copy by a pointer once the mesh stays constant
  std::vector<Element> const& elements = i_mesh.getElements();
  for( unsigned int l_cell = 0; l_cell < elements.size(); l_cell++ ) {
    m_cells.push_back( elements[l_cell] );
  }
  std::vector<Fault> const& fault = i_mesh.getFault();
  for( unsigned int i = 0; i < fault.size(); ++i ) {
    m_fault.push_back( fault[i] );
  }

  m_cellTimeStepWidths = new double[       m_cells.size() ];
  m_cellClusterIds     = new unsigned int[ m_cells.size() ];

  // initialize with invalid values
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    m_cellTimeStepWidths[l_cell] = std::numeric_limits<double>::min();
    m_cellClusterIds[l_cell]     = std::numeric_limits<unsigned int>::max();
  }
}

void seissol::initializers::time_stepping::LtsLayout::setTimeStepWidth( unsigned int i_cellId,
                                                                        double       i_timeStepWidth ) {
  if( i_cellId >= m_cells.size() ) logError() << "cell id >= mesh size: " << i_cellId << m_cells.size() << "aborting";

  // set time step width
  m_cellTimeStepWidths[i_cellId] = i_timeStepWidth;
}

faceType seissol::initializers::time_stepping::LtsLayout::getFaceType( int i_meshFaceType ) {
  if(      i_meshFaceType == 0 ) return regular;
  else if( i_meshFaceType == 1 ) return freeSurface;
  else if( i_meshFaceType == 3 ) return dynamicRupture;
  else if( i_meshFaceType == 5 ) return outflow;
  else if( i_meshFaceType == 6 ) return periodic;
  else logError() << "face type" << i_meshFaceType << "not supported.";
  return regular;
}

void seissol::initializers::time_stepping::LtsLayout::derivePlainCopyInterior() {
	const int rank = seissol::MPI::mpi.rank();

  // unique set of neighboring ranks
  std::set< int > l_neighboringRanks;

  // derive neighboring ranks
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      if(  m_cells[l_cell].neighborRanks[l_face] != rank ) {
        l_neighboringRanks.insert( m_cells[l_cell].neighborRanks[l_face] ) ;
      }
    }
  }

  // convert set to vector
  m_plainNeighboringRanks = std::vector< int > ( l_neighboringRanks.begin(), l_neighboringRanks.end() );

  // allocate data structure for the copy and ghost layer
  m_plainCopyRegions  = new std::vector< unsigned int >[ m_plainNeighboringRanks.size() ];

  // derive copy regions (split by ranks alone) and interior
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    bool l_copyCell = false;
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      if(  m_cells[l_cell].neighborRanks[l_face] != rank ) {
        // derive id of the local region
        int l_region = getPlainRegion( m_cells[l_cell].neighborRanks[l_face] );

        // assert cell ids are increasing
        assert(  m_plainCopyRegions[l_region].size() == 0 ||
               *(m_plainCopyRegions[l_region].end()-1) <= l_cell );

        // put the cell to the right copy region (if not present already)
        if(  m_plainCopyRegions[l_region].size() == 0 ||
           *(m_plainCopyRegions[l_region].end()-1) != l_cell ) {
          m_plainCopyRegions[l_region].push_back( l_cell );
        }

        l_copyCell = true;
      }
    }

    if( !l_copyCell ) m_plainInterior.push_back( l_cell );
  }
}

void seissol::initializers::time_stepping::LtsLayout::derivePlainGhost() {
  /*
   * Get sizes of ghost regions.
   */
  // number of copy cells
  unsigned int *l_numberOfCopyCells = new unsigned int[ m_plainNeighboringRanks.size() ];

  // number of ghost cells
  m_numberOfPlainGhostCells = new unsigned int[ m_plainNeighboringRanks.size() ];

#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  MPI_Request *l_requests = new MPI_Request[ m_plainNeighboringRanks.size() * 2 ];

  for( unsigned int l_neighbor = 0; l_neighbor < m_plainNeighboringRanks.size(); l_neighbor++ ) {
    // get number of copy cells
    l_numberOfCopyCells[l_neighbor] = m_plainCopyRegions[l_neighbor].size();

    // send this ranks copy layer size
    MPI_Isend( l_numberOfCopyCells+l_neighbor,                        // buffer
               1,                                                     // size
               MPI_UNSIGNED,                                          // data type
               m_plainNeighboringRanks[l_neighbor],                   // destination
               deriveGhostPlain,                                      // message tag
               seissol::MPI::mpi.comm(),                              // communicator
               l_requests + l_neighbor );                             // mpi request

    // receive neighboring ranks copy layer size
    MPI_Irecv( m_numberOfPlainGhostCells+l_neighbor,                       // buffer
               1,                                                          // size
               MPI_UNSIGNED,                                               // data type
               m_plainNeighboringRanks[l_neighbor],                        // source
               deriveGhostPlain,                                           // message tag
               seissol::MPI::mpi.comm(),                                   // communicator
               l_requests + l_neighbor + m_plainNeighboringRanks.size() ); // mpi request
  }

  // wait for sends/receives
  MPI_Waitall( m_plainNeighboringRanks.size()*2, // size
               l_requests,                       // array of requests
               MPI_STATUS_IGNORE );              // mpi status

  delete[] l_requests;
#endif // USE_MPI

  // free memory
  delete[] l_numberOfCopyCells;
}

void seissol::initializers::time_stepping::LtsLayout::deriveDynamicRupturePlainCopyInterior()
{
  m_dynamicRupturePlainInterior.resize( m_localClusters.size() );
  m_dynamicRupturePlainCopy.resize(     m_localClusters.size() );
  for (unsigned face = 0; face < m_fault.size(); ++face) {
    int meshId = (m_fault[face].element >= 0) ? m_fault[face].element : m_fault[face].neighborElement;
    unsigned localCluster = getLocalClusterId( m_cellClusterIds[meshId] );
    
    assert(localCluster < m_localClusters.size());
    
    // Local dynamic rupture face
    if (m_fault[face].element >= 0 && m_fault[face].neighborElement >= 0) {
      m_dynamicRupturePlainInterior[localCluster].push_back(face);
    // Dynamic rupture face with one neighbour in the ghost layer
    } else {
      m_dynamicRupturePlainCopy[localCluster].push_back(face);
    }
  }

  int* localClusterHistogram = new int[m_numberOfGlobalClusters];
  for (unsigned gc = 0; gc < m_numberOfGlobalClusters; ++gc) {
    localClusterHistogram[gc] = 0;
  }
  for (unsigned cluster = 0; cluster < m_localClusters.size(); ++cluster) {
    unsigned gc = m_localClusters[cluster];
    localClusterHistogram[gc] = m_dynamicRupturePlainInterior[cluster].size() + m_dynamicRupturePlainCopy[cluster].size();
  }

  const int rank = seissol::MPI::mpi.rank();
  int* globalClusterHistogram = NULL;
#ifdef USE_MPI
  if (rank == 0) {
    globalClusterHistogram = new int[m_numberOfGlobalClusters];
  }
  MPI_Reduce(localClusterHistogram, globalClusterHistogram, m_numberOfGlobalClusters, MPI_INT, MPI_SUM, 0, seissol::MPI::mpi.comm());
#else
  globalClusterHistogram = localClusterHistogram;
#endif
  if (rank == 0) {
    logInfo(rank) << "Number of elements in dynamic rupture time clusters:";
    for (unsigned cluster = 0; cluster < m_numberOfGlobalClusters; ++cluster) {
      logInfo(rank) << utils::nospace << cluster << " (dr):" << utils::space << globalClusterHistogram[cluster];
    }
#ifdef USE_MPI
    delete[] globalClusterHistogram;
#endif
  }
  delete[] localClusterHistogram;
}

void seissol::initializers::time_stepping::LtsLayout::normalizeMpiIndices() {
	const int rank = seissol::MPI::mpi.rank();

  /*
   * Derive local mappings
   */
  // mapping from local mpi-faces to cell ids
  std::vector< std::vector< unsigned int > > l_faceToCellIdMappings;
  l_faceToCellIdMappings.resize( m_plainNeighboringRanks.size() );

  // iterate over mesh and derive mapping
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      if(  m_cells[l_cell].neighborRanks[l_face] != rank ) {
        // asert this is a regular, dynamic rupture or periodic face
        assert( getFaceType( m_cells[l_cell].boundaries[l_face] ) == regular        ||
                getFaceType( m_cells[l_cell].boundaries[l_face] ) == dynamicRupture ||
                getFaceType( m_cells[l_cell].boundaries[l_face] ) == periodic   );

        // derive id of the local region
        int l_region = getPlainRegion( m_cells[l_cell].neighborRanks[l_face] );

        // get unique mpi face id
        int l_mpiIndex = m_cells[l_cell].mpiIndices[l_face];

        // resize mapping array to hold this index if required
        if( l_faceToCellIdMappings[l_region].size() <= static_cast<unsigned int>(l_mpiIndex) ) {
          l_faceToCellIdMappings[l_region].resize( l_mpiIndex+1 );
        }

        // add the face to the mapping of this region
        l_faceToCellIdMappings[l_region][l_mpiIndex] = l_cell;
      }
    }
  }

  /*
   * Exchange and check sizes (debugging only)
   */
  // sizes of the mappings for the local and the neighboring ranks
  std::vector< unsigned int > l_localMappingSizes, l_remoteMappingSizes;
  l_localMappingSizes.resize(  m_plainNeighboringRanks.size() );
  l_remoteMappingSizes.resize( m_plainNeighboringRanks.size() );

#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  // exchange information about the size of the mappings (for debugging purposes only)
  MPI_Request *l_requests = new MPI_Request[ m_plainNeighboringRanks.size() * 2 ];
  for( unsigned int l_neighbor = 0; l_neighbor < m_plainNeighboringRanks.size(); l_neighbor++ ) {
    l_localMappingSizes[l_neighbor] = l_faceToCellIdMappings[l_neighbor].size();

    // send this ranks copy layer size
    MPI_Isend( &l_localMappingSizes[l_neighbor],     // buffer
                1,                                   // size
                MPI_UNSIGNED,                        // data type
                m_plainNeighboringRanks[l_neighbor], // destination
                normalizeIndices,                    // message tag
                seissol::MPI::mpi.comm(),            // communicator
                l_requests + l_neighbor );           // mpi request

    // receive neighboring ranks copy layer size
    MPI_Irecv( &l_remoteMappingSizes[l_neighbor],                            // buffer
                 1,                                                          // size
                 MPI_UNSIGNED,                                               // data type
                 m_plainNeighboringRanks[l_neighbor],                        // source
                 normalizeIndices,                                           // message tag
                 seissol::MPI::mpi.comm(),                                   // communicator
                 l_requests + l_neighbor + m_plainNeighboringRanks.size() ); // mpi request
  }

  // wait for sends/receives
  MPI_Waitall( m_plainNeighboringRanks.size()*2, // size
               l_requests,                       // array of requests
               MPI_STATUS_IGNORE );              // mpi status
#endif

  // make sure the sizes of the local mapping and neighboring mapping match
  for( unsigned int l_region = 0; l_region < m_plainNeighboringRanks.size(); l_region++ ) {
    if( l_faceToCellIdMappings[l_region].size() != l_remoteMappingSizes[l_region] ) {
      logError() << "mapping sizes don't match" << l_faceToCellIdMappings[l_region].size() << l_remoteMappingSizes[l_region];
    }
  }

  /*
   * Exchange the mappings
   */
  // remote mappings
  std::vector< std::vector< unsigned int > > l_remoteFaceToCellIdMappings;
  l_remoteFaceToCellIdMappings.resize( m_plainNeighboringRanks.size() );
  for( unsigned int l_region = 0; l_region < m_plainNeighboringRanks.size(); l_region++ ) {
    l_remoteFaceToCellIdMappings[l_region].resize( l_faceToCellIdMappings[l_region].size() );
  }

#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  // exchange the mappings
  for( unsigned int l_neighbor = 0; l_neighbor < m_plainNeighboringRanks.size(); l_neighbor++ ) {
    // send this ranks copy layer size
    MPI_Isend( &l_faceToCellIdMappings[l_neighbor][0],     // buffer
                l_faceToCellIdMappings[l_neighbor].size(), // size
                MPI_UNSIGNED,                              // data type
                m_plainNeighboringRanks[l_neighbor],       // destination
                normalizeIndices,                          // message tag
                seissol::MPI::mpi.comm(),                  // communicator
                l_requests + l_neighbor );                 // mpi request

    // receive neighboring ranks copy layer size
    MPI_Irecv( &l_remoteFaceToCellIdMappings[l_neighbor][0],                // buffer
                l_remoteFaceToCellIdMappings[l_neighbor].size(),            // size
                MPI_UNSIGNED,                                               // data type
                m_plainNeighboringRanks[l_neighbor],                        // source
                normalizeIndices,                                           // message tag
                seissol::MPI::mpi.comm(),                                   // communicator
                l_requests + l_neighbor + m_plainNeighboringRanks.size() ); // mpi request
  }

  // wait for sends/receives
  MPI_Waitall( m_plainNeighboringRanks.size()*2, // size
               l_requests,                       // array of requests
               MPI_STATUS_IGNORE );              // mpi status
#endif // USE_MPI

  /*
   * Replace the useless mpi-indices by the neighboring cell id
   */
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      if( m_cells[l_cell].neighborRanks[l_face] != rank ) {
        // derive id of the local region
        int l_region = getPlainRegion( m_cells[l_cell].neighborRanks[l_face] );

        // assert we have a corresponding mappiong
        assert( m_cells[l_cell].mpiIndices[l_face] < static_cast<int>(l_remoteFaceToCellIdMappings[l_region].size()) );

        // replace mpi index by the cell id
        m_cells[l_cell].mpiIndices[l_face] = l_remoteFaceToCellIdMappings[l_region][ m_cells[l_cell].mpiIndices[l_face] ];
      }
    }
  }

  /*
   * Convert the neighboring mappings to unique and sorted lists of the neighbors
   */
  for( unsigned int l_region = 0; l_region < m_plainNeighboringRanks.size(); l_region++ ) {
    // sort
    std::sort( l_remoteFaceToCellIdMappings[l_region].begin(), l_remoteFaceToCellIdMappings[l_region].end() );
    // unique
    std::vector< unsigned int >::iterator l_overhead = std::unique( l_remoteFaceToCellIdMappings[l_region].begin(), l_remoteFaceToCellIdMappings[l_region].end() );
    l_remoteFaceToCellIdMappings[l_region].erase( l_overhead, l_remoteFaceToCellIdMappings[l_region].end() );

    // store the results
    m_plainGhostCellIds.push_back( l_remoteFaceToCellIdMappings[l_region] );
  }

  /*
   * Replace neighboring cell id mpi indices by plain ghost region indices.
   */
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      if( m_cells[l_cell].neighborRanks[l_face] != rank ) {
        // derive id of the local region
        int l_region = getPlainRegion( m_cells[l_cell].neighborRanks[l_face] );

        // get ghost index
        std::vector<unsigned int>::iterator l_ghostIdIterator = std::find( m_plainGhostCellIds[l_region].begin(), m_plainGhostCellIds[l_region].end(), (unsigned int) m_cells[l_cell].mpiIndices[l_face] );
        unsigned int l_ghostId = std::distance( m_plainGhostCellIds[l_region].begin(), l_ghostIdIterator );

        // assert a match
        assert( l_ghostId < m_plainGhostCellIds[l_region].size() );

        // replace cell id with ghost index
        m_cells[l_cell].mpiIndices[l_face] = l_ghostId;
      }
    }
  }
}

void seissol::initializers::time_stepping::LtsLayout::synchronizePlainGhostData(unsigned* cellData, unsigned** plainGhostData) {
  // buffer for copy cell cluster ids
  std::vector< std::vector< unsigned int > > l_copyBuffer;

  // fill copy buffer
  l_copyBuffer.resize( m_plainNeighboringRanks.size() );
  for( unsigned int l_region = 0; l_region < m_plainNeighboringRanks.size(); l_region++ ) {
    for( unsigned int l_copyCell = 0; l_copyCell < m_plainCopyRegions[l_region].size(); l_copyCell++ ) {
      l_copyBuffer[l_region].push_back( cellData[ m_plainCopyRegions[l_region][l_copyCell] ] );
    }
  }

#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  // mpi requests
  MPI_Request *l_requests = new MPI_Request[ m_plainNeighboringRanks.size() * 2 ];

  // exchange copy/ghost cluster ids
  for( unsigned int l_region = 0; l_region < m_plainNeighboringRanks.size(); l_region++ ) {
    // send this ranks copy layer size
    MPI_Isend( &l_copyBuffer[l_region][0],         // buffer
                l_copyBuffer[l_region].size(),     // size
                MPI_UNSIGNED,                      // data type
                m_plainNeighboringRanks[l_region], // destination
                synchronizeClusters,               // message tag
                seissol::MPI::mpi.comm(),          // communicator
                l_requests + l_region );           // mpi request

    // receive neighboring ranks copy layer size
    MPI_Irecv( plainGhostData[l_region],                     // buffer
               m_numberOfPlainGhostCells[l_region],                      // size
               MPI_UNSIGNED,                                             // data type
               m_plainNeighboringRanks[l_region],                        // source
               synchronizeClusters,                                      // message tag
               seissol::MPI::mpi.comm(),                                 // communicator
               l_requests + l_region + m_plainNeighboringRanks.size() ); // mpi request
  }

  // wait for sends/receives
  MPI_Waitall( m_plainNeighboringRanks.size()*2, // size
               l_requests,                       // array of requests
               MPI_STATUS_IGNORE );              // mpi status

  delete[] l_requests;
#endif // USE_MPI
}

void seissol::initializers::time_stepping::LtsLayout::synchronizePlainGhostClusterIds() {
  synchronizePlainGhostData(m_cellClusterIds, m_plainGhostCellClusterIds);
}

unsigned seissol::initializers::time_stepping::LtsLayout::enforceDynamicRuptureGTS() {
  const int rank = seissol::MPI::mpi.rank();
  unsigned reductions = 0;
  
  for( std::vector<Fault>::const_iterator fault = m_fault.begin(); fault < m_fault.end(); ++fault ) {
    int meshId, face;
    if (fault->element >= 0) {
      meshId = fault->element;
      face = fault->side;
    } else {
      meshId = fault->neighborElement;
      face = fault->neighborSide;
    }
    if (m_cells[meshId].neighborRanks[face] == rank ) {
      unsigned neighborId = m_cells[meshId].neighbors[face];
      if (m_cellClusterIds[meshId] != m_cellClusterIds[neighborId]) {
        unsigned minCluster = std::min(m_cellClusterIds[meshId], m_cellClusterIds[neighborId]);
        m_cellClusterIds[meshId]     = minCluster;
        m_cellClusterIds[neighborId] = minCluster;
        ++reductions;
      }
    } else {
      unsigned region = getPlainRegion( m_cells[meshId].neighborRanks[face] );
      unsigned localGhostCell = m_cells[meshId].mpiIndices[face];
      assert( localGhostCell < m_numberOfPlainGhostCells[region] );
      if (m_cellClusterIds[meshId] > m_plainGhostCellClusterIds[region][localGhostCell]) {
        m_cellClusterIds[meshId] = m_plainGhostCellClusterIds[region][localGhostCell];
        ++reductions;
      }
    }
  }
  
  return reductions;
}

unsigned int seissol::initializers::time_stepping::LtsLayout::enforceMaximumDifference( unsigned int i_difference ) {
	const int rank = seissol::MPI::mpi.rank();

  // get up-to-date cluster ids of the ghost layer before starting
  synchronizePlainGhostClusterIds();

  // number of reductions per iteration
  unsigned int l_numberOfReductions      = 1;

  // total number of reductions
  unsigned int l_totalNumberOfReductions = 0;

  while( l_numberOfReductions != 0 ) {
    l_numberOfReductions = 0;

    // iterate over mesh
    for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
      unsigned int l_minimumNeighborId = std::numeric_limits<unsigned int>::max();

      // get the ids
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        faceType l_faceType = getFaceType( m_cells[l_cell].boundaries[l_face] );
        if( l_faceType == regular  ||
            l_faceType == periodic ||
            l_faceType == dynamicRupture ) {
          // neighbor cell is part of copy layer/interior
          if( m_cells[l_cell].neighborRanks[l_face] == rank ) {
            unsigned int l_neighborId = m_cells[l_cell].neighbors[l_face];
            l_minimumNeighborId = std::min( l_minimumNeighborId, m_cellClusterIds[l_neighborId] );
          }
          // cell is part of the ghost layer
          else {
            // derive id of the ghost region
            unsigned int l_region = getPlainRegion( m_cells[l_cell].neighborRanks[l_face] );

            // local id in the ghost region
            unsigned int l_localGhostCell = m_cells[l_cell].mpiIndices[l_face];

            assert( l_localGhostCell < m_numberOfPlainGhostCells[l_region] );

            l_minimumNeighborId = std::min( l_minimumNeighborId, m_plainGhostCellClusterIds[l_region][l_localGhostCell] );
          }
        }
      }

      // assert we have a face neighbor
      assert( l_minimumNeighborId != std::numeric_limits<unsigned int>::max() );

      // lower id of the cell if required
      if( m_cellClusterIds[l_cell] > l_minimumNeighborId + i_difference ) {
        m_cellClusterIds[l_cell] = l_minimumNeighborId + i_difference;
        l_numberOfReductions++;
      }

    }

    l_totalNumberOfReductions += l_numberOfReductions;
  }

  return l_totalNumberOfReductions;
}

unsigned int seissol::initializers::time_stepping::LtsLayout::enforceSingleBuffer() {
  // TODO: Implementation required
  return 0;
}

void seissol::initializers::time_stepping::LtsLayout::normalizeClustering() {
  const int rank = seissol::MPI::mpi.rank();
  // allocate memory for the cluster ids of the ghost layer
  m_plainGhostCellClusterIds = new unsigned int*[ m_plainNeighboringRanks.size() ];
  for( unsigned int l_neighbor = 0; l_neighbor < m_plainNeighboringRanks.size(); l_neighbor++ ) {
    m_plainGhostCellClusterIds[l_neighbor] = new unsigned int[ m_numberOfPlainGhostCells[l_neighbor] ];
  }

  // enforce requirements until mesh is valid
  unsigned int l_maximumDifference      = 0;
  unsigned int l_dynamicRupture         = 0;
  unsigned int l_singleBuffer           = 0;
  unsigned int l_totalMaximumDifference = 0;
  unsigned int l_totalDynamicRupture    = 0;
  unsigned int l_totalSingleBuffer      = 0;

  int l_globalContinue = 1;

  // continue until all ranks converged to a normalized mesh
  while( l_globalContinue ) {
    // get up-to-date cluster ids of the ghost layer before starting
    synchronizePlainGhostClusterIds();

    l_dynamicRupture = enforceDynamicRuptureGTS();

    // enforce maximum difference of cluster ids
    if( m_clusteringStrategy == single ) {
      l_maximumDifference = enforceMaximumDifference( 0 );
    }
    else if( m_clusteringStrategy == multiRate ) {
      l_maximumDifference = enforceMaximumDifference( 1 );
    }
    else logError() << "clustering stategy not supported";

    // TODO: missing implementation (works only for max. difference 0 or 1)
    l_singleBuffer = enforceSingleBuffer();

    l_totalMaximumDifference += l_maximumDifference;
    l_totalDynamicRupture    += l_dynamicRupture;
    l_totalSingleBuffer      += l_singleBuffer;

    // check if this rank requires another iteration
    int l_localContinue = l_maximumDifference + l_dynamicRupture + l_singleBuffer;

#ifdef USE_MPI
    // continue if any rank is required to continue
    MPI_Allreduce( &l_localContinue, &l_globalContinue, 1, MPI_INT, MPI_MAX, seissol::MPI::mpi.comm() );
#else
    l_globalContinue = l_localContinue;
#endif
  }

  //logInfo() << "Performed a total of" << l_totalMaximumDifference << "reductions (max. diff.) for" << m_cells.size() << "cells," << l_totalDynamicRupture << "reductions (dyn. rup.) for" << m_fault.size() << "faces.";
  
  int* localClusterHistogram = new int[m_numberOfGlobalClusters];
  for (unsigned cluster = 0; cluster < m_numberOfGlobalClusters; ++cluster) {
    localClusterHistogram[cluster] = 0;
  }
  for (unsigned cell = 0; cell < m_cells.size(); ++cell) {
    ++localClusterHistogram[ m_cellClusterIds[cell] ];
  }

  int* globalClusterHistogram = NULL;
#ifdef USE_MPI
  if (rank == 0) {
    globalClusterHistogram = new int[m_numberOfGlobalClusters];
  }
  MPI_Reduce(localClusterHistogram, globalClusterHistogram, m_numberOfGlobalClusters, MPI_INT, MPI_SUM, 0, seissol::MPI::mpi.comm());
#else
  globalClusterHistogram = localClusterHistogram;
#endif
  if (rank == 0) {
    logInfo(rank) << "Number of elements in time clusters:";
    for (unsigned cluster = 0; cluster < m_numberOfGlobalClusters; ++cluster) {
      logInfo(rank) << utils::nospace << cluster << ":" << utils::space << globalClusterHistogram[cluster];
    }
#ifdef USE_MPI
    delete[] globalClusterHistogram;
#endif
  }
  delete[] localClusterHistogram;
}

void seissol::initializers::time_stepping::LtsLayout::getTheoreticalSpeedup( double &o_perCellTimeStepWidths,
                                                                             double &o_clustering  ) {
  // use khan sum
  // 0: true sum
  // 1: compensation
  double l_localPerCellSpeedup[2];
  l_localPerCellSpeedup[0] = l_localPerCellSpeedup[1] = 0;
  double l_localClusteringSpeedup[2];
  l_localClusteringSpeedup[0] = l_localClusteringSpeedup[1] = 0;

  double l_t, l_y;

  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    l_y = (m_globalTimeStepWidths[m_numberOfGlobalClusters-1] /  m_cellTimeStepWidths[l_cell]) - l_localPerCellSpeedup[1];
    l_t = l_localPerCellSpeedup[0] + l_y;
    l_localPerCellSpeedup[1] = (l_t - l_localPerCellSpeedup[0] ) - l_y;
    l_localPerCellSpeedup[0] = l_t;

    l_y = (m_globalTimeStepWidths[m_numberOfGlobalClusters-1] / m_globalTimeStepWidths[ m_cellClusterIds[l_cell] ] ) - l_localClusteringSpeedup[1];
    l_t = l_localClusteringSpeedup[0] + l_y;
    l_localClusteringSpeedup[1] = (l_t - l_localClusteringSpeedup[0]) - l_y;
    l_localClusteringSpeedup[0] = l_t;
  }

  unsigned int l_localNumberOfCells = m_cells.size();
  unsigned int l_globalNumberOfCells = 0;

  // derive global number of cells
#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  MPI_Allreduce( &l_localNumberOfCells, &l_globalNumberOfCells, 1, MPI_UNSIGNED, MPI_SUM, seissol::MPI::mpi.comm() );
#else // USE_MPI
  l_globalNumberOfCells = l_localNumberOfCells;
#endif // USE_MPI

  // derive global "speedup"
#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  MPI_Allreduce( l_localPerCellSpeedup,    &o_perCellTimeStepWidths, 1, MPI_DOUBLE, MPI_SUM, seissol::MPI::mpi.comm() );
  MPI_Allreduce( l_localClusteringSpeedup, &o_clustering,            1, MPI_DOUBLE, MPI_SUM, seissol::MPI::mpi.comm() );
#else // USE_MPI
  o_perCellTimeStepWidths = l_localPerCellSpeedup[0];
  o_clustering = l_localClusteringSpeedup[0];
#endif // USE_MPI

  o_perCellTimeStepWidths = (l_globalNumberOfCells * ( m_globalTimeStepWidths[m_numberOfGlobalClusters-1] / m_globalTimeStepWidths[0] ) ) / o_perCellTimeStepWidths;
  o_clustering            = (l_globalNumberOfCells * ( m_globalTimeStepWidths[m_numberOfGlobalClusters-1] / m_globalTimeStepWidths[0] ) ) / o_clustering;
}

void seissol::initializers::time_stepping::LtsLayout::addClusteredCopyCell( unsigned int i_cellId,
                                                                            unsigned int i_globalClusterId,
                                                                            unsigned int i_neighboringRank,
                                                                            unsigned int i_neighboringGlobalClusterId ) {
  // get local cluster id
  unsigned int l_localClusterId = getLocalClusterId( i_globalClusterId );

  // create first copy region if not present by now
  if( m_clusteredCopy[l_localClusterId].size() == 0 ) {
    m_clusteredCopy[l_localClusterId].push_back( clusterCopyRegion() );
    m_clusteredCopy[l_localClusterId][0].first[0] = i_neighboringRank;
    m_clusteredCopy[l_localClusterId][0].first[1] = i_neighboringGlobalClusterId;
    m_clusteredCopy[l_localClusterId][0].second.push_back( i_cellId );
    return;
  }

  // iterate over copy regions
  for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_localClusterId].size(); l_region++ ) {
    // copy region already present: add cell
    if( m_clusteredCopy[l_localClusterId][l_region].first[0]  == i_neighboringRank &&
        m_clusteredCopy[l_localClusterId][l_region].first[1]  == i_neighboringGlobalClusterId ) {
      // assert ordering is preserved
      assert( *(m_clusteredCopy[l_localClusterId][l_region].second.end()-1) <= i_cellId );

      // only add a cell if not present already
      if( *(m_clusteredCopy[l_localClusterId][l_region].second.end()-1) != i_cellId ) {
        m_clusteredCopy[l_localClusterId][l_region].second.push_back( i_cellId );
      }
      break;
    }
    // copy region with higher rank or same rank and higher neighboring cluster present: insert before
    else if( m_clusteredCopy[l_localClusterId][l_region].first[0] > i_neighboringRank ||
             (  m_clusteredCopy[l_localClusterId][l_region].first[0] == i_neighboringRank
             && m_clusteredCopy[l_localClusterId][l_region].first[1] > i_neighboringGlobalClusterId ) ) {
      m_clusteredCopy[l_localClusterId].insert( m_clusteredCopy[l_localClusterId].begin()+l_region, clusterCopyRegion() );
      m_clusteredCopy[l_localClusterId][l_region].first[0] = i_neighboringRank;
      m_clusteredCopy[l_localClusterId][l_region].first[1] = i_neighboringGlobalClusterId;
      m_clusteredCopy[l_localClusterId][l_region].second.push_back( i_cellId );
      break;
    }
    // no matches: this cell comes into a new copy region at the very end
    else if( l_region == m_clusteredCopy[l_localClusterId].size() - 1 ) {
      m_clusteredCopy[l_localClusterId].push_back( clusterCopyRegion() );
      m_clusteredCopy[l_localClusterId][l_region+1].first[0] = i_neighboringRank;
      m_clusteredCopy[l_localClusterId][l_region+1].first[1] = i_neighboringGlobalClusterId;
      m_clusteredCopy[l_localClusterId][l_region+1].second.push_back( i_cellId );
    }
  }
}

void seissol::initializers::time_stepping::LtsLayout::sortClusteredCopyGts( clusterCopyRegion &io_copyRegion ) {
	const int rank = seissol::MPI::mpi.rank();

  // buffers holding cells sending either derivatives or buffers
  std::vector< unsigned int > l_derivatives;
  std::vector< unsigned int > l_buffers;

  for( unsigned int l_copyCell = 0; l_copyCell < io_copyRegion.second.size(); l_copyCell++ ) {
    // per default no reordering required
    bool l_reorder = false;

    // get associated mesh id
    unsigned int l_meshId = io_copyRegion.second[ l_copyCell ];

    // assert this is GTS
    assert( io_copyRegion.first[1] == m_cellClusterIds[l_meshId] );

    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      // check if the cell qualifies for reordering because of dynamic rupture
      if( // m_cells[l_meshId].neighborRanks[l_face] == io_copyRegion.first[0] && // Dynamic rupture cells in the copy layer communicate derivatives only
                                                                                  // TODO: Here's some minor potential for optimizations
          getFaceType( m_cells[l_meshId].boundaries[l_face] ) == dynamicRupture ) {
        l_reorder = true;
      }

      // check if the cells qualifies for reordering because of "GTS on der"
      if( m_cells[l_meshId].neighborRanks[l_face] == rank &&
          ( getFaceType( m_cells[l_meshId].boundaries[l_face] ) == regular         ||
            getFaceType( m_cells[l_meshId].boundaries[l_face] ) == periodic ) ) {
        // get neighboring mesh id
        unsigned int l_neighboringMeshId = m_cells[l_meshId].neighbors[l_face];

        // check for true GTS buffer
        if( m_cellClusterIds[l_neighboringMeshId] > m_cellClusterIds[l_meshId] ) {
          l_reorder = true;
        }
      }
      // reorder also when the "GTS on der"-case is related to another copy region
      // TODO: Strictly speaking this reordering is not required as the cell could operate
      //       on the buffer for GTS; However the additional overhead is minimal
      else if( m_cells[l_meshId].neighborRanks[l_face] != rank ) {
        unsigned int l_region = getPlainRegion( m_cells[l_meshId].neighborRanks[l_face] );
        unsigned int l_ghostId = m_cells[l_meshId].mpiIndices[l_face];
        unsigned int l_ghostClusterId = m_plainGhostCellClusterIds[l_region][l_ghostId];

        if( l_ghostClusterId > io_copyRegion.first[1] ) {
          l_reorder = true;
        }
      }
    }

    // add cell either to derivatives or buffers
    if( l_reorder ) l_derivatives.push_back( l_meshId );
    else            l_buffers.push_back(     l_meshId );
  }

  // save the number of derivatives
  io_copyRegion.first[2] = l_derivatives.size();

  // perform the reordering
  for( unsigned int l_cell = 0; l_cell < l_derivatives.size(); l_cell++ ) {
    io_copyRegion.second[l_cell] = l_derivatives[l_cell];
  }
  for( unsigned int l_cell = 0; l_cell < l_buffers.size();     l_cell++ ) {
    io_copyRegion.second[l_derivatives.size()+l_cell] = l_buffers[l_cell];
  }
}

void seissol::initializers::time_stepping::LtsLayout::deriveClusteredCopyInterior() {
	const int rank = seissol::MPI::mpi.rank();

  /*
   * get local clusters
   */
  // unique set of local cluster
  std::set< unsigned int > l_localClusters;

  // derive local clusters
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    l_localClusters.insert( m_cellClusterIds[l_cell] );
  }

  // convert set to vector
  m_localClusters = std::vector< unsigned int > ( l_localClusters.begin(), l_localClusters.end() );

  /*
   * Add cells to clustered copy layers
   */
  // resize interior and copy layer clusters to cover all local clusters
  m_clusteredInterior.resize( m_localClusters.size() );
  m_clusteredCopy.resize(     m_localClusters.size() );

  // iterate over all cells and add the respective layers
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    bool l_copyCell = false;

    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
      // copy cell
      if( m_cells[l_cell].neighborRanks[l_face] != rank ) {
        l_copyCell = true;

        // plain region of the ghost cell
        unsigned int l_plainRegion = getPlainRegion( m_cells[l_cell].neighborRanks[l_face] );

        // local id in the ghost region
        unsigned int l_localGhostCell = m_cells[l_cell].mpiIndices[l_face];
        assert( l_localGhostCell < m_numberOfPlainGhostCells[l_plainRegion] );

        // neighboring cluster id
        unsigned int l_neighboringClusterId = m_plainGhostCellClusterIds[l_plainRegion][l_localGhostCell];

        // add the cell to the respective copy region
        addClusteredCopyCell( l_cell,
                              m_cellClusterIds[l_cell],
                              m_cells[l_cell].neighborRanks[l_face],
                              l_neighboringClusterId );
      }
    }

    if( !l_copyCell ) {
      // get local cluster id
      unsigned int l_localClusterId = getLocalClusterId( m_cellClusterIds[l_cell] );

      // add cell to interior cluster
      m_clusteredInterior[l_localClusterId].push_back( l_cell );
    }
  }

  /*
   * Sort GTS regions: DR and "GTS on der" comes first.
   */
  for( unsigned int l_cluster = 0; l_cluster < m_localClusters.size(); l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_cluster].size(); l_region++ ) {
      // check for GTS and perform sorting in case
      if( m_localClusters[l_cluster] == m_clusteredCopy[l_cluster][l_region].first[1] ) {
        sortClusteredCopyGts( m_clusteredCopy[l_cluster][l_region] );
      }
    }
  }

  /*
   * Set number of derivatives for non-GTS regions.
   */
  for( unsigned int l_cluster = 0; l_cluster < m_localClusters.size(); l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_cluster].size(); l_region++ ) {
      // LTS "<"
      if( m_localClusters[l_cluster] < m_clusteredCopy[l_cluster][l_region].first[1] ) {
        m_clusteredCopy[l_cluster][l_region].first[2] = 0;
      }
      // LTS ">"
      if( m_localClusters[l_cluster] > m_clusteredCopy[l_cluster][l_region].first[1] ) {
        m_clusteredCopy[l_cluster][l_region].first[2] = m_clusteredCopy[l_cluster][l_region].second.size();
      }
    }
  }

  /*
   * make sure we have all interior cells and at least the plain copy layer cells clustered layout
   */
  unsigned int l_numberOfClusteredInteriorCells = 0;
  unsigned int l_numberOfClusteredCopyCells     = 0;
  for( unsigned int l_localClusterId = 0; l_localClusterId < m_localClusters.size(); l_localClusterId++ ) {
    l_numberOfClusteredInteriorCells += m_clusteredInterior[l_localClusterId].size();

    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_localClusterId].size(); l_region++ ) {
      l_numberOfClusteredCopyCells += m_clusteredCopy[l_localClusterId][l_region].second.size();
    }
  }

  unsigned int l_numberOfPlainInteriorCells = m_plainInterior.size();
  unsigned int l_numberOfPlainCopyCells     = 0;
  for( unsigned int l_region = 0; l_region < m_plainNeighboringRanks.size(); l_region++ ) {
    l_numberOfPlainCopyCells += m_plainCopyRegions[l_region].size();
  }

  if( l_numberOfPlainInteriorCells != l_numberOfClusteredInteriorCells ) {
    logError() << "mismatch in the number of interior cells"
               << l_numberOfPlainInteriorCells << l_numberOfClusteredInteriorCells;
  }

  if( l_numberOfPlainCopyCells > l_numberOfClusteredCopyCells ) {
    logError() << "clustered copy layer incomplete"
               << l_numberOfPlainCopyCells << l_numberOfClusteredCopyCells;
  }
}

void seissol::initializers::time_stepping::LtsLayout::deriveClusteredGhost() {
  /*
   * Get sizes of the ghost regions
   */
  // get raw sizes of the copy and ghost regions
  std::vector< std::vector < unsigned int > > l_copyBuffer;
  std::vector< std::vector < unsigned int > > l_clusteredGhostSizes;

  // total number of one sided mpi requests
  unsigned int l_numberOfMpiRequests = 0;

  // resize buffer and target structure to hold all clusters
  l_copyBuffer.resize(          m_clusteredCopy.size() );
  l_clusteredGhostSizes.resize( m_clusteredCopy.size() );

   // assemble copy buffer and target structure
  for( unsigned int l_localCluster = 0; l_localCluster < m_clusteredCopy.size(); l_localCluster++ ) {
    // resize buffer and target structure to hold the number of cells for all per-cluster communication regions
    l_copyBuffer[l_localCluster].resize(          m_clusteredCopy[l_localCluster].size() );
    l_clusteredGhostSizes[l_localCluster].resize( m_clusteredCopy[l_localCluster].size() );

    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_localCluster].size(); l_region++ ) {
      // "neighboring" number of elements
      l_copyBuffer[l_localCluster][l_region] = m_clusteredCopy[l_localCluster][l_region].second.size();
    }

    l_numberOfMpiRequests += m_clusteredCopy[l_localCluster].size();
  }

#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  // mpi requests spawned
  MPI_Request *l_requests = new MPI_Request[ l_numberOfMpiRequests*2 ];

  unsigned int l_request = 0;

  /*
   * communicate the sizes
   */
  for( unsigned int l_localCluster = 0; l_localCluster < m_clusteredCopy.size(); l_localCluster++ ) {
    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_localCluster].size(); l_region++ ) {
      assert( l_request < l_numberOfMpiRequests );

      // generate a unqiue identifier for the message: per rank each cluster is allowed to send to / receive from every other time cluster
      int l_sendIdentifier =    clusteredGhost                                                                // identifier of clustered ghost derivation
                                + m_localClusters[l_localCluster]                   *m_numberOfGlobalClusters // encoding of the local cluster id
                                + m_clusteredCopy[l_localCluster][l_region].first[1];                         // id of the neighboring cluster

      // receive identifier is the send identifier as computed in the neighbor
      int l_receiveIdentifier = clusteredGhost                                                                // identifier of clustered ghost derivation
                                + m_clusteredCopy[l_localCluster][l_region].first[1]*m_numberOfGlobalClusters // encoding of the "local" cluster id
                                + m_localClusters[l_localCluster];                                            // id of the "neighboring" cluster

      // send the size of this copy region
      MPI_Isend( &l_copyBuffer[l_localCluster][l_region],             // buffer
                  1,                                                  // size
                  MPI_UNSIGNED,                                       // data type
                  m_clusteredCopy[l_localCluster][l_region].first[0], // destination
                  l_sendIdentifier,                                   // message tag
                  seissol::MPI::mpi.comm(),                           // communicator
                  l_requests+l_request );                             // mpi request

      // receive the size of the corresponding ghost region
      MPI_Irecv( &l_clusteredGhostSizes[l_localCluster][l_region],    // buffer
                  1,                                                  // size
                  MPI_UNSIGNED,                                       // data type
                  m_clusteredCopy[l_localCluster][l_region].first[0], // source
                  l_receiveIdentifier,                                // message tag
                  seissol::MPI::mpi.comm(),                           // communicator
                  l_requests+l_numberOfMpiRequests+l_request );       // mpi request

      l_request++;
    }
  }

  // wait for sends/receives
  MPI_Waitall( l_numberOfMpiRequests*2, // size
               l_requests,              // array of requests
               MPI_STATUS_IGNORE );     // mpi status
#endif // USE_MPI

  /*
   * Get cell ids of the ghost regions.
   */

  // setup target data structure
  m_clusteredGhost.resize( m_clusteredCopy.size() );
  for( unsigned int l_cluster = 0; l_cluster < m_clusteredGhost.size(); l_cluster++ ) {
    m_clusteredGhost[l_cluster].resize( m_clusteredCopy[l_cluster].size() );

    for( unsigned int l_region = 0; l_region < m_clusteredGhost[l_cluster].size(); l_region++ ) {
      m_clusteredGhost[l_cluster][l_region].second.resize( l_clusteredGhostSizes[l_cluster][l_region] );
    }
  }

#ifdef USE_MPI
  // TODO please check if this ifdef is correct

  // reset number of requests
  l_request = 0;

  // communicate cell ids
  for( unsigned int l_cluster = 0; l_cluster < m_clusteredGhost.size(); l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_clusteredGhost[l_cluster].size(); l_region++ ) {
      assert( l_request < l_numberOfMpiRequests );

      // generate a unqiue identifier for the message: per rank each cluster is allowed to send to / receive from every other time cluster
      int l_sendIdentifier =    clusteredGhost                                                           // identifier of clustered ghost derivation
                                + m_localClusters[l_cluster]                   *m_numberOfGlobalClusters // encoding of the local cluster id
                                + m_clusteredCopy[l_cluster][l_region].first[1];                         // id of the neighboring cluster

      // receive identifier is the send identifier as computed in the neighbor
      int l_receiveIdentifier = clusteredGhost                                                           // identifier of clustered ghost derivation
                                + m_clusteredCopy[l_cluster][l_region].first[1]*m_numberOfGlobalClusters // encoding of the "local" cluster id
                                + m_localClusters[l_cluster];                                            // id of the "neighboring" cluster

      // send the size of this copy region
      MPI_Isend( &m_clusteredCopy[l_cluster][l_region].second[0],     // buffer
                  m_clusteredCopy[l_cluster][l_region].second.size(), // size
                  MPI_UNSIGNED,                                       // data type
                  m_clusteredCopy[l_cluster][l_region].first[0],      // destination
                  l_sendIdentifier,                                   // message tag
                  seissol::MPI::mpi.comm(),                           // communicator
                  l_requests+l_request );                             // mpi request

      // receive the size of the corresponding ghost region
      MPI_Irecv( &m_clusteredGhost[l_cluster][l_region].second[0], // buffer
                  l_clusteredGhostSizes[l_cluster][l_region],      // size
                  MPI_UNSIGNED,                                    // data type
                  m_clusteredCopy[l_cluster][l_region].first[0],   // source
                  l_receiveIdentifier,                             // message tag
                  seissol::MPI::mpi.comm(),                        // communicator
                  l_requests+l_numberOfMpiRequests+l_request );    // mpi request

      l_request++;
    }
  }

  // wait for sends/receives
  MPI_Waitall( l_numberOfMpiRequests*2, // size
               l_requests,              // array of requests
               MPI_STATUS_IGNORE );     // mpi status

  /*
   * communicate number of derivatives
   */
  l_request = 0;
  for( unsigned int l_localCluster = 0; l_localCluster < m_clusteredCopy.size(); l_localCluster++ ) {
    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_localCluster].size(); l_region++ ) {
      assert( l_request < l_numberOfMpiRequests );

      l_copyBuffer[l_localCluster][l_region] = m_clusteredCopy[l_localCluster][l_region].first[2];

      // generate a unqiue identifier for the message: per rank each cluster is allowed to send to / receive from every other time cluster
      int l_sendIdentifier =    clusteredGhost                                                                // identifier of clustered ghost derivation
                                + m_localClusters[l_localCluster]                   *m_numberOfGlobalClusters // encoding of the local cluster id
                                + m_clusteredCopy[l_localCluster][l_region].first[1];                         // id of the neighboring cluster

      // receive identifier is the send identifier as computed in the neighbor
      int l_receiveIdentifier = clusteredGhost                                                                // identifier of clustered ghost derivation
                                + m_clusteredCopy[l_localCluster][l_region].first[1]*m_numberOfGlobalClusters // encoding of the "local" cluster id
                                + m_localClusters[l_localCluster];                                            // id of the "neighboring" cluster

      // send the size of this copy region
      MPI_Isend( &l_copyBuffer[l_localCluster][l_region],             // buffer
                  1,                                                  // size
                  MPI_UNSIGNED,                                       // data type
                  m_clusteredCopy[l_localCluster][l_region].first[0], // destination
                  l_sendIdentifier,                                   // message tag
                  seissol::MPI::mpi.comm(),                           // communicator
                  l_requests+l_request );                             // mpi request

      // receive the size of the corresponding ghost region
      MPI_Irecv( &m_clusteredGhost[l_localCluster][l_region].first,   // buffer
                  1,                                                  // size
                  MPI_UNSIGNED,                                       // data type
                  m_clusteredCopy[l_localCluster][l_region].first[0], // source
                  l_receiveIdentifier,                                // message tag
                  seissol::MPI::mpi.comm(),                           // communicator
                  l_requests+l_numberOfMpiRequests+l_request );       // mpi request

      l_request++;
    }
  }

  // wait for sends/receives
  MPI_Waitall( l_numberOfMpiRequests*2, // size
               l_requests,              // array of requests
               MPI_STATUS_IGNORE );     // mpi status

  // free memory
  delete[] l_requests;
#endif // USE_MPI
}

void seissol::initializers::time_stepping::LtsLayout::deriveLayout( enum TimeClustering i_timeClustering,
                                                                    unsigned int        i_clusterRate ) {
	const int rank = seissol::MPI::mpi.rank();

  m_clusteringStrategy = i_timeClustering;

  // derive time stepping clusters and per-cell cluster ids (w/o normalizations)
  if( m_clusteringStrategy == single ) {
    MultiRate::deriveClusterIds( m_cells.size(),
                                 std::numeric_limits<unsigned int>::max(),
                                 m_cellTimeStepWidths,
                                 m_cellClusterIds,
                                 m_numberOfGlobalClusters,
                                 m_globalTimeStepWidths,
                                 m_globalTimeStepRates  );
  }
  else if ( m_clusteringStrategy == multiRate ) {
    MultiRate::deriveClusterIds( m_cells.size(),
                                 i_clusterRate,
                                 m_cellTimeStepWidths,
                                 m_cellClusterIds,
                                 m_numberOfGlobalClusters,
                                 m_globalTimeStepWidths,
                                 m_globalTimeStepRates );
  }

  // derive plain copy and the interior
  derivePlainCopyInterior();

  // derive plain ghost regions
  derivePlainGhost();

  // normalize mpi indices
  normalizeMpiIndices();

  // normalize clustering
  normalizeClustering();

  // get maximum speedups compared to GTS
  double l_perCellSpeedup, l_clusteringSpeedup;
  getTheoreticalSpeedup( l_perCellSpeedup, l_clusteringSpeedup );

  // get maximum speedup
  logInfo(rank) << "maximum theoretical speedup (compared to GTS):"
                  << l_perCellSpeedup << "per cell LTS," << l_clusteringSpeedup << "with the used clustering.";

  // derive clustered copy and interior layout
  deriveClusteredCopyInterior();

  // derive the region sizes of the ghost layer
  deriveClusteredGhost();
  
  // derive dynamic rupture layers
  deriveDynamicRupturePlainCopyInterior();
}

void seissol::initializers::time_stepping::LtsLayout::getCrossClusterTimeStepping( struct TimeStepping &o_timeStepping ) {
  // set number of global clusters
  o_timeStepping.numberOfGlobalClusters = m_numberOfGlobalClusters;

  o_timeStepping.globalTimeStepRates     = new unsigned int[ o_timeStepping.numberOfGlobalClusters ];
  o_timeStepping.globalCflTimeStepWidths = new double[ o_timeStepping.numberOfGlobalClusters ];

  // set global time step rates
  for( unsigned int l_cluster = 0; l_cluster < o_timeStepping.numberOfGlobalClusters; l_cluster++ ) {
    o_timeStepping.globalTimeStepRates[l_cluster]     = m_globalTimeStepRates[l_cluster];
    o_timeStepping.globalCflTimeStepWidths[l_cluster] = m_globalTimeStepWidths[l_cluster];
  }

  // set synchronization time invalid
  o_timeStepping.synchronizationTime = std::numeric_limits<double>::min();

  // set number of local clusters
  o_timeStepping.numberOfLocalClusters = m_localClusters.size();

  o_timeStepping.clusterIds = new unsigned int[ o_timeStepping.numberOfLocalClusters ];

  for( unsigned int l_cluster = 0; l_cluster < o_timeStepping.numberOfLocalClusters; l_cluster++ ) {
    o_timeStepping.clusterIds[l_cluster] = m_localClusters[l_cluster];
  }
}

void seissol::initializers::time_stepping::LtsLayout::getCellInformation( CellLocalInformation* io_cellLocalInformation,
                                                                          unsigned int         *&o_ltsToMesh,
                                                                          unsigned int          &o_numberOfMeshCells ) {
	const int rank = seissol::MPI::mpi.rank();

  // total sizes of the communication layers covering all clusters
  unsigned int l_totalGhostLayerSize, l_totalCopyLayerSize, l_totalInteriorSize;
  l_totalGhostLayerSize = l_totalCopyLayerSize = l_totalInteriorSize = 0;

  // ghost region offsets offsets relative to the absolute start of the lts setup
  std::vector< std::vector< unsigned int > > l_ghostOffsets;
  l_ghostOffsets.resize( m_clusteredCopy.size() );

  // copy region offsets offsets relative to the absolute start of the lts setup
  std::vector< std::vector< unsigned int > > l_copyOffsets;
  l_copyOffsets.resize( m_clusteredGhost.size() );

  // interior offsets offsets relative to the absolute start of the lts setup
  std::vector< unsigned int >  l_interiorOffsets;
  l_interiorOffsets.resize( m_clusteredInterior.size() );

  // last offset
  unsigned int l_lastOffset = 0;

  // derive offsets
  for( unsigned int l_cluster = 0; l_cluster < m_clusteredInterior.size(); l_cluster++ ) {
    l_ghostOffsets[l_cluster].resize( m_clusteredGhost[l_cluster].size() );
    // iterate over ghost regions
    for( unsigned int l_ghostRegion = 0; l_ghostRegion < m_clusteredGhost[l_cluster].size(); l_ghostRegion++ ) {
      l_ghostOffsets[l_cluster][l_ghostRegion] = l_lastOffset;
      l_lastOffset += m_clusteredGhost[l_cluster][l_ghostRegion].second.size();
    }

    l_copyOffsets[l_cluster].resize( m_clusteredCopy[l_cluster].size() );
    // iterate over copy regions
    for( unsigned int l_copyRegion = 0; l_copyRegion < m_clusteredCopy[l_cluster].size(); l_copyRegion++ ) {
      l_copyOffsets[l_cluster][l_copyRegion] = l_lastOffset;
      l_lastOffset += m_clusteredCopy[l_cluster][l_copyRegion].second.size();
    }

    l_interiorOffsets[l_cluster] = l_lastOffset;
    l_lastOffset += m_clusteredInterior[l_cluster].size();
  }

  // derive size of the ghost and copy layer
  for( unsigned int l_cluster = 0; l_cluster < m_clusteredCopy.size(); l_cluster++ ) {
    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_cluster].size(); l_region++ ) {
      l_totalGhostLayerSize += m_clusteredGhost[l_cluster][l_region].second.size();
      l_totalCopyLayerSize  += m_clusteredCopy[l_cluster][l_region].second.size();
    }
  }

  // derive size of the interior
  for( unsigned int l_cluster = 0; l_cluster < m_clusteredInterior.size(); l_cluster++ ) {
    l_totalInteriorSize += m_clusteredInterior[l_cluster].size();
  }

  // set number of mesh and lts cells
  o_numberOfMeshCells         = m_cells.size();
  unsigned numberOfLtsCells   = l_totalGhostLayerSize + l_totalCopyLayerSize + l_totalInteriorSize;


  // allocate memory
  // TODO: free sometime somewhere
  o_ltsToMesh            = new unsigned int[ numberOfLtsCells ];

  // current lts cell
  unsigned int l_ltsCell = 0;

  // iterate over the setup an derive a linear layout
  for( unsigned int l_cluster = 0; l_cluster < m_clusteredInterior.size(); l_cluster++ ) {
    /*
     * iterate over ghost layer
     */
    for( unsigned int l_region = 0; l_region < m_clusteredGhost[l_cluster].size(); l_region++ ) {
      for( unsigned int l_ghostCell = 0; l_ghostCell < m_clusteredGhost[l_cluster][l_region].second.size(); l_ghostCell++ ) {
        // get values
        unsigned int l_clusterId = m_clusteredCopy[l_cluster][l_region].first[1];

        // store values
        io_cellLocalInformation[l_ltsCell].clusterId = l_clusterId;

        // set mapping invalid
        o_ltsToMesh[l_ltsCell] = std::numeric_limits<unsigned int>::max();

        l_ltsCell++;
      }
    }

    /*
     * iterate over copy layer
     */
    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_cluster].size(); l_region++ ) {
      for( unsigned int l_copyCell = 0; l_copyCell < m_clusteredCopy[l_cluster][l_region].second.size(); l_copyCell++ ) {
        // get values
        unsigned int l_clusterId  = m_localClusters[l_cluster];
        unsigned int l_meshId     = m_clusteredCopy[l_cluster][l_region].second[l_copyCell];

        // store face independent information
        io_cellLocalInformation[l_ltsCell].clusterId = l_clusterId;

        // set mappings
        o_ltsToMesh[l_ltsCell]                   = l_meshId;

        // iterate over faces
        for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
          // store mpi-independent values
          io_cellLocalInformation[l_ltsCell].faceTypes[l_face]        = getFaceType( m_cells[l_meshId].boundaries[l_face] );
          io_cellLocalInformation[l_ltsCell].faceRelations[l_face][0] = m_cells[l_meshId].neighborSides[l_face];
          io_cellLocalInformation[l_ltsCell].faceRelations[l_face][1] = m_cells[l_meshId].sideOrientations[l_face];

          // neighboring cell is part of the corresponding ghost region
          if( m_cells[l_meshId].neighborRanks[l_face] != rank ) {
            // remark: it is not sufficient to just search in the corresponding ghost region as copy cells can have more than one mpi-neighbor

            // global neighboring cluster id
            unsigned int l_plainRegion = getPlainRegion( m_cells[l_meshId].neighborRanks[l_face] );
            unsigned int l_globalNeighboringCluster = m_plainGhostCellClusterIds[l_plainRegion][ m_cells[l_meshId].mpiIndices[l_face] ];

            // find local neighboring region
            unsigned int l_localNeighboringRegion;
            for( l_localNeighboringRegion = 0; l_localNeighboringRegion < m_clusteredCopy[l_cluster].size(); l_localNeighboringRegion++ ) {
              // match if the ghost regions matches the neighboring rank and cluster
              if( m_clusteredCopy[l_cluster][l_localNeighboringRegion].first[0]
						== static_cast<unsigned int>(m_cells[l_meshId].neighborRanks[l_face]) &&
                  m_clusteredCopy[l_cluster][l_localNeighboringRegion].first[1] == l_globalNeighboringCluster ) break;

              // make sure there's a match
              if( l_localNeighboringRegion == m_clusteredCopy[l_cluster].size() -1 ) logError() << "no matching neighboring ghost region";
            }

            // get mesh id in the neighboring domain
            unsigned int l_neighboringMeshId = m_plainGhostCellIds[l_plainRegion][ m_cells[l_meshId].mpiIndices[l_face] ];

            unsigned int l_localGhostId = searchClusteredGhostCell( l_neighboringMeshId,
                                                                    l_cluster,
                                                                    l_localNeighboringRegion );

            // store value
            io_cellLocalInformation[l_ltsCell].faceNeighborIds[l_face] = l_ghostOffsets[l_cluster][l_localNeighboringRegion] + l_localGhostId;
          }
          // else neighboring cell is part of the interior or copy layer
          else if( io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == regular  ||
                   io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == periodic ||
                   io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == dynamicRupture ) {
            // neighboring mesh id
            unsigned int l_neighboringMeshId = m_cells[l_meshId].neighbors[l_face];

            // neighboring cell is part of the copy layers
            if( m_cells[l_neighboringMeshId].neighborRanks[0] != rank ||
                m_cells[l_neighboringMeshId].neighborRanks[1] != rank ||
                m_cells[l_neighboringMeshId].neighborRanks[2] != rank ||
                m_cells[l_neighboringMeshId].neighborRanks[3] != rank    ) {

              // get neighboring copy positions
              unsigned int l_localNeighboringClusterId;
              // TODO check if these values are really initialized by searchClusterCopyCell()
			  unsigned int l_neighboringCopyRegion = std::numeric_limits<unsigned int>::max();
			  unsigned int l_localNeighboringCellId = std::numeric_limits<unsigned int>::max();

              searchClusteredCopyCell( l_neighboringMeshId,
                                       l_localNeighboringClusterId,
                                       l_neighboringCopyRegion,
                                       l_localNeighboringCellId );

              io_cellLocalInformation[l_ltsCell].faceNeighborIds[l_face] = l_copyOffsets[l_localNeighboringClusterId][l_neighboringCopyRegion] + l_localNeighboringCellId;
            }
            // neighboring cell is part of the interior: search for the position of the neighboring cell in the interior of the cluster
            else {
              // get neighboring interior cell positions
              unsigned int l_localNeighboringClusterId, l_localNeighboringCellId;

              searchClusteredInteriorCell( l_neighboringMeshId,
                                           l_localNeighboringClusterId,
                                           l_localNeighboringCellId );

              // store values
              io_cellLocalInformation[l_ltsCell].faceNeighborIds[l_face] = l_interiorOffsets[l_localNeighboringClusterId] + l_localNeighboringCellId;
            }
          }
        }

        l_ltsCell++;
      }
    }

    /*
     * interate over interior
     */
    for( unsigned int l_interiorCell = 0; l_interiorCell <  m_clusteredInterior[l_cluster].size(); l_interiorCell++ ) {
      // get values
      unsigned int l_clusterId  = m_localClusters[l_cluster];
      unsigned int l_meshId     = m_clusteredInterior[l_cluster][l_interiorCell];

      // store face independent information
      io_cellLocalInformation[l_ltsCell].clusterId = l_clusterId;

      // set mappings
      o_ltsToMesh[l_ltsCell]                   = l_meshId;

      // iterate over faces
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // store mpi-independent values
        io_cellLocalInformation[l_ltsCell].faceTypes[l_face]        = getFaceType( m_cells[l_meshId].boundaries[l_face] );
        io_cellLocalInformation[l_ltsCell].faceRelations[l_face][0] = m_cells[l_meshId].neighborSides[l_face];
        io_cellLocalInformation[l_ltsCell].faceRelations[l_face][1] = m_cells[l_meshId].sideOrientations[l_face];

        if( io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == regular  ||
            io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == periodic ||
            io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == dynamicRupture ) {
          // neighboring mesh id
          unsigned int l_neighboringMeshId = m_cells[l_meshId].neighbors[l_face];

          // neighboring cell is part of the copy layers
          if( m_cells[l_neighboringMeshId].neighborRanks[0] != rank ||
              m_cells[l_neighboringMeshId].neighborRanks[1] != rank ||
              m_cells[l_neighboringMeshId].neighborRanks[2] != rank ||
              m_cells[l_neighboringMeshId].neighborRanks[3] != rank    ) {

            // get neighboring copy positions
            unsigned int l_localNeighboringClusterId;
            // TODO check if these values are really initialized by searchClusterCopyCell()
			unsigned int l_neighboringCopyRegion = std::numeric_limits<unsigned int>::max();
			unsigned int l_localNeighboringCellId = std::numeric_limits<unsigned int>::max();

            searchClusteredCopyCell( l_neighboringMeshId,
                                     l_localNeighboringClusterId,
                                     l_neighboringCopyRegion,
                                     l_localNeighboringCellId );

            io_cellLocalInformation[l_ltsCell].faceNeighborIds[l_face] = l_copyOffsets[l_localNeighboringClusterId][l_neighboringCopyRegion] + l_localNeighboringCellId;
          }
          // neighboring cell is part of the interior: search for the position of the neighboring cell in the interior of the cluster
          else {
            // get neighboring interior cell positions
            unsigned int l_localNeighboringClusterId, l_localNeighboringCellId;

            searchClusteredInteriorCell( l_neighboringMeshId,
                                         l_localNeighboringClusterId,
                                         l_localNeighboringCellId );

            // store values
            io_cellLocalInformation[l_ltsCell].faceNeighborIds[l_face] = l_interiorOffsets[l_localNeighboringClusterId] + l_localNeighboringCellId;
          }
        }
      }

      l_ltsCell++;
    }
  }
}

void seissol::initializers::time_stepping::LtsLayout::getDynamicRuptureInformation( unsigned*&  ltsToFace,
                                                                                    unsigned*&   numberOfDRCopyFaces,
                                                                                    unsigned*&   numberOfDRInteriorFaces )
{
  assert( m_dynamicRupturePlainCopy.size() == m_dynamicRupturePlainInterior.size() );
  
  numberOfDRCopyFaces     = new unsigned[ m_dynamicRupturePlainCopy.size()     ];
  numberOfDRInteriorFaces = new unsigned[ m_dynamicRupturePlainInterior.size() ];
  
  unsigned numberOfDRFaces = 0;
  for (unsigned cluster = 0; cluster < m_dynamicRupturePlainCopy.size(); ++cluster) {
    numberOfDRCopyFaces[cluster]      = m_dynamicRupturePlainCopy[cluster].size();
    numberOfDRInteriorFaces[cluster]  = m_dynamicRupturePlainInterior[cluster].size();
    numberOfDRFaces += numberOfDRCopyFaces[cluster] + numberOfDRInteriorFaces[cluster];
  }
  
  ltsToFace = new unsigned[numberOfDRFaces];
  
  unsigned ltsId = 0;
  for (unsigned cluster = 0; cluster < m_dynamicRupturePlainCopy.size(); ++cluster) {
    for (std::vector<int>::const_iterator it = m_dynamicRupturePlainCopy[cluster].begin(); it != m_dynamicRupturePlainCopy[cluster].end(); ++it) {
      ltsToFace[ltsId++] = *it;
    }  
    for (std::vector<int>::const_iterator it = m_dynamicRupturePlainInterior[cluster].begin(); it != m_dynamicRupturePlainInterior[cluster].end(); ++it) {
      ltsToFace[ltsId++] = *it;
    }
  }
}

void seissol::initializers::time_stepping::LtsLayout::getMeshStructure( MeshStructure *&o_meshStructure ) {
  // allocate data for per cluster mesh structure
  o_meshStructure = new MeshStructure[ m_localClusters.size() ];

  // iterate over clusters
  for( unsigned int l_cluster = 0; l_cluster < m_localClusters.size(); l_cluster++ ) {
    // set number of regions
    o_meshStructure[l_cluster].numberOfRegions = m_clusteredCopy[l_cluster].size();

    o_meshStructure[l_cluster].neighboringClusters = (int (*) [2]) malloc( o_meshStructure[l_cluster].numberOfRegions * 2 * sizeof(int) );

    o_meshStructure[l_cluster].numberOfGhostRegionCells                   = new unsigned int[ o_meshStructure[l_cluster].numberOfRegions ];
    o_meshStructure[l_cluster].numberOfGhostRegionDerivatives             = new unsigned int[ o_meshStructure[l_cluster].numberOfRegions ];
    o_meshStructure[l_cluster].ghostRegions                               = new real*[        o_meshStructure[l_cluster].numberOfRegions ];
    o_meshStructure[l_cluster].ghostRegionSizes                           = new unsigned int[ o_meshStructure[l_cluster].numberOfRegions ];

    o_meshStructure[l_cluster].numberOfCopyRegionCells                    = new unsigned int[ o_meshStructure[l_cluster].numberOfRegions ];
    o_meshStructure[l_cluster].numberOfCommunicatedCopyRegionDerivatives  = new unsigned int[ o_meshStructure[l_cluster].numberOfRegions ];
    o_meshStructure[l_cluster].copyRegions                                = new real*[        o_meshStructure[l_cluster].numberOfRegions ];
    o_meshStructure[l_cluster].copyRegionSizes                            = new unsigned int[ o_meshStructure[l_cluster].numberOfRegions ];

    o_meshStructure[l_cluster].numberOfGhostCells = 0;
    o_meshStructure[l_cluster].numberOfCopyCells = 0;

    o_meshStructure[l_cluster].sendIdentifiers    = new int[ o_meshStructure[l_cluster].numberOfRegions ];
    o_meshStructure[l_cluster].receiveIdentifiers = new int[ o_meshStructure[l_cluster].numberOfRegions ];

    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_cluster].size(); l_region++ ) {
      // set neighboring cluster information
      o_meshStructure[l_cluster].neighboringClusters[l_region][0] = m_clusteredCopy[l_cluster][l_region].first[0];
      o_meshStructure[l_cluster].neighboringClusters[l_region][1] = m_clusteredCopy[l_cluster][l_region].first[1];

      // set number of ghost region cells
      o_meshStructure[l_cluster].numberOfGhostRegionCells[l_region] = m_clusteredGhost[l_cluster][l_region].second.size();

      // add ghost region to ghost cells
      o_meshStructure[l_cluster].numberOfGhostCells += o_meshStructure[l_cluster].numberOfGhostRegionCells[l_region];

      // set number of ghost region derivatives
      o_meshStructure[l_cluster].numberOfGhostRegionDerivatives[l_region] = m_clusteredGhost[l_cluster][l_region].first;

      // set number of copy region cells
      o_meshStructure[l_cluster].numberOfCopyRegionCells[l_region] = m_clusteredCopy[l_cluster][l_region].second.size();

      // set number of copy region derivatives
      o_meshStructure[l_cluster].numberOfCommunicatedCopyRegionDerivatives[l_region] = m_clusteredCopy[l_cluster][l_region].first[2];

      // add copy region to copy cells
      o_meshStructure[l_cluster].numberOfCopyCells += o_meshStructure[l_cluster].numberOfCopyRegionCells[l_region];

      // set unique send identifier
      o_meshStructure[l_cluster].sendIdentifiers[l_region]    = m_localClusters[l_cluster]                    * m_numberOfGlobalClusters +
                                                                m_clusteredCopy[l_cluster][l_region].first[1];

      // set unique receive identifier
      o_meshStructure[l_cluster].receiveIdentifiers[l_region] = m_clusteredCopy[l_cluster][l_region].first[1] * m_numberOfGlobalClusters +
                                                                m_localClusters[l_cluster];
    }

    o_meshStructure[l_cluster].numberOfInteriorCells = m_clusteredInterior[l_cluster].size();

#ifdef USE_MPI
    o_meshStructure[l_cluster].sendRequests    = new MPI_Request[ m_clusteredCopy[l_cluster].size() ];
    o_meshStructure[l_cluster].receiveRequests = new MPI_Request[ m_clusteredCopy[l_cluster].size() ];
#endif
  }
}
