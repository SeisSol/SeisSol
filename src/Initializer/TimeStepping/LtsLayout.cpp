// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Rettenberger

#include "Parallel/MPI.h"

#include "utils/logger.h"

#include "LtsLayout.h"
#include "GlobalTimestep.h"
#include <Initializer/BasicTypedefs.h>
#include <Initializer/CellLocalInformation.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/TimeStepping/Halo.h>
#include <Monitoring/Unit.h>
#include <math.h>
#include <math.h>
#include <iterator>

#include "Initializer/ParameterDB.h"

#include <iomanip>
#include <limits>

seissol::initializer::time_stepping::LtsLayout::LtsLayout(const seissol::initializer::parameters::SeisSolParameters& parameters):
 m_plainCopyRegions(         NULL ),
 m_numberOfPlainGhostCells(  NULL ),
 seissolParams(parameters) {}

seissol::initializer::time_stepping::LtsLayout::~LtsLayout() {
  // free memory of member variables
  delete[] m_numberOfPlainGhostCells;
  delete[] m_plainCopyRegions;
}

void seissol::initializer::time_stepping::LtsLayout::setMesh( const seissol::geometry::MeshReader &i_mesh ) {
  // TODO: remove the copy by a pointer once the mesh stays constant
  m_cells = i_mesh.getElements();
  m_fault = i_mesh.getFault();

  m_mesh = &i_mesh;

  m_cellClusterIds.resize(m_cells.size());
  m_cellTimeStepWidths.resize(m_cells.size());

  m_numberOfGlobalClusters = 0;
  for (std::size_t i = 0; i < m_cells.size(); ++i) {
    m_cellClusterIds[i] = m_cells[i].clusterId;
    m_cellTimeStepWidths[i] = m_cells[i].timestep;

    m_numberOfGlobalClusters = std::max(m_numberOfGlobalClusters, m_cellClusterIds[i] + 1);
  }

  MPI_Allreduce( MPI_IN_PLACE, &m_numberOfGlobalClusters, 1, MPI_INT, MPI_MAX, seissol::MPI::mpi.comm() );
}

seissol::FaceType seissol::initializer::time_stepping::LtsLayout::getFaceType(int i_meshFaceType) {
  if (i_meshFaceType < 0 || i_meshFaceType > 7) {
    logError() << "face type" << i_meshFaceType << "not supported.";
  }
  return static_cast<FaceType>(i_meshFaceType);
}

void seissol::initializer::time_stepping::LtsLayout::derivePlainCopyInterior() {
	const int rank = seissol::MPI::mpi.rank();

  // unique set of neighboring ranks
  std::set< int > l_neighboringRanks;

  // derive neighboring ranks
  for( unsigned int l_cell = 0; l_cell < m_cells.size(); l_cell++ ) {
    for( unsigned int l_face = 0; l_face < Cell::NumFaces; l_face++ ) {
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
    for( unsigned int l_face = 0; l_face < Cell::NumFaces; l_face++ ) {
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

void seissol::initializer::time_stepping::LtsLayout::derivePlainGhost() {
  /*
   * Get sizes of ghost regions.
   */
  // number of copy cells
  unsigned int *l_numberOfCopyCells = new unsigned int[ m_plainNeighboringRanks.size() ];

  // number of ghost cells
  m_numberOfPlainGhostCells = new unsigned int[ m_plainNeighboringRanks.size() ];

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

  // free memory
  delete[] l_numberOfCopyCells;
}

void seissol::initializer::time_stepping::LtsLayout::deriveDynamicRupturePlainCopyInterior()
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
    // Dynamic rupture face with one neighbor in the ghost layer
    } else {
      m_dynamicRupturePlainCopy[localCluster].push_back(face);
    }
  }
}

void seissol::initializer::time_stepping::LtsLayout::normalizeMpiIndices() {
  /*
   * Convert the neighboring mappings to unique and sorted lists of the neighbors
   */
  for( unsigned int l_region = 0; l_region < m_plainNeighboringRanks.size(); l_region++ ) {
    std::vector<unsigned int> l_remoteFaceToCellIdMappings(m_mesh->getGhostlayerMetadata().at(m_plainNeighboringRanks[l_region]).size());
    for (std::size_t i = 0; i < l_remoteFaceToCellIdMappings.size(); ++i) {
      l_remoteFaceToCellIdMappings[i] = m_mesh->getGhostlayerMetadata().at(m_plainNeighboringRanks[l_region])[i].localId;
    }
    // sort
    std::sort( l_remoteFaceToCellIdMappings.begin(), l_remoteFaceToCellIdMappings.end() );
    // unique
    auto l_overhead = std::unique( l_remoteFaceToCellIdMappings.begin(), l_remoteFaceToCellIdMappings.end() );
    l_remoteFaceToCellIdMappings.erase( l_overhead, l_remoteFaceToCellIdMappings.end() );

    // store the results
    m_plainGhostCellIds.push_back( l_remoteFaceToCellIdMappings );
  }
}

void seissol::initializer::time_stepping::LtsLayout::addClusteredCopyCell( unsigned int i_cellId,
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

void seissol::initializer::time_stepping::LtsLayout::sortClusteredCopyGts( clusterCopyRegion &io_copyRegion ) {
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

    for( unsigned int l_face = 0; l_face < Cell::NumFaces; l_face++ ) {
      // check if the cell qualifies for reordering because of dynamic rupture
      if( // m_cells[l_meshId].neighborRanks[l_face] == io_copyRegion.first[0] && // Dynamic rupture cells in the copy layer communicate derivatives only
                                                                                  // TODO: Here's some minor potential for optimizations
	 getFaceType( m_cells[l_meshId].boundaries[l_face] ) == FaceType::DynamicRupture ) {
        l_reorder = true;
      }

      // check if the cells qualifies for reordering because of "GTS on der"
      if(m_cells[l_meshId].neighborRanks[l_face] == rank &&
         ( getFaceType( m_cells[l_meshId].boundaries[l_face] ) == FaceType::Regular ||
           getFaceType( m_cells[l_meshId].boundaries[l_face] ) == FaceType::Periodic )) {
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
        unsigned int l_ghostId = m_cells[l_meshId].mpiIndices[l_face];
        unsigned int l_ghostClusterId = m_mesh->getGhostlayerMetadata().at(m_cells[l_meshId].neighborRanks[l_face])[l_ghostId].clusterId;

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

void seissol::initializer::time_stepping::LtsLayout::deriveClusteredCopyInterior() {
	const int rank = seissol::MPI::mpi.rank();

  /*
   * get local clusters
   */
  // unique set of local cluster
  std::set< unsigned int > l_localClusters;

  // derive local clusters / HACK: no, always take all clusters
  for( unsigned int i = 0; i < m_numberOfGlobalClusters; i++ ) {
    l_localClusters.insert( i );
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

    for( unsigned int l_face = 0; l_face < Cell::NumFaces; l_face++ ) {
      // copy cell
      if( m_cells[l_cell].neighborRanks[l_face] != rank ) {
        l_copyCell = true;

        // local id in the ghost region
        unsigned int l_localGhostCell = m_cells[l_cell].mpiIndices[l_face];
        assert( l_localGhostCell < m_mesh->getGhostlayerMetadata().at(m_cells[l_cell].neighborRanks[l_face]).size() );

        // neighboring cluster id
        unsigned int l_neighboringClusterId = m_mesh->getGhostlayerMetadata().at(m_cells[l_cell].neighborRanks[l_face])[l_localGhostCell].clusterId;

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

void seissol::initializer::time_stepping::LtsLayout::deriveClusteredGhost() {
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
}

void seissol::initializer::time_stepping::LtsLayout::deriveLayout(  ) {
  // derive plain copy and the interior
  derivePlainCopyInterior();

  // derive plain ghost regions
  derivePlainGhost();

  // normalize mpi indices (not anymore; just set up the ghost layers)
  normalizeMpiIndices();

  // derive clustered copy and interior layout
  deriveClusteredCopyInterior();

  // derive the region sizes of the ghost layer
  deriveClusteredGhost();
  
  // derive dynamic rupture layers
  deriveDynamicRupturePlainCopyInterior();
}

seissol::initializer::ClusterLayout seissol::initializer::time_stepping::LtsLayout::clusterLayout() const {
  return ClusterLayout({m_globalTimeStepRates[0]}, m_globalTimeStepWidths[0], m_numberOfGlobalClusters);
}

void seissol::initializer::time_stepping::LtsLayout::getDynamicRuptureInformation( unsigned*&  ltsToFace )
{
  assert( m_dynamicRupturePlainCopy.size() == m_dynamicRupturePlainInterior.size() );
  
  unsigned numberOfDRFaces = 0;
  for (unsigned cluster = 0; cluster < m_dynamicRupturePlainCopy.size(); ++cluster) {
    numberOfDRFaces += m_dynamicRupturePlainCopy[cluster].size() + m_dynamicRupturePlainInterior[cluster].size();
  }
  
  ltsToFace = new unsigned[numberOfDRFaces];
  
  unsigned ltsId = 0;
  for (unsigned cluster = 0; cluster < m_dynamicRupturePlainCopy.size(); ++cluster) {
    for (auto it = m_dynamicRupturePlainCopy[cluster].begin(); it != m_dynamicRupturePlainCopy[cluster].end(); ++it) {
      ltsToFace[ltsId++] = *it;
    }  
    for (auto it = m_dynamicRupturePlainInterior[cluster].begin(); it != m_dynamicRupturePlainInterior[cluster].end(); ++it) {
      ltsToFace[ltsId++] = *it;
    }
  }
}

std::vector<std::size_t> seissol::initializer::time_stepping::LtsLayout::volumeSizes() const {
  std::vector<std::size_t> sizes(3 * m_localClusters.size());
  for( std::size_t i = 0; i < m_localClusters.size(); i++ ) {
    const std::size_t interior = m_clusteredInterior[i].size();
    std::size_t copy = 0;
    std::size_t ghost = 0;

    for( std::size_t j = 0; j < m_clusteredCopy[i].size(); j++ ) {
      copy += m_clusteredCopy[i][j].second.size();
      ghost += m_clusteredGhost[i][j].second.size();
    }

    sizes[0 + i * 3] = ghost;
    sizes[1 + i * 3] = copy;
    sizes[2 + i * 3] = interior;
  }
  return sizes;
}
std::vector<std::size_t> seissol::initializer::time_stepping::LtsLayout::drSizes() const {
  std::vector<std::size_t> sizes(3 * m_localClusters.size());
  for( std::size_t i = 0; i < m_localClusters.size(); i++ ) {
    const std::size_t interior = m_dynamicRupturePlainInterior[i].size();
    const std::size_t copy = m_dynamicRupturePlainCopy[i].size();
    const std::size_t ghost = 0;

    sizes[0 + i * 3] = ghost;
    sizes[1 + i * 3] = copy;
    sizes[2 + i * 3] = interior;
  }
  return sizes;
}

