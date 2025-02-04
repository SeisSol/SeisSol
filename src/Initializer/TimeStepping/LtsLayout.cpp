// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
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

#ifdef USE_MPI
  MPI_Allreduce( MPI_IN_PLACE, &m_numberOfGlobalClusters, 1, MPI_INT, MPI_MAX, seissol::MPI::mpi.comm() );
#endif

  m_globalTimeStepWidths.resize(m_numberOfGlobalClusters);
  m_globalTimeStepWidths[0] = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < m_cells.size(); ++i) {
    m_globalTimeStepWidths[0] = std::min(m_globalTimeStepWidths[0], m_cellTimeStepWidths[i]);
  }
#ifdef USE_MPI
  MPI_Allreduce( MPI_IN_PLACE, m_globalTimeStepWidths.data(), 1, MPI_DOUBLE, MPI_MIN, seissol::MPI::mpi.comm() );
#endif

  const auto wiggle = seissolParams.timeStepping.lts.getWiggleFactor();
  if (wiggle == 1) {
    logInfo() << "Minimum timestep:" << seissol::UnitTime.formatPrefix(m_globalTimeStepWidths[0]).c_str();
  }
  else {
    logInfo() << "Minimum timestep (pre-wiggle):" << seissol::UnitTime.formatPrefix(m_globalTimeStepWidths[0]).c_str();

    // apply wiggle here
    m_globalTimeStepWidths[0] *= wiggle;
    logInfo() << "Minimum timestep (with wiggle" << wiggle << "):" << seissol::UnitTime.formatPrefix(m_globalTimeStepWidths[0]).c_str();
  }
  logInfo() << "Global cluster count:" << m_numberOfGlobalClusters;
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

void seissol::initializer::time_stepping::LtsLayout::derivePlainGhost() {
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
    logInfo() << "Number of elements in dynamic rupture time clusters:";
    for (unsigned cluster = 0; cluster < m_numberOfGlobalClusters; ++cluster) {
      logInfo() << utils::nospace << cluster << " (dr):" << utils::space << globalClusterHistogram[cluster];
    }
#ifdef USE_MPI
    delete[] globalClusterHistogram;
#endif
  }
  delete[] localClusterHistogram;
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

void seissol::initializer::time_stepping::LtsLayout::normalizeClustering() {
  const int rank = seissol::MPI::mpi.rank();
  
  std::vector<int> localClusterHistogram(m_numberOfGlobalClusters);
  for (unsigned cluster = 0; cluster < m_numberOfGlobalClusters; ++cluster) {
    localClusterHistogram[cluster] = 0;
  }
  for (unsigned cell = 0; cell < m_cells.size(); ++cell) {
    ++localClusterHistogram[ m_cellClusterIds[cell] ];
  }

  std::vector<int> globalClusterHistogram;
#ifdef USE_MPI
  globalClusterHistogram.resize(m_numberOfGlobalClusters);
  MPI_Reduce(localClusterHistogram.data(), globalClusterHistogram.data(), m_numberOfGlobalClusters, MPI_INT, MPI_SUM, 0, seissol::MPI::mpi.comm());
#else
  globalClusterHistogram = localClusterHistogram;
#endif
  if (rank == 0) {
    logInfo() << "Number of elements in time clusters:";
    for (unsigned cluster = 0; cluster < m_numberOfGlobalClusters; ++cluster) {
      logInfo() << utils::nospace << cluster << ":" << utils::space << globalClusterHistogram[cluster];
    }
  }
}

void seissol::initializer::time_stepping::LtsLayout::getTheoreticalSpeedup( double &o_perCellTimeStepWidths,
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

    for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
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
        unsigned int l_region = getPlainRegion( m_cells[l_meshId].neighborRanks[l_face] );
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

void seissol::initializer::time_stepping::LtsLayout::deriveLayout( TimeClustering i_timeClustering,
                                                                    unsigned int        i_clusterRate ) {
	const int rank = seissol::MPI::mpi.rank();

  m_globalTimeStepRates.resize(1);
  m_globalTimeStepRates[0] = i_clusterRate;
  
  for (std::size_t i = 1; i < m_numberOfGlobalClusters; ++i) {
    m_globalTimeStepWidths[i] = m_globalTimeStepWidths[i - 1] * i_clusterRate;
  }

  // derive plain copy and the interior
  derivePlainCopyInterior();

  // derive plain ghost regions
  derivePlainGhost();

  // normalize mpi indices (not anymore; just set up the ghost layers)
  normalizeMpiIndices();

  // normalize clustering (not anymore, just prints some info; done in LtsWeights/somewhere else)
  normalizeClustering();

  // get maximum speedups compared to GTS
  double perCellSpeedup = 0;
  double clusteringSpeedup = 0;
  getTheoreticalSpeedup( perCellSpeedup, clusteringSpeedup );

  // The speedups above are computed without considering the wiggle factor
  const auto wiggleFactor = seissolParams.timeStepping.lts.getWiggleFactor();
  perCellSpeedup *= wiggleFactor;
  clusteringSpeedup *= wiggleFactor;

  // get maximum speedup
  logInfo() << "Theoretical speedup to GTS:" << perCellSpeedup << "elementwise LTS;" << clusteringSpeedup << "clustered LTS (current setup)";

  // derive clustered copy and interior layout
  deriveClusteredCopyInterior();

  // derive the region sizes of the ghost layer
  deriveClusteredGhost();
  
  // derive dynamic rupture layers
  deriveDynamicRupturePlainCopyInterior();
}

void seissol::initializer::time_stepping::LtsLayout::getCrossClusterTimeStepping( struct TimeStepping &o_timeStepping ) {
  // set number of global clusters
  o_timeStepping.numberOfGlobalClusters = m_numberOfGlobalClusters;

  o_timeStepping.globalTimeStepRates     = new unsigned int[ 1 ];
  o_timeStepping.globalCflTimeStepWidths = new double[ o_timeStepping.numberOfGlobalClusters ];

  o_timeStepping.globalTimeStepRates[0]     = m_globalTimeStepRates[0];

  // set global time step rates
  for( unsigned int l_cluster = 0; l_cluster < o_timeStepping.numberOfGlobalClusters; l_cluster++ ) {
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

void seissol::initializer::time_stepping::LtsLayout::getCellInformation( CellLocalInformation* io_cellLocalInformation,
                                                                          SecondaryCellLocalInformation* secondaryInformation,
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

  std::vector<int> duplicate(m_cells.size());


  // allocate memory
  // TODO: free sometime somewhere
  o_ltsToMesh            = new unsigned int[ numberOfLtsCells ];

  // current lts cell
  unsigned int l_ltsCell = 0;
  unsigned int layerId = 0;

  // iterate over the setup an derive a linear layout
  for( unsigned int l_cluster = 0; l_cluster < m_clusteredInterior.size(); l_cluster++ ) {
    /*
     * iterate over ghost layer
     */
    for( unsigned int l_region = 0; l_region < m_clusteredGhost[l_cluster].size(); l_region++ ) {
      for( unsigned int l_ghostCell = 0; l_ghostCell < m_clusteredGhost[l_cluster][l_region].second.size(); l_ghostCell++ ) {
        // get values
        unsigned int l_clusterId = m_clusteredCopy[l_cluster][l_region].first[1];

        // set mapping invalid
        o_ltsToMesh[l_ltsCell] = std::numeric_limits<unsigned int>::max();

        secondaryInformation[l_ltsCell].clusterId = l_clusterId;
        secondaryInformation[l_ltsCell].meshId = l_ghostCell;
        secondaryInformation[l_ltsCell].duplicate = 0;
        secondaryInformation[l_ltsCell].halo = HaloType::Ghost;

        secondaryInformation[l_ltsCell].globalId = m_mesh->getGhostlayerMetadata().at(m_clusteredCopy[l_cluster][l_region].first[0])[l_ghostCell].globalId; // TODO: check
        secondaryInformation[l_ltsCell].rank = m_clusteredCopy[l_cluster][l_region].first[0];
        secondaryInformation[l_ltsCell].layerId = layerId;
        secondaryInformation[l_ltsCell].configId = 0;

        l_ltsCell++;
        ++layerId;
      }
    }
    layerId = 0;

    /*
     * iterate over copy layer
     */
    for( unsigned int l_region = 0; l_region < m_clusteredCopy[l_cluster].size(); l_region++ ) {
      for( unsigned int l_copyCell = 0; l_copyCell < m_clusteredCopy[l_cluster][l_region].second.size(); l_copyCell++ ) {
        // get values
        unsigned int l_clusterId  = m_localClusters[l_cluster];
        unsigned int l_meshId     = m_clusteredCopy[l_cluster][l_region].second[l_copyCell];

        // store face independent information

        // set mappings
        o_ltsToMesh[l_ltsCell]                   = l_meshId;

        secondaryInformation[l_ltsCell].clusterId = l_clusterId;
        secondaryInformation[l_ltsCell].meshId = l_meshId;
        secondaryInformation[l_ltsCell].duplicate = duplicate[l_meshId];
        ++duplicate[l_meshId];

        secondaryInformation[l_ltsCell].halo = HaloType::Copy;
        secondaryInformation[l_ltsCell].globalId = m_cells[l_meshId].globalId;
        secondaryInformation[l_ltsCell].rank = rank;
        secondaryInformation[l_ltsCell].layerId = layerId;
        secondaryInformation[l_ltsCell].configId = 0;
        std::memcpy(secondaryInformation[l_ltsCell].neighborRanks, m_cells[l_meshId].neighborRanks, sizeof(int[4]));

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
            unsigned int l_globalNeighboringCluster = m_mesh->getGhostlayerMetadata().at(m_cells[l_meshId].neighborRanks[l_face])[ m_cells[l_meshId].mpiIndices[l_face]].clusterId;

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
            unsigned int l_neighboringMeshId = m_mesh->getGhostlayerMetadata().at(m_cells[l_meshId].neighborRanks[l_face])[ m_cells[l_meshId].mpiIndices[l_face]].localId;

            unsigned int l_localGhostId = searchClusteredGhostCell( l_neighboringMeshId,
                                                                    l_cluster,
                                                                    l_localNeighboringRegion );

            // store value
            secondaryInformation[l_ltsCell].faceNeighborIds[l_face] = l_ghostOffsets[l_cluster][l_localNeighboringRegion] + l_localGhostId;
          }
          // else neighboring cell is part of the interior or copy layer
          else if (io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == FaceType::Regular ||
                   io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == FaceType::Periodic ||
                   io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == FaceType::DynamicRupture) {
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

              secondaryInformation[l_ltsCell].faceNeighborIds[l_face] = l_copyOffsets[l_localNeighboringClusterId][l_neighboringCopyRegion] + l_localNeighboringCellId;
            }
            // neighboring cell is part of the interior: search for the position of the neighboring cell in the interior of the cluster
            else {
              // get neighboring interior cell positions
              unsigned int l_localNeighboringClusterId, l_localNeighboringCellId;

              searchClusteredInteriorCell( l_neighboringMeshId,
                                           l_localNeighboringClusterId,
                                           l_localNeighboringCellId );

              // store values
              secondaryInformation[l_ltsCell].faceNeighborIds[l_face] = l_interiorOffsets[l_localNeighboringClusterId] + l_localNeighboringCellId;
            }
          }
        }

        l_ltsCell++;
        ++layerId;
      }
    }
    layerId = 0;

    /*
     * interate over interior
     */
    for( unsigned int l_interiorCell = 0; l_interiorCell <  m_clusteredInterior[l_cluster].size(); l_interiorCell++ ) {
      // get values
      unsigned int l_clusterId  = m_localClusters[l_cluster];
      unsigned int l_meshId     = m_clusteredInterior[l_cluster][l_interiorCell];

      // store face independent information
      secondaryInformation[l_ltsCell].clusterId = l_clusterId;
      secondaryInformation[l_ltsCell].meshId = l_meshId;
      secondaryInformation[l_ltsCell].duplicate = 0;
      secondaryInformation[l_ltsCell].halo = HaloType::Interior;
      secondaryInformation[l_ltsCell].globalId = m_cells[l_meshId].globalId;
      secondaryInformation[l_ltsCell].rank = rank;
      secondaryInformation[l_ltsCell].layerId = layerId;
      secondaryInformation[l_ltsCell].configId = 0;
      std::memcpy(secondaryInformation[l_ltsCell].neighborRanks, m_cells[l_meshId].neighborRanks, sizeof(int[4]));

      // set mappings
      o_ltsToMesh[l_ltsCell]                   = l_meshId;

      // iterate over faces
      for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
        // store mpi-independent values
        io_cellLocalInformation[l_ltsCell].faceTypes[l_face]        = getFaceType( m_cells[l_meshId].boundaries[l_face] );
        io_cellLocalInformation[l_ltsCell].faceRelations[l_face][0] = m_cells[l_meshId].neighborSides[l_face];
        io_cellLocalInformation[l_ltsCell].faceRelations[l_face][1] = m_cells[l_meshId].sideOrientations[l_face];

        if(io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == FaceType::Regular ||
           io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == FaceType::Periodic ||
           io_cellLocalInformation[l_ltsCell].faceTypes[l_face] == FaceType::DynamicRupture) {
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

            secondaryInformation[l_ltsCell].faceNeighborIds[l_face] = l_copyOffsets[l_localNeighboringClusterId][l_neighboringCopyRegion] + l_localNeighboringCellId;
          }
          // neighboring cell is part of the interior: search for the position of the neighboring cell in the interior of the cluster
          else {
            // get neighboring interior cell positions
            unsigned int l_localNeighboringClusterId, l_localNeighboringCellId;

            searchClusteredInteriorCell( l_neighboringMeshId,
                                         l_localNeighboringClusterId,
                                         l_localNeighboringCellId );

            // store values
            secondaryInformation[l_ltsCell].faceNeighborIds[l_face] = l_interiorOffsets[l_localNeighboringClusterId] + l_localNeighboringCellId;
          }
        }
      }

      l_ltsCell++;
      ++layerId;
    }
    layerId = 0;
  }
}

void seissol::initializer::time_stepping::LtsLayout::getDynamicRuptureInformation( unsigned*&  ltsToFace,
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


void seissol::initializer::time_stepping::LtsLayout::getMeshStructure( MeshStructure *&o_meshStructure ) {
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

