// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSLAYOUT_H_

#include "Initializer/Typedefs.h"

#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"

#include "Initializer/Parameters/SeisSolParameters.h"

#include <array>
#include <cassert>
#include <limits>


  
namespace seissol::initializer::time_stepping {
  class LtsLayout;
} // namespace seissol::initializer::time_stepping
  

class seissol::initializer::time_stepping::LtsLayout {
  //private:
    const seissol::geometry::MeshReader* m_mesh;

    //! cells in the local domain
    std::vector<Element> m_cells;

    //! fault in the local domain
    std::vector<Fault> m_fault;

    //! time step widths of the cells (cfl)
    std::vector<double>       m_cellTimeStepWidths;

    //! cluster ids of the cells
    std::vector<unsigned int> m_cellClusterIds;

    //! number of clusters in the global domain
    unsigned int  m_numberOfGlobalClusters;

    //! time step widths of all clusters
    std::vector<double> m_globalTimeStepWidths;

    //! time step rates of all clusters
    std::vector<unsigned int> m_globalTimeStepRates;

    //! mpi tags used for communication
    enum mpiTag {
      deriveGhostPlain    = 0,
      normalizeIndices    = 1,
      synchronizeClusters = 2,
      clusteredGhost      = 3
    };

    /*
     * Plain characteristics: Used for internal setup only.
     * Only communication related issues (as in GTS) are relevant.
     * No setup of a per-cluster layer.
     * Ordering is:
     *   1) neighboring rank (if applicable)
     *   2) cell id
     */
    //! plain neighboring ranks
    std::vector< int > m_plainNeighboringRanks;

    //! cells in the iterior of the local domain
    std::vector< unsigned int > m_plainInterior;

    //! plain copy regions (duplicated entries only for multiple mpi neighbors)
    std::vector< unsigned int > *m_plainCopyRegions;

    //! cell ids of the plain ghost cells in the neighoring domain
    std::vector< std::vector< unsigned int > > m_plainGhostCellIds;

    //! number of ghost cells per rank
    unsigned int *m_numberOfPlainGhostCells;

    //! cluster ids of the cells in the ghost layer
    unsigned int **m_plainGhostCellClusterIds;

    //! face ids of interior dr faces
    std::vector< std::vector<int> > m_dynamicRupturePlainInterior;

    //! face ids of copy dr faces
    std::vector< std::vector<int> > m_dynamicRupturePlainCopy;

    /*
     * Cluster dependent characteristics.
     * This is the final layout used in the computational scheme.
     * Ordering:
     *  1) local cluster
     *  2) neighboring rank (if applicable)
     *  3) neighboring cluster (if applicable)
     *  4) cell id (reordering for communication possible)
     */
    //! clusters present in the local computational domain
    std::vector< unsigned int > m_localClusters;

    /**
     * per cluster-neighboring rank
     * [*][ ] : cluster
     * [ ][*] : mpi-rank of the local communication region
     **/
    std::vector< std::vector< int > > m_clusterNeighboringRanks;

    /**
     * cell of a time stepping cluster
     **/
    typedef unsigned int clusterCell;

    /**
     * clusters in the interior
     * [*][ ]        : cluster
     * [ ][*]        : cluster local cell
     **/
    std::vector< std::vector< clusterCell > > m_clusteredInterior;

    /**
     * copy region of a time stepping cluster.
     * first[0]: mpi rank of the neighboring cluster
     * first[1]: global cluster id of the neighboring cluster
     * first[2]: number of derivatives cells communicated
     * second  : cluster cells in the copy region
     **/
    typedef std::pair< std::array<unsigned int, 3>, std::vector< clusterCell > > clusterCopyRegion;

    /**
     * per cluster copy layer.
     **/
    typedef std::vector< clusterCopyRegion > clusterCopyLayer;

    /**
     * copy layer consisting of all  per cluster copy layers.
     **/
    std::vector< clusterCopyLayer > m_clusteredCopy;

    /**
     * cells in the clustered ghost regions.
     * [*][ ][ ]: local cluster
     * [ ][*][ ]: communication region
     * [ ][ ].first: number of derivative cells communicated
     * [ ][ ].second[*]: cell id in the neighboring domain
     **/
    std::vector< std::vector< std::pair< unsigned int, std::vector< unsigned int > > > > m_clusteredGhost;

    const seissol::initializer::parameters::SeisSolParameters& seissolParams;

    /**
     * Gets the associated plain local region of the given mpi rank.
     *
     * @param i_mpiRank rank for which the local region is requested.
     **/
    unsigned int getPlainRegion( int i_mpiRank ) {
      std::vector<int>::iterator l_regionIterator = std::find( m_plainNeighboringRanks.begin(), m_plainNeighboringRanks.end(), i_mpiRank );
      unsigned int l_region = std::distance( m_plainNeighboringRanks.begin(), l_regionIterator );

      assert( l_region < m_plainNeighboringRanks.size() );
      return l_region;
    }

    /**
     * Get the id of the global cluster id in the local cluster setting
     *
     * @param i_clusterId global cluster id.
     * @return cluster id in the local setup.
     **/
    unsigned int getLocalClusterId( unsigned int i_clusterId ) {
      std::vector<unsigned int>::iterator l_clusterIterator = std::find( m_localClusters.begin(), m_localClusters.end(), i_clusterId );
      unsigned int l_clusterId = std::distance( m_localClusters.begin(), l_clusterIterator );

      assert( l_clusterId < m_localClusters.size() );
      return l_clusterId;
    }

    /**
     * Gets the enum face type from a mesh face type id.
     *
     * @param i_meshFaceType face type as stored in the mesh.
     **/
    FaceType getFaceType(int i_meshFaceType);

    /**
     * Derives plain copy regions and the interior.
     **/
    void derivePlainCopyInterior();

    /**
     * Derives plain ghost regions.
     **/
    void derivePlainGhost();

    /**
     * Derives plain copy and interior regions for dynamic rupture.
     **/
    void deriveDynamicRupturePlainCopyInterior();

    /**
     * Overwrite the given setting of mpi indices with required information.
     *   Default: Mpi index represents a per-region unique id of the mpi-face between two face neighbors.
     *   Overwrite: Local cell-id in the plain ghost region; additional a corresponding
     *   vector containing the cell ids in the neighboring computational domain is derived.
     **/
    void normalizeMpiIndices();

    /**
     * Normalizes the clustering.
     **/
    void normalizeClustering();

    /**
     * Gets the maximum possible speedups.
     *
     * @param o_perCellTimeStepWidths maximum possible speedup, when every cell is allowed to do an individual time step width.
     * @param o_clustering maximum possible speedup with the derived clustering strategy.
     **/
    void getTheoreticalSpeedup( double &o_perCellTimeStepWidths,
                                double &o_clustering );

    /**
     * Sorts a clustered copy region neighboring to a copy region in GTS fashion.
     * Copy cells send either buffers or derivatives to neighboring cells, never both.
     * For a copy region with clusterd id $l$ neighboring to a cluster with clusterd $n$ the following holds:
     *  l < n: Send buffers for all copy cells in the region.
     *  l > n: Send derivatives for all copy cells.
     *  l = n: Global time stepping, in general a copy cell sends derivatives; exceptions are:
     *         dynamic rupture mpi-face, which required derivatives
     *         face-neighbors, which provide their buffers in a true LTS fashion to cells with larger time steps -> GTS uses derivatives.
     *  As we need a linear ordering for MPI in the final scheme we sort the mesh: "Dynamic rupture" and "GTS on derivative"-cells come first.
     *  By additionally separating buffers and derivatives in the final scheme, a linear setting is derived:
     * <verb>
     *   __________________________________________________________________________________________________________________________________
     *  |                                        |                        ||                                  |                            |
     *  | Buffers: "DR"- and "GTS on der."-cells | Buffers: reg. elements || Derivatives: "DR"- "GTS on der." | Derivatives: reg. elements |
     *  |________________________________________|________________________||__________________________________|____________________________|
     *                                           |<--------------------- MPI Message ------------------------>|
     *                                           |____________________________________________________________|
     *
     * </verb>
     */
    void sortClusteredCopyGts( clusterCopyRegion &io_copyRegion);

    /**
     * Adds a specific cell with given cluster id, neighboring rank and neighboring cluster id to the respective copy region (if not present already).
     *
     * @param i_cellId id of the cell in the mesh.
     * @param i_globalClusterId global cluster id of the cell.
     * @param i_neighboringRank rank of the neighboring cell.
     * @param i_neighboringGlobalClusterId global cluster id of the neighboring cell.
     **/
    void addClusteredCopyCell( unsigned int i_cellId,
                               unsigned int i_globalClusterId,
                               unsigned int i_neighboringRank,
                               unsigned int i_neighboringGlobalClusterId );

    /**
     * Derives the layout of the copy regions and interior regions with respect to the clustering.
     **/
    void deriveClusteredCopyInterior();

    /**
     * Derives the clustered ghost region (cell ids in then neighboring domain).
     **/
    void deriveClusteredGhost();

    /**
     * Searches for the position of  cell in the specified ghost region.
     *
     * @param i_meshId mesh id of the cell searched for.
     * @param i_cluster copy layer of which cluster.
     * @param i_region copy region in the copy layer where to search.
     **/
    unsigned int searchClusteredGhostCell( unsigned int i_meshId,
                                           unsigned int i_cluster,
                                           unsigned int i_region ) {
      unsigned int l_localGhostId = 0;
      std::vector< unsigned int >::iterator l_searchResult;

      // non-gts neighbors have a linear ordering
      if( m_clusteredCopy[i_cluster][i_region].first[1] != m_localClusters[i_cluster] ) {
        // search for the right cell in the ghost region (exploits sorting by mesh ids)
        l_searchResult = std::lower_bound( m_clusteredGhost[i_cluster][i_region].second.begin(), // start of the search
                                           m_clusteredGhost[i_cluster][i_region].second.end(),   // end of the search
                                           i_meshId );                                           // value to search for
      }
      // gts neighbors not necessarily
      else {
        l_searchResult = std::find( m_clusteredGhost[i_cluster][i_region].second.begin(), // start of the search
                                    m_clusteredGhost[i_cluster][i_region].second.end(),   // end of the search
                                    i_meshId );                                           // value to search for
      }

      l_localGhostId = l_searchResult - m_clusteredGhost[i_cluster][i_region].second.begin();

      // ensure there's a valid result
      if(  l_localGhostId > m_clusteredGhost[i_cluster][i_region].second.size() - 1 ||
          *l_searchResult != i_meshId ) logError() << "no matching neighboring ghost region cell";

      return l_localGhostId;
    }

    /**
     * Searches for the position of a cell in the specified copy layer.
     *
     * @param i_meshId mesh id of the cell.
     * @param o_localClusterId set to the local id of the time cluster.
     * @param o_copyRegion set to the local id of the copy region containing the cell.
     * @param o_localCellId set to the local id of the cell in the copy region.
     **/
    void searchClusteredCopyCell( unsigned int  i_meshId,
                                  unsigned int &o_localClusterId,
                                  unsigned int &o_copyRegion,
                                  unsigned int &o_localCellId ) {
      // local cluster id
      o_localClusterId = m_cellClusterIds[ i_meshId ];
      o_localClusterId = getLocalClusterId( o_localClusterId );

      // search cell in the possible copy regions
      for( unsigned int l_region = 0; l_region < m_clusteredCopy[o_localClusterId].size(); l_region++ ) {
        std::vector< unsigned int >::iterator l_searchResult;
        unsigned int l_localCellId;

        // non-gts neighbors have a linear ordering
        if( m_clusteredCopy[o_localClusterId][l_region].first[1] != m_cellClusterIds[i_meshId] ) {
          l_searchResult = std::lower_bound( m_clusteredCopy[o_localClusterId][l_region].second.begin(), // start of the search
                                             m_clusteredCopy[o_localClusterId][l_region].second.end(),   // end of the search
                                             i_meshId );                                                 // value to search for
          l_localCellId = l_searchResult - m_clusteredCopy[o_localClusterId][l_region].second.begin();
        }
        // gts neighbors not necessarily
        else {
          l_searchResult = std::find( m_clusteredCopy[o_localClusterId][l_region].second.begin(), // start of the search
                                      m_clusteredCopy[o_localClusterId][l_region].second.end(),   // end of the search
                                      i_meshId );                                                 // value to search for
          l_localCellId = l_searchResult - m_clusteredCopy[o_localClusterId][l_region].second.begin();
        }

        // not found continue
        if( l_localCellId > m_clusteredCopy[o_localClusterId][l_region].second.size() - 1 ||
            *l_searchResult != i_meshId ) {
          assert( l_region != m_clusteredCopy[o_localClusterId].size() - 1 );
        }
        // success: store value and abort search over regions
        else {
          o_copyRegion  = l_region;
          o_localCellId = l_localCellId;
          return;
        }
      }
    }

    /**
     * Searches for the position of a cell in the interior.
     *
     * @param i_meshId mesh id of the cell.
     * @param o_localClusterId set to the local id of the time cluster.
     * @param o_localCellId set to the local id of the cell in the copy region.
     **/
    void searchClusteredInteriorCell( unsigned int  i_meshId,
                                      unsigned int &o_localClusterId,
                                      unsigned int &o_localCellId ) {
      // local cluster id
      o_localClusterId = m_cellClusterIds[ i_meshId ];
      o_localClusterId = getLocalClusterId( o_localClusterId );

      std::vector< unsigned int >::iterator l_searchResult = std::lower_bound( m_clusteredInterior[o_localClusterId].begin(), // start of the search
                                                                               m_clusteredInterior[o_localClusterId].end(),   // end of the search
                                                                               i_meshId );                                    // value to search for
      o_localCellId = l_searchResult - m_clusteredInterior[o_localClusterId].begin();

      // ensure a valid value
      if( o_localCellId > m_clusteredInterior[o_localClusterId].size() - 1 ||
          *l_searchResult != i_meshId ) logError() << "no matching neighboring interior cell";
    }

  public:
    /**
     * Constructor which initializes all pointers to NULL.
     **/
    LtsLayout(const seissol::initializer::parameters::SeisSolParameters& parameters);

    /**
     * Destructor which frees all dynamically allocated memory of the class members.
     **/
    ~LtsLayout();

    /**
     * Sets the mesh and mesh-dependent default values.
     *
     * @param i_mesh mesh.
     **/
    void setMesh( const seissol::geometry::MeshReader &i_mesh );

    /**
     * Derives the layout of the LTS scheme.
     *
     * @param i_timeClustering clustering strategy.
     * @param i_clusterRate cluster rate in the case of a multi-rate scheme.
     **/
    void deriveLayout( enum TimeClustering i_timeClustering,
                       unsigned int        i_clusterRate = std::numeric_limits<unsigned int>::max() );

    /**
     * Gets the cross-cluster time stepping information.
     *
     * @param o_timeStepping cross-cluster time stepping information.
     **/
    void getCrossClusterTimeStepping( struct TimeStepping &o_timeStepping );

    /**
     * Initializes the data structures required for computation.
     *
     * The ordering of the cell local information is:
     *  1) local cluster.
     *  2) ghost, copy, interior.
     *  3) neighboring rank (ghost and copy), neighboring cluster (ghost and copy).
     *  4) cell id in the mesh (reordering for communicatio possible).
     *
     * @param io_cellLocalInformation set to: cell local information of all computational cells.
     * @param o_ltsToMesh mapping from the global (accross all clusters and layers) lts id to the mesh id.
     * @param o_numberOfMeshCells number of cells in the mesh.
     **/
    void getCellInformation( CellLocalInformation* io_cellLocalInformation,
                             SecondaryCellLocalInformation* secondaryInformation,
                             unsigned int         *&o_ltsToMesh,
                             unsigned int          &o_numberOfMeshCells );

    void getDynamicRuptureInformation(  unsigned*&  ltsToFace,
                                        unsigned*&  numberOfDRCopyFaces,
                                        unsigned*&  numberOfDRInteriorFaces );

    /**
     * Get the per cluster mesh structure.
     *
     * @param mesh structure.
     **/
    void getMeshStructure( MeshStructure *&o_meshStructure );
};


#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSLAYOUT_H_

