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
 * Layout the LTS schemes compute on.
 **/

#ifndef LTSLAYOUT_H_
#define LTSLAYOUT_H_

#include <Initializer/typedefs.hpp>

#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshReader.h>

#include <array>
#include <limits>
#include <cassert>

namespace seissol {
  namespace initializers {
    namespace time_stepping {
      class LtsLayout;
    }
  }
}

/**
 * Layout used by the LTS schemes for computation.
 **/
class seissol::initializers::time_stepping::LtsLayout {
  //private:
    //! used clustering strategy
    enum TimeClustering m_clusteringStrategy;

    //! cells in the local domain
    std::vector<Element> m_cells;

    //! fault in the local domain
    std::vector<Fault> m_fault;

    //! time step widths of the cells (cfl)
    double       *m_cellTimeStepWidths;

    //! cluster ids of the cells
    unsigned int *m_cellClusterIds;

    //! number of clusters in the global domain
    unsigned int  m_numberOfGlobalClusters;

    //! time step widths of all clusters
    double       *m_globalTimeStepWidths;

    //! time step rates of all clusters
    unsigned int *m_globalTimeStepRates;

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
    faceType getFaceType( int i_meshFaceType );

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
     * Helper function for synchronization.
     */
    void synchronizePlainGhostData(unsigned* cellData, unsigned** plainGhostData);

    /**
     * Synchronizes the cluster ids of the cells in the plain ghost layer.
     **/
    void synchronizePlainGhostClusterIds();

    /**
     * Enforces the same time step for all cells with dynamic rupture faces.
     */
    unsigned enforceDynamicRuptureGTS();

    /**
     * Enforces a maximum cluster difference between all cells.
     * 0: GTS (no difference of cluster ids allowed).
     * 1: Only a single difference in the cluster id is allowed, for example 2 is allowed to neighbor 1,2,3 but not 0 or 4.
     * [...]
     *
     * @param i_difference maximum allowed difference.
     * @return number of performed per-cell adjustments.
     **/
    unsigned int enforceMaximumDifference( unsigned int i_difference );

    /**
     * Lowers the time step width of all cells, which would require neighboring
     * cells to have more than one time buffer.
     *  Example:
     *   A cell has time step 2*dt with neighboring time steps dt, 2*dt, 4*dt and 8*dt.
     *   The scheme only provides one derivative (-> 2*dt - neighbor) and one buffer (-> 4*dt neighbor).
     *   Therefore we have to lower the time step of the 8*dt neighbor to 4*dt.
     *
     * TODO: This function is not implemented but required for maximum differences other than 0 or 1.
     *
     * @return number of performed normalizations.
     **/
    unsigned int enforceSingleBuffer();

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
    LtsLayout();

    /**
     * Destructor which frees all dynamically allocated memory of the class members.
     **/
    ~LtsLayout();

    /**
     * Sets the mesh and mesh-dependent default values.
     *
     * @param i_mesh mesh.
     **/
    void setMesh( const MeshReader &i_mesh );

    /**
     * Sets the time step width of a specidic cell.
     *
     * @param i_cellId id of the cell.
     * @param i_timeStepWidth time step width of the cell.
     **/
    void setTimeStepWidth( unsigned int i_cellId,
                           double       i_timeStepWidth );

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

#endif
