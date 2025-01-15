// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_MULTIRATE_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_MULTIRATE_H_

#include "Common.h"
#include <algorithm>
#include <limits>

namespace seissol {
  class SeisSol;
  namespace initializer {
    namespace time_stepping {
      class MultiRate;
    }
  }
}

/**
 * Multi-rate scheme for local time stepping.
 **/
class seissol::initializer::time_stepping::MultiRate {
    /**
     * Gets the information of the cluster belonging to a cell with the given time step width in a multi-rate scheme.
     *
     * @param i_timeStepWidth time step width.
     * @param i_minimumTimeStepWidth global minimum time step width.
     * @param i_multiRate rate of the multi-rate scheme.
     * @param o_clusterTimeStepWidth time step width of the corresponding cluster.
     * @param o_clusterId global id of the corrensponding cluster.
     **/
  static void getMultiRateInfo(double i_timeStepWidth,
                               double i_minimumTimeStepWidth,
                               unsigned int i_multiRate,
                               unsigned int maxClusterId,
                               double& o_clusterTimeStepWidth,
                               unsigned int& o_clusterId) {
      // first multi-rate interval
      double l_lower = i_minimumTimeStepWidth;
      double l_upper = i_multiRate * l_lower;

      for( unsigned int l_id = 0; ; l_id++ ) {
        // the first cluster with an upper bound above the time step width is our
        // limit cluster id to maximum
        if (l_id >= maxClusterId || l_upper > i_timeStepWidth) {
          o_clusterTimeStepWidth = l_lower;
          o_clusterId = l_id;
          return;
        }

        // update interval and continue searching
        l_lower = l_upper;
        l_upper = i_multiRate * l_lower;
      }
    }

    /**
     * Derives the global time step with limits present in the overall computational domain.
     * @param i_numberOfCells local number of cells.
     * @param i_cellTimeSteppWidths time step widths of the cells.
     * @param o_minimumTimeStepWidth global minimum time step width.
     * @param o_maximumTimeStepWidth global maximum time step width.
     **/
    static void deriveTimeStepWidthLimits( unsigned int  i_numberOfCells,
                                           const double *i_cellTimeStepWidths,
                                                 double &o_minimumTimeStepWidth,
                                                 double &o_maximumTimeStepWidth ) {
      // derive local minimum and maximum time step widths
      double l_minimumTimeStepWidth = std::numeric_limits<double>::max();
      double l_maximumTimeStepWidth = std::numeric_limits<double>::min();

      for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
        l_minimumTimeStepWidth = std::min( l_minimumTimeStepWidth, i_cellTimeStepWidths[l_cell] );
        l_maximumTimeStepWidth = std::max( l_maximumTimeStepWidth, i_cellTimeStepWidths[l_cell] );
      }

      assert(l_minimumTimeStepWidth != std::numeric_limits<double>::max());
      assert(l_maximumTimeStepWidth != std::numeric_limits<double>::min());

#ifdef USE_MPI
      // derive global minimum and maximum time step widths
      MPI_Allreduce( &l_minimumTimeStepWidth, &o_minimumTimeStepWidth, 1, MPI_DOUBLE, MPI_MIN, seissol::MPI::mpi.comm() );
      MPI_Allreduce( &l_maximumTimeStepWidth, &o_maximumTimeStepWidth, 1, MPI_DOUBLE, MPI_MAX, seissol::MPI::mpi.comm() );
#else
      o_minimumTimeStepWidth = l_minimumTimeStepWidth;
      o_maximumTimeStepWidth = l_maximumTimeStepWidth;
#endif
    }

    /**
     * Derives the cluster ids the cells belong to.
     *
     * @param i_minimumTimeStepWidth global minimum time step width.
     * @param i_multiRate multi-rate configuration.
     * @param i_numberOfCells number of cells in the local domain.
     * @param i_cellTimeStepWidth time step widths of the cells.
     * @param o_cellClusterIds set to global cluster ids of the cells.
     **/
    static void deriveCellIds(       double        i_minimumTimeStepWidth,
                                     unsigned int  i_multiRate,
                                     unsigned int  i_numberOfCells,
                                     double        i_wiggleFactor,
                                     double        i_maxNumberOfClusters,
                               const double       *i_cellTimeStepWidths,
                                     unsigned int *o_cellClusterIds) {
      logInfo(seissol::MPI::mpi.rank()) << "Deriving clusters ids for min. time step width / multiRate:" << i_minimumTimeStepWidth << "/"
                                                                                 << i_multiRate;
      logInfo(seissol::MPI::mpi.rank())
          << "Due to wiggle factor of" << i_wiggleFactor << "the minimum timestep size is reduced to"
          << i_wiggleFactor * i_minimumTimeStepWidth;
      const auto maxClusterId = i_maxNumberOfClusters - 1;

      i_minimumTimeStepWidth *= i_wiggleFactor;
      // iterate over all cells
      for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
        double l_clusterTimeStepWidth;

        // get the id for this cell
        getMultiRateInfo(i_cellTimeStepWidths[l_cell],
                         i_minimumTimeStepWidth,
                         i_multiRate,
                         maxClusterId,
                         l_clusterTimeStepWidth,
                         o_cellClusterIds[l_cell]);
      }
    }

    /**
     * Derives the global setup -- number of clusters, ids, time step widths -- of the time clustering.
     *
     * @param i_minimumTimeStepWidth global minimum time step width.
     * @param i_maximumTimeStepWidth global maximum time step width.
     * @param i_multiRate rate of the multi-rate scheme.
     * @param o_numberOfCluster number of clusters in the global setting.
     * @param o_clusterTimeStepWidths time step widths of the clusters; memory will be allocated.
     * @param o_timeStepRates time step rate with respect to the next cluster
     **/
    static void deriveGlobalClusters( double         i_minimumTimeStepWidth,
                                      double         i_maximumTimeStepWidth,
                                      unsigned int   i_multiRate,
                                      double         i_wiggleFactor,
                                      unsigned int   i_maxNumberOfClusters,
                                      unsigned int  &o_numberOfClusters,
                                      double       *&o_clusterTimeStepWidths,
                                      unsigned int *&o_timeStepRates ) {
       double l_currentMaximumTime = i_wiggleFactor * i_minimumTimeStepWidth;
       l_currentMaximumTime *= i_multiRate;

       for( o_numberOfClusters = 1; l_currentMaximumTime <= i_maximumTimeStepWidth; o_numberOfClusters++ ) {
         l_currentMaximumTime *= i_multiRate;
       }

       o_numberOfClusters = std::min(
               o_numberOfClusters,
               static_cast<unsigned int>(i_maxNumberOfClusters)
       );
       assert(o_numberOfClusters > 0);
       // allocate memory for the time step widths and rate; TODO: Free
       o_clusterTimeStepWidths = new double[        o_numberOfClusters ];
       o_timeStepRates         = new unsigned int [ o_numberOfClusters ];

       // set the time step widths
       o_clusterTimeStepWidths[0] = i_wiggleFactor * i_minimumTimeStepWidth;
       for( unsigned int l_cluster = 1; l_cluster < o_numberOfClusters; l_cluster++ ) {
         o_clusterTimeStepWidths[l_cluster] = o_clusterTimeStepWidths[l_cluster-1] * i_multiRate;
       }

       // set time step rates
       o_timeStepRates[o_numberOfClusters-1] = 1; // last cluster has a rate of 1 (resets buffers in every time step)
       for( unsigned int l_cluster = 0; l_cluster < o_numberOfClusters-1; l_cluster++ ) {
         o_timeStepRates[l_cluster] = i_multiRate; // constant time step rate for the rest
       }
     }


  public:
    /**
     * Derives the cluster ids of the cells.
     *
     * @param i_numberOfCells number of cells in the local domain.
     * @param i_multiRate multi rate configuration.
     * @param i_timeStepWidths time step widths of the cells.
     * @param o_cellClusterIds set to the cluster ids of the cells.
     * @param o_numberOfGlobalClusters set to number of global clusters.
     * @param o_globalTimeStepWidths set to time step widths of the global clusters.
     * @param o_globalTimeStepRates set to the time step rates of the global clusters.
     **/
    static void deriveClusterIds( unsigned int   i_numberOfCells,
                                  unsigned int   i_multiRate,
                                  double         i_wiggleFactor,
                                  unsigned int   i_maxNumberOfClusters,
                                  double        *i_timeStepWidths,
                                  unsigned int  *o_cellClusterIds,
                                  unsigned int  &o_numberOfGlobalClusters,
                                  double       *&o_globalTimeStepWidths,
                                  unsigned int *&o_globalTimeStepRates ) {
      double l_minimumTimeStepWidth = 0, l_maximumTimeStepWidth = 0; // TODO check if 0 is current in serial version

      // derive global minimum and maximum
      deriveTimeStepWidthLimits( i_numberOfCells,
                                 i_timeStepWidths,
                                 l_minimumTimeStepWidth,
                                 l_maximumTimeStepWidth );
      // Note: l_minimumTimeStepWidth and l_maximumTimeStepWidth do not consider the wiggle factor!

      // derive the number and time step widths of the global clusters
      deriveGlobalClusters( l_minimumTimeStepWidth,
                            l_maximumTimeStepWidth,
                            i_multiRate,
                            i_wiggleFactor,
                            i_maxNumberOfClusters,
                            o_numberOfGlobalClusters,
                            o_globalTimeStepWidths,
                            o_globalTimeStepRates );

      // derive cluster ids of the cells
      deriveCellIds( l_minimumTimeStepWidth,
                     i_multiRate,
                     i_numberOfCells,
                     i_wiggleFactor,
                     i_maxNumberOfClusters,
                     i_timeStepWidths,
                     o_cellClusterIds );
    }
};


#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_MULTIRATE_H_

