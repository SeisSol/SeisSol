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
 * Multi-rate scheme.
 **/

#include "Parallel/MPI.h"

#include "common.hpp"
#include <limits>

#ifndef MULTIRATE_HPP
#define MULTIRATE_HPP

namespace seissol {
  namespace initializers {
    namespace time_stepping {
      class MultiRate;
    }
  }
}

/**
 * Multi-rate scheme for local time stepping.
 **/
class seissol::initializers::time_stepping::MultiRate {
  //private
    /**
     * Gets the information of the cluster belonging to a cell with the given time step width in a multi-rate scheme.
     *
     * @param i_timeStepWidth time step width.
     * @param i_minimumTimeStepWidth global minimum time step width.
     * @param i_multiRate rate of the multi-rate scheme.
     * @param o_clusterTimeStepWidth time step width of the corresponding cluster.
     * @param o_clusterId global id of the corrensponding cluster.
     **/
    static void getMultiRateInfo( double        i_timeStepWidth,
                                  double        i_minimumTimeStepWidth,
                                  unsigned int  i_multiRate,
                                  double       &o_clusterTimeStepWidth,
                                  unsigned int &o_clusterId ) {
      // first multi-rate interval
      double l_lower = i_minimumTimeStepWidth;
      double l_upper = i_multiRate*l_lower;

      for( unsigned int l_id = 0; ; l_id++ ) {
        // the first cluster with an upper bound above the time step width is our
        if( l_upper > i_timeStepWidth ) {
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
                               const double       *i_cellTimeStepWidths,
                                     unsigned int *o_cellClusterIds ) {
      logInfo(seissol::MPI::mpi.rank()) << "Deriving clusters ids for min. time step width / multiRate:" << i_minimumTimeStepWidth << "/"
                                                                                 << i_multiRate;
      // iterate over all cells
      for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
        double l_clusterTimeStepWidth;

        // get the id for this cell
        getMultiRateInfo( i_cellTimeStepWidths[l_cell],
                          i_minimumTimeStepWidth,
                          i_multiRate,
                          l_clusterTimeStepWidth,
                          o_cellClusterIds[l_cell] );
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
                                      unsigned int  &o_numberOfClusters,
                                      double       *&o_clusterTimeStepWidths,
                                      unsigned int *&o_timeStepRates ) {
       double l_currentMaximumTime = i_minimumTimeStepWidth;
       l_currentMaximumTime *= i_multiRate;

       for( o_numberOfClusters = 1; l_currentMaximumTime <= i_maximumTimeStepWidth; o_numberOfClusters++ ) {
         l_currentMaximumTime *= i_multiRate;
       }

       // allocate memory for the time step widths and rate; TODO: Free
       o_clusterTimeStepWidths = new double[        o_numberOfClusters ];
       o_timeStepRates         = new unsigned int [ o_numberOfClusters ];

       // set the time step widths
       o_clusterTimeStepWidths[0] = i_minimumTimeStepWidth;
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

      // derive the number and time step widths of the global clusters
      deriveGlobalClusters( l_minimumTimeStepWidth,
                            l_maximumTimeStepWidth,
                            i_multiRate,
                            o_numberOfGlobalClusters,
                            o_globalTimeStepWidths,
                            o_globalTimeStepRates );

      // derive cluster ids of the cells
      deriveCellIds( l_minimumTimeStepWidth,
                     i_multiRate,
                     i_numberOfCells,
                     i_timeStepWidths,
                     o_cellClusterIds );
    }
};

#endif
