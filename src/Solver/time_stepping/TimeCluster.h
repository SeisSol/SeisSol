/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * 
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
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
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * LTS cluster in SeisSol.
 **/

#ifndef TIMECLUSTER_H_
#define TIMECLUSTER_H_

#ifdef USE_MPI
#include <mpi.h>
#include <list>
#endif

#include <Initializer/typedefs.hpp>
#include <SourceTerm/typedefs.hpp>
#include <utils/logger.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>

#include <Kernels/Time.h>

#ifdef REQUIRE_SOURCE_MATRIX
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#else
#include <Kernels/Volume.h>
#include <Kernels/Boundary.h>
#endif

// some check for correct functionality
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
#ifndef _OPENMP
#error NUMBER_OF_THREADS_PER_GLOBALDATA_COPY requires OpenMP to be enabled
#endif
#endif

namespace seissol {
  namespace time_stepping {
    class TimeCluster;
  }
}

/**
 * Time cluster, which represents a collection of elements havign the same time step width.
 **/
class seissol::time_stepping::TimeCluster
{
public:
    //! cluster id on this rank
    const unsigned int m_clusterId;

    //! global cluster cluster id
    const unsigned int m_globalClusterId;

private:
    //! number of time steps
    unsigned long m_numberOfTimeSteps;

    //! sampling of the receivers
    double m_receiverSampling;

    /*
     * integrators
     */
    //! time kernel
    kernels::Time     &m_timeKernel;

    //! local kernel
    kernels::Local   &m_localKernel;

    //! neighbor kernel
    kernels::Neighbor &m_neighborKernel;

    /*
     * mesh structure
     */
    struct MeshStructure *m_meshStructure;

    /*
     * global data
     */
     //! global data structures
    struct GlobalData *m_globalData;
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
#if NUMBER_OF_THREADS_PER_GLOBALDATA_COPY > 0
    struct GlobalData *m_globalDataCopies;
#else
#error NUMBER_OF_THREADS_PER_GLOBALDATA_COPY needs to be larger than 0 if defined
#endif
#endif

    /*
     * element data and mpi queues
     */     
#ifdef USE_MPI
    //! pending copy region sends
    std::list< MPI_Request* > m_sendQueue;

    //! pending ghost region receives
    std::list< MPI_Request* > m_receiveQueue;
#endif    
    seissol::initializers::TimeCluster* m_clusterData;
    seissol::initializers::LTS*         m_lts;

    //! receivers
    std::vector< int > m_receivers;
    
    //! Mapping of cells to point sources
    sourceterm::CellToPointSourcesMapping const* m_cellToPointSources;

    //! Number of mapping of cells to point sources
    unsigned m_numberOfCellToPointSourcesMappings;

    //! Point sources
    sourceterm::PointSources const* m_pointSources;

    //! true if dynamic rupture faces are present
    bool m_dynamicRuptureFaces;
    
    enum ComputePart {
      LocalInterior = 0,
      NeighborInterior,
#ifdef USE_MPI
      LocalCopy,
      NeighborCopy,
#endif
      NUM_COMPUTE_PARTS
    };
    
    long long m_flops_nonZero[NUM_COMPUTE_PARTS];
    long long m_flops_hardware[NUM_COMPUTE_PARTS];

#ifdef USE_MPI
    /**
     * Receives the copy layer data from relevant neighboring MPI clusters.
     **/
    void receiveGhostLayer();

    /**
     * Sends the associated regions of the copy layer to relevant neighboring MPI clusters
     **/
    void sendCopyLayer();

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
    /**
     * Inits Receives the copy layer data from relevant neighboring MPI clusters, active when using communication thread
     **/
    void initReceiveGhostLayer();

    /**
     * Inits Sends the associated regions of the copy layer to relevant neighboring MPI clusters, active when using communication thread
     **/
    void initSendCopyLayer();

    /**
     * Waits until the initialization of the communication is finished.
     **/
    void waitForInits();
#endif

    /**
     * Tests for pending ghost layer communication.
     **/
    bool testForGhostLayerReceives();

    /**
     * Tests for pending copy layer communication.
     **/
    bool testForCopyLayerSends();
#endif

    /**
     * Writes the receiver output if applicable (receivers present, receivers have to be written).
     **/
    void writeReceivers();

    /**
     * Computes the source terms if applicable.
     **/
    void computeSources();

    /**
     * Computes dynamic rupture.
     **/
    void computeDynamicRupture();

    /**
     * Computes all cell local integration.
     *
     * This are:
     *  * time integration
     *  * volume integration
     *  * local boundary integration
     *
     * Remark: After this step the DOFs are only updated half with the boundary contribution
     *         of the neighborings cells missing.
     *
     * @param i_numberOfCells number of cells.
     * @param i_cellInformation cell local information.
     * @param i_cellData cell data.
     * @param io_buffers time integration buffers.
     * @param io_derivatives time derivatives.
     * @param io_dofs degrees of freedom.
     **/
    void computeLocalIntegration( seissol::initializers::Layer&  i_layerData );

    /**
     * Computes the contribution of the neighboring cells to the boundary integral.
     *
     * Remark: After this step (in combination with the local integration) the DOFs are at the next time step.
     * TODO: This excludes dynamic rupture contribution.
     *
     * @param i_numberOfCells number of cells.
     * @param i_cellInformation cell local information.
     * @param i_cellData cell data.
     * @param i_faceNeighbors pointers to neighboring time buffers or derivatives.
     * @param io_dofs degrees of freedom.
     **/
    void computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData );

    void computeLocalIntegrationFlops(  unsigned                    numberOfCells,
                                        CellLocalInformation const* cellInformation,
                                        long long&                  nonZeroFlops,
                                        long long&                  hardwareFlops  );

    void computeNeighborIntegrationFlops( unsigned                    numberOfCells,
                                          CellLocalInformation const* cellInformation,
                                          long long&                  nonZeroFlops,
                                          long long&                  hardwareFlops);
    void computeFlops();

  public:
    //! flags identifiying if the respective part is allowed to be updated
    struct {
      bool localCopy;
      bool neighboringCopy;
      bool localInterior;
      bool neighboringInterior;
    } m_updatable;

#ifdef USE_MPI
    //! send true LTS buffers
    volatile bool m_sendLtsBuffers;
#endif

    //! reset lts buffers before performing time predictions
    volatile bool m_resetLtsBuffers;

    //! time step width of the performed time step.
    double m_timeStepWidth;

    /* Sub start time of width respect to the next cluster; use 0 if not relevant, for example in GTS.
     * LTS requires to evaluate a partial time integration of the derivatives. The point zero in time refers to the derivation of the surrounding time derivatives, which
     * coincides with the last completed time step of the next cluster. The start/end of the time step is the start/end of this clusters time step relative to the zero point.
     *   Example:
     *   <verb>
     *                                              5 dt
     *   |-----------------------------------------------------------------------------------------| <<< Time stepping of the next cluster (Cn) (5x larger than the current).
     *   |                 |                 |                 |                 |                 |
     *   |*****************|*****************|+++++++++++++++++|                 |                 | <<< Status of the current cluster.
     *   |                 |                 |                 |                 |                 |
     *   |-----------------|-----------------|-----------------|-----------------|-----------------| <<< Time stepping of the current cluster (Cc).
     *   0                 dt               2dt               3dt               4dt               5dt
     *   </verb>
     *
     *   In the example above two clusters are illustrated: Cc and Cn. Cc is the current cluster under consideration and Cn the next cluster with respect to LTS terminology.
     *   Cn is currently at time 0 and provided Cc with derivatives valid until 5dt. Cc updated already twice and did its last full update to reach 2dt (== subTimeStart). Next
     *   computeNeighboringCopy is called to accomplish the next full update to reach 3dt (+++). Besides working on the buffers of own buffers and those of previous clusters,
     *   Cc needs to evaluate the time prediction of Cn in the interval [2dt, 3dt].
     */
    double m_subTimeStart;

    //! number of full updates the cluster has performed since the last synchronization
    unsigned int m_numberOfFullUpdates;

    //! simulation time of the last full update (this is a complete volume and boundary integration)
    double m_fullUpdateTime;

    //! final time of the prediction (derivatives and time integrated DOFs).
    double m_predictionTime;

    //! time of the next receiver output
    double m_receiverTime;

    /**
     * Constructs a new LTS cluster.
     *
     * @param i_clusterId id of this cluster with respect to the current rank.
     * @param i_globalClusterId global id of this cluster.
     * @param i_timeKernel time integration kernel.
     * @param i_volumeKernel volume integration kernel.
     * @param i_boundaryKernel boundary integration kernel.
     * @param i_meshStructure mesh structure of this cluster.
     * @param i_copyCellInformation cell information in the copy layer.
     * @param i_interiorCellInformation cell information in the interior.
     * @param i_globalData global data.
     * @param i_copyCellData cell data in the copy layer.
     * @param i_interiorCellData cell data in the interior.
     * @param i_cells degrees of freedom, time buffers, time derivatives.
     **/
    TimeCluster( unsigned int                   i_clusterId,
                 unsigned int                   i_globalClusterId,
                 kernels::Time                 &i_timeKernel,
                 kernels::Local                &i_localKernel,
                 kernels::Neighbor             &i_neighborKernel,
                 struct MeshStructure          *i_meshStructure,
                 struct GlobalData             *i_globalData,
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
                 struct GlobalData             *i_globalDataCopies,
#endif
                 seissol::initializers::TimeCluster* i_clusterData,
                 seissol::initializers::LTS*         i_lts );

    /**
     * Destructor of a LTS cluster.
     * TODO: Currently prints only statistics in debug mode.
     **/
    ~TimeCluster();

    /**
     * Adds a source to the cluster.
     *
     * @param i_meshId mesh id of the point of interest.
     **/
    void addSource( unsigned int i_meshId );
    
    /**
     * Sets the pointer to the cluster's point sources
     * 
     * @param i_cellToPointSources Contains mappings of 1 cell offset to m point sources
     * @param i_numberOfCellToPointSourcesMappings Size of i_cellToPointSources
     * @param i_pointSources pointer to all point sources used on this cluster
     */
    void setPointSources( sourceterm::CellToPointSourcesMapping const* i_cellToPointSources,
                          unsigned i_numberOfCellToPointSourcesMappings,
                          sourceterm::PointSources const* i_pointSources );

    /**
     * Adds a receiver to the cluster.
     *
     * @param i_receiverId id of the receiver as used in Fortran.
     * @param i_meshId mesh id.
     **/
    void addReceiver( unsigned int i_receiverId,
                      unsigned int i_meshId );

    /**
     * Sets the receiver sampling for this cluster.
     *
     * @param i_receiverSampling receiver sampling.
     **/
    void setReceiverSampling( double i_receiverSampling );

    /**
     * Enables dynamic rupture call-backs in every time step.
     **/
    void enableDynamicRupture();

#ifdef USE_MPI
    /**
     * Computes cell local integration of all cells in the copy layer and initiates the corresponding communication.
     * LTS buffers (updated more than once in general) are reset to zero up on request; GTS-Buffers are reset independently of the request.
     *
     * Cell local integration is:
     *  * time integration
     *  * volume integration
     *  * local boundary integration
     *
     * @return true if the update (incl. communication requests), false if the update failed due to unfinshed sends of copy data to MPI neighbors.
     **/
    bool computeLocalCopy();
#endif

    /**
     * Computes cell local integration of all cells in the interior.
     * LTS buffers (updated more than once in general) are reset to zero up on request; GTS-Buffers are reset independently of the request.
     *
     * Cell local integration is:
     *  * time integration
     *  * volume integration
     *  * local boundary integration
     **/
    void computeLocalInterior();

#ifdef USE_MPI
    /**
     * Computes the neighboring contribution to the boundary integral for all cells in the copy layer.
     *
     * @return true if the update (incl. communication requests), false if the update failed due to missing data from neighboring ranks.
     **/
    bool computeNeighboringCopy();
#endif

    /**
     * Computes the neighboring contribution to the boundary integral for all cells in the interior.
     **/
    void computeNeighboringInterior();

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
    /**
     * Tests for pending ghost layer communication, active when using communication thread 
     **/
    void pollForGhostLayerReceives();

    /**
     * Polls for pending copy layer communication, active when using communication thread 
     **/
    void pollForCopyLayerSends();

    /**
     * Start Receives the copy layer data from relevant neighboring MPI clusters, active when using communication thread
     **/
    void startReceiveGhostLayer();

    /**
     * start Sends the associated regions of the copy layer to relevant neighboring MPI clusters, active when using communication thread
     **/
    void startSendCopyLayer();
#endif
};

#endif
