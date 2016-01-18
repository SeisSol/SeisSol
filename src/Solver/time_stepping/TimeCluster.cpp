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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "TimeCluster.h"
#include <Solver/Interoperability.h>
#include <SourceTerm/PointSource.h>
#include <Kernels/TimeCommon.h>

#ifndef NDEBUG
extern long long g_SeisSolNonZeroFlopsLocal;
extern long long g_SeisSolHardwareFlopsLocal;
extern long long g_SeisSolNonZeroFlopsNeighbor;
extern long long g_SeisSolHardwareFlopsNeighbor;
#endif

#include <cstring>

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
extern volatile unsigned int* volatile g_handleRecvs;
extern volatile unsigned int* volatile g_handleSends;
#endif

//! fortran interoperability
extern seissol::Interoperability e_interoperability;

seissol::time_stepping::TimeCluster::TimeCluster( unsigned int                   i_clusterId,
                                                  unsigned int                   i_globalClusterId,
                                                  kernels::Time                 &i_timeKernel,
                                                  kernels::Volume               &i_volumeKernel,
                                                  kernels::Boundary             &i_boundaryKernel,
                                                  struct MeshStructure          *i_meshStructure,
#ifdef USE_MPI
                                                  struct CellLocalInformation   *i_copyCellInformation,
#endif
                                                  struct CellLocalInformation   *i_interiorCellInformation,
                                                  struct GlobalData             *i_globalData,
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
                                                  struct GlobalData             *i_globalDataCopies,
#endif
#ifdef USE_MPI
                                                  struct CellData               *i_copyCellData,
#endif
                                                  struct CellData               *i_interiorCellData,
                                                  struct Cells                  *i_cells ):
 // cluster ids
 m_clusterId(               i_clusterId                ),
 m_globalClusterId(         i_globalClusterId          ),
 // kernels
 m_timeKernel(              i_timeKernel               ),
 m_volumeKernel(            i_volumeKernel             ),
 m_boundaryKernel(          i_boundaryKernel           ),
 // mesh structure
 m_meshStructure(           i_meshStructure            ),
 // global data
 m_globalData(              i_globalData               ),
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
 // global data copies
 m_globalDataCopies(        i_globalDataCopies         ),
#endif
 // cell info
#ifdef USE_MPI
 m_copyCellInformation(     i_copyCellInformation      ),
#endif
 m_interiorCellInformation( i_interiorCellInformation  ),
 // cell data
#ifdef USE_MPI
 m_copyCellData(            i_copyCellData             ),
#endif
 m_interiorCellData(        i_interiorCellData         ),
 // cells
 m_cells(                   i_cells                    ),
 m_cellToPointSources(      NULL                       ),
 m_numberOfCellToPointSourcesMappings(0                ),
 m_pointSources(            NULL                       )
{
    // assert all pointers are valid
    assert( m_meshStructure                            != NULL );
#ifdef USE_MPI
    assert( m_copyCellInformation                      != NULL );
#endif
    assert( m_interiorCellInformation                  != NULL );
    assert( m_globalData                               != NULL );
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
    assert( m_globalDataCopies                         != NULL );
#endif
#ifdef USE_MPI
    assert( m_copyCellData                             != NULL );
    assert( m_copyCellData->localIntegration           != NULL );
    assert( m_copyCellData->neighboringIntegration     != NULL );
#endif
    assert( m_interiorCellData                         != NULL );
    assert( m_interiorCellData->localIntegration       != NULL );
    assert( m_interiorCellData->neighboringIntegration != NULL );
    assert( m_cells                                    != NULL );

  // default: no updates are allowed
  m_updatable.localCopy           = false;
  m_updatable.neighboringCopy     = false;
  m_updatable.localInterior       = false;
  m_updatable.neighboringInterior = false;
#ifdef USE_MPI
  m_sendLtsBuffers                = false;
#endif
  m_resetLtsBuffers               = false;
  // set timings to zero
  m_numberOfTimeSteps             = 0;
  m_receiverTime                  = 0;
  m_receiverSampling              = 0;
  m_timeStepWidth                 = 0;
  m_subTimeStart                  = 0;
  m_numberOfFullUpdates           = 0;
  m_fullUpdateTime                = 0;
  m_predictionTime                = 0;

  // disable dynamic rupture by default
  m_dynamicRuptureFaces = false;
}

seissol::time_stepping::TimeCluster::~TimeCluster() {  
#ifndef NDEBUG
  logInfo() << "#(time steps):" << m_numberOfTimeSteps;
#endif
}

void seissol::time_stepping::TimeCluster::setPointSources( sourceterm::CellToPointSourcesMapping const* i_cellToPointSources,
                                                           unsigned i_numberOfCellToPointSourcesMappings,
                                                           sourceterm::PointSources const* i_pointSources )
{
  m_cellToPointSources = i_cellToPointSources;
  m_numberOfCellToPointSourcesMappings = i_numberOfCellToPointSourcesMappings;
  m_pointSources = i_pointSources;
}

void seissol::time_stepping::TimeCluster::addReceiver( unsigned int i_receiverId,
                                                       unsigned int i_meshId ) {
  logInfo() << "cluster" << m_clusterId << "adding local receiver" << i_receiverId << "located at cell" << i_meshId;
  m_receivers.push_back( i_receiverId );
}

void seissol::time_stepping::TimeCluster::setReceiverSampling( double i_receiverSampling ) {
  m_receiverSampling = i_receiverSampling;
}

void seissol::time_stepping::TimeCluster::writeReceivers() {
  SCOREP_USER_REGION( "writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION )

  if( m_receivers.size() > 0 && m_fullUpdateTime + m_timeStepWidth > m_receiverTime ) {
    logDebug() << "cluster" << m_clusterId << "is writing a total of" << m_receivers.size() << "receivers at time" << m_receiverTime;

    e_interoperability.writeReceivers( m_fullUpdateTime,
                                       m_timeStepWidth,
                                       m_receiverTime,
                                       m_receivers );

    // increase receiver time until larger than next update time
    while( m_fullUpdateTime + m_timeStepWidth > m_receiverTime ) {
      m_receiverTime += m_receiverSampling;
    }
  }
}

void seissol::time_stepping::TimeCluster::enableDynamicRupture() {
  m_dynamicRuptureFaces = true;
}

void seissol::time_stepping::TimeCluster::computeSources() {
  SCOREP_USER_REGION( "computeSources", SCOREP_USER_REGION_TYPE_FUNCTION )
  
  // Return when point sources not initialised. This might happen if there
  // are no point sources on this rank.
  if (m_numberOfCellToPointSourcesMappings != 0) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned mapping = 0; mapping < m_numberOfCellToPointSourcesMappings; ++mapping) {
      real dofUpdate[NUMBER_OF_ALIGNED_DOFS];
      memset(dofUpdate, 0, NUMBER_OF_ALIGNED_DOFS * sizeof(real));

      unsigned startSource = m_cellToPointSources[mapping].pointSourcesOffset;
      unsigned endSource = m_cellToPointSources[mapping].pointSourcesOffset + m_cellToPointSources[mapping].numberOfPointSources;
      if (m_pointSources->mode == sourceterm::PointSources::NRF) {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceNRF( m_pointSources->mInvJInvPhisAtSources[source],
                                                       m_pointSources->tensor[source],
                                                       m_pointSources->muA[source],
                                                       m_pointSources->lambdaA[source],
                                                       m_pointSources->slipRates[source],
                                                       m_fullUpdateTime,
                                                       m_fullUpdateTime + m_timeStepWidth,
                                                       dofUpdate );
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceFSRM( m_pointSources->mInvJInvPhisAtSources[source],
                                                        m_pointSources->tensor[source],
                                                        &m_pointSources->slipRates[source][0],
                                                        m_fullUpdateTime,
                                                        m_fullUpdateTime + m_timeStepWidth,
                                                        dofUpdate );
        }
      }
      
      unsigned offset = m_cellToPointSources[mapping].copyInteriorOffset;
      for (unsigned dof = 0; dof < NUMBER_OF_ALIGNED_DOFS; ++dof) {
#ifdef USE_MPI
        m_cells->copyDofs[offset][dof] += dofUpdate[dof];
#else
        m_cells->interiorDofs[offset][dof] += dofUpdate[dof];
#endif
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeDynamicRupture() {
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  if( m_dynamicRuptureFaces == true ) {
    e_interoperability.computeDynamicRupture( m_fullUpdateTime,
                                              m_timeStepWidth );

#ifdef USE_MPI
    // TODO: This is not optimal as we are syncing all copy layers
    e_interoperability.synchronizeCopyLayerDofs();
#endif
  }
}

#ifdef USE_MPI
/*
 * MPI-Communication during the simulation; exchange of DOFs.
 */
void seissol::time_stepping::TimeCluster::receiveGhostLayer(){
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )

  /*
   * Receive data of the ghost regions
   */
  for( unsigned int l_region = 0; l_region < m_meshStructure->numberOfRegions; l_region++ ) {
    // continue only if the cluster qualifies for communication
    if( m_resetLtsBuffers || m_meshStructure->neighboringClusters[l_region][1] <= static_cast<int>(m_globalClusterId) ) {
      // post receive request
      MPI_Irecv(   m_meshStructure->ghostRegions[l_region],                // initial address
                   m_meshStructure->ghostRegionSizes[l_region],            // number of elements in the receive buffer
                   MPI_C_REAL,                                               // datatype of each receive buffer element
                   m_meshStructure->neighboringClusters[l_region][0],      // rank of source
                   timeData+m_meshStructure->receiveIdentifiers[l_region], // message tag
                   MPI_COMM_WORLD,                                         // communicator
                   m_meshStructure->receiveRequests + l_region             // communication request
               );

      // add receive request to list of receives
      m_receiveQueue.push_back( m_meshStructure->receiveRequests + l_region );
    }
  }
}

void seissol::time_stepping::TimeCluster::sendCopyLayer(){
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )

  /*
   * Send data of the copy regions
   */
  for( unsigned int l_region = 0; l_region < m_meshStructure->numberOfRegions; l_region++ ) {
    if( m_sendLtsBuffers || m_meshStructure->neighboringClusters[l_region][1] <= static_cast<int>(m_globalClusterId) ) {
      // post send request
      MPI_Isend(   m_meshStructure->copyRegions[l_region],              // initial address
                   m_meshStructure->copyRegionSizes[l_region],          // number of elements in the send buffer
                   MPI_C_REAL,                                            // datatype of each send buffer element
                   m_meshStructure->neighboringClusters[l_region][0],   // rank of destination
                   timeData+m_meshStructure->sendIdentifiers[l_region], // message tag
                   MPI_COMM_WORLD,                                      // communicator
                   m_meshStructure->sendRequests + l_region             // communication request
               );

      // add send request to list of sends
      m_sendQueue.push_back(m_meshStructure->sendRequests + l_region );
    }
  }
}

bool seissol::time_stepping::TimeCluster::testForGhostLayerReceives(){
  SCOREP_USER_REGION( "testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION )

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  bool l_return;
  if (g_handleRecvs[m_clusterId] == 0) {
    l_return = true;
  } else {
    l_return = false;
  }
  return l_return;
#else
  // iterate over all pending receives
  for( std::list<MPI_Request*>::iterator l_receive = m_receiveQueue.begin(); l_receive != m_receiveQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the receive is complete
    MPI_Test( *l_receive, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending receives if completed
    if( l_mpiStatus == 1 )   l_receive = m_receiveQueue.erase( l_receive );
    // continue otherwise
    else                   ++l_receive;
  }

  // return true if the communication is finished
  return m_receiveQueue.empty();
#endif
}

bool seissol::time_stepping::TimeCluster::testForCopyLayerSends(){
  SCOREP_USER_REGION( "testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION )

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  bool l_return;
  if (g_handleSends[m_clusterId] == 0) {
    l_return = true;
  } else {
    l_return = false;
  }
  return l_return;
#else
  for( std::list<MPI_Request*>::iterator l_send = m_sendQueue.begin(); l_send != m_sendQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the send is complete
    MPI_Test( *l_send, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending sends if completed
    if( l_mpiStatus == 1 )   l_send = m_sendQueue.erase( l_send );
    // continue otherwise
    else                   ++l_send;
  }

  // return true if the communication is finished
  return m_sendQueue.empty();
#endif
}

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
void seissol::time_stepping::TimeCluster::initReceiveGhostLayer(){
  g_handleRecvs[m_clusterId] = 1;
}

void seissol::time_stepping::TimeCluster::initSendCopyLayer(){
  g_handleSends[m_clusterId] = 1;
}

void seissol::time_stepping::TimeCluster::waitForInits() {
  while( g_handleRecvs[m_clusterId] == 1 || g_handleSends[m_clusterId] == 1 );
}
#endif

#endif

void seissol::time_stepping::TimeCluster::computeLocalIntegration( unsigned int           i_numberOfCells,
                                                                   CellLocalInformation  *i_cellInformation,
                                                                   CellData              *i_cellData,
                                                                   real                 **io_buffers,
                                                                   real                 **io_derivatives,
                                                                   real                 (*io_dofs)[NUMBER_OF_ALIGNED_DOFS] ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  // local integration buffer
  real l_integrationBuffer[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));

  // pointer for the call of the ADER-function
  real *l_bufferPointer;

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer) schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
    // overwrite cell buffer
    // TODO: Integrate this step into the kernel

    bool l_buffersProvided = (i_cellInformation[l_cell].ltsSetup >> 8)%2 == 1; // buffers are provided
    bool l_resetBuffers = l_buffersProvided && ( (i_cellInformation[l_cell].ltsSetup >> 10) %2 == 0 || m_resetLtsBuffers ); // they should be reset

    if( l_resetBuffers ) {
      // assert presence of the buffer
      assert( io_buffers[l_cell] != NULL );

      l_bufferPointer = io_buffers[l_cell];
    }
    // work on local buffer
    else {
      l_bufferPointer = l_integrationBuffer;
    }

    m_timeKernel.computeAder(              m_timeStepWidth,
                                           m_globalData->stiffnessMatricesTransposed,
                                           io_dofs[l_cell],
                                           i_cellData->localIntegration[l_cell].starMatrices,
#ifdef REQUIRE_SOURCE_MATRIX
                                           i_cellData->localIntegration[l_cell].sourceMatrix,
#endif
                                           l_bufferPointer,
                                           io_derivatives[l_cell] );

    m_volumeKernel.computeIntegral(        m_globalData->stiffnessMatrices,
                                           l_bufferPointer,
                                           i_cellData->localIntegration[l_cell].starMatrices,
#ifdef REQUIRE_SOURCE_MATRIX
                                           i_cellData->localIntegration[l_cell].sourceMatrix,
#endif
                                           io_dofs[l_cell] );

    m_boundaryKernel.computeLocalIntegral( i_cellInformation[l_cell].faceTypes,
                                           m_globalData->fluxMatrices,
                                           l_bufferPointer,
                                           i_cellData->localIntegration[l_cell].nApNm1,
#ifdef ENABLE_STREAM_MATRIX_PREFETCH
                                           io_dofs[l_cell],
                                           io_buffers[l_cell+1],
                                           io_dofs[l_cell+1] );
#else
                                           io_dofs[l_cell] );
#endif

#ifndef NDEBUG
    unsigned int l_tempHardwareFlops = 0;
    unsigned int l_tempNonZeroFlops = 0;
    m_timeKernel.flopsAder(              l_tempNonZeroFlops,
                                         l_tempHardwareFlops);
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolNonZeroFlopsLocal += (long long)l_tempNonZeroFlops;
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolHardwareFlopsLocal += (long long)l_tempHardwareFlops;
    
    m_volumeKernel.flopsIntegral(        l_tempNonZeroFlops,
                                         l_tempHardwareFlops);
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolNonZeroFlopsLocal += (long long)l_tempNonZeroFlops;
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolHardwareFlopsLocal += (long long)l_tempHardwareFlops;

    m_boundaryKernel.flopsLocalIntegral( i_cellInformation[l_cell].faceTypes, 
                                         l_tempNonZeroFlops,
                                         l_tempHardwareFlops);
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolNonZeroFlopsLocal += (long long)l_tempNonZeroFlops;
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolHardwareFlopsLocal += (long long)l_tempHardwareFlops;
#endif

    // update lts buffers if required
    // TODO: Integrate this step into the kernel
    if( !l_resetBuffers && l_buffersProvided ) {
      // assert presence of the buffer
      assert( io_buffers[l_cell] != NULL );

      for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof++ ) {
        io_buffers[l_cell][l_dof] += l_integrationBuffer[l_dof];
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( unsigned int            i_numberOfCells,
                                                                         CellLocalInformation   *i_cellInformation,
                                                                         CellData               *i_cellData,
                                                                         real                 *(*i_faceNeighbors)[4],
                                                                         real                  (*io_dofs)[NUMBER_OF_ALIGNED_DOFS],
																		 real                  (*io_pstrain)[7] ) {
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  real *l_timeIntegrated[4];
#ifdef ENABLE_MATRIX_PREFETCH
  real *l_faceNeighbors_prefetch[4];
  real *l_fluxMatricies_prefetch[4];
#endif

#ifdef _OPENMP
#ifdef ENABLE_MATRIX_PREFETCH
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch, l_fluxMatricies_prefetch)
#else
  #pragma omp parallel for schedule(static) private(l_timeIntegrated)
#endif
#endif
  for( unsigned int l_cell = 0; l_cell < i_numberOfCells; l_cell++ ) {
    seissol::kernels::TimeCommon::computeIntegrals( m_timeKernel,
                                                    i_cellInformation[l_cell].ltsSetup,
                                                    i_cellInformation[l_cell].faceTypes,
                                                    m_subTimeStart,
                                                    m_timeStepWidth,
                                                    m_globalData,
                                                    i_cellData->neighboringIntegration[l_cell].timeIntegration,
                                                    i_faceNeighbors[l_cell],
#ifdef _OPENMP
                                                    *reinterpret_cast<real (*)[4][NUMBER_OF_ALIGNED_DOFS]>(&(m_globalData->integrationBufferLTS[omp_get_thread_num()*4*NUMBER_OF_ALIGNED_DOFS])),
#else
                                                    *reinterpret_cast<real (*)[4][NUMBER_OF_ALIGNED_DOFS]>(m_globalData->integrationBufferLTS),
#endif
                                                    l_timeIntegrated );
#ifdef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY
    GlobalData* l_globalData = &(m_globalDataCopies[omp_get_thread_num()/NUMBER_OF_THREADS_PER_GLOBALDATA_COPY]);
#else
    GlobalData* l_globalData = m_globalData;
#endif

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
    // first face's prefetches
    int l_face = 1;
    l_faceNeighbors_prefetch[0] = i_faceNeighbors[l_cell][l_face];
    l_fluxMatricies_prefetch[0] = l_globalData->fluxMatrices[4+(l_face*12)
                                                             +(i_cellInformation[l_cell].faceRelations[l_face][0]*3)
                                                             +(i_cellInformation[l_cell].faceRelations[l_face][1])];
    // second face's prefetches
    l_face = 2;
    l_faceNeighbors_prefetch[1] = i_faceNeighbors[l_cell][l_face];
    l_fluxMatricies_prefetch[1] = l_globalData->fluxMatrices[4+(l_face*12)
                                                             +(i_cellInformation[l_cell].faceRelations[l_face][0]*3)
                                                             +(i_cellInformation[l_cell].faceRelations[l_face][1])];
    // third face's prefetches
    l_face = 3;
    l_faceNeighbors_prefetch[2] = i_faceNeighbors[l_cell][l_face];
    l_fluxMatricies_prefetch[2] = l_globalData->fluxMatrices[4+(l_face*12)
                                                             +(i_cellInformation[l_cell].faceRelations[l_face][0]*3)
                                                             +(i_cellInformation[l_cell].faceRelations[l_face][1])];
    // fourth face's prefetches
    if (l_cell < (i_numberOfCells-1) ) {
      l_face = 0;
      l_faceNeighbors_prefetch[3] = i_faceNeighbors[l_cell+1][l_face];
      l_fluxMatricies_prefetch[3] = l_globalData->fluxMatrices[4+(l_face*12)
                                                               +(i_cellInformation[l_cell+1].faceRelations[l_face][0]*3)
                                                               +(i_cellInformation[l_cell+1].faceRelations[l_face][1])];
    } else {
      l_faceNeighbors_prefetch[3] = i_faceNeighbors[l_cell][3];
      l_fluxMatricies_prefetch[3] = l_globalData->fluxMatrices[4+(3*12)
                                                               +(i_cellInformation[l_cell].faceRelations[l_face][0]*3)
                                                               +(i_cellInformation[l_cell].faceRelations[l_face][1])];
    }
#endif

    // @TODO in case of multiple global data copies, choose a distribution which
    //       cannot generate a 0-id copy reference in the end as remainder handling
    m_boundaryKernel.computeNeighborsIntegral( i_cellInformation[l_cell].faceTypes,
                                               i_cellInformation[l_cell].faceRelations,
                                               l_globalData->fluxMatrices,
                                               l_timeIntegrated,
                                               i_cellData->neighboringIntegration[l_cell].nAmNm1,
#ifdef ENABLE_MATRIX_PREFETCH
                                               io_dofs[l_cell],
                                               l_faceNeighbors_prefetch,
                                               l_fluxMatricies_prefetch );
#else
                                               io_dofs[l_cell] );
#endif

#ifdef USE_PLASTICITY
  e_interoperability.computePlasticity(  m_timeStepWidth,
                                         i_cellData->plasticity[l_cell].initialLoading,
                                         io_dofs[l_cell],
										 io_pstrain[l_cell] );
#endif

#ifndef NDEBUG
    unsigned int l_tempHardwareFlops = 0;
    unsigned int l_tempNonZeroFlops = 0;
    m_boundaryKernel.flopsNeighborsIntegral( i_cellInformation[l_cell].faceTypes,
                                             i_cellInformation[l_cell].faceRelations,
                                             l_tempNonZeroFlops,
                                             l_tempHardwareFlops);
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolNonZeroFlopsNeighbor += (long long)l_tempNonZeroFlops;
#ifdef _OPENMP
    #pragma omp atomic
#endif
    g_SeisSolHardwareFlopsNeighbor += (long long)l_tempHardwareFlops;
#endif
  }
}

#ifdef USE_MPI
bool seissol::time_stepping::TimeCluster::computeLocalCopy(){
  SCOREP_USER_REGION( "computeLocalCopy", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.localCopy ) {
    logError() << "Invalid call of computeLocalCopy, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // continue only if copy layer sends are complete
  if( !testForCopyLayerSends() ) return false;

  // post receive requests
#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  initReceiveGhostLayer();
#else
  receiveGhostLayer();
#endif

  // MPI checks for receiver writes receivers either in the copy layer or interior
  if( m_updatable.localInterior ) writeReceivers();

  // integrate copy layer locally
  computeLocalIntegration( m_meshStructure->numberOfCopyCells,
                           m_copyCellInformation,
                           m_copyCellData,
                           m_cells->copyBuffers,
                           m_cells->copyDerivatives,
                           m_cells->copyDofs );

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  initSendCopyLayer();
#else
  sendCopyLayer();
#endif

#ifndef USE_COMM_THREAD
  // continue with communication
  testForGhostLayerReceives();
#endif

  // compute sources, update simulation time
  if( !m_updatable.localInterior ) {
    computeSources();
    m_predictionTime += m_timeStepWidth;
  }

  // update finished
  m_updatable.localCopy  = false;

  // wait until communication thread finished initializing the receives
#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  waitForInits();
#endif

  return true;
}
#endif

void seissol::time_stepping::TimeCluster::computeLocalInterior(){
  SCOREP_USER_REGION( "computeLocalInterior", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.localInterior ) {
    logError() << "Invalid call of computeLocalInterior, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // MPI checks for receiver writes receivers either in the copy layer or interior
#ifdef USE_MPI
  if( m_updatable.localCopy ) writeReceivers();
#else
  // non-MPI checks for write in the interior
  writeReceivers();
#endif

  // integrate interior cells locally
  computeLocalIntegration( m_meshStructure->numberOfInteriorCells,
                           m_interiorCellInformation,
                           m_interiorCellData,
                           m_cells->interiorBuffers,
                           m_cells->interiorDerivatives,
                           m_cells->interiorDofs );

#ifdef USE_MPI
#ifndef USE_COMM_THREAD
  // continue with communication
  testForGhostLayerReceives();
  testForCopyLayerSends();
#endif
#endif

  // compute sources, update simulation time
  if( !m_updatable.localCopy ) {
    computeSources();
    m_predictionTime += m_timeStepWidth;
  }

  // update finished
  m_updatable.localInterior = false;
}

#ifdef USE_MPI
bool seissol::time_stepping::TimeCluster::computeNeighboringCopy() {
  SCOREP_USER_REGION( "computeNeighboringCopy", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.neighboringCopy ) {
    logError() << "Invalid call of computeNeighboringCopy aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // continue only of ghost layer receives are complete
  if( !testForGhostLayerReceives() ) return false;

#ifndef USE_COMM_THREAD
  // continue with communication
  testForCopyLayerSends();
#endif

  computeNeighboringIntegration( m_meshStructure->numberOfCopyCells,
                                 m_copyCellInformation,
                                 m_copyCellData,
                                 m_cells->copyFaceNeighbors,
                                 m_cells->copyDofs,
								 m_cells->copyPstrain);

#ifndef USE_COMM_THREAD
  // continue with communication
  testForCopyLayerSends();
#endif

  // compute dynamic rupture, update simulation time and statistics
  if( !m_updatable.neighboringInterior ) {
    computeDynamicRupture();
    m_fullUpdateTime      += m_timeStepWidth;
    m_subTimeStart        += m_timeStepWidth;
    m_numberOfFullUpdates += 1;
    m_numberOfTimeSteps   += 1;
  }

  // update finished
  m_updatable.neighboringCopy = false;

  return true;
}
#endif

void seissol::time_stepping::TimeCluster::computeNeighboringInterior() {
  SCOREP_USER_REGION( "computeNeighboringInterior", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.neighboringInterior ) {
    logError() << "Invalid call of computeNeighboringInterior, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // Update all cells in the interior with the neighboring boundary contribution.
  computeNeighboringIntegration( m_meshStructure->numberOfInteriorCells,
                                 m_interiorCellInformation,
                                 m_interiorCellData,
                                 m_cells->interiorFaceNeighbors,
                                 m_cells->interiorDofs,
								 m_cells->interiorPstrain );

  // compute dynamic rupture, update simulation time and statistics
  if( !m_updatable.neighboringCopy ) {
    computeDynamicRupture();
    m_fullUpdateTime      += m_timeStepWidth;
    m_subTimeStart        += m_timeStepWidth;
    m_numberOfFullUpdates += 1;
    m_numberOfTimeSteps   += 1;
  }

  // update finished
  m_updatable.neighboringInterior = false;
}

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
void seissol::time_stepping::TimeCluster::pollForCopyLayerSends(){
  for( std::list<MPI_Request*>::iterator l_send = m_sendQueue.begin(); l_send != m_sendQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the send is complete
    MPI_Test( *l_send, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending sends if completed
    if( l_mpiStatus == 1 )   l_send = m_sendQueue.erase( l_send );
    // continue otherwise
    else                   ++l_send;
  }

  if (m_sendQueue.empty()) {
    g_handleSends[m_clusterId] = 0;
  }
}

void seissol::time_stepping::TimeCluster::pollForGhostLayerReceives(){
  // iterate over all pending receives
  for( std::list<MPI_Request*>::iterator l_receive = m_receiveQueue.begin(); l_receive != m_receiveQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the receive is complete
    MPI_Test( *l_receive, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending receives if completed
    if( l_mpiStatus == 1 )   l_receive = m_receiveQueue.erase( l_receive );
    // continue otherwise
    else                   ++l_receive;
  }

  if (m_receiveQueue.empty()) {
    g_handleRecvs[m_clusterId] = 0;
  }
}

void seissol::time_stepping::TimeCluster::startReceiveGhostLayer() {
  receiveGhostLayer();
}

void seissol::time_stepping::TimeCluster::startSendCopyLayer() {
  sendCopyLayer();
}
#endif

