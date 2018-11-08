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
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2013-2017, SeisSol Group
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

#include "Parallel/MPI.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SeisSol.h"
#include "TimeCluster.h"
#include <Solver/Interoperability.h>
#include <SourceTerm/PointSource.h>
#include <Kernels/TimeCommon.h>
#include <Kernels/DynamicRupture.h>
#include <Monitoring/FlopCounter.hpp>

#include <cassert>
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
                                                  kernels::Local                &i_localKernel,
                                                  kernels::Neighbor             &i_neighborKernel,
                                                  struct MeshStructure          *i_meshStructure,
                                                  struct GlobalData             *i_globalData,
                                                  seissol::initializers::TimeCluster* i_clusterData,
                                                  seissol::initializers::TimeCluster* i_dynRupClusterData,
                                                  seissol::initializers::LTS*         i_lts,
                                                  seissol::initializers::DynamicRupture* i_dynRup,
                                                  LoopStatistics*                        i_loopStatistics,
                                                  writer::ReceiverWriter*                receiverWriter ):
 // cluster ids
 m_clusterId(               i_clusterId                ),
 m_globalClusterId(         i_globalClusterId          ),
 // kernels
 m_timeKernel(              i_timeKernel               ),
 m_localKernel(             i_localKernel              ),
 m_neighborKernel(          i_neighborKernel           ),
 // mesh structure
 m_meshStructure(           i_meshStructure            ),
 // global data
 m_globalData(              i_globalData               ),
 m_clusterData(             i_clusterData              ),
 m_dynRupClusterData(       i_dynRupClusterData        ),
 m_lts(                     i_lts                      ),
 m_dynRup(                  i_dynRup                   ),
 // cells
 m_cellToPointSources(      NULL                       ),
 m_numberOfCellToPointSourcesMappings(0                ),
 m_pointSources(            NULL                       ),

 m_loopStatistics(          i_loopStatistics           ),
 m_receiverWriter(          receiverWriter             )
{
    // assert all pointers are valid
    assert( m_meshStructure                            != NULL );
    assert( m_globalData                               != NULL );
    assert( m_clusterData                              != NULL );

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
  m_receiverSampling              = std::numeric_limits<double>::max();
  m_timeStepWidth                 = 0;
  m_subTimeStart                  = 0;
  m_numberOfFullUpdates           = 0;
  m_fullUpdateTime                = 0;
  m_predictionTime                = 0;

  m_dynamicRuptureFaces = (i_dynRupClusterData->child<Ghost>().getNumberOfCells() > 0)
	|| (i_dynRupClusterData->child<Copy>().getNumberOfCells() > 0)
	|| (i_dynRupClusterData->child<Interior>().getNumberOfCells() > 0);

  computeFlops();

  m_regionComputeLocalIntegration = m_loopStatistics->getRegion("computeLocalIntegration");
  m_regionComputeNeighboringIntegration = m_loopStatistics->getRegion("computeNeighboringIntegration");
  m_regionComputeDynamicRupture = m_loopStatistics->getRegion("computeDynamicRupture");
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

void seissol::time_stepping::TimeCluster::setReceiverSampling( double i_receiverSampling ) {
  m_receiverSampling = i_receiverSampling;
}

void seissol::time_stepping::TimeCluster::writeReceivers() {
  SCOREP_USER_REGION( "writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_receiverTime = m_receiverWriter->writeReceivers(m_clusterId, m_receiverTime, m_fullUpdateTime, m_timeStepWidth, m_receiverSampling);
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
                                                       *m_cellToPointSources[mapping].dofs );
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceFSRM( m_pointSources->mInvJInvPhisAtSources[source],
                                                        m_pointSources->tensor[source],
                                                        &m_pointSources->slipRates[source][0],
                                                        m_fullUpdateTime,
                                                        m_fullUpdateTime + m_timeStepWidth,
                                                        *m_cellToPointSources[mapping].dofs );
        }
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData ) {
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  DRFaceInformation*                    faceInformation                                                   = layerData.var(m_dynRup->faceInformation);
  DRGodunovData*                        godunovData                                                       = layerData.var(m_dynRup->godunovData);
  real**                                timeDerivativePlus                                                = layerData.var(m_dynRup->timeDerivativePlus);
  real**                                timeDerivativeMinus                                               = layerData.var(m_dynRup->timeDerivativeMinus);
  real                                (*godunov)[CONVERGENCE_ORDER][seissol::model::godunovState::reals]  = layerData.var(m_dynRup->godunov);
  real                                (*imposedStatePlus)[seissol::model::godunovState::reals]            = layerData.var(m_dynRup->imposedStatePlus);
  real                                (*imposedStateMinus)[seissol::model::godunovState::reals]           = layerData.var(m_dynRup->imposedStateMinus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                                    = layerData.var(m_dynRup->waveSpeedsPlus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                                   = layerData.var(m_dynRup->waveSpeedsMinus);
  real                                (*absoluteSlip)[seissol::model::godunovState::rows]                 = layerData.var(m_dynRup->absoluteSlip);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynamicRuptureKernel.computeGodunovState( faceInformation[face],
                                                m_globalData,
                                               &godunovData[face],
                                                timeDerivativePlus[face],
                                                timeDerivativeMinus[face],
                                                godunov[face],
                                                timeDerivativePlus[prefetchFace],
                                                timeDerivativeMinus[prefetchFace],
                                                absoluteSlip[face] );

    e_interoperability.evaluateFrictionLaw( static_cast<int>(faceInformation[face].meshFace),
                                            godunov[face],
                                            imposedStatePlus[face],
                                            imposedStateMinus[face],
                                            absoluteSlip[face],
                                            m_fullUpdateTime,
                                            m_dynamicRuptureKernel.timePoints,
                                            m_dynamicRuptureKernel.timeWeights,
                                            waveSpeedsPlus[face],
                                            waveSpeedsMinus[face] );
  }

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells());
}


void seissol::time_stepping::TimeCluster::computeDynamicRuptureFlops( seissol::initializers::Layer& layerData,
                                                                      long long&                    nonZeroFlops,
                                                                      long long&                    hardwareFlops )
{
  nonZeroFlops = 0;
  hardwareFlops = 0;

  DRFaceInformation* faceInformation = layerData.var(m_dynRup->faceInformation);

  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    long long faceNonZeroFlops, faceHardwareFlops;
    m_dynamicRuptureKernel.flopsGodunovState( faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
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
                   seissol::MPI::mpi.comm(),                               // communicator
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
                   seissol::MPI::mpi.comm(),                            // communicator
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

void seissol::time_stepping::TimeCluster::computeLocalIntegration( seissol::initializers::Layer&  i_layerData ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // local integration buffer
  real l_integrationBuffer[NUMBER_OF_ALIGNED_DOFS] __attribute__((aligned(ALIGNMENT)));

  // pointer for the call of the ADER-function
  real *l_bufferPointer;

  real                (*dofs)[NUMBER_OF_ALIGNED_DOFS] = i_layerData.var(m_lts->dofs);
  real**                buffers                       = i_layerData.var(m_lts->buffers);
  real**                derivatives                   = i_layerData.var(m_lts->derivatives);
  LocalIntegrationData* localIntegration              = i_layerData.var(m_lts->localIntegration);
  CellLocalInformation* cellInformation               = i_layerData.var(m_lts->cellInformation);
  real**                displacements                 = i_layerData.var(m_lts->displacements);

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer) schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
    // overwrite cell buffer
    // TODO: Integrate this step into the kernel

    bool l_buffersProvided = (cellInformation[l_cell].ltsSetup >> 8)%2 == 1; // buffers are provided
    bool l_resetBuffers = l_buffersProvided && ( (cellInformation[l_cell].ltsSetup >> 10) %2 == 0 || m_resetLtsBuffers ); // they should be reset

    if( l_resetBuffers ) {
      // assert presence of the buffer
      assert( buffers[l_cell] != NULL );

      l_bufferPointer = buffers[l_cell];
    }
    // work on local buffer
    else {
      l_bufferPointer = l_integrationBuffer;
    }

    m_timeKernel.computeAder( m_timeStepWidth,
                              m_globalData,
                              &localIntegration[l_cell],
                              dofs[l_cell],
                              l_bufferPointer,
                              derivatives[l_cell] );
    m_localKernel.computeIntegral( cellInformation[l_cell].faceTypes,
                                   m_globalData,
                                   &localIntegration[l_cell],
                                   l_bufferPointer,
                                   dofs[l_cell] );

    if (displacements[l_cell] != NULL) {
      for (unsigned dof = 0; dof < NUMBER_OF_ALIGNED_VELOCITY_DOFS; ++dof) {
        displacements[l_cell][dof] += l_bufferPointer[NUMBER_OF_ALIGNED_STRESS_DOFS + dof];
      }
    }

    // update lts buffers if required
    // TODO: Integrate this step into the kernel
    if( !l_resetBuffers && l_buffersProvided ) {
      // assert presence of the buffer
      assert( buffers[l_cell] != NULL );

      for( unsigned int l_dof = 0; l_dof < NUMBER_OF_ALIGNED_DOFS; l_dof++ ) {
        buffers[l_cell][l_dof] += l_integrationBuffer[l_dof];
      }
    }
  }

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells());
}

void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData ) {
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  real                      (*dofs)[NUMBER_OF_ALIGNED_DOFS] = i_layerData.var(m_lts->dofs);
  real*                     (*faceNeighbors)[4]             = i_layerData.var(m_lts->faceNeighbors);
  CellDRMapping             (*drMapping)[4]                 = i_layerData.var(m_lts->drMapping);
  NeighboringIntegrationData* neighboringIntegration        = i_layerData.var(m_lts->neighboringIntegration);
  CellLocalInformation*       cellInformation               = i_layerData.var(m_lts->cellInformation);
#ifdef USE_PLASTICITY
  PlasticityData*             plasticity                    = i_layerData.var(m_lts->plasticity);
  real                      (*energy)[3]                    = i_layerData.var(m_lts->energy);
  real                      (*pstrain)[7]                   = i_layerData.var(m_lts->pstrain);
  unsigned                   numberOTetsWithPlasticYielding = 0;
#endif

  real *l_timeIntegrated[4];
  real *l_faceNeighbors_prefetch[4];

#ifdef _OPENMP
#ifdef USE_PLASTICITY
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch) reduction(+:numberOTetsWithPlasticYielding)
#else
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch)
#endif
#endif
  for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
    seissol::kernels::TimeCommon::computeIntegrals( m_timeKernel,
                                                    cellInformation[l_cell].ltsSetup,
                                                    cellInformation[l_cell].faceTypes,
                                                    m_subTimeStart,
                                                    m_timeStepWidth,
                                                    faceNeighbors[l_cell],
#ifdef _OPENMP
                                                    *reinterpret_cast<real (*)[4][NUMBER_OF_ALIGNED_DOFS]>(&(m_globalData->integrationBufferLTS[omp_get_thread_num()*4*NUMBER_OF_ALIGNED_DOFS])),
#else
                                                    *reinterpret_cast<real (*)[4][NUMBER_OF_ALIGNED_DOFS]>(m_globalData->integrationBufferLTS),
#endif
                                                    l_timeIntegrated );

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
    l_faceNeighbors_prefetch[0] = (cellInformation[l_cell].faceTypes[1] != dynamicRupture) ? faceNeighbors[l_cell][1] : drMapping[l_cell][1].godunov;
    l_faceNeighbors_prefetch[1] = (cellInformation[l_cell].faceTypes[2] != dynamicRupture) ? faceNeighbors[l_cell][2] : drMapping[l_cell][2].godunov;
    l_faceNeighbors_prefetch[2] = (cellInformation[l_cell].faceTypes[3] != dynamicRupture) ? faceNeighbors[l_cell][3] : drMapping[l_cell][3].godunov;

    // fourth face's prefetches
    if (l_cell < (i_layerData.getNumberOfCells()-1) ) {
      l_faceNeighbors_prefetch[3] = (cellInformation[l_cell+1].faceTypes[0] != dynamicRupture) ? faceNeighbors[l_cell+1][0] : drMapping[l_cell+1][0].godunov;
    } else {
      l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
    }
#endif

    m_neighborKernel.computeNeighborsIntegral( cellInformation[l_cell].faceTypes,
                                               cellInformation[l_cell].faceRelations,
                                               drMapping[l_cell],
                                               m_globalData,
                                               &neighboringIntegration[l_cell],
                                               l_timeIntegrated,
#ifdef ENABLE_MATRIX_PREFETCH
                                               l_faceNeighbors_prefetch,
#endif
                                               dofs[l_cell] );

#ifdef USE_PLASTICITY
  numberOTetsWithPlasticYielding += seissol::kernels::Plasticity::computePlasticity( m_relaxTime,
                                                                                     m_timeStepWidth,
                                                                                     m_globalData,
                                                                                     &plasticity[l_cell],
                                                                                     dofs[l_cell],
                                                                                     pstrain[l_cell] );
#endif
#ifdef INTEGRATE_QUANTITIES
  seissol::SeisSol::main.postProcessor().integrateQuantities( m_timeStepWidth,
                                                              i_layerData,
                                                              l_cell,
                                                              dofs[l_cell] );
#endif // INTEGRATE_QUANTITIES
  }

  #ifdef USE_PLASTICITY
  g_SeisSolNonZeroFlopsPlasticity += i_layerData.getNumberOfCells() * m_flops_nonZero[PlasticityCheck] + numberOTetsWithPlasticYielding * m_flops_nonZero[PlasticityYield];
  g_SeisSolHardwareFlopsPlasticity += i_layerData.getNumberOfCells() * m_flops_hardware[PlasticityCheck] + numberOTetsWithPlasticYielding * m_flops_hardware[PlasticityYield];
  #endif

  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells());
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
  if( m_updatable.localInterior ) {
    writeReceivers();
  }

  // integrate copy layer locally
  computeLocalIntegration( m_clusterData->child<Copy>() );

  g_SeisSolNonZeroFlopsLocal += m_flops_nonZero[LocalCopy];
  g_SeisSolHardwareFlopsLocal += m_flops_hardware[LocalCopy];

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
  if( m_updatable.localCopy ) {
    writeReceivers();
  }
#else
  // non-MPI checks for write in the interior
  writeReceivers();
#endif

  // integrate interior cells locally
  computeLocalIntegration( m_clusterData->child<Interior>() );

  g_SeisSolNonZeroFlopsLocal += m_flops_nonZero[LocalInterior];
  g_SeisSolHardwareFlopsLocal += m_flops_hardware[LocalInterior];

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

  if (m_dynamicRuptureFaces == true) {
    if (m_updatable.neighboringInterior) {
      computeDynamicRupture(m_dynRupClusterData->child<Interior>());
      g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawInterior];
      g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawInterior];
    }

    computeDynamicRupture(m_dynRupClusterData->child<Copy>());
    g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawCopy];
    g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawCopy];
  }

  computeNeighboringIntegration( m_clusterData->child<Copy>() );

  g_SeisSolNonZeroFlopsNeighbor += m_flops_nonZero[NeighborCopy];
  g_SeisSolHardwareFlopsNeighbor += m_flops_hardware[NeighborCopy];
  g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRNeighborCopy];
  g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRNeighborCopy];

#ifndef USE_COMM_THREAD
  // continue with communication
  testForCopyLayerSends();
#endif

  // compute dynamic rupture, update simulation time and statistics
  if( !m_updatable.neighboringInterior ) {
    // First cluster calls fault receiver output
    // TODO: Change from iteration based to time based
    if (m_clusterId == 0) {
      e_interoperability.faultOutput( m_fullUpdateTime, m_timeStepWidth );
    }

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

  if (m_dynamicRuptureFaces == true && m_updatable.neighboringCopy == true) {
    computeDynamicRupture(m_dynRupClusterData->child<Interior>());
    g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawInterior];
    g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawInterior];
  }

  // Update all cells in the interior with the neighboring boundary contribution.
  computeNeighboringIntegration( m_clusterData->child<Interior>() );

  g_SeisSolNonZeroFlopsNeighbor += m_flops_nonZero[NeighborInterior];
  g_SeisSolHardwareFlopsNeighbor += m_flops_hardware[NeighborInterior];
  g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRNeighborInterior];
  g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRNeighborInterior];

  // compute dynamic rupture, update simulation time and statistics
  if( !m_updatable.neighboringCopy ) {
    // First cluster calls fault receiver output
    // TODO: Change from iteration based to time based
    if (m_clusterId == 0) {
      e_interoperability.faultOutput( m_fullUpdateTime, m_timeStepWidth );
    }

    m_fullUpdateTime      += m_timeStepWidth;
    m_subTimeStart        += m_timeStepWidth;
    m_numberOfFullUpdates += 1;
    m_numberOfTimeSteps   += 1;
  }

  // update finished
  m_updatable.neighboringInterior = false;
}

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops( unsigned                    numberOfCells,
                                                                        CellLocalInformation const* cellInformation,
                                                                        long long&                  nonZeroFlops,
                                                                        long long&                  hardwareFlops  )
{
  nonZeroFlops = 0;
  hardwareFlops = 0;

  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    unsigned cellNonZero, cellHardware;
    m_timeKernel.flopsAder(cellNonZero, cellHardware);
    nonZeroFlops += cellNonZero;
    hardwareFlops += cellHardware;
    m_localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
    nonZeroFlops += cellNonZero;
    hardwareFlops += cellHardware;
  }
}

void seissol::time_stepping::TimeCluster::computeNeighborIntegrationFlops(  unsigned                    numberOfCells,
                                                                            CellLocalInformation const* cellInformation,
                                                                            long long&                  nonZeroFlops,
                                                                            long long&                  hardwareFlops,
                                                                            long long&                  drNonZeroFlops,
                                                                            long long&                  drHardwareFlops )
{
  nonZeroFlops = 0;
  hardwareFlops = 0;
  drNonZeroFlops = 0;
  drHardwareFlops = 0;

  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    unsigned cellNonZero, cellHardware;
    long long cellDRNonZero, cellDRHardware;
    m_neighborKernel.flopsNeighborsIntegral(  cellInformation[cell].faceTypes,
                                              cellInformation[cell].faceRelations,
                                              cellNonZero,
                                              cellHardware,
                                              cellDRNonZero,
                                              cellDRHardware );
    nonZeroFlops += cellNonZero;
    hardwareFlops += cellHardware;
    drNonZeroFlops += cellDRNonZero;
    drHardwareFlops += cellDRHardware;

    /// \todo add lts time integration
    /// \todo add plasticity
  }
}

void seissol::time_stepping::TimeCluster::computeFlops()
{
#ifdef USE_MPI
  computeLocalIntegrationFlops( m_meshStructure->numberOfCopyCells,
                                m_clusterData->child<Copy>().var(m_lts->cellInformation),
                                m_flops_nonZero[LocalCopy],
                                m_flops_hardware[LocalCopy] );
#endif

  computeLocalIntegrationFlops( m_meshStructure->numberOfInteriorCells,
                                m_clusterData->child<Interior>().var(m_lts->cellInformation),
                                m_flops_nonZero[LocalInterior],
                                m_flops_hardware[LocalInterior] );

#ifdef USE_MPI
  computeNeighborIntegrationFlops(  m_meshStructure->numberOfCopyCells,
                                    m_clusterData->child<Copy>().var(m_lts->cellInformation),
                                    m_flops_nonZero[NeighborCopy],
                                    m_flops_hardware[NeighborCopy],
                                    m_flops_nonZero[DRNeighborCopy],
                                    m_flops_hardware[DRNeighborCopy] );
#endif

  computeNeighborIntegrationFlops(  m_meshStructure->numberOfInteriorCells,
                                    m_clusterData->child<Interior>().var(m_lts->cellInformation),
                                    m_flops_nonZero[NeighborInterior],
                                    m_flops_hardware[NeighborInterior],
                                    m_flops_nonZero[DRNeighborInterior],
                                    m_flops_hardware[DRNeighborInterior] );

  computeDynamicRuptureFlops( m_dynRupClusterData->child<Copy>(), m_flops_nonZero[DRFrictionLawCopy], m_flops_hardware[DRFrictionLawCopy] );
  computeDynamicRuptureFlops( m_dynRupClusterData->child<Interior>(), m_flops_nonZero[DRFrictionLawInterior], m_flops_hardware[DRFrictionLawInterior] );

  seissol::kernels::Plasticity::flopsPlasticity(  m_flops_nonZero[PlasticityCheck],
                                                  m_flops_hardware[PlasticityCheck],
                                                  m_flops_nonZero[PlasticityYield],
                                                  m_flops_hardware[PlasticityYield] );
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

