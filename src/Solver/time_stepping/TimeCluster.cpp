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
#include <Kernels/Receiver.h>
#include <Monitoring/FlopCounter.hpp>

#include <cassert>
#include <cstring>

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
extern volatile unsigned int* volatile g_handleRecvs;
extern volatile unsigned int* volatile g_handleSends;
#endif

//! fortran interoperability
extern seissol::Interoperability e_interoperability;

seissol::time_stepping::TimeCluster::TimeCluster(unsigned int i_clusterId,
                                                 unsigned int i_globalClusterId,
                                                 double maxTimeStepSize,
                                                 int timeStepRate,
                                                 double timeTolerance,
                                                 struct GlobalData *i_globalData,
                                                 seissol::initializers::Layer *i_clusterData,
                                                 seissol::initializers::Layer *i_dynRupClusterData,
                                                 seissol::initializers::LTS *i_lts,
                                                 seissol::initializers::DynamicRupture *i_dynRup,
                                                 LoopStatistics *i_loopStatistics) :
    AbstractTimeCluster(maxTimeStepSize, timeTolerance, timeStepRate),
    // cluster ids
    m_clusterId(i_clusterId),
    m_globalClusterId(i_globalClusterId),
    // global data
    m_globalData(i_globalData),
    m_clusterData(i_clusterData),
    m_dynRupClusterData(i_dynRupClusterData),
    m_lts(i_lts),
    m_dynRup(i_dynRup),
    // cells
    m_cellToPointSources(nullptr),
    m_numberOfCellToPointSourcesMappings(0),
    m_pointSources(nullptr),

    m_loopStatistics(i_loopStatistics),
    m_receiverCluster(nullptr)
{
    assert( m_globalData                               != nullptr);
    assert( m_clusterData                              != nullptr);

  // default: no updates are allowed
#ifdef USE_MPI
  //m_sendLtsBuffers                = false;
#endif
  // set timings to zero
  m_receiverTime                  = 0;

  m_dynamicRuptureFaces = i_dynRupClusterData->getNumberOfCells() > 0;

  m_timeKernel.setGlobalData(m_globalData);
  m_localKernel.setGlobalData(m_globalData);
  m_localKernel.setInitConds(&e_interoperability.getInitialConditions());
  m_neighborKernel.setGlobalData(m_globalData);
  m_dynamicRuptureKernel.setGlobalData(m_globalData);

  computeFlops();

  m_regionComputeLocalIntegration = m_loopStatistics->getRegion("computeLocalIntegration");
  m_regionComputeNeighboringIntegration = m_loopStatistics->getRegion("computeNeighboringIntegration");
  m_regionComputeDynamicRupture = m_loopStatistics->getRegion("computeDynamicRupture");
}

seissol::time_stepping::TimeCluster::~TimeCluster() {
#ifndef NDEBUG
  logInfo() << "#(time steps):" << numberOfTimeSteps;
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

void seissol::time_stepping::TimeCluster::writeReceivers() {
  SCOREP_USER_REGION( "writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION )

  if (m_receiverCluster != nullptr) {
    m_receiverTime = m_receiverCluster->calcReceivers(m_receiverTime, ct.correctionTime, timeStepSize());
  }
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
                                                       ct.correctionTime,
                                                       ct.correctionTime + timeStepSize(),
                                                       *m_cellToPointSources[mapping].dofs );
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceFSRM( m_pointSources->mInvJInvPhisAtSources[source],
                                                        m_pointSources->tensor[source],
                                                        &m_pointSources->slipRates[source][0],
                                                        ct.correctionTime,
                                                        ct.correctionTime + timeStepSize(),
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
  real                                (*imposedStatePlus)[tensor::QInterpolated::size()]                  = layerData.var(m_dynRup->imposedStatePlus);
  real                                (*imposedStateMinus)[tensor::QInterpolated::size()]                 = layerData.var(m_dynRup->imposedStateMinus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                                    = layerData.var(m_dynRup->waveSpeedsPlus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                                   = layerData.var(m_dynRup->waveSpeedsMinus);

  alignas(ALIGNMENT) real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  alignas(ALIGNMENT) real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];

  m_dynamicRuptureKernel.setTimeStepWidth(timeStepSize());
#ifdef _OPENMP
  #pragma omp parallel for schedule(static) private(QInterpolatedPlus,QInterpolatedMinus)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynamicRuptureKernel.spaceTimeInterpolation(  faceInformation[face],
                                                    m_globalData,
                                                   &godunovData[face],
                                                    timeDerivativePlus[face],
                                                    timeDerivativeMinus[face],
                                                    QInterpolatedPlus,
                                                    QInterpolatedMinus,
                                                    timeDerivativePlus[prefetchFace],
                                                    timeDerivativeMinus[prefetchFace] );

    e_interoperability.evaluateFrictionLaw( static_cast<int>(faceInformation[face].meshFace),
                                            QInterpolatedPlus,
                                            QInterpolatedMinus,
                                            imposedStatePlus[face],
                                            imposedStateMinus[face],
                                            ct.correctionTime,
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
  /*
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )

  //
  // Receive data of the ghost regions
  //
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
   */
}

void seissol::time_stepping::TimeCluster::sendCopyLayer(){
  /*
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )

  //
  //Send data of the copy regions
  //
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
  */
}

bool seissol::time_stepping::TimeCluster::testForGhostLayerReceives(){
  /*
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
   */
  return false;
}

bool seissol::time_stepping::TimeCluster::testForCopyLayerSends(){
  /*
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
   */
  return false;
}

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
void seissol::time_stepping::TimeCluster::initReceiveGhostLayer(){
  //g_handleRecvs[m_clusterId] = 1;
}

void seissol::time_stepping::TimeCluster::initSendCopyLayer(){
  //g_handleSends[m_clusterId] = 1;
}

void seissol::time_stepping::TimeCluster::waitForInits() {
  //while( g_handleRecvs[m_clusterId] == 1 || g_handleSends[m_clusterId] == 1 );
}
#endif

#endif

void seissol::time_stepping::TimeCluster::computeLocalIntegration( seissol::initializers::Layer&  i_layerData, bool resetBuffers ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // local integration buffer
  real l_integrationBuffer[tensor::I::size()] __attribute__((aligned(ALIGNMENT)));

  // pointer for the call of the ADER-function
  real* l_bufferPointer;

  real** buffers = i_layerData.var(m_lts->buffers);
  real** derivatives = i_layerData.var(m_lts->derivatives);
  real** displacements = i_layerData.var(m_lts->displacements);
  CellMaterialData* materialData = i_layerData.var(m_lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp{};

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer, tmp) schedule(static)
#endif
  for (unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++) {
    auto data = loader.entry(l_cell);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided = (data.cellInformation.ltsSetup >> 8) % 2 == 1; // buffers are provided
    const bool resetMyBuffers = buffersProvided && ( (data.cellInformation.ltsSetup >> 10) %2 == 0 || resetBuffers ); // they should be reset

    if (resetMyBuffers) {
      // assert presence of the buffer
      assert(buffers[l_cell] != nullptr);

      l_bufferPointer = buffers[l_cell];
    } else {
      // work on local buffer
      l_bufferPointer = l_integrationBuffer;
    }

    m_timeKernel.computeAder(timeStepSize(),
                             data,
                             tmp,
                             l_bufferPointer,
                             derivatives[l_cell]);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping (*boundaryMapping)[4] = i_layerData.var(m_lts->boundaryMapping);
    m_localKernel.computeIntegral(l_bufferPointer,
                                  data,
                                  tmp,
                                  &materialData[l_cell],
                                  &boundaryMapping[l_cell],
                                  ct.correctionTime,
                                  timeStepSize()
    );
    
    // Update displacement
    if (displacements[l_cell] != nullptr) {
      kernel::addVelocity krnl;
      krnl.I = l_bufferPointer;
      krnl.selectVelocity = init::selectVelocity::Values;
      krnl.displacement = displacements[l_cell];
      krnl.execute();
    }

    // TODO: Integrate this step into the kernel
    // We've used a temporary buffer -> need to accumulate update in
    // shared buffer.
    if (!resetMyBuffers && buffersProvided) {
      assert(buffers[l_cell] != nullptr);

      for (unsigned int l_dof = 0; l_dof < tensor::I::size(); ++l_dof) {
        buffers[l_cell][l_dof] += l_integrationBuffer[l_dof];
      }
    }
  }

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells());
}

void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData, double subTimeStart ) {
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
  CellDRMapping (*drMapping)[4] = i_layerData.var(m_lts->drMapping);
  CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);
#ifdef USE_PLASTICITY
  updateRelaxTime();
  PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
  real (*pstrain)[7] = i_layerData.var(m_lts->pstrain);
  unsigned numberOTetsWithPlasticYielding = 0;
#endif

  kernels::NeighborData::Loader loader;
  loader.load(*m_lts, i_layerData);

  real *l_timeIntegrated[4];
  real *l_faceNeighbors_prefetch[4];

  //std::cout << "dt_max=" << ct.maxTimeStepSize << " dt=" << timeStepSize() << " t_sub=" << subTimeStart << std::endl;

#ifdef _OPENMP
#ifdef USE_PLASTICITY
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch) reduction(+:numberOTetsWithPlasticYielding)
#else
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch)
#endif
#endif
  for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
    auto data = loader.entry(l_cell);
    seissol::kernels::TimeCommon::computeIntegrals(m_timeKernel,
                                                   data.cellInformation.ltsSetup,
                                                   data.cellInformation.faceTypes,
                                                   subTimeStart,
                                                   timeStepSize(),
                                                   faceNeighbors[l_cell],
#ifdef _OPENMP
                                                   *reinterpret_cast<real (*)[4][tensor::I::size()]>(&(m_globalData->integrationBufferLTS[omp_get_thread_num()*4*tensor::I::size()])),
#else
                                                   *reinterpret_cast<real (*)[4][tensor::I::size()]>(m_globalData->integrationBufferLTS),
#endif
                                                   l_timeIntegrated);

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
    l_faceNeighbors_prefetch[0] = (cellInformation[l_cell].faceTypes[1] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][1] :
      drMapping[l_cell][1].godunov;
    l_faceNeighbors_prefetch[1] = (cellInformation[l_cell].faceTypes[2] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][2] :
      drMapping[l_cell][2].godunov;
    l_faceNeighbors_prefetch[2] = (cellInformation[l_cell].faceTypes[3] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][3] :
      drMapping[l_cell][3].godunov;

    // fourth face's prefetches
    if (l_cell < (i_layerData.getNumberOfCells()-1) ) {
      l_faceNeighbors_prefetch[3] = (cellInformation[l_cell+1].faceTypes[0] != FaceType::dynamicRupture) ?
	faceNeighbors[l_cell+1][0] :
	drMapping[l_cell+1][0].godunov;
    } else {
      l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
    }
#endif

    m_neighborKernel.computeNeighborsIntegral( data,
                                               drMapping[l_cell],
#ifdef ENABLE_MATRIX_PREFETCH
                                               l_timeIntegrated, l_faceNeighbors_prefetch
#else
                                               l_timeIntegrated
#endif
                                               );

#ifdef USE_PLASTICITY
  numberOTetsWithPlasticYielding += seissol::kernels::Plasticity::computePlasticity( m_relaxTime,
                                                                                     timeStepSize(),
                                                                                     m_globalData,
                                                                                     &plasticity[l_cell],
                                                                                     data.dofs,
                                                                                     pstrain[l_cell] );
#endif
#ifdef INTEGRATE_QUANTITIES
  seissol::SeisSol::main.postProcessor().integrateQuantities( timeStepSize(),
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
  /*
  SCOREP_USER_REGION( "computeLocalCopy", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.localCopy ) {
    logError() << "Invalid call of computeLocalCopy, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << timeStepSize()   << m_subTimeStart      << m_resetLtsBuffers;
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
  computeLocalIntegration(*m_clusterData);

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
    m_predictionTime += timeStepSize();
  }

  // update finished
  m_updatable.localCopy  = false;

  // wait until communication thread finished initializing the receives
#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  waitForInits();
#endif

  return true;
  */
  return false;
}
#endif

void seissol::time_stepping::TimeCluster::computeLocalInterior(){
  /*
  SCOREP_USER_REGION( "computeLocalInterior", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.localInterior ) {
    logError() << "Invalid call of computeLocalInterior, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << timeStepSize()   << m_subTimeStart      << m_resetLtsBuffers;
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
    m_predictionTime += timeStepSize();
  }

  // update finished
  m_updatable.localInterior = false;
  */
}

#ifdef USE_MPI
bool seissol::time_stepping::TimeCluster::computeNeighboringCopy() {
  /*
  SCOREP_USER_REGION( "computeNeighboringCopy", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.neighboringCopy ) {
    logError() << "Invalid call of computeNeighboringCopy aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << timeStepSize()   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // continue only of ghost layer receives are complete
  if( !testForGhostLayerReceives() ) return false;

#ifndef USE_COMM_THREAD
  // continue with communication
  testForCopyLayerSends();
#endif

  if (m_dynamicRuptureFaces) {
    if (m_updatable.neighboringInterior) {
      computeDynamicRupture(*m_dynRupClusterData);
      g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawInterior];
      g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawInterior];
    }

    computeDynamicRupture(*m_dynRupClusterData);
    g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawCopy];
    g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawCopy];
  }

  //computeNeighboringIntegration(*m_clusterData);

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
      e_interoperability.faultOutput( m_fullUpdateTime, timeStepSize() );
    }

    m_fullUpdateTime      += timeStepSize();
    m_subTimeStart        += timeStepSize();
    m_numberOfFullUpdates += 1;
    m_numberOfTimeSteps   += 1;
  }

  // update finished
  m_updatable.neighboringCopy = false;

  return true;
   */
  return false;
}
#endif

void seissol::time_stepping::TimeCluster::computeNeighboringInterior() {
  /*
  SCOREP_USER_REGION( "computeNeighboringInterior", SCOREP_USER_REGION_TYPE_FUNCTION )
  // ensure a valid call
  if( !m_updatable.neighboringInterior ) {
    logError() << "Invalid call of computeNeighboringInterior, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << timeStepSize()   << m_subTimeStart      << m_resetLtsBuffers;
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
      e_interoperability.faultOutput( m_fullUpdateTime, timeStepSize() );
    }

    m_fullUpdateTime      += timeStepSize();
    m_subTimeStart        += timeStepSize();
    m_numberOfFullUpdates += 1;
    m_numberOfTimeSteps   += 1;
  }

  // update finished
  m_updatable.neighboringInterior = false;
   */
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
    // TODO(Lukas) Maybe include avg. displacement computation here at some point.
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
                                                                            CellDRMapping const       (*drMapping)[4],
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
                                              drMapping[cell],
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
  /*
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
                                    m_clusterData->child<Copy>().var(m_lts->drMapping),
                                    m_flops_nonZero[NeighborCopy],
                                    m_flops_hardware[NeighborCopy],
                                    m_flops_nonZero[DRNeighborCopy],
                                    m_flops_hardware[DRNeighborCopy] );
#endif

  computeNeighborIntegrationFlops(  m_meshStructure->numberOfInteriorCells,
                                    m_clusterData->child<Interior>().var(m_lts->cellInformation),
                                    m_clusterData->child<Interior>().var(m_lts->drMapping),
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
*/
   }

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
/*
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
 */
#endif

namespace seissol::time_stepping {
void TimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  if (neighborCluster.ct.maxTimeStepSize > ct.maxTimeStepSize) {
    lastSubTime = neighborCluster.ct.correctionTime;
  }
}
void TimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster&) {
  // Doesn't do anything
}
void TimeCluster::predict() {
  assert(state == ActorState::Corrected);
  bool resetBuffers = true;
  for (auto& neighbor : neighbors) {
    if (neighbor.ct.maxTimeStepSize > ct.maxTimeStepSize &&
        // ct.correctionTime > neighbor.ct.correctionTime
        ct.correctionTime - neighbor.ct.correctionTime > timeTolerance) {
      resetBuffers = false;
      /*
      std::cout
      << "Resetting predict buffers, our rate > their rate = "
      << int{ct.timeStepRate > neighbor.ct.timeStepRate}
      << " timeStepFractions = " << rateDiff
      << " our steps / rateFractions" << ct.stepsSinceLastSync / rateDiff
      << " their steps = " << neighbor.ct.stepsSinceLastSync
      << std::endl;
       */
      break;
    }
  }
  /*
  if( clusters[l_cluster]->m_numberOfFullUpdates % m_timeStepping.globalTimeStepRates[l_globalClusterId] == 0 ) {
   */
  resetBuffers = resetBuffersOld; // TODO(Lukas) :(
  assert(resetBuffers == resetBuffersOld);
  if (numberOfTimeSteps == 0) {
    assert(resetBuffers);
  }

  writeReceivers(); // TODO(Lukas) Is this correct?
  computeLocalIntegration(*m_clusterData, resetBuffers);
  computeSources();
  std::cout << m_globalClusterId << ": Predicted to t=" << ct.correctionTime + timeStepSize() << std::endl;
}
void TimeCluster::correct() {
  assert(state == ActorState::Predicted);

  double subTimeStart = ct.correctionTime - lastSubTime;
  // Note, the following is likely wrong
  // TODO(Lukas) Check scheduling if this is correct
  if (m_dynamicRuptureFaces > 0) {
    computeDynamicRupture(*m_dynRupClusterData);
  }
  computeNeighboringIntegration(*m_clusterData, subTimeStart);

  // First cluster calls fault receiver output
  // TODO: Change from iteration based to time based
  if (m_clusterId == 0) {
    //e_interoperability.faultOutput(ct.correctionTime + timeStepSize(), timeStepSize());
  }
  std::cout << m_globalClusterId << ": Corrected to t=" << ct.correctionTime + timeStepSize() << std::endl;
}
}
