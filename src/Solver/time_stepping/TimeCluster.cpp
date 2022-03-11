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
                                                  bool usePlasticity,
                                                  MeshStructure                 *i_meshStructure,
                                                  CompoundGlobalData             i_globalData,
                                                  seissol::initializers::TimeCluster* i_clusterData,
                                                  seissol::initializers::TimeCluster* i_dynRupClusterData,
                                                  seissol::initializers::LTS*         i_lts,
                                                  seissol::initializers::DynamicRupture* i_dynRup,
                                                  LoopStatistics*                        i_loopStatistics ):
 // cluster ids
 m_clusterId(               i_clusterId                ),
 m_globalClusterId(         i_globalClusterId          ),
 usePlasticity(usePlasticity),
 // mesh structure
 m_meshStructure(           i_meshStructure            ),
 // global data
 m_globalDataOnHost( i_globalData.onHost ),
 m_globalDataOnDevice(i_globalData.onDevice ),
 m_clusterData(             i_clusterData              ),
 m_dynRupClusterData(       i_dynRupClusterData        ),
 m_lts(                     i_lts                      ),
 m_dynRup(                  i_dynRup                   ),
 // cells
 m_cellToPointSources(      NULL                       ),
 m_numberOfCellToPointSourcesMappings(0                ),
 m_pointSources(            NULL                       ),

 m_loopStatistics(          i_loopStatistics           ),
 m_receiverCluster(          nullptr                   )
{
    // assert all pointers are valid
    assert( m_meshStructure                            != nullptr );
    assert( m_clusterData                              != NULL );
    assert( m_globalDataOnHost                         != nullptr );
    if constexpr (seissol::isDeviceOn()) {
        assert( m_globalDataOnDevice                   != nullptr );
    }

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
  m_timeStepWidth                 = 0;
  m_subTimeStart                  = 0;
  m_numberOfFullUpdates           = 0;
  m_fullUpdateTime                = 0;
  m_predictionTime                = 0;

  m_dynamicRuptureFaces = (i_dynRupClusterData->child<Ghost>().getNumberOfCells() > 0)
	|| (i_dynRupClusterData->child<Copy>().getNumberOfCells() > 0)
	|| (i_dynRupClusterData->child<Interior>().getNumberOfCells() > 0);
  
  m_timeKernel.setGlobalData(i_globalData);
  m_localKernel.setGlobalData(i_globalData);
  m_localKernel.setInitConds(&e_interoperability.getInitialConditions());
  m_neighborKernel.setGlobalData(i_globalData);
  m_dynamicRuptureKernel.setGlobalData(i_globalData);

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

void seissol::time_stepping::TimeCluster::writeReceivers() {
  SCOREP_USER_REGION( "writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION )

  if (m_receiverCluster != nullptr) {
    m_receiverTime = m_receiverCluster->calcReceivers(m_receiverTime, m_fullUpdateTime, m_timeStepWidth);
  }
}

void seissol::time_stepping::TimeCluster::computeSources() {
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
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
                                                       m_pointSources->A[source],
                                                       m_pointSources->stiffnessTensor[source],
                                                       m_pointSources->slipRates[source],
                                                       m_fullUpdateTime,
                                                       m_fullUpdateTime + m_timeStepWidth,
                                                       *m_cellToPointSources[mapping].dofs );
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
  
          sourceterm::addTimeIntegratedPointSourceFSRM( m_pointSources->mInvJInvPhisAtSources[source],
                                                        m_pointSources->tensor[source],
                                                        m_pointSources->slipRates[source][0],
                                                        m_fullUpdateTime,
                                                        m_fullUpdateTime + m_timeStepWidth,
                                                        *m_cellToPointSources[mapping].dofs );
        }
      }
    }
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
}

#ifndef ACL_DEVICE
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

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) private(QInterpolatedPlus,QInterpolatedMinus)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynamicRuptureKernel.spaceTimeInterpolation(  faceInformation[face],
                                                    m_globalDataOnHost,
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
                                            m_fullUpdateTime,
                                            m_dynamicRuptureKernel.timePoints,
                                            m_dynamicRuptureKernel.timeWeights,
                                            waveSpeedsPlus[face],
                                            waveSpeedsMinus[face] );
  }

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells());
}
#else

void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData ) {
  device.api->putProfilingMark("computeDynamicRupture", device::ProfilingColors::Cyan);
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  if (layerData.getNumberOfCells() > 0) {
    // compute space time interpolation part
    ConditionalBatchTableT &table = layerData.getCondBatchTable();
    m_dynamicRuptureKernel.batchedSpaceTimeInterpolation(table);

    // compute friction part
    using namespace dr::pipeline;
    DrContext context;
    context.QInterpolatedPlusOnDevice = static_cast<real *>(layerData.getScratchpadMemory(m_dynRup->QInterpolatedPlusOnDevice));
    context.QInterpolatedMinusOnDevice = static_cast<real *>(layerData.getScratchpadMemory(m_dynRup->QInterpolatedMinusOnDevice));
    context.QInterpolatedPlusOnHost = static_cast<DrContext::QInterpolatedPtrT>(layerData.getScratchpadMemory(m_dynRup->QInterpolatedPlusOnHost));
    context.QInterpolatedMinusOnHost = static_cast<DrContext::QInterpolatedPtrT>(layerData.getScratchpadMemory(m_dynRup->QInterpolatedMinusOnHost));

    context.faceInformation = layerData.var(m_dynRup->faceInformation);
    context.devImposedStatePlus = layerData.var(m_dynRup->imposedStatePlus);
    context.devImposedStateMinus = layerData.var(m_dynRup->imposedStateMinus);
    context.waveSpeedsPlus = layerData.var(m_dynRup->waveSpeedsPlus);
    context.waveSpeedsMinus = layerData.var(m_dynRup->waveSpeedsMinus);

    context.imposedStatePlusOnHost = static_cast<DrContext::imposedStatePlusT>(layerData.getScratchpadMemory(m_dynRup->imposedStatePlusOnHost));
    context.imposedStateMinusOnHost = static_cast<DrContext::imposedStatePlusT>(layerData.getScratchpadMemory(m_dynRup->imposedStateMinusOnHost));

    struct AsyncCopyFrom : public DrBaseCallBack {
      explicit AsyncCopyFrom(DrContext userContext, TimeCluster *cluster) : DrBaseCallBack(userContext, cluster) {
        assert(context.QInterpolatedPlusOnDevice != nullptr);
        assert(context.QInterpolatedMinusOnDevice != nullptr);
        assert(context.QInterpolatedPlusOnHost != nullptr);
        assert(context.QInterpolatedMinusOnHost != nullptr);
        copyFromStream = device.api->getNextCircularStream();
      }

      void operator()(size_t begin, size_t batchSize, size_t callCounter) override {
        constexpr size_t QInterpolatedSize = CONVERGENCE_ORDER * tensor::QInterpolated::size();
        const size_t upperStageOffset = (callCounter % DrPipeline::TailSize) * DrPipeline::DefaultBatchSize;

        device.api->syncStreamFromCircularBuffer(copyFromStream);
        device.api->copyFromAsync(reinterpret_cast<real *>(&context.QInterpolatedPlusOnHost[upperStageOffset]),
                                  reinterpret_cast<real *>(&context.QInterpolatedPlusOnDevice[begin * QInterpolatedSize]),
                                  batchSize * QInterpolatedSize * sizeof(real),
                                  copyFromStream);

        device.api->copyFromAsync(reinterpret_cast<real *>(&context.QInterpolatedMinusOnHost[upperStageOffset]),
                                  reinterpret_cast<real *>(&context.QInterpolatedMinusOnDevice[begin * QInterpolatedSize]),
                                  batchSize * QInterpolatedSize * sizeof(real),
                                  copyFromStream);
      }

      void finalize() override {
        device.api->syncStreamFromCircularBuffer(copyFromStream);
      }
    private:
      void *copyFromStream{nullptr};
    } asyncCopyFrom(context, this);

    struct ComputeFriction : public DrBaseCallBack {
      explicit ComputeFriction(DrContext userContext, TimeCluster *cluster) : DrBaseCallBack(userContext, cluster) {
        assert(context.faceInformation != nullptr);
        assert(context.imposedStatePlusOnHost != nullptr);
        assert(context.imposedStateMinusOnHost != nullptr);
        assert(context.waveSpeedsPlus != nullptr);
        assert(context.waveSpeedsMinus != nullptr);
      }

      void operator()(size_t begin, size_t batchSize, size_t callCounter) override {
        const size_t upperStageOffset = (callCounter % DrPipeline::TailSize) * DrPipeline::DefaultBatchSize;
        const size_t lowerStageOffset = (callCounter % DrPipeline::NumStages) * DrPipeline::DefaultBatchSize;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (unsigned face = 0; face < batchSize; ++face) {
          e_interoperability.evaluateFrictionLaw(static_cast<int>(context.faceInformation[begin + face].meshFace),
                                                 context.QInterpolatedPlusOnHost[upperStageOffset + face],
                                                 context.QInterpolatedMinusOnHost[upperStageOffset + face],
                                                 context.imposedStatePlusOnHost[lowerStageOffset + face],
                                                 context.imposedStateMinusOnHost[lowerStageOffset + face],
                                                 cluster->m_fullUpdateTime,
                                                 cluster->m_dynamicRuptureKernel.timePoints,
                                                 cluster->m_dynamicRuptureKernel.timeWeights,
                                                 context.waveSpeedsPlus[begin + face],
                                                 context.waveSpeedsMinus[begin + face]);
        }
      }
      void finalize() override {}
    } computeFriction(context, this);


    struct AsyncCopyBack : public DrBaseCallBack {
      explicit AsyncCopyBack(DrContext userContext, TimeCluster *cluster) : DrBaseCallBack(userContext, cluster) {
        assert(context.devImposedStatePlus != nullptr);
        assert(context.devImposedStateMinus != nullptr);
        copyBackStream = device.api->getNextCircularStream();
      }

      void operator()(size_t begin, size_t batchSize, size_t callCounter) override {
        const size_t lowerStageOffset = (callCounter % DrPipeline::NumStages) * DrPipeline::DefaultBatchSize;
        const size_t imposedStateSize = tensor::QInterpolated::size() * batchSize * sizeof(real);

        device.api->syncStreamFromCircularBuffer(copyBackStream);
        device.api->copyToAsync(reinterpret_cast<real *>(&context.devImposedStatePlus[begin]),
                                reinterpret_cast<real *>(&context.imposedStatePlusOnHost[lowerStageOffset]),
                                imposedStateSize,
                                copyBackStream);

        device.api->copyToAsync(reinterpret_cast<real *>(&context.devImposedStateMinus[begin]),
                                reinterpret_cast<real *>(&context.imposedStateMinusOnHost[lowerStageOffset]),
                                imposedStateSize,
                                copyBackStream);
      }

      void finalize() override {
        device.api->syncStreamFromCircularBuffer(copyBackStream);
      }
    private:
      void *copyBackStream{nullptr};
    } asyncCopyBack(context, this);

    drPipeline.registerCallBack(0, &asyncCopyFrom);
    drPipeline.registerCallBack(1, &computeFriction);
    drPipeline.registerCallBack(2, &asyncCopyBack);
    drPipeline.run(layerData.getNumberOfCells());

    device.api->resetCircularStreamCounter();
  }
  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells());
  device.api->popLastProfilingMark();
}
#endif


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

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegration( seissol::initializers::Layer&  i_layerData ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // local integration buffer
  real l_integrationBuffer[tensor::I::size()] __attribute__((aligned(ALIGNMENT)));

  // pointer for the call of the ADER-function
  real* l_bufferPointer;

  real** buffers = i_layerData.var(m_lts->buffers);
  real** derivatives = i_layerData.var(m_lts->derivatives);
  CellMaterialData* materialData = i_layerData.var(m_lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp;

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer, tmp) schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
    auto data = loader.entry(l_cell);
    // overwrite cell buffer
    // TODO: Integrate this step into the kernel

    bool l_buffersProvided = (data.cellInformation.ltsSetup >> 8)%2 == 1; // buffers are provided
    bool l_resetBuffers = l_buffersProvided && ( (data.cellInformation.ltsSetup >> 10) %2 == 0 || m_resetLtsBuffers ); // they should be reset

    if (l_resetBuffers) {
      // assert presence of the buffer
      assert(buffers[l_cell] != nullptr);

      l_bufferPointer = buffers[l_cell];
    } else {
      // work on local buffer
      l_bufferPointer = l_integrationBuffer;
    }

    m_timeKernel.computeAder(m_timeStepWidth,
                             data,
                             tmp,
                             l_bufferPointer,
                             derivatives[l_cell],
                             m_fullUpdateTime,
                             true);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping (*boundaryMapping)[4] = i_layerData.var(m_lts->boundaryMapping);
    m_localKernel.computeIntegral(l_bufferPointer,
                                  data,
                                  tmp,
                                  &materialData[l_cell],
                                  &boundaryMapping[l_cell],
                                  m_fullUpdateTime,
                                  m_timeStepWidth
    );

    for (unsigned face = 0; face < 4; ++face) {
      auto& curFaceDisplacements = data.faceDisplacements[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr
          && data.cellInformation.faceTypes[face] != FaceType::freeSurfaceGravity) {
        kernel::addVelocity addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = m_globalDataOnHost->V3mTo2nFace;
        addVelocityKrnl.selectVelocity = init::selectVelocity::Values;
        addVelocityKrnl.faceDisplacement = data.faceDisplacements[face];
        addVelocityKrnl.I = l_bufferPointer;
        addVelocityKrnl.execute(face);
      }
    }

    // update lts buffers if required
    // TODO: Integrate this step into the kernel
    if (!l_resetBuffers && l_buffersProvided) {
      assert (buffers[l_cell] != nullptr);

      for (unsigned int l_dof = 0; l_dof < tensor::I::size(); ++l_dof) {
        buffers[l_cell][l_dof] += l_integrationBuffer[l_dof];
      }
    }
  }

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells());
}
#else // ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegration( seissol::initializers::Layer&  i_layerData ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  ConditionalBatchTableT& table = i_layerData.getCondBatchTable();
  kernels::LocalTmp tmp;

  m_timeKernel.computeBatchedAder(m_timeStepWidth, tmp, table);
  m_localKernel.computeBatchedIntegral(table, tmp);
  auto defaultStream = device.api->getDefaultStream();

  for (unsigned face = 0; face < 4; ++face) {
    ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
    if (table.find(key) != table.end()) {
      BatchTable &entry = table[key];
      // NOTE: integrated velocities have been computed implicitly, i.e
      // it is 6th, 7the and 8th columns of integrated dofs

      kernel::gpu_addVelocity displacementKrnl;
      displacementKrnl.faceDisplacement = entry.content[*EntityId::FaceDisplacement]->getPointers();
      displacementKrnl.integratedVelocities = const_cast<real const**>(entry.content[*EntityId::Ivelocities]->getPointers());
      displacementKrnl.V3mTo2nFace = m_globalDataOnDevice->V3mTo2nFace;

      // Note: this kernel doesn't require tmp. memory
      displacementKrnl.numElements = entry.content[*EntityId::FaceDisplacement]->getSize();
      displacementKrnl.streamPtr = device.api->getDefaultStream();
      displacementKrnl.execute(face);
    }
  }

  ConditionalKey key = ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
  if (table.find(key) != table.end()) {
    BatchTable &entry = table[key];

    if (m_resetLtsBuffers) {
      device.algorithms.streamBatchedData((entry.content[*EntityId::Idofs])->getPointers(),
                                          (entry.content[*EntityId::Buffers])->getPointers(),
                                          tensor::I::Size,
                                          (entry.content[*EntityId::Idofs])->getSize(),
                                          defaultStream);
    }
    else {
      device.algorithms.accumulateBatchedData((entry.content[*EntityId::Idofs])->getPointers(),
                                              (entry.content[*EntityId::Buffers])->getPointers(),
                                              tensor::I::Size,
                                              (entry.content[*EntityId::Idofs])->getSize(),
                                              defaultStream);
    }
  }

  device.api->synchDevice();
  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells());
  device.api->popLastProfilingMark();
}
#endif // ACL_DEVICE

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData ) {
  if (usePlasticity) {
    computeNeighboringIntegrationImplementation<true>(i_layerData);
  } else {
    computeNeighboringIntegrationImplementation<false>(i_layerData);
  }
}
#else // ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData ) {
  device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  ConditionalBatchTableT &table = i_layerData.getCondBatchTable();

  seissol::kernels::TimeCommon::computeBatchedIntegrals(m_timeKernel,
                                                        m_subTimeStart,
                                                        m_timeStepWidth,
                                                        table);
  m_neighborKernel.computeBatchedNeighborsIntegral(table);

  if (usePlasticity) {
    PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
    unsigned numAdjustedDofs = seissol::kernels::Plasticity::computePlasticityBatched(m_oneMinusIntegratingFactor,
                                                                                      m_timeStepWidth,
                                                                                      m_tv,
                                                                                      m_globalDataOnDevice,
                                                                                      table,
                                                                                      plasticity);

    g_SeisSolNonZeroFlopsPlasticity += i_layerData.getNumberOfCells() * m_flops_nonZero[PlasticityCheck] + numAdjustedDofs * m_flops_nonZero[PlasticityYield];
    g_SeisSolHardwareFlopsPlasticity += i_layerData.getNumberOfCells() * m_flops_hardware[PlasticityCheck] + numAdjustedDofs * m_flops_hardware[PlasticityYield];
  }

  device.api->synchDevice();
  device.api->popLastProfilingMark();
  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells());
}
#endif // ACL_DEVICE

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

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops(
    unsigned numberOfCells,
    CellLocalInformation const* cellInformation,
    long long& nonZeroFlops,
    long long& hardwareFlops  )
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
    // Contribution from displacement/integrated displacement
    for (unsigned face = 0; face < 4; ++face) {
      if (cellInformation->faceTypes[face] == FaceType::freeSurfaceGravity) {
        const auto [nonZeroFlopsDisplacement, hardwareFlopsDisplacement] =
        GravitationalFreeSurfaceBc::getFlopsDisplacementFace(face,
                                                             cellInformation[cell].faceTypes[face]);
        nonZeroFlops += nonZeroFlopsDisplacement;
        hardwareFlops += hardwareFlopsDisplacement;
      }
    }
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
}

long seissol::time_stepping::TimeCluster::getNumberOfCells() const {
  return m_clusterData->child<Copy>().getNumberOfCells() +
         m_clusterData->child<Interior>().getNumberOfCells();
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

