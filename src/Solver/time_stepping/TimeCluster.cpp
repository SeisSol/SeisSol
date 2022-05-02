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

//! fortran interoperability
extern seissol::Interoperability e_interoperability;

seissol::time_stepping::TimeCluster::TimeCluster(unsigned int i_clusterId, unsigned int i_globalClusterId,
                                                 bool usePlasticity,
                                                 LayerType layerType, double maxTimeStepSize,
                                                 long timeStepRate, bool printProgress,
                                                 DynamicRuptureScheduler *dynamicRuptureScheduler,
                                                 CompoundGlobalData i_globalData,
                                                 seissol::initializers::Layer *i_clusterData,
                                                 seissol::initializers::Layer *dynRupInteriorData,
                                                 seissol::initializers::Layer *dynRupCopyData,
                                                 seissol::initializers::LTS *i_lts,
                                                 seissol::initializers::DynamicRupture *i_dynRup,
                                                 LoopStatistics *i_loopStatistics,
                                                 ActorStateStatistics* actorStateStatistics) :
    AbstractTimeCluster(maxTimeStepSize, timeStepRate),
    // cluster ids
    usePlasticity(usePlasticity),
    m_globalDataOnHost( i_globalData.onHost ),
    m_globalDataOnDevice(i_globalData.onDevice ),
    m_clusterData(i_clusterData),
    // global data
    dynRupInteriorData(dynRupInteriorData),
    dynRupCopyData(dynRupCopyData),
    m_lts(i_lts),
    m_dynRup(i_dynRup),
    m_cellToPointSources(nullptr),
    m_numberOfCellToPointSourcesMappings(0),
    m_pointSources(nullptr),
    // cells
    m_loopStatistics(i_loopStatistics),
    actorStateStatistics(actorStateStatistics),
    m_receiverCluster(nullptr),
    layerType(layerType),
    printProgress(printProgress),
    m_clusterId(i_clusterId),
    m_globalClusterId(i_globalClusterId),
    dynamicRuptureScheduler(dynamicRuptureScheduler)
{
    // assert all pointers are valid
    assert( m_clusterData                              != nullptr );
    assert( m_globalDataOnHost                         != nullptr );
    if constexpr (seissol::isDeviceOn()) {
        assert( m_globalDataOnDevice                   != nullptr );
    }

  // set timings to zero
  m_receiverTime                  = 0;

  m_timeKernel.setGlobalData(i_globalData);
  m_localKernel.setGlobalData(i_globalData);
  m_localKernel.setInitConds(&e_interoperability.getInitialConditions());
  m_neighborKernel.setGlobalData(i_globalData);
  m_dynamicRuptureKernel.setGlobalData(i_globalData);

  computeFlops();

  m_regionComputeLocalIntegration = m_loopStatistics->getRegion("computeLocalIntegration");
  m_regionComputeNeighboringIntegration = m_loopStatistics->getRegion("computeNeighboringIntegration");
  m_regionComputeDynamicRupture = m_loopStatistics->getRegion("computeDynamicRupture");
};

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
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (m_receiverCluster != nullptr) {
    m_receiverTime = m_receiverCluster->calcReceivers(m_receiverTime, ct.correctionTime, timeStepSize());
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
      unsigned endSource =
          m_cellToPointSources[mapping].pointSourcesOffset + m_cellToPointSources[mapping].numberOfPointSources;
      if (m_pointSources->mode == sourceterm::PointSources::NRF) {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceNRF(m_pointSources->mInvJInvPhisAtSources[source],
                                                      m_pointSources->tensor[source],
                                                      m_pointSources->A[source],
                                                      m_pointSources->stiffnessTensor[source],
                                                      m_pointSources->slipRates[source],
                                                      ct.correctionTime,
                                                      ct.correctionTime + timeStepSize(),
                                                      *m_cellToPointSources[mapping].dofs);
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceFSRM(m_pointSources->mInvJInvPhisAtSources[source],
                                                       m_pointSources->tensor[source],
                                                       m_pointSources->slipRates[source][0],
                                                       ct.correctionTime,
                                                       ct.correctionTime + timeStepSize(),
                                                       *m_cellToPointSources[mapping].dofs);
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
  DROutput     *                        drOutput                                                          = layerData.var(m_dynRup->drOutput);
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
                                                    m_globalDataOnHost,
                                                   &godunovData[face],
                                                   &drOutput[face],
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

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_globalClusterId);
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
                                                 cluster->ct.correctionTime,
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
  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_globalClusterId);

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
    m_dynamicRuptureKernel.flopsGodunovState(faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegration(seissol::initializers::Layer& i_layerData, bool resetBuffers ) {
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
                             derivatives[l_cell],
                             ct.correctionTime,
                             true);

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

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);
}
#else // ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegration(seissol::initializers::Layer& i_layerData, bool resetBuffers ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
  ConditionalBatchTableT& table = i_layerData.getCondBatchTable();
  kernels::LocalTmp tmp;

  m_timeKernel.computeBatchedAder(timeStepSize(), tmp, table);
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

    if (resetBuffers) {
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
  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);

  device.api->popLastProfilingMark();
}
#endif // ACL_DEVICE

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeNeighboringIntegration(seissol::initializers::Layer& i_layerData,
                                                                        double subTimeStart) {
  if (usePlasticity) {
    computeNeighboringIntegrationImplementation<true>(i_layerData, subTimeStart);
  } else {
    computeNeighboringIntegrationImplementation<false>(i_layerData, subTimeStart);
  }
}
#else // ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData,
                                                                         double subTimeStart) {
  device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  ConditionalBatchTableT &table = i_layerData.getCondBatchTable();

  seissol::kernels::TimeCommon::computeBatchedIntegrals(m_timeKernel,
                                                        subTimeStart,
                                                        timeStepSize(),
                                                        table);
  m_neighborKernel.computeBatchedNeighborsIntegral(table);

  if (usePlasticity) {
    PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
    unsigned numAdjustedDofs = seissol::kernels::Plasticity::computePlasticityBatched(m_oneMinusIntegratingFactor,
                                                                                      timeStepSize(),
                                                                                      m_tv,
                                                                                      m_globalDataOnDevice,
                                                                                      table,
                                                                                      plasticity);

    g_SeisSolNonZeroFlopsPlasticity +=
        i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)];
    g_SeisSolHardwareFlopsPlasticity +=
        i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)];
  }

  device.api->synchDevice();
  device.api->popLastProfilingMark();
  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);
}
#endif // ACL_DEVICE

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Local)];
  flopsNonZero = 0;
  flopsHardware = 0;

  auto* cellInformation = layerData.var(m_lts->cellInformation);
  for (unsigned cell = 0; cell < layerData.getNumberOfCells(); ++cell) {
    unsigned cellNonZero, cellHardware;
    m_timeKernel.flopsAder(cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    m_localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    // Contribution from displacement/integrated displacement
    for (unsigned face = 0; face < 4; ++face) {
      if (cellInformation->faceTypes[face] == FaceType::freeSurfaceGravity) {
        const auto [nonZeroFlopsDisplacement, hardwareFlopsDisplacement] =
        GravitationalFreeSurfaceBc::getFlopsDisplacementFace(face,
                                                             cellInformation[cell].faceTypes[face]);
        flopsNonZero += nonZeroFlopsDisplacement;
        flopsHardware += hardwareFlopsDisplacement;
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeNeighborIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)];
  auto& drFlopsHardware = m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)];
  flopsNonZero = 0;
  flopsHardware = 0;
  drFlopsNonZero = 0;
  drFlopsHardware = 0;

  auto* cellInformation = layerData.var(m_lts->cellInformation);
  auto* drMapping = layerData.var(m_lts->drMapping);
  for (unsigned cell = 0; cell < layerData.getNumberOfCells(); ++cell) {
    unsigned cellNonZero, cellHardware;
    long long cellDRNonZero, cellDRHardware;
    m_neighborKernel.flopsNeighborsIntegral(cellInformation[cell].faceTypes,
                                            cellInformation[cell].faceRelations,
                                            drMapping[cell],
                                            cellNonZero,
                                            cellHardware,
                                            cellDRNonZero,
                                            cellDRHardware );
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    drFlopsNonZero += cellDRNonZero;
    drFlopsHardware += cellDRHardware;

    /// \todo add lts time integration
    /// \todo add plasticity
  }
}

void seissol::time_stepping::TimeCluster::computeFlops() {
  computeLocalIntegrationFlops(*m_clusterData);
  computeNeighborIntegrationFlops(*m_clusterData);
  computeDynamicRuptureFlops(*dynRupInteriorData,
                             m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)],
                             m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
  computeDynamicRuptureFlops(*dynRupCopyData,
                             m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)],
                             m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
  seissol::kernels::Plasticity::flopsPlasticity(
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)],
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)],
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)],
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]
          );
}

namespace seissol::time_stepping {
ActResult TimeCluster::act() {
  actorStateStatistics->enter(state);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics->enter(state);
  return result;
}

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
      if (neighbor.ct.timeStepRate > ct.timeStepRate
          && ct.stepsSinceLastSync > neighbor.ct.stepsSinceLastSync) {
          resetBuffers = false;
        }
  }
  if (ct.stepsSinceLastSync == 0) {
    resetBuffers = true;
  }

  // These methods compute the receivers/sources for both interior and copy cluster
  // and are called in actors for both copy AND interior.
  writeReceivers();
  computeLocalIntegration(*m_clusterData, resetBuffers);
  computeSources();

  g_SeisSolNonZeroFlopsLocal += m_flops_nonZero[static_cast<int>(ComputePart::Local)];
  g_SeisSolHardwareFlopsLocal += m_flops_hardware[static_cast<int>(ComputePart::Local)];
}
void TimeCluster::correct() {
  assert(state == ActorState::Predicted);

  /* Sub start time of width respect to the next cluster; use 0 if not relevant, for example in GTS.
   * LTS requires to evaluate a partial time integration of the derivatives. The point zero in time refers to the derivation of the surrounding time derivatives, which
   * coincides with the last completed time step of the next cluster. The start/end of the time step is the start/end of this clusters time step relative to the zero point.
   *   Example:
   *                                              5 dt
   *   |-----------------------------------------------------------------------------------------| <<< Time stepping of the next cluster (Cn) (5x larger than the current).
   *   |                 |                 |                 |                 |                 |
   *   |*****************|*****************|+++++++++++++++++|                 |                 | <<< Status of the current cluster.
   *   |                 |                 |                 |                 |                 |
   *   |-----------------|-----------------|-----------------|-----------------|-----------------| <<< Time stepping of the current cluster (Cc).
   *   0                 dt               2dt               3dt               4dt               5dt
   *
   *   In the example above two clusters are illustrated: Cc and Cn. Cc is the current cluster under consideration and Cn the next cluster with respect to LTS terminology.
   *   Cn is currently at time 0 and provided Cc with derivatives valid until 5dt. Cc updated already twice and did its last full update to reach 2dt (== subTimeStart). Next
   *   computeNeighboringCopy is called to accomplish the next full update to reach 3dt (+++). Besides working on the buffers of own buffers and those of previous clusters,
   *   Cc needs to evaluate the time prediction of Cn in the interval [2dt, 3dt].
   */
  double subTimeStart = ct.correctionTime - lastSubTime;

  // Note, if this is a copy layer actor, we need the FL_Copy and the FL_Int.
  // Otherwise, this is an interior layer actor, and we need only the FL_Int.
  // We need to avoid computing it twice.
  if (dynamicRuptureScheduler->hasDynamicRuptureFaces()) {
    if (dynamicRuptureScheduler->mayComputeInterior(ct.stepsSinceStart)) {
      computeDynamicRupture(*dynRupInteriorData);
      g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)];
      g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)];
      dynamicRuptureScheduler->setLastCorrectionStepsInterior(ct.stepsSinceStart);
    }
    if (layerType == Copy) {
      computeDynamicRupture(*dynRupCopyData);
      g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)];
      g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)];
      dynamicRuptureScheduler->setLastCorrectionStepsCopy((ct.stepsSinceStart));
    }

  }
  computeNeighboringIntegration(*m_clusterData, subTimeStart);

  g_SeisSolNonZeroFlopsNeighbor += m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)];
  g_SeisSolHardwareFlopsNeighbor += m_flops_hardware[static_cast<int>(ComputePart::Neighbor)];
  g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)];
  g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)];

  // First cluster calls fault receiver output
  // Call fault output only if both interior and copy parts of DR were computed
  // TODO: Change from iteration based to time based
  if (m_clusterId == 0
      && dynamicRuptureScheduler->mayComputeFaultOutput(ct.stepsSinceStart)) {
    e_interoperability.faultOutput(ct.correctionTime + timeStepSize(), timeStepSize());
    dynamicRuptureScheduler->setLastFaultOutput(ct.stepsSinceStart);
  }



  // TODO(Lukas) Adjust with time step rate? Relevant is maximum cluster is not on this node
  const auto nextCorrectionSteps = ct.nextCorrectionSteps();
  if constexpr (USE_MPI) {
    if (printProgress && (((nextCorrectionSteps / timeStepRate) % 100) == 0)) {
      const int rank = MPI::mpi.rank();
      logInfo(rank) << "#max-updates since sync: " << nextCorrectionSteps
                    << " @ " << ct.nextCorrectionTime(syncTime);

      }
  }

}

void TimeCluster::reset() {
    AbstractTimeCluster::reset();
}

void TimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  const auto rank = MPI::mpi.rank();
  logWarning(rank)
  << "No update since " << timeSinceLastUpdate.count()
  << "[s] for global cluster " << m_globalClusterId
  << " with local cluster id " << m_clusterId
  << " at state " << actorStateToString(state)
  << " predTime = " << ct.predictionTime
  << " predictionsSinceSync = " << ct.predictionsSinceLastSync
  << " corrTime = " << ct.correctionTime
  << " correctionsSinceSync = " << ct.stepsSinceLastSync
  << " stepsTillSync = " << ct.stepsUntilSync
  << " mayPredict = " << mayPredict()
  << " mayCorrect = " << mayCorrect()
  << " maySync = " << maySync();
  for (auto& neighbor : neighbors) {
    logWarning(rank)
    << "Neighbor with rate = " << neighbor.ct.timeStepRate
    << "PredTime = " << neighbor.ct.predictionTime
    << "CorrTime = " << neighbor.ct.correctionTime
    << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
    << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
  }

}

unsigned int TimeCluster::getClusterId() const {
  return m_clusterId;
}

unsigned int TimeCluster::getGlobalClusterId() const {
  return m_globalClusterId;
}

LayerType TimeCluster::getLayerType() const {
  return layerType;
}
void TimeCluster::setReceiverTime(double receiverTime) {
  m_receiverTime = receiverTime;
}

}
