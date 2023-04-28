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
#include <Monitoring/instrumentation.fpp>

#include <cassert>
#include <cstring>

#include <generated_code/kernel.h>

//! fortran interoperability
extern seissol::Interoperability e_interoperability;

seissol::time_stepping::TimeCluster::TimeCluster(unsigned int i_clusterId, unsigned int i_globalClusterId,
                                                 unsigned int profilingId,
                                                 bool usePlasticity,
                                                 LayerType layerType, double maxTimeStepSize,
                                                 long timeStepRate, bool printProgress,
                                                 DynamicRuptureScheduler *dynamicRuptureScheduler,
                                                 CompoundGlobalData i_globalData,
                                                 seissol::initializers::Layer *i_clusterData,
                                                 seissol::initializers::Layer *dynRupInteriorData,
                                                 seissol::initializers::Layer *dynRupCopyData,
                                                 seissol::initializers::LTS *i_lts,
                                                 seissol::initializers::DynamicRupture* i_dynRup,
                                                 seissol::dr::friction_law::FrictionSolver* i_FrictionSolver,
                                                 dr::output::OutputManager* i_faultOutputManager,
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
    frictionSolver(i_FrictionSolver),
    faultOutputManager(i_faultOutputManager),
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
    m_profilingId(profilingId),
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

  frictionSolverLocalInterior = std::move(frictionSolver->clone());
  frictionSolverLocalCopy = std::move(frictionSolver->clone());
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
  if (layerData.getNumberOfCells() == 0) return;
  SCOREP_USER_REGION_DEFINE(myRegionHandle)
  SCOREP_USER_REGION_BEGIN(myRegionHandle, "computeDynamicRuptureSpaceTimeInterpolation", SCOREP_USER_REGION_TYPE_COMMON )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  DRFaceInformation* faceInformation = layerData.var(m_dynRup->faceInformation);
  DRGodunovData* godunovData = layerData.var(m_dynRup->godunovData);
  DREnergyOutput* drEnergyOutput = layerData.var(m_dynRup->drEnergyOutput);
  real** timeDerivativePlus = layerData.var(m_dynRup->timeDerivativePlus);
  real** timeDerivativeMinus = layerData.var(m_dynRup->timeDerivativeMinus);
  auto* qInterpolatedPlus = layerData.var(m_dynRup->qInterpolatedPlus);
  auto* qInterpolatedMinus = layerData.var(m_dynRup->qInterpolatedMinus);

  m_dynamicRuptureKernel.setTimeStepWidth(timeStepSize());
  frictionSolver->computeDeltaT(m_dynamicRuptureKernel.timePoints);

#pragma omp parallel 
  {
  LIKWID_MARKER_START("computeDynamicRuptureSpaceTimeInterpolation");
  }
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynamicRuptureKernel.spaceTimeInterpolation(faceInformation[face],
                                                  m_globalDataOnHost,
                                                  &godunovData[face],
                                                  &drEnergyOutput[face],
                                                  timeDerivativePlus[face],
                                                  timeDerivativeMinus[face],
                                                  qInterpolatedPlus[face],
                                                  qInterpolatedMinus[face],
                                                  timeDerivativePlus[prefetchFace],
                                                  timeDerivativeMinus[prefetchFace]);
  }
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel 
  {
  LIKWID_MARKER_STOP("computeDynamicRuptureSpaceTimeInterpolation");
  LIKWID_MARKER_START("computeDynamicRuptureFrictionLaw");
  }

  SCOREP_USER_REGION_BEGIN(myRegionHandle, "computeDynamicRuptureFrictionLaw", SCOREP_USER_REGION_TYPE_COMMON )
  frictionSolver->evaluate(layerData,
                           m_dynRup,
                           ct.correctionTime,
                           m_dynamicRuptureKernel.timeWeights);
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel 
  {
  LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_profilingId);
}
#else

void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData ) {
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  if (layerData.getNumberOfCells() > 0) {
    // compute space time interpolation part

    const double stepSizeWidth = timeStepSize();
    ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    auto computeGraphKey = initializers::GraphKey(graphType, stepSizeWidth);
    auto computeGraphHandle = layerData.getDeviceComputeGraphHandle(computeGraphKey);

    auto& table = layerData.getConditionalTable<inner_keys::Dr>();
    m_dynamicRuptureKernel.setTimeStepWidth(stepSizeWidth);

    device.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    if (!computeGraphHandle) {
      device.api->streamBeginCapture();

      m_dynamicRuptureKernel.batchedSpaceTimeInterpolation(table);
      assert(device.api->isCircularStreamsJoinedWithDefault() &&
             "circular streams must be joined with the default stream");

      device.api->streamEndCapture();

      computeGraphHandle = device.api->getLastGraphHandle();
      layerData.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
      device.api->syncDefaultStreamWithHost();
    }

    if (computeGraphHandle.isInitialized()) {
      device.api->launchGraph(computeGraphHandle);
      device.api->syncGraph(computeGraphHandle);
    }
    device.api->popLastProfilingMark();

    device.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
    frictionSolver->computeDeltaT(m_dynamicRuptureKernel.timePoints);
    frictionSolver->evaluate(layerData,
                             m_dynRup,
                             ct.correctionTime,
                             m_dynamicRuptureKernel.timeWeights);
    device.api->popLastProfilingMark();
  }
  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_profilingId);
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

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_profilingId);
}
#else // ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegration(
  seissol::initializers::Layer& i_layerData,
  bool resetBuffers) {

  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
  auto& dataTable = i_layerData.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = i_layerData.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = i_layerData.getConditionalTable<inner_keys::Indices>();

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp;

  const double timeStepWidth = timeStepSize();

  ComputeGraphType graphType{ComputeGraphType::LocalIntegral};
  auto computeGraphKey = initializers::GraphKey(graphType, timeStepWidth, true);
  auto computeGraphHandle = i_layerData.getDeviceComputeGraphHandle(computeGraphKey);

  if (!computeGraphHandle) {
    device.api->streamBeginCapture();

    m_timeKernel.computeBatchedAder(timeStepWidth,
                                    tmp,
                                    dataTable,
                                    materialTable,
                                    true);
    assert(device.api->isCircularStreamsJoinedWithDefault() &&
           "circular streams must be joined with the default stream");

    m_localKernel.computeBatchedIntegral(dataTable,
                                         materialTable,
                                         indicesTable,
                                         loader,
                                         tmp,
                                         timeStepWidth);
    assert(device.api->isCircularStreamsJoinedWithDefault() &&
           "circular streams must be joined with the default stream");

    device.api->streamEndCapture();

    computeGraphHandle = device.api->getLastGraphHandle();
    i_layerData.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    device.api->syncDefaultStreamWithHost();
  }

  if (computeGraphHandle.isInitialized()) {
    device.api->launchGraph(computeGraphHandle);
    device.api->syncGraph(computeGraphHandle);
  }

  m_localKernel.evaluateBatchedTimeDependentBc(dataTable,
                                               indicesTable,
                                               loader,
                                               ct.correctionTime,
                                               timeStepWidth);

  graphType = resetBuffers ? ComputeGraphType::AccumulatedVelocities : ComputeGraphType::StreamedVelocities;
  computeGraphKey = initializers::GraphKey(graphType);
  computeGraphHandle = i_layerData.getDeviceComputeGraphHandle(computeGraphKey);

  if (!computeGraphHandle) {
    device.api->streamBeginCapture();

    auto defaultStream = device.api->getDefaultStream();

    for (unsigned face = 0; face < 4; ++face) {
      ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
      if (dataTable.find(key) != dataTable.end()) {
        auto& entry = dataTable[key];
        // NOTE: integrated velocities have been computed implicitly, i.e
        // it is 6th, 7the and 8th columns of integrated dofs

        kernel::gpu_addVelocity displacementKrnl;
        displacementKrnl.faceDisplacement = entry.get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
        displacementKrnl.integratedVelocities = const_cast<real const**>(entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
        displacementKrnl.V3mTo2nFace = m_globalDataOnDevice->V3mTo2nFace;

        // Note: this kernel doesn't require tmp. memory
        displacementKrnl.numElements = entry.get(inner_keys::Wp::Id::FaceDisplacement)->getSize();
        displacementKrnl.streamPtr = defaultStream;
        displacementKrnl.execute(face);
      }
    }

    ConditionalKey key = ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
    if (dataTable.find(key) != dataTable.end()) {
      auto& entry = dataTable[key];

      if (resetBuffers) {
        device.algorithms.streamBatchedData(
            (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
            (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
            tensor::I::Size,
            (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
            defaultStream);
      } else {
        device.algorithms.accumulateBatchedData(
            (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
            (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
            tensor::I::Size,
            (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
            defaultStream);
      }
    }

    device.api->streamEndCapture();

    computeGraphHandle = device.api->getLastGraphHandle();
    i_layerData.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    device.api->syncDefaultStreamWithHost();
  }

  if (computeGraphHandle.isInitialized()) {
    device.api->launchGraph(computeGraphHandle);
    device.api->syncGraph(computeGraphHandle);
  }

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_profilingId);
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

  const double timeStepWidth = timeStepSize();
  auto& table = i_layerData.getConditionalTable<inner_keys::Wp>();

  seissol::kernels::TimeCommon::computeBatchedIntegrals(m_timeKernel,
                                                        subTimeStart,
                                                        timeStepWidth,
                                                        table);
  assert(device.api->isCircularStreamsJoinedWithDefault() &&
         "circular streams must be joined with the default stream");

  ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializers::GraphKey(graphType);
  auto computeGraphHandle = i_layerData.getDeviceComputeGraphHandle(computeGraphKey);

  if (!computeGraphHandle) {
    device.api->streamBeginCapture();

    m_neighborKernel.computeBatchedNeighborsIntegral(table);
    assert(device.api->isCircularStreamsJoinedWithDefault() &&
           "circular streams must be joined with the default stream");

    device.api->streamEndCapture();

    computeGraphHandle = device.api->getLastGraphHandle();
    i_layerData.updateDeviceComputeGraphHandle(computeGraphKey, computeGraphHandle);
    device.api->syncDefaultStreamWithHost();
  }

  if (computeGraphHandle.isInitialized()) {
    // Note: graph stream needs to wait the default stream
    // (used in `computeBatchedIntegrals`)
    device.api->syncDefaultStreamWithHost();
    device.api->launchGraph(computeGraphHandle);
    device.api->syncGraph(computeGraphHandle);
  }

  if (usePlasticity) {
    PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
    unsigned numAdjustedDofs = seissol::kernels::Plasticity::computePlasticityBatched(m_oneMinusIntegratingFactor,
                                                                                      timeStepWidth,
                                                                                      m_tv,
                                                                                      m_globalDataOnDevice,
                                                                                      table,
                                                                                      plasticity);

    seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsPlasticity(
        i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)]);
    seissol::SeisSol::main.flopCounter().incrementHardwareFlopsPlasticity(
        i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]);
  }

  device.api->syncDefaultStreamWithHost();
  device.api->popLastProfilingMark();
  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_profilingId);
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
  if (m_clusterData->getNumberOfCells() == 0) return;

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

  writeReceivers();
  computeLocalIntegration(*m_clusterData, resetBuffers);
  computeSources();

  seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsLocal(m_flops_nonZero[static_cast<int>(ComputePart::Local)]);
  seissol::SeisSol::main.flopCounter().incrementHardwareFlopsLocal(m_flops_hardware[static_cast<int>(ComputePart::Local)]);
}
void TimeCluster::correct() {
  assert(state == ActorState::Predicted);

  /* Sub start time of width respect to the next cluster; use 0 if not relevant, for example in GTS.
   * LTS requires to evaluate a partial time integration of the derivatives. The point zero in time
   * refers to the derivation of the surrounding time derivatives, which coincides with the last
   * completed time step of the next cluster. The start/end of the time step is the start/end of
   * this clusters time step relative to the zero point.
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
      seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      dynamicRuptureScheduler->setLastCorrectionStepsInterior(ct.stepsSinceStart);
    }
    if (layerType == Copy) {
      computeDynamicRupture(*dynRupCopyData);
      seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      dynamicRuptureScheduler->setLastCorrectionStepsCopy((ct.stepsSinceStart));
    }

  }
  computeNeighboringIntegration(*m_clusterData, subTimeStart);

  seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsNeighbor(m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)]);
  seissol::SeisSol::main.flopCounter().incrementHardwareFlopsNeighbor(m_flops_hardware[static_cast<int>(ComputePart::Neighbor)]);
  seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)]);
  seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)]);

  // First cluster calls fault receiver output
  // Call fault output only if both interior and copy parts of DR were computed
  // TODO: Change from iteration based to time based
  if (dynamicRuptureScheduler->isFirstClusterWithDynamicRuptureFaces()
      && dynamicRuptureScheduler->mayComputeFaultOutput(ct.stepsSinceStart)) {
    faultOutputManager->writePickpointOutput(ct.correctionTime + timeStepSize(), timeStepSize());
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

void TimeCluster::collectiveComputeLocalIntegration(const std::vector<TimeCluster*>& clusters) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  for (TimeCluster* cluster : clusters) {
    cluster->m_loopStatistics->begin(cluster->m_regionComputeLocalIntegration);
  }

  std::vector<bool> resetBuffers(clusters.size());
  std::vector<kernels::LocalData::Loader> loaders(clusters.size());

  std::vector<size_t> clusterOffset(clusters.size() + 1);

  for (size_t i = 0; i < clusters.size(); ++i) {
    resetBuffers[i] = true;
    for (const auto& neighbor : clusters[i]->neighbors) {
        if (neighbor.ct.timeStepRate > clusters[i]->ct.timeStepRate
            && clusters[i]->ct.stepsSinceLastSync > neighbor.ct.stepsSinceLastSync) {
            resetBuffers[i] = false;
          }
    }
    if (clusters[i]->ct.stepsSinceLastSync == 0) {
      resetBuffers[i] = true;
    }

    loaders[i].load(*clusters[i]->m_lts, *clusters[i]->m_clusterData);

    clusterOffset[i+1] = clusterOffset[i] + clusters[i]->m_clusterData->getNumberOfCells();
  }
  

  // local integration buffer
  real l_integrationBuffer[tensor::I::size()] alignas(ALIGNMENT);
  // pointer for the call of the ADER-function
  real* l_bufferPointer;
  kernels::LocalTmp tmp{};

  auto totalCellCount = *clusterOffset.rbegin();

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer, tmp) schedule(static)
#endif
  //for (unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++) {
  for (size_t i = 0; i < totalCellCount; ++i) {
    // https://en.cppreference.com/w/cpp/algorithm/upper_bound
    auto idElement = std::upper_bound(clusterOffset.begin(), clusterOffset.end(), i);
    size_t clusterId = std::distance(clusterOffset.begin(), idElement) - 1;

    assert(i > clusterOffset[clusterId]);
    unsigned l_cell = i - clusterOffset[clusterId];

    auto data = loaders[clusterId].entry(l_cell);
    real** buffers = clusters[clusterId]->m_clusterData->var(clusters[clusterId]->m_lts->buffers);
    real** derivatives = clusters[clusterId]->m_clusterData->var(clusters[clusterId]->m_lts->derivatives);
    CellMaterialData* materialData = clusters[clusterId]->m_clusterData->var(clusters[clusterId]->m_lts->material);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided = (data.cellInformation.ltsSetup >> 8) % 2 == 1; // buffers are provided
    const bool resetMyBuffers = buffersProvided && ( (data.cellInformation.ltsSetup >> 10) %2 == 0 || resetBuffers[clusterId] ); // they should be reset

    if (resetMyBuffers) {
      // assert presence of the buffer
      assert(buffers[l_cell] != nullptr);

      l_bufferPointer = buffers[l_cell];
    } else {
      // work on local buffer
      l_bufferPointer = l_integrationBuffer;
    }

    clusters[clusterId]->m_timeKernel.computeAder(clusters[clusterId]->timeStepSize(),
                             data,
                             tmp,
                             l_bufferPointer,
                             derivatives[l_cell],
                             true);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping (*boundaryMapping)[4] = clusters[clusterId]->m_clusterData->var(clusters[clusterId]->m_lts->boundaryMapping);
    clusters[clusterId]->m_localKernel.computeIntegral(l_bufferPointer,
                                  data,
                                  tmp,
                                  &materialData[l_cell],
                                  &boundaryMapping[l_cell],
                                  clusters[clusterId]->ct.correctionTime,
                                  clusters[clusterId]->timeStepSize()
    );

    for (unsigned face = 0; face < 4; ++face) {
      auto& curFaceDisplacements = data.faceDisplacements[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr
          && data.cellInformation.faceTypes[face] != FaceType::freeSurfaceGravity) {
        kernel::addVelocity addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = clusters[clusterId]->m_globalDataOnHost->V3mTo2nFace;
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

  for (TimeCluster* cluster : clusters) {
    cluster->m_loopStatistics->end(cluster->m_regionComputeLocalIntegration, cluster->m_clusterData->getNumberOfCells(), cluster->m_profilingId);
  }
}

void TimeCluster::collectiveComputeSources(const std::vector<TimeCluster*>& clusters) {
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
  SCOREP_USER_REGION( "computeSources", SCOREP_USER_REGION_TYPE_FUNCTION )

  std::vector<size_t> clusterOffset(clusters.size() + 1);
  for (size_t i = 0; i < clusters.size(); ++i) {
    clusterOffset[i+1] = clusterOffset[i] + clusters[i]->m_numberOfCellToPointSourcesMappings;
  }

  auto sourceCount = *clusterOffset.rbegin();

  // Return when point sources not initialised. This might happen if there
  // are no point sources on this rank.
  if (sourceCount > 0) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned i = 0; i < sourceCount; ++i) {
      auto idElement = std::upper_bound(clusterOffset.begin(), clusterOffset.end(), i);
      size_t clusterId = std::distance(clusterOffset.begin(), idElement) - 1;

      assert(i > clusterOffset[clusterId]);
      unsigned mapping = i - clusterOffset[clusterId];

      unsigned startSource = clusters[clusterId]->m_cellToPointSources[mapping].pointSourcesOffset;
      unsigned endSource = startSource + clusters[clusterId]->m_cellToPointSources[mapping].numberOfPointSources;
      if (clusters[clusterId]->m_pointSources->mode == sourceterm::PointSources::NRF) {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceNRF(clusters[clusterId]->m_pointSources->mInvJInvPhisAtSources[source],
                                                      clusters[clusterId]->m_pointSources->tensor[source],
                                                      clusters[clusterId]->m_pointSources->A[source],
                                                      clusters[clusterId]->m_pointSources->stiffnessTensor[source],
                                                      clusters[clusterId]->m_pointSources->slipRates[source],
                                                      clusters[clusterId]->ct.correctionTime,
                                                      clusters[clusterId]->ct.correctionTime + clusters[clusterId]->timeStepSize(),
                                                      *clusters[clusterId]->m_cellToPointSources[mapping].dofs);
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceFSRM(clusters[clusterId]->m_pointSources->mInvJInvPhisAtSources[source],
                                                       clusters[clusterId]->m_pointSources->tensor[source],
                                                       clusters[clusterId]->m_pointSources->slipRates[source][0],
                                                       clusters[clusterId]->ct.correctionTime,
                                                       clusters[clusterId]->ct.correctionTime + clusters[clusterId]->timeStepSize(),
                                                       *clusters[clusterId]->m_cellToPointSources[mapping].dofs);
        }
      }
    }
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
}

void TimeCluster::collectivePredict(const std::vector<TimeCluster*>& clusters) {
  for (TimeCluster* cluster : clusters) {
    assert(cluster->state == ActorState::Corrected);
  }

  /*logInfo(seissol::MPI::mpi.rank()) << "Collective prediction. " << clusters.size() << "clusters.";
  for (TimeCluster* cluster : clusters) {
    logInfo(seissol::MPI::mpi.rank()) << "Cluster ID:" << cluster->m_clusterId << ". Layer type: " << static_cast<int>(cluster->layerType);
  }*/

  if (clusters.empty()) {
    return;
  }

  // writing receivers is still serial
  for (TimeCluster* cluster : clusters) {
    cluster->writeReceivers();
  }
  collectiveComputeLocalIntegration(clusters);
  collectiveComputeSources(clusters);

  // TODO(David): need parallelization here? Most likely not. But just in case... That's why this TODO's here.
  for (TimeCluster* cluster : clusters) {
    seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsLocal(cluster->m_flops_nonZero[static_cast<int>(ComputePart::Local)]);
    seissol::SeisSol::main.flopCounter().incrementHardwareFlopsLocal(cluster->m_flops_hardware[static_cast<int>(ComputePart::Local)]);
  }

  // TODO(David): is it really worth it to put a parallel for here? Maybe we should also directly update all neighbors which are in the time cluster list...
  #pragma omp parallel for
  for (TimeCluster* cluster : clusters) {
    cluster->postPredict();
  }
}

std::vector<TimeCluster*> TimeCluster::collectMayPredict(const std::vector<TimeCluster*>& highPriority, const std::vector<TimeCluster*>& lowPriority, size_t threshold) {
    std::vector<TimeCluster*> clusters;

    // take all high-priority clusters which may predict
    size_t cellCount = 0;
    for (TimeCluster* cluster : highPriority) {
      auto numCells = cluster->m_clusterData->getNumberOfCells();
      if (cluster->getNextLegalAction() == ActorAction::Predict) {
        clusters.push_back(cluster);
        cellCount += numCells;
      }
    }

    // for the low-priority clusters, take all which still fit in (do this sort greedily, for now at least)
    // or, if we have not found a high-priority cluster, take the first-best low-priority cluster at least
    if (cellCount < threshold) {
      size_t cellCount = 0;
      for (TimeCluster* cluster : lowPriority) {
        auto numCells = cluster->m_clusterData->getNumberOfCells();
        if (cluster->getNextLegalAction() == ActorAction::Predict && (cellCount + numCells < threshold || cellCount == 0)) {
          clusters.push_back(cluster);
          cellCount += numCells;
        }
      }
    }

    return clusters;
  }

void TimeCluster::collectiveCorrect(const std::vector<TimeCluster*>& clusters) {
  for (TimeCluster* cluster : clusters) {
    assert(cluster->state == ActorState::Predicted);
  }

  /*logInfo(seissol::MPI::mpi.rank()) << "Collective correction. " << clusters.size() << "clusters.";
  for (TimeCluster* cluster : clusters) {
    logInfo(seissol::MPI::mpi.rank()) << "Cluster ID:" << cluster->m_clusterId << ". Layer type: " << static_cast<int>(cluster->layerType) << "" << cluster->ct.predictionTime << "" << cluster->ct.correctionTime;
  }*/

  // Note, if this is a copy layer actor, we need the FL_Copy and the FL_Int.
  // Otherwise, this is an interior layer actor, and we need only the FL_Int.
  // We need to avoid computing it twice.

  if (clusters.empty()) {
    return;
  }

  /*for (TimeCluster* cluster : clusters) {
    if (cluster->dynamicRuptureScheduler->hasDynamicRuptureFaces()) {
      if (cluster->dynamicRuptureScheduler->mayComputeInterior(cluster->ct.stepsSinceStart)) {
        cluster->computeDynamicRupture(*cluster->dynRupInteriorData);
        seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(cluster->m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
        seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(cluster->m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
        cluster->dynamicRuptureScheduler->setLastCorrectionStepsInterior(cluster->ct.stepsSinceStart);
      }
      if (cluster->layerType == Copy) {
        cluster->computeDynamicRupture(*cluster->dynRupCopyData);
        seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(cluster->m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
        seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(cluster->m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
        cluster->dynamicRuptureScheduler->setLastCorrectionStepsCopy((cluster->ct.stepsSinceStart));
      }

    }
  }*/

  collectiveComputeDynamicRupture(clusters);

  // TODO(David): adjust plasticity check
  if (clusters.size() > 0 && clusters[0]->usePlasticity) {
    collectiveComputeNeighboringIntegrationImplementation<true>(clusters);
  }
  else {
    collectiveComputeNeighboringIntegrationImplementation<false>(clusters);
  }

  for (TimeCluster* cluster : clusters) {
    seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(cluster->m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
    seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(cluster->m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);

    seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(cluster->m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
    seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(cluster->m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);

    seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsNeighbor(cluster->m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)]);
    seissol::SeisSol::main.flopCounter().incrementHardwareFlopsNeighbor(cluster->m_flops_hardware[static_cast<int>(ComputePart::Neighbor)]);
    seissol::SeisSol::main.flopCounter().incrementNonZeroFlopsDynamicRupture(cluster->m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)]);
    seissol::SeisSol::main.flopCounter().incrementHardwareFlopsDynamicRupture(cluster->m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)]);

    // First cluster calls fault receiver output
    // Call fault output only if both interior and copy parts of DR were computed
    // TODO: Change from iteration based to time based

    if (cluster->dynamicRuptureScheduler->isFirstClusterWithDynamicRuptureFaces()
        && cluster->dynamicRuptureScheduler->mayComputeFaultOutput(cluster->ct.stepsSinceStart)) {
      cluster->faultOutputManager->writePickpointOutput(cluster->ct.correctionTime + cluster->timeStepSize(), cluster->timeStepSize());
      cluster->dynamicRuptureScheduler->setLastFaultOutput(cluster->ct.stepsSinceStart);
    }

    // TODO(Lukas) Adjust with time step rate? Relevant is maximum cluster is not on this node
    const auto nextCorrectionSteps = cluster->ct.nextCorrectionSteps();
    if constexpr (USE_MPI) {
      if (cluster->printProgress && (((nextCorrectionSteps / cluster->timeStepRate) % 100) == 0)) {
        const int rank = MPI::mpi.rank();
        logInfo(rank) << "#max-updates since sync: " << nextCorrectionSteps
                      << " @ " << cluster->ct.nextCorrectionTime(cluster->syncTime);

        }
    }
  }

  // TODO(David): is it really worth it to put a parallel for here? Maybe we should also directly update all neighbors which are in the time cluster list...
  #pragma omp parallel for
  for (TimeCluster* cluster : clusters) {
    cluster->postCorrect();
  }
}

void TimeCluster::collectiveComputeDynamicRupture(const std::vector<TimeCluster*>& clusters) {
  for (TimeCluster* cluster : clusters) {
    cluster->m_loopStatistics->begin(cluster->m_regionComputeDynamicRupture);
  }

  std::set<DynamicRuptureScheduler*> schedulers;
  std::vector<size_t> clusterOffset(2*clusters.size() + 1);
  for (size_t i = 0; i < clusters.size(); ++i) {
    size_t numCells
      = clusters[i]->dynamicRuptureScheduler->hasDynamicRuptureFaces()
      && clusters[i]->dynamicRuptureScheduler->mayComputeInterior(clusters[i]->ct.stepsSinceStart)
      && clusters[i]->layerType == Interior // that's new.
      ? clusters[i]->dynRupInteriorData->getNumberOfCells() : 0;
    clusterOffset[i+1] = clusterOffset[i] + numCells;
  }
  for (size_t i = 0; i < clusters.size(); ++i) {
    size_t numCells
      = clusters[i]->dynamicRuptureScheduler->hasDynamicRuptureFaces()
      && clusters[i]->layerType == Copy
      ? clusters[i]->dynRupCopyData->getNumberOfCells() : 0;
    clusterOffset[clusters.size() + i+1] = clusterOffset[clusters.size() + i] + numCells;
  }

  for (TimeCluster* cluster : clusters) {
    cluster->m_dynamicRuptureKernel.setTimeStepWidth(cluster->timeStepSize());

    cluster->frictionSolverLocalInterior->computeDeltaT(cluster->m_dynamicRuptureKernel.timePoints);
    cluster->frictionSolverLocalCopy->computeDeltaT(cluster->m_dynamicRuptureKernel.timePoints);

    cluster->frictionSolverLocalInterior->copyData(*cluster->dynRupInteriorData, cluster->m_dynRup, cluster->ct.correctionTime);
    cluster->frictionSolverLocalCopy->copyData(*cluster->dynRupCopyData, cluster->m_dynRup, cluster->ct.correctionTime);
  }

  auto totalCount = *clusterOffset.rbegin();  

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < totalCount; ++i) {

    auto idElement = std::upper_bound(clusterOffset.begin(), clusterOffset.end(), i);
    size_t clusterOffsetId = std::distance(clusterOffset.begin(), idElement) - 1;

    assert(i > clusterOffset[clusterOffsetId]);
    unsigned face = i - clusterOffset[clusterOffsetId];

    size_t clusterId = clusterOffsetId % clusters.size();
    auto& layerData = clusterId == clusterOffsetId ? *clusters[clusterId]->dynRupInteriorData : *clusters[clusterId]->dynRupCopyData;
    auto* cluster = clusters[clusterId];

    auto& layerDataPrefetch = clusterId == clusterOffsetId ? *clusters[clusterId]->dynRupInteriorData : *clusters[clusterId]->dynRupCopyData;
    auto* clusterPrefetch = clusterId < clusterOffset.size() - 2 ? clusters[(clusterId+1) % clusters.size()] : clusters[clusterId];

    DRFaceInformation* faceInformation = layerData.var(cluster->m_dynRup->faceInformation);
    DRGodunovData* godunovData = layerData.var(cluster->m_dynRup->godunovData);
    DREnergyOutput* drEnergyOutput = layerData.var(cluster->m_dynRup->drEnergyOutput);
    real** timeDerivativePlus = layerData.var(cluster->m_dynRup->timeDerivativePlus);
    real** timeDerivativeMinus = layerData.var(cluster->m_dynRup->timeDerivativeMinus);
    auto* qInterpolatedPlus = layerData.var(cluster->m_dynRup->qInterpolatedPlus);
    auto* qInterpolatedMinus = layerData.var(cluster->m_dynRup->qInterpolatedMinus);

    real** timeDerivativePlusPrefetch = layerDataPrefetch.var(clusterPrefetch->m_dynRup->timeDerivativePlus);
    real** timeDerivativeMinusPrefetch = layerDataPrefetch.var(clusterPrefetch->m_dynRup->timeDerivativeMinus);

    // right now, this ignores the very last cell (for which it may be better just to reload the same face again)
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : 0;

    cluster->m_dynamicRuptureKernel.spaceTimeInterpolation(faceInformation[face],
                                                  cluster->m_globalDataOnHost,
                                                  &godunovData[face],
                                                  &drEnergyOutput[face],
                                                  timeDerivativePlus[face],
                                                  timeDerivativeMinus[face],
                                                  qInterpolatedPlus[face],
                                                  qInterpolatedMinus[face],
                                                  timeDerivativePlusPrefetch[prefetchFace],
                                                  timeDerivativeMinusPrefetch[prefetchFace]);

    // TODO(David): check if this part should be moved into some loop of its own
    auto* frictionSolver = clusterId == clusterOffsetId ? clusters[clusterId]->frictionSolverLocalInterior.get() : clusters[clusterId]->frictionSolverLocalCopy.get();

    frictionSolver->evaluateSingle(face, cluster->m_dynamicRuptureKernel.timeWeights);
  }

  for (size_t i = 0; i < clusters.size(); ++i) {
    TimeCluster* cluster = clusters[i];
    size_t numCells = clusterOffset[i+1] - clusterOffset[i] + clusterOffset[i+1+clusters.size()] - clusterOffset[i+clusters.size()];

    if (clusterOffset[i+1] > clusterOffset[i]) { cluster->dynamicRuptureScheduler->setLastCorrectionStepsInterior(cluster->ct.stepsSinceStart); }
    if (clusterOffset[i+1+clusters.size()] > clusterOffset[i+clusters.size()]) { cluster->dynamicRuptureScheduler->setLastCorrectionStepsCopy(cluster->ct.stepsSinceStart); }

    cluster->m_loopStatistics->end(cluster->m_regionComputeDynamicRupture, numCells, cluster->m_profilingId);
  }
}

std::vector<TimeCluster*> TimeCluster::collectMayCorrect(const std::vector<TimeCluster*>& highPriority, const std::vector<TimeCluster*>& lowPriority, size_t threshold) {
    std::vector<TimeCluster*> clusters;

    // take all high-priority clusters which may predict
    size_t cellCount = 0;
    for (TimeCluster* cluster : highPriority) {
      auto numCells = cluster->m_clusterData->getNumberOfCells();
      if (cluster->getNextLegalAction() == ActorAction::Correct) {
        clusters.push_back(cluster);
        cellCount += numCells;
      }
    }

    // for the low-priority clusters, take all which still fit in (do this sort greedily, for now at least)
    // or, if we have not found a high-priority cluster, take the first-best low-priority cluster at least
    if (cellCount < threshold) {
      size_t cellCount = 0;
      for (TimeCluster* cluster : lowPriority) {
        auto numCells = cluster->m_clusterData->getNumberOfCells();
        if (cluster->getNextLegalAction() == ActorAction::Correct && (cellCount + numCells < threshold || cellCount == 0)) {
          clusters.push_back(cluster);
          cellCount += numCells;
        }
      }
    }

    return clusters;
  }

}

