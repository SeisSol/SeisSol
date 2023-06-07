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
#include <SourceTerm/PointSource.h>
#include <Kernels/DynamicRupture.h>
#include <Kernels/Receiver.h>
#include <Monitoring/FlopCounter.hpp>
#include <Monitoring/instrumentation.hpp>
#include "Common/cellconfigconv.hpp"
#include "Kernels/Plasticity.h"

#include <cassert>
#include <cstring>

#include <generated_code/kernel.h>

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
    dynamicRuptureScheduler(dynamicRuptureScheduler),
    waveProp(std::move(waveprop::getDispatcher(*i_lts, *i_clusterData, defaultConfig(usePlasticity))))
{
    // assert all pointers are valid
    assert( m_clusterData                              != nullptr );
    assert( m_globalDataOnHost                         != nullptr );
    if constexpr (seissol::isDeviceOn()) {
        assert( m_globalDataOnDevice                   != nullptr );
    }

  // set timings to zero
  m_receiverTime                  = 0;

  m_dynamicRuptureKernel.setGlobalData(i_globalData);

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

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Local)];
  
  waveProp->computePredictFlops(flopsNonZero, flopsHardware);
}

void seissol::time_stepping::TimeCluster::computeNeighborIntegrationFlops(seissol::initializers::Layer& layerData) {
  auto& flopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = m_flops_hardware[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)];
  auto& drFlopsHardware = m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)];
  
  waveProp->computeCorrectFlops(flopsNonZero, flopsHardware, drFlopsNonZero, drFlopsHardware);
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
  m_loopStatistics->begin(m_regionComputeLocalIntegration);
  waveProp->dispatchPredict(timeStepSize(), ct.correctionTime, resetBuffers);
  m_loopStatistics->end(m_regionComputeLocalIntegration, m_clusterData->getNumberOfCells(), m_profilingId);
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

  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);
  waveProp->dispatchNeighborCorrect(timeStepSize(), subTimeStart);
  m_loopStatistics->end(m_regionComputeNeighboringIntegration, m_clusterData->getNumberOfCells(), m_profilingId);

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

}

