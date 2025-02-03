// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

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
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 **/

#include "Parallel/MPI.h"
#include <AbstractAPI.h>
#include <Common/Executor.h>
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/LTS.h>
#include <Initializer/Tree/Layer.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/GravitationalFreeSurfaceBC.h>
#include <Kernels/Interface.h>
#include <Kernels/Plasticity.h>
#include <Kernels/PointSourceCluster.h>
#include <Kernels/Precision.h>
#include <Monitoring/ActorStateStatistics.h>
#include <Monitoring/LoopStatistics.h>
#include <Parallel/Helper.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <chrono>
#include <init.h>
#include <tensor.h>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <xdmfwriter/scorep_wrapper.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Kernels/Receiver.h"
#include "Kernels/TimeCommon.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "SeisSol.h"
#include "TimeCluster.h"

#include <cassert>
#include <cstring>

#include "generated_code/kernel.h"

namespace seissol::solver::clustering::computation {

TimeCluster::TimeCluster(unsigned int clusterId,
                         unsigned int globalClusterId,
                         unsigned int profilingId,
                         bool usePlasticity,
                         double maxTimeStepSize,
                         long timeStepRate,
                         bool printProgress,
                         CompoundGlobalData globalData,
                         seissol::initializer::Layer* clusterData,
                         seissol::initializer::LTS* lts,
                         seissol::SeisSol& seissolInstance,
                         LoopStatistics* loopStatistics,
                         ActorStateStatistics* actorStateStatistics,
                         const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                         double priority)
    : CellCluster(maxTimeStepSize,
                  timeStepRate,
#ifdef ACL_DEVICE
                  clusterData->getNumberOfCells() >= deviceHostSwitch() ? Executor::Device
                                                                        : Executor::Host,
#else
                  Executor::Host,
#endif
                  cpuExecutor,
                  priority),
      // cluster ids
      usePlasticity(usePlasticity), seissolInstance(seissolInstance),
      globalDataOnHost(globalData.onHost), globalDataOnDevice(globalData.onDevice),
      layer(clusterData),
      // global data
      lts(lts), sourceCluster(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
      // cells
      loopStatistics(loopStatistics), actorStateStatistics(actorStateStatistics),
      receiverCluster(nullptr), printProgress(printProgress), clusterId(clusterId),
      globalClusterId(globalClusterId), profilingId(profilingId) {
  // assert all pointers are valid
  assert(layer != nullptr);
  assert(globalDataOnHost != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(globalDataOnDevice != nullptr);
  }

  // set timings to zero
  receiverTime = 0;

  timeKernel.setGlobalData(globalData);
  localKernel.setGlobalData(globalData);
  localKernel.setInitConds(&seissolInstance.getMemoryManager().getInitialConditions());
  localKernel.setGravitationalAcceleration(seissolInstance.getGravitationSetup().acceleration);
  neighborKernel.setGlobalData(globalData);

  computeFlops();

  regionComputeLocalIntegration = loopStatistics->getRegion("computeLocalIntegration");
  regionComputeNeighboringIntegration = loopStatistics->getRegion("computeNeighboringIntegration");
  regionComputePointSources = loopStatistics->getRegion("computePointSources");
}

TimeCluster::~TimeCluster() {
#ifndef NDEBUG
  logInfo() << "#(time steps):" << numberOfTimeSteps;
#endif
}

void TimeCluster::setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) {
  this->sourceCluster = std::move(sourceCluster);
}

void TimeCluster::writeReceivers() {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (receiverCluster != nullptr) {
    receiverTime = receiverCluster->calcReceivers(
        receiverTime, ct.time.at(ComputeStep::Correct), timeStepSize(), executor, streamRuntime);
  }
}

std::vector<NeighborCluster>* TimeCluster::getNeighborClusters() { return &neighbors; }

void TimeCluster::computeSources() {
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
  SCOREP_USER_REGION("computeSources", SCOREP_USER_REGION_TYPE_FUNCTION)

  // Return when point sources not initialized. This might happen if there
  // are no point sources on this rank.
  auto* pointSourceCluster = [&]() -> kernels::PointSourceCluster* {
#ifdef ACL_DEVICE
    if (executor == Executor::Device) {
      return sourceCluster.device.get();
    } else {
      return sourceCluster.host.get();
    }
#else
    return sourceCluster.host.get();
#endif
  }();

  if (pointSourceCluster) {
    loopStatistics->begin(regionComputePointSources);
    auto timeStepSizeLocal = timeStepSize();
    pointSourceCluster->addTimeIntegratedPointSources(ct.time.at(ComputeStep::Correct),
                                                      ct.time.at(ComputeStep::Correct) +
                                                          timeStepSizeLocal,
                                                      streamRuntime);
    loopStatistics->end(regionComputePointSources, pointSourceCluster->size(), profilingId);
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
}

void TimeCluster::computeLocalIntegration(bool resetBuffers) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeLocalIntegration);
  // pointer for the call of the ADER-function
  real** buffers = layer->var(lts->buffers);
  real** derivatives = layer->var(lts->derivatives);
  CellMaterialData* materialData = layer->var(lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*lts, *layer);
  const auto acceleration = seissolInstance.getGravitationSetup().acceleration;

  streamRuntime.enqueueOmpFor(layer->getNumberOfCells(), [=](std::size_t cell) {
    auto data = loader.entry(cell);

    // local integration buffer
    alignas(Alignment) real integrationBuffer[tensor::I::size()];
    kernels::LocalTmp tmp(acceleration);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided =
        (data.cellInformation().ltsSetup >> 8) % 2 == 1; // buffers are provided
    const bool resetMyBuffers =
        buffersProvided &&
        ((data.cellInformation().ltsSetup >> 10) % 2 == 0 || resetBuffers); // they should be reset

    real* bufferPointer;
    if (resetMyBuffers) {
      // assert presence of the buffer
      assert(buffers[cell] != nullptr);

      bufferPointer = buffers[cell];
    } else {
      // work on local buffer
      bufferPointer = integrationBuffer;
    }

    timeKernel.computeAder(timeStepSize(), data, tmp, bufferPointer, derivatives[cell], true);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping(*boundaryMapping)[4] = layer->var(lts->boundaryMapping);
    localKernel.computeIntegral(bufferPointer,
                                data,
                                tmp,
                                &materialData[cell],
                                &boundaryMapping[cell],
                                ct.time[ComputeStep::Correct],
                                timeStepSize());

    for (unsigned face = 0; face < 4; ++face) {
      auto& curFaceDisplacements = data.faceDisplacements()[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr &&
          data.cellInformation().faceTypes[face] != FaceType::FreeSurfaceGravity) {
        kernel::addVelocity addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = globalDataOnHost->V3mTo2nFace;
        addVelocityKrnl.selectVelocity = init::selectVelocity::Values;
        addVelocityKrnl.faceDisplacement = data.faceDisplacements()[face];
        addVelocityKrnl.I = bufferPointer;
        addVelocityKrnl.execute(face);
      }
    }

    // TODO: Integrate this step into the kernel
    // We've used a temporary buffer -> need to accumulate update in
    // shared buffer.
    if (!resetMyBuffers && buffersProvided) {
      assert(buffers[cell] != nullptr);

      for (unsigned int dof = 0; dof < tensor::I::size(); ++dof) {
        buffers[cell][dof] += integrationBuffer[dof];
      }
    }
  });

  loopStatistics->end(regionComputeLocalIntegration, layer->getNumberOfCells(), profilingId);
}

#ifdef ACL_DEVICE
void TimeCluster::computeLocalIntegrationDevice(bool resetBuffers) {

  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  loopStatistics->begin(regionComputeLocalIntegration);

  auto& dataTable = layer->getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer->getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer->getConditionalTable<inner_keys::Indices>();

  kernels::LocalData::Loader loader;
  loader.load(*lts, *layer);
  kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);

  const double timeStepWidth = timeStepSize();

  ComputeGraphType graphType{ComputeGraphType::LocalIntegral};
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, true);
  streamRuntime.runGraph(
      computeGraphKey, *layer, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        timeKernel.computeBatchedAder(
            timeStepWidth, tmp, dataTable, materialTable, true, streamRuntime);

        localKernel.computeBatchedIntegral(
            dataTable, materialTable, indicesTable, loader, tmp, timeStepWidth, streamRuntime);
      });

  localKernel.evaluateBatchedTimeDependentBc(dataTable,
                                             indicesTable,
                                             loader,
                                             *layer,
                                             *lts,
                                             ct.time.at(ComputeStep::Correct),
                                             timeStepWidth,
                                             streamRuntime);

  graphType =
      resetBuffers ? ComputeGraphType::AccumulatedVelocities : ComputeGraphType::StreamedVelocities;
  computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime.runGraph(
      computeGraphKey, *layer, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        for (unsigned face = 0; face < 4; ++face) {
          const ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
          if (dataTable.find(key) != dataTable.end()) {
            auto& entry = dataTable[key];
            // NOTE: integrated velocities have been computed implicitly, i.e
            // it is 6th, 7the and 8th columns of integrated dofs

            kernel::gpu_addVelocity displacementKrnl;
            displacementKrnl.faceDisplacement =
                entry.get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
            displacementKrnl.integratedVelocities = const_cast<const real**>(
                entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
            displacementKrnl.V3mTo2nFace = globalDataOnDevice->V3mTo2nFace;

            // Note: this kernel doesn't require tmp. memory
            displacementKrnl.numElements =
                entry.get(inner_keys::Wp::Id::FaceDisplacement)->getSize();
            displacementKrnl.streamPtr = streamRuntime.stream();
            displacementKrnl.execute(face);
          }
        }

        const ConditionalKey key =
            ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
        if (dataTable.find(key) != dataTable.end()) {
          auto& entry = dataTable[key];

          if (resetBuffers) {
            device.algorithms.streamBatchedData(
                (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
                (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime.stream());
          } else {
            device.algorithms.accumulateBatchedData(
                (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
                (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime.stream());
          }
        }
      });

  loopStatistics->end(regionComputeLocalIntegration, layer->getNumberOfCells(), profilingId);
  device.api->popLastProfilingMark();
}
#endif // ACL_DEVICE

void TimeCluster::computeNeighboringIntegration(double subTimeStart) {
  if (usePlasticity) {
    computeNeighboringIntegrationImplementation<true>(subTimeStart);
  } else {
    computeNeighboringIntegrationImplementation<false>(subTimeStart);
  }
}
#ifdef ACL_DEVICE
void TimeCluster::computeNeighboringIntegrationDevice(double subTimeStart) {
  device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  loopStatistics->begin(regionComputeNeighboringIntegration);

  const double timeStepWidth = timeStepSize();
  auto& table = layer->getConditionalTable<inner_keys::Wp>();

  seissol::kernels::TimeCommon::computeBatchedIntegrals(
      timeKernel, subTimeStart, timeStepWidth, table, streamRuntime);

  const ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime.runGraph(
      computeGraphKey, *layer, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        neighborKernel.computeBatchedNeighborsIntegral(table, streamRuntime);
      });

  if (usePlasticity) {
    auto oneMinusIntegratingFactor = getRelaxTime();
    auto* plasticity = layer->var(lts->plasticity, seissol::initializer::AllocationPlace::Device);
    const unsigned numAdjustedDofs =
        seissol::kernels::Plasticity::computePlasticityBatched(oneMinusIntegratingFactor,
                                                               timeStepWidth,
                                                               tv,
                                                               globalDataOnDevice,
                                                               table,
                                                               plasticity,
                                                               streamRuntime);

    seissolInstance.flopCounter().incrementNonZeroFlopsPlasticity(
        layer->getNumberOfCells() * flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)] +
        numAdjustedDofs * flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)]);
    seissolInstance.flopCounter().incrementHardwareFlopsPlasticity(
        layer->getNumberOfCells() * flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)] +
        numAdjustedDofs * flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]);
  }

  device.api->popLastProfilingMark();
  loopStatistics->end(regionComputeNeighboringIntegration, layer->getNumberOfCells(), profilingId);
}
#endif // ACL_DEVICE

void TimeCluster::computeLocalIntegrationFlops() {
  auto& flopsNonZero = flops_nonZero[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = flops_hardware[static_cast<int>(ComputePart::Local)];
  flopsNonZero = 0;
  flopsHardware = 0;

  auto* cellInformation = layer->var(lts->cellInformation);
  for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
    unsigned cellNonZero, cellHardware;
    timeKernel.flopsAder(cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    // Contribution from displacement/integrated displacement
    for (unsigned face = 0; face < 4; ++face) {
      if (cellInformation->faceTypes[face] == FaceType::FreeSurfaceGravity) {
        const auto [nonZeroFlopsDisplacement, hardwareFlopsDisplacement] =
            GravitationalFreeSurfaceBc::getFlopsDisplacementFace(
                face, cellInformation[cell].faceTypes[face]);
        flopsNonZero += nonZeroFlopsDisplacement;
        flopsHardware += hardwareFlopsDisplacement;
      }
    }
  }
}

void TimeCluster::computeNeighborIntegrationFlops() {
  auto& flopsNonZero = flops_nonZero[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = flops_hardware[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)];
  auto& drFlopsHardware = flops_hardware[static_cast<int>(ComputePart::DRNeighbor)];
  flopsNonZero = 0;
  flopsHardware = 0;
  drFlopsNonZero = 0;
  drFlopsHardware = 0;

  auto* cellInformation = layer->var(lts->cellInformation);
  auto* drMapping = layer->var(lts->drMapping);
  for (unsigned cell = 0; cell < layer->getNumberOfCells(); ++cell) {
    unsigned cellNonZero, cellHardware;
    long long cellDRNonZero, cellDRHardware;
    neighborKernel.flopsNeighborsIntegral(cellInformation[cell].faceTypes,
                                          cellInformation[cell].faceRelations,
                                          drMapping[cell],
                                          cellNonZero,
                                          cellHardware,
                                          cellDRNonZero,
                                          cellDRHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    drFlopsNonZero += cellDRNonZero;
    drFlopsHardware += cellDRHardware;

    /// \todo add lts time integration
    /// \todo add plasticity
  }
}

void TimeCluster::computeFlops() {
  computeLocalIntegrationFlops();
  computeNeighborIntegrationFlops();
  seissol::kernels::Plasticity::flopsPlasticity(
      flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)],
      flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)],
      flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)],
      flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]);
}

ActResult TimeCluster::act() {
  actorStateStatistics->enter(state.step);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics->enter(state.step);
  return result;
}

void TimeCluster::handleAdvancedComputeTimeMessage(ComputeStep step,
                                                   const NeighborCluster& neighborCluster) {
  if (step == ComputeStep::Interact && neighborCluster.ct.maxTimeStepSize > ct.maxTimeStepSize) {
    if (neighborCluster.ct.time.find(ComputeStep::Correct) != neighborCluster.ct.time.end()) {
      lastSubTime = neighborCluster.ct.time.at(ComputeStep::Correct);
    } else {
      lastSubTime = 0;
    }
  }
}
void TimeCluster::predict() {
  if (layer->getNumberOfCells() == 0)
    return;

  bool resetBuffers = true;
  for (auto& neighbor : neighbors) {
    if (neighbor.ct.timeStepRate > ct.timeStepRate &&
        ct.computeSinceLastSync.at(ComputeStep::Correct) >
            neighbor.ct.computeSinceLastSync.at(ComputeStep::Correct)) {
      resetBuffers = false;
    }
  }
  if (ct.computeSinceLastSync.at(ComputeStep::Correct) == 0) {
    resetBuffers = true;
  }

  writeReceivers();
#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeLocalIntegrationDevice(resetBuffers);
  } else {
    computeLocalIntegration(resetBuffers);
  }
#else
  computeLocalIntegration(resetBuffers);
#endif
  computeSources();

  seissolInstance.flopCounter().incrementNonZeroFlopsLocal(
      flops_nonZero[static_cast<int>(ComputePart::Local)]);
  seissolInstance.flopCounter().incrementHardwareFlopsLocal(
      flops_hardware[static_cast<int>(ComputePart::Local)]);
#ifdef ACL_DEVICE
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                              : seissol::initializer::AllocationPlace::Device;
    layer->bucketSynchronizeTo(lts->buffersDerivatives, other, streamRuntime.stream());
  }
#endif
}

void TimeCluster::correct() {
  /* Sub start time of width respect to the next cluster; use 0 if not relevant, for example in GTS.
   * LTS requires to evaluate a partial time integration of the derivatives. The point zero in time
   * refers to the derivation of the surrounding time derivatives, which coincides with the last
   * completed time step of the next cluster. The start/end of the time step is the start/end of
   * this clusters time step relative to the zero point.
   *   Example:
   *                                              5 dt
   *   |-----------------------------------------------------------------------------------------|
   * <<< Time stepping of the next cluster (Cn) (5x larger than the current). |                 | |
   * |                 |                 |
   *   |*****************|*****************|+++++++++++++++++|                 |                 |
   * <<< Status of the current cluster. |                 |                 |                 | | |
   *   |-----------------|-----------------|-----------------|-----------------|-----------------|
   * <<< Time stepping of the current cluster (Cc). 0                 dt               2dt 3dt 4dt
   * 5dt
   *
   *   In the example above two clusters are illustrated: Cc and Cn. Cc is the current cluster under
   * consideration and Cn the next cluster with respect to LTS terminology. Cn is currently at time
   * 0 and provided Cc with derivatives valid until 5dt. Cc updated already twice and did its last
   * full update to reach 2dt (== subTimeStart). Next computeNeighboringCopy is called to accomplish
   * the next full update to reach 3dt (+++). Besides working on the buffers of own buffers and
   * those of previous clusters, Cc needs to evaluate the time prediction of Cn in the interval
   * [2dt, 3dt].
   */
  const double subTimeStart = ct.time.at(ComputeStep::Correct) - lastSubTime;

#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeNeighboringIntegrationDevice(subTimeStart);
  } else {
    computeNeighboringIntegration(subTimeStart);
  }
#else
  computeNeighboringIntegration(subTimeStart);
#endif

  seissolInstance.flopCounter().incrementNonZeroFlopsNeighbor(
      flops_nonZero[static_cast<int>(ComputePart::Neighbor)]);
  seissolInstance.flopCounter().incrementHardwareFlopsNeighbor(
      flops_hardware[static_cast<int>(ComputePart::Neighbor)]);

  // TODO(Lukas) Adjust with time step rate? Relevant is maximum cluster is not on this node
  const auto nextCorrectionSteps = ct.nextSteps();
  if constexpr (USE_MPI) {
    if (printProgress && (((nextCorrectionSteps / timeStepRate) % 100) == 0)) {
      logInfo() << "#max-updates since sync: " << nextCorrectionSteps << " @ "
                << ct.nextComputeTime(ComputeStep::Correct, syncTime);
    }
  }
}

void TimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logWarning(true) << "No update since " << timeSinceLastUpdate.count() << "[s] for global cluster "
                   << globalClusterId << " with local cluster id " << clusterId << " at state "
                   << actorStateToString(state)
                   << " predTime = " << ct.time.at(ComputeStep::Predict)
                   << " predictionsSinceSync = " << ct.computeSinceLastSync.at(ComputeStep::Predict)
                   << " corrTime = " << ct.time.at(ComputeStep::Correct)
                   << " correctionsSinceSync = " << ct.computeSinceLastSync.at(ComputeStep::Correct)
                   << " stepsTillSync = " << ct.stepsUntilSync << " maySync = " << maySynchronize();
  for (auto& neighbor : neighbors) {
    logWarning(true) << "Neighbor with rate = " << neighbor.ct.timeStepRate
                     << "PredTime = " << neighbor.ct.time.at(ComputeStep::Predict)
                     << "CorrTime = " << neighbor.ct.time.at(ComputeStep::Correct)
                     << "predictionsSinceSync = "
                     << neighbor.ct.computeSinceLastSync.at(ComputeStep::Predict)
                     << "correctionsSinceSync = "
                     << neighbor.ct.computeSinceLastSync.at(ComputeStep::Correct);
  }
}

unsigned int TimeCluster::getClusterId() const { return clusterId; }

unsigned int TimeCluster::getGlobalClusterId() const { return globalClusterId; }

void TimeCluster::setReceiverTime(double receiverTime) { this->receiverTime = receiverTime; }

void TimeCluster::finalize() {
  // sourceCluster.host.reset(nullptr);
  // sourceCluster.device.reset(nullptr);
  streamRuntime.dispose();
}

template <bool UsePlasticity>
std::pair<long, long>
    TimeCluster::computeNeighboringIntegrationImplementation(double subTimeStart) {
  if (layer->getNumberOfCells() == 0)
    return {0, 0};
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  auto tv = this->tv;
  auto oneMinusIntegratingFactor = getRelaxTime();
  auto timeStep = timeStepSize();
  const std::size_t numberOfCells = layer->getNumberOfCells();

  loopStatistics->begin(regionComputeNeighboringIntegration);

  auto* faceNeighbors = layer->var(lts->faceNeighbors);
  auto* drMapping = layer->var(lts->drMapping);
  auto* cellInformation = layer->var(lts->cellInformation);
  auto* plasticity = layer->var(lts->plasticity);
  auto* pstrain = layer->var(lts->pstrain);

  kernels::NeighborData::Loader loader;
  loader.load(*lts, *layer);

  real* timeIntegrated[4];
  real* faceNeighborsPrefetch[4];

  auto globalDataOnHost = this->globalDataOnHost;
  auto* self = this;

  streamRuntime.enqueueOmpFor(numberOfCells, [=](std::size_t cell) {
    real* timeIntegrated[4];
    real* faceNeighborsPrefetch[4];
    auto data = loader.entry(cell);
    seissol::kernels::TimeCommon::computeIntegrals(
        self->timeKernel,
        data.cellInformation().ltsSetup,
        data.cellInformation().faceTypes,
        subTimeStart,
        timeStep,
        faceNeighbors[cell],
#ifdef _OPENMP
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(&(
            globalDataOnHost->integrationBufferLTS[omp_get_thread_num() * 4 * tensor::I::size()])),
#else
          *reinterpret_cast<real(*)[4][tensor::I::size()]>(globalDataOnHost->integrationBufferLTS),
#endif
        timeIntegrated);

    faceNeighborsPrefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
                                   ? faceNeighbors[cell][1]
                                   : drMapping[cell][1].godunov;
    faceNeighborsPrefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
                                   ? faceNeighbors[cell][2]
                                   : drMapping[cell][2].godunov;
    faceNeighborsPrefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
                                   ? faceNeighbors[cell][3]
                                   : drMapping[cell][3].godunov;

    // fourth face's prefetches
    if (cell < (numberOfCells - 1)) {
      faceNeighborsPrefetch[3] =
          (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
              ? faceNeighbors[cell + 1][0]
              : drMapping[cell + 1][0].godunov;
    } else {
      faceNeighborsPrefetch[3] = faceNeighbors[cell][3];
    }

    self->neighborKernel.computeNeighborsIntegral(
        data, drMapping[cell], timeIntegrated, faceNeighborsPrefetch);

    if constexpr (UsePlasticity) {
      seissol::kernels::Plasticity::computePlasticity(oneMinusIntegratingFactor,
                                                      timeStep,
                                                      tv,
                                                      globalDataOnHost,
                                                      &plasticity[cell],
                                                      data.dofs(),
                                                      pstrain[cell]);
    }
#ifdef INTEGRATE_QUANTITIES
    seissolInstance.postProcessor().integrateQuantities(timeStepWidth, *layer, cell, dofs[cell]);
#endif // INTEGRATE_QUANTITIES
  });

  const auto numberOTetsWithPlasticYielding = 0;
  const long long nonZeroFlopsPlasticity =
      numberOfCells * flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)] +
      numberOTetsWithPlasticYielding *
          flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)];
  const long long hardwareFlopsPlasticity =
      numberOfCells * flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)] +
      numberOTetsWithPlasticYielding *
          flops_hardware[static_cast<int>(ComputePart::PlasticityYield)];

  loopStatistics->end(regionComputeNeighboringIntegration, numberOfCells, profilingId);

  return {0, 0};
}

void TimeCluster::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
#ifdef ACL_DEVICE
  if ((place == initializer::AllocationPlace::Host && executor == Executor::Device) ||
      (place == initializer::AllocationPlace::Device && executor == Executor::Host)) {
    layer->synchronizeTo(place, stream);
  }
#endif
}

void TimeCluster::runCompute(ComputeStep step) {
  if (step == ComputeStep::Predict) {
    predict();
  }
  if (step == ComputeStep::Correct) {
    correct();
  }
}

} // namespace seissol::solver::clustering::computation
