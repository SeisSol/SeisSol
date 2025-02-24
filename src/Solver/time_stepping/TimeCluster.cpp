// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)
// SPDX-FileContributor: Sebastian Rettenberger

#include "Parallel/MPI.h"
#include <Common/Executor.h>
#include <Memory/Tree/Layer.h>
#include <Kernels/PointSourceCluster.h>
#include <SourceTerm/Manager.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SeisSol.h"
#include "TimeCluster.h"
#include "SourceTerm/PointSource.h"
#include "Kernels/TimeCommon.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Receiver.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"

#include <cassert>
#include <cstring>

#include "generated_code/kernel.h"

seissol::time_stepping::TimeCluster::TimeCluster(unsigned int i_clusterId, unsigned int i_globalClusterId,
                                                 unsigned int profilingId,
                                                 bool usePlasticity,
                                                 LayerType layerType, double maxTimeStepSize,
                                                 long timeStepRate, bool printProgress,
                                                 DynamicRuptureScheduler *dynamicRuptureScheduler,
                                                 CompoundGlobalData i_globalData,
                                                 seissol::initializer::Layer *i_clusterData,
                                                 seissol::initializer::Layer *dynRupInteriorData,
                                                 seissol::initializer::Layer *dynRupCopyData,
                                                 seissol::initializer::LTS *i_lts,
                                                 seissol::initializer::DynamicRupture* i_dynRup,
                                                 seissol::dr::friction_law::FrictionSolver* i_FrictionSolver,
                                                 seissol::dr::friction_law::FrictionSolver* i_FrictionSolverDevice,
                                                 dr::output::OutputManager* i_faultOutputManager,
                                                 seissol::SeisSol& seissolInstance,
                                                 LoopStatistics *i_loopStatistics,
                                                 ActorStateStatistics* actorStateStatistics) :
    AbstractTimeCluster(maxTimeStepSize, timeStepRate,
#ifdef ACL_DEVICE
      seissolInstance.executionPlace(i_clusterData->getNumberOfCells())
#else
      Executor::Host
#endif
    ),
    // cluster ids
    usePlasticity(usePlasticity),
    seissolInstance(seissolInstance),
    m_globalDataOnHost( i_globalData.onHost ),
    m_globalDataOnDevice(i_globalData.onDevice ),
    m_clusterData(i_clusterData),
    // global data
    dynRupInteriorData(dynRupInteriorData),
    dynRupCopyData(dynRupCopyData),
    m_lts(i_lts),
    m_dynRup(i_dynRup),
    frictionSolver(i_FrictionSolver),
    frictionSolverDevice(i_FrictionSolverDevice),
    faultOutputManager(i_faultOutputManager),
    m_sourceCluster(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
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
  m_localKernel.setInitConds(&seissolInstance.getMemoryManager().getInitialConditions());
  m_localKernel.setGravitationalAcceleration(seissolInstance.getGravitationSetup().acceleration);
  m_neighborKernel.setGlobalData(i_globalData);
  m_dynamicRuptureKernel.setGlobalData(i_globalData);

  computeFlops();

  m_regionComputeLocalIntegration = m_loopStatistics->getRegion("computeLocalIntegration");
  m_regionComputeNeighboringIntegration = m_loopStatistics->getRegion("computeNeighboringIntegration");
  m_regionComputeDynamicRupture = m_loopStatistics->getRegion("computeDynamicRupture");
  m_regionComputePointSources = m_loopStatistics->getRegion("computePointSources");
}

seissol::time_stepping::TimeCluster::~TimeCluster() {
#ifndef NDEBUG
  logInfo() << "#(time steps):" << numberOfTimeSteps;
#endif
}

void seissol::time_stepping::TimeCluster::setPointSources(
    seissol::kernels::PointSourceClusterPair sourceCluster) {
  m_sourceCluster = std::move(sourceCluster);
}

void seissol::time_stepping::TimeCluster::writeReceivers() {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (m_receiverCluster != nullptr) {
    m_receiverTime = m_receiverCluster->calcReceivers(m_receiverTime, ct.correctionTime, timeStepSize(), executor, nullptr);
  }
}

std::vector<seissol::time_stepping::NeighborCluster>*
    seissol::time_stepping::TimeCluster::getNeighborClusters() {
  return &neighbors;
}

void seissol::time_stepping::TimeCluster::computeSources() {
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
  SCOREP_USER_REGION( "computeSources", SCOREP_USER_REGION_TYPE_FUNCTION )

  // Return when point sources not initialized. This might happen if there
  // are no point sources on this rank.
  auto* pointSourceCluster = [&]() -> kernels::PointSourceCluster* {
#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    return m_sourceCluster.device.get();
  }
  else {
    return m_sourceCluster.host.get();
  }
#else
  return m_sourceCluster.host.get();
#endif
  }();

  if (pointSourceCluster) {
    m_loopStatistics->begin(m_regionComputePointSources);
    auto timeStepSizeLocal = timeStepSize();
    pointSourceCluster->addTimeIntegratedPointSources(ct.correctionTime, ct.correctionTime + timeStepSizeLocal, streamRuntime);
    m_loopStatistics->end(m_regionComputePointSources, pointSourceCluster->size(), m_profilingId);
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
}

void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializer::Layer&  layerData ) {
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
                           m_dynamicRuptureKernel.timeWeights,
                           streamRuntime);
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel 
  {
  LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_profilingId);
}

#ifdef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeDynamicRuptureDevice( seissol::initializer::Layer&  layerData ) {
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  if (layerData.getNumberOfCells() > 0) {
    // compute space time interpolation part

    const double stepSizeWidth = timeStepSize();
    ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    auto computeGraphKey = initializer::GraphKey(graphType, stepSizeWidth);
    auto& table = layerData.getConditionalTable<inner_keys::Dr>();
    m_dynamicRuptureKernel.setTimeStepWidth(stepSizeWidth);
    streamRuntime.runGraph(computeGraphKey, layerData, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
      m_dynamicRuptureKernel.batchedSpaceTimeInterpolation(table, streamRuntime);
    });
    device.api->popLastProfilingMark();
    if (frictionSolverDevice->allocationPlace() == initializer::AllocationPlace::Host) {
      layerData.varSynchronizeTo(m_dynRup->qInterpolatedPlus, initializer::AllocationPlace::Host, streamRuntime.stream());
      layerData.varSynchronizeTo(m_dynRup->qInterpolatedMinus, initializer::AllocationPlace::Host, streamRuntime.stream());
      streamRuntime.wait();
    }

    device.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
    frictionSolverDevice->computeDeltaT(m_dynamicRuptureKernel.timePoints);
    frictionSolverDevice->evaluate(layerData,
                             m_dynRup,
                             ct.correctionTime,
                             m_dynamicRuptureKernel.timeWeights,
                             streamRuntime);
    device.api->popLastProfilingMark();
    if (frictionSolverDevice->allocationPlace() == initializer::AllocationPlace::Host) {
      layerData.varSynchronizeTo(m_dynRup->fluxSolverMinus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo(m_dynRup->fluxSolverPlus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo(m_dynRup->imposedStateMinus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo(m_dynRup->imposedStatePlus, initializer::AllocationPlace::Device, streamRuntime.stream());
    }
    streamRuntime.wait();
  }
  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells(), m_profilingId);
}
#endif


void seissol::time_stepping::TimeCluster::computeDynamicRuptureFlops( seissol::initializer::Layer& layerData,
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

void seissol::time_stepping::TimeCluster::computeLocalIntegration(seissol::initializer::Layer& i_layerData, bool resetBuffers ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // local integration buffer
  alignas(Alignment) real l_integrationBuffer[tensor::I::size()];

  // pointer for the call of the ADER-function
  real* l_bufferPointer;

  real** buffers = i_layerData.var(m_lts->buffers);
  real** derivatives = i_layerData.var(m_lts->derivatives);
  CellMaterialData* materialData = i_layerData.var(m_lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer), firstprivate(tmp) schedule(static)
#endif
  for (unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++) {
    auto data = loader.entry(l_cell);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided = (data.cellInformation().ltsSetup >> 8) % 2 == 1; // buffers are provided
    const bool resetMyBuffers = buffersProvided && ( (data.cellInformation().ltsSetup >> 10) %2 == 0 || resetBuffers ); // they should be reset

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
      auto& curFaceDisplacements = data.faceDisplacements()[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr
          && data.cellInformation().faceTypes[face] != FaceType::FreeSurfaceGravity) {
        kernel::addVelocity addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = m_globalDataOnHost->V3mTo2nFace;
        addVelocityKrnl.selectVelocity = init::selectVelocity::Values;
        addVelocityKrnl.faceDisplacement = data.faceDisplacements()[face];
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
#ifdef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeLocalIntegrationDevice(
  seissol::initializer::Layer& i_layerData,
  bool resetBuffers) {

  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  auto& dataTable = i_layerData.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = i_layerData.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = i_layerData.getConditionalTable<inner_keys::Indices>();

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);

  const double timeStepWidth = timeStepSize();

  ComputeGraphType graphType = resetBuffers ? ComputeGraphType::AccumulatedVelocities : ComputeGraphType::StreamedVelocities;
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, true);
  streamRuntime.runGraph(computeGraphKey, i_layerData, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
    m_timeKernel.computeBatchedAder(timeStepWidth,
                                    tmp,
                                    dataTable,
                                    materialTable,
                                    true,
                                    streamRuntime);

    m_localKernel.computeBatchedIntegral(dataTable,
                                         materialTable,
                                         indicesTable,
                                         loader,
                                         tmp,
                                         timeStepWidth,
                                         streamRuntime);

    m_localKernel.evaluateBatchedTimeDependentBc(dataTable,
                                               indicesTable,
                                               loader,
                                               i_layerData,
                                               *m_lts,
                                               ct.correctionTime,
                                               timeStepWidth,
                                               streamRuntime);

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
        displacementKrnl.streamPtr = streamRuntime.stream();
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

  streamRuntime.wait();

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells(), m_profilingId);
  device.api->popLastProfilingMark();
}
#endif // ACL_DEVICE

void seissol::time_stepping::TimeCluster::computeNeighboringIntegration(seissol::initializer::Layer& i_layerData,
                                                                        double subTimeStart) {
  if (usePlasticity) {
    computeNeighboringIntegrationImplementation<true>(i_layerData, subTimeStart);
  } else {
    computeNeighboringIntegrationImplementation<false>(i_layerData, subTimeStart);
  }
}
#ifdef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeNeighboringIntegrationDevice( seissol::initializer::Layer&  i_layerData,
                                                                         double subTimeStart) {
  device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )
  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  const double timeStepWidth = timeStepSize();
  auto& table = i_layerData.getConditionalTable<inner_keys::Wp>();

  seissol::kernels::TimeCommon::computeBatchedIntegrals(m_timeKernel,
                                                        subTimeStart,
                                                        timeStepWidth,
                                                        table,
                                                        streamRuntime);

  ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime.runGraph(computeGraphKey, i_layerData, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
    m_neighborKernel.computeBatchedNeighborsIntegral(table, streamRuntime);
  });

  if (usePlasticity) {
    updateRelaxTime();
    auto* plasticity = i_layerData.var(m_lts->plasticity, seissol::initializer::AllocationPlace::Device);
    unsigned numAdjustedDofs = seissol::kernels::Plasticity::computePlasticityBatched(m_oneMinusIntegratingFactor,
                                                                                      timeStepWidth,
                                                                                      m_tv,
                                                                                      m_globalDataOnDevice,
                                                                                      table,
                                                                                      plasticity,
                                                                                      streamRuntime);

    seissolInstance.flopCounter().incrementNonZeroFlopsPlasticity(
        i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)]);
    seissolInstance.flopCounter().incrementHardwareFlopsPlasticity(
        i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)]
        + numAdjustedDofs * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)]);
  }

  device.api->popLastProfilingMark();
  streamRuntime.wait();
  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_profilingId);
}
#endif // ACL_DEVICE

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops(seissol::initializer::Layer& layerData) {
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
      if (cellInformation->faceTypes[face] == FaceType::FreeSurfaceGravity) {
        const auto [nonZeroFlopsDisplacement, hardwareFlopsDisplacement] =
        GravitationalFreeSurfaceBc::getFlopsDisplacementFace(face,
                                                             cellInformation[cell].faceTypes[face]);
        flopsNonZero += nonZeroFlopsDisplacement;
        flopsHardware += hardwareFlopsDisplacement;
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeNeighborIntegrationFlops(
    seissol::initializer::Layer& layerData) {
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
#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeLocalIntegrationDevice(*m_clusterData, resetBuffers);
  }
  else {
    computeLocalIntegration(*m_clusterData, resetBuffers);
  }
#else
  computeLocalIntegration(*m_clusterData, resetBuffers);
#endif
  computeSources();

  seissolInstance.flopCounter().incrementNonZeroFlopsLocal(m_flops_nonZero[static_cast<int>(ComputePart::Local)]);
  seissolInstance.flopCounter().incrementHardwareFlopsLocal(m_flops_hardware[static_cast<int>(ComputePart::Local)]);
#ifdef ACL_DEVICE
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host : seissol::initializer::AllocationPlace::Device;
    m_clusterData->bucketSynchronizeTo(m_lts->buffersDerivatives, other, streamRuntime.stream());
    streamRuntime.wait();
  }
#endif
}

void TimeCluster::handleDynamicRupture(initializer::Layer& layerData) {
#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeDynamicRuptureDevice(layerData);
  }
  else {
    computeDynamicRupture(layerData);
  }

  // TODO(David): restrict to copy/interior of same cluster type
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host : seissol::initializer::AllocationPlace::Device;
    layerData.varSynchronizeTo(m_dynRup->fluxSolverMinus, other, streamRuntime.stream());
    layerData.varSynchronizeTo(m_dynRup->fluxSolverPlus, other, streamRuntime.stream());
    layerData.varSynchronizeTo(m_dynRup->imposedStateMinus, other, streamRuntime.stream());
    layerData.varSynchronizeTo(m_dynRup->imposedStatePlus, other, streamRuntime.stream());
    streamRuntime.wait();
  }
#else
  computeDynamicRupture(layerData);
#endif
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
      handleDynamicRupture(*dynRupInteriorData);
      seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      dynamicRuptureScheduler->setLastCorrectionStepsInterior(ct.stepsSinceStart);
    }
    if (layerType == Copy) {
      handleDynamicRupture(*dynRupCopyData);
      seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      dynamicRuptureScheduler->setLastCorrectionStepsCopy((ct.stepsSinceStart));
    }

  }

#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeNeighboringIntegrationDevice(*m_clusterData, subTimeStart);
  }
  else {
    computeNeighboringIntegration(*m_clusterData, subTimeStart);
  }
#else
  computeNeighboringIntegration(*m_clusterData, subTimeStart);
#endif

  seissolInstance.flopCounter().incrementNonZeroFlopsNeighbor(m_flops_nonZero[static_cast<int>(ComputePart::Neighbor)]);
  seissolInstance.flopCounter().incrementHardwareFlopsNeighbor(m_flops_hardware[static_cast<int>(ComputePart::Neighbor)]);
  seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRNeighbor)]);
  seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRNeighbor)]);

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
      logInfo() << "#max-updates since sync: " << nextCorrectionSteps
                    << " @ " << ct.nextCorrectionTime(syncTime);

      }
  }
}

void TimeCluster::reset() {
    AbstractTimeCluster::reset();
}

void TimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logWarning(true)
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
    logWarning(true)
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

void TimeCluster::finalize() {
  streamRuntime.dispose();
}

template<bool usePlasticity>
    std::pair<long, long> TimeCluster::computeNeighboringIntegrationImplementation(seissol::initializer::Layer& i_layerData,
                                                                      double subTimeStart) {
      if (i_layerData.getNumberOfCells() == 0) return {0,0};
      SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

      m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

      real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
      CellDRMapping (*drMapping)[4] = i_layerData.var(m_lts->drMapping);
      CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);
      auto* plasticity = i_layerData.var(m_lts->plasticity);
      auto* pstrain = i_layerData.var(m_lts->pstrain);
      unsigned numberOTetsWithPlasticYielding = 0;

      kernels::NeighborData::Loader loader;
      loader.load(*m_lts, i_layerData);

      real *l_timeIntegrated[4];
      real *l_faceNeighbors_prefetch[4];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) private(l_timeIntegrated, l_faceNeighbors_prefetch) shared(cellInformation, loader, faceNeighbors, pstrain, i_layerData, plasticity, drMapping, subTimeStart) reduction(+:numberOTetsWithPlasticYielding)
#endif
      for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
        auto data = loader.entry(l_cell);
        seissol::kernels::TimeCommon::computeIntegrals(m_timeKernel,
                                                       data.cellInformation().ltsSetup,
                                                       data.cellInformation().faceTypes,
                                                       subTimeStart,
                                                       timeStepSize(),
                                                       faceNeighbors[l_cell],
#ifdef _OPENMP
                                                       *reinterpret_cast<real (*)[4][tensor::I::size()]>(&(m_globalDataOnHost->integrationBufferLTS[omp_get_thread_num()*4*tensor::I::size()])),
#else
            *reinterpret_cast<real (*)[4][tensor::I::size()]>(m_globalDataOnHost->integrationBufferLTS),
#endif
                                                       l_timeIntegrated);

        l_faceNeighbors_prefetch[0] = (cellInformation[l_cell].faceTypes[1] != FaceType::DynamicRupture) ?
                                      faceNeighbors[l_cell][1] :
                                      drMapping[l_cell][1].godunov;
        l_faceNeighbors_prefetch[1] = (cellInformation[l_cell].faceTypes[2] != FaceType::DynamicRupture) ?
                                      faceNeighbors[l_cell][2] :
                                      drMapping[l_cell][2].godunov;
        l_faceNeighbors_prefetch[2] = (cellInformation[l_cell].faceTypes[3] != FaceType::DynamicRupture) ?
                                      faceNeighbors[l_cell][3] :
                                      drMapping[l_cell][3].godunov;

        // fourth face's prefetches
        if (l_cell < (i_layerData.getNumberOfCells()-1) ) {
          l_faceNeighbors_prefetch[3] = (cellInformation[l_cell+1].faceTypes[0] != FaceType::DynamicRupture) ?
                                        faceNeighbors[l_cell+1][0] :
                                        drMapping[l_cell+1][0].godunov;
        } else {
          l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
        }

        m_neighborKernel.computeNeighborsIntegral( data,
                                                   drMapping[l_cell],
                                                   l_timeIntegrated, l_faceNeighbors_prefetch
        );

        if constexpr (usePlasticity) {
          updateRelaxTime();
          numberOTetsWithPlasticYielding += seissol::kernels::Plasticity::computePlasticity( m_oneMinusIntegratingFactor,
                                                                                             timeStepSize(),
                                                                                             m_tv,
                                                                                             m_globalDataOnHost,
                                                                                             &plasticity[l_cell],
                                                                                             data.dofs(),
                                                                                             pstrain[l_cell] );
        }
#ifdef INTEGRATE_QUANTITIES
        seissolInstance.postProcessor().integrateQuantities( m_timeStepWidth,
                                                              i_layerData,
                                                              l_cell,
                                                              dofs[l_cell] );
#endif // INTEGRATE_QUANTITIES
      }

      const long long nonZeroFlopsPlasticity =
          i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)] +
          numberOTetsWithPlasticYielding * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)];
      const long long hardwareFlopsPlasticity =
          i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)] +
          numberOTetsWithPlasticYielding * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)];

      m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_profilingId);

      return {nonZeroFlopsPlasticity, hardwareFlopsPlasticity};
    }

void TimeCluster::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
#ifdef ACL_DEVICE
  if ((place == initializer::AllocationPlace::Host && executor == Executor::Device) || (place == initializer::AllocationPlace::Device && executor == Executor::Host)) {
    m_clusterData->synchronizeTo(place, stream);
    if (layerType == Interior) {
      dynRupInteriorData->synchronizeTo(place, stream);
    }
    if (layerType == Copy) {
      dynRupCopyData->synchronizeTo(place, stream);
    }
  }
#endif
}

} // namespace seissol::time_stepping

