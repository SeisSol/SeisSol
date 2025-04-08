// SPDX-FileCopyrightText: 2013 SeisSol Group
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

#include "Kernels/Receiver.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"
#include "TimeCluster.h"

#include <cassert>
#include <cstring>

#include "generated_code/kernel.h"

seissol::time_stepping::TimeCluster::TimeCluster(
    unsigned int clusterId,
    unsigned int globalClusterId,
    unsigned int profilingId,
    bool usePlasticity,
    LayerType layerType,
    double maxTimeStepSize,
    long timeStepRate,
    bool printProgress,
    DynamicRuptureScheduler*
        dynamicRuptureScheduler, // Need only one scheduler and multiple data structures
    CompoundGlobalData globalData,
    seissol::initializer::Layer* clusterData,
    std::array<seissol::initializer::Layer*, seissol::multisim::NumSimulations> dynRupInteriorData,
    std::array<seissol::initializer::Layer*, seissol::multisim::NumSimulations> dynRupCopyData,
    seissol::initializer::LTS* lts,
    std::array<std::shared_ptr<seissol::initializer::DynamicRupture>, seissol::multisim::NumSimulations>
        dynRup,
    std::array<std::shared_ptr<seissol::dr::friction_law::FrictionSolver>, seissol::multisim::NumSimulations>
        frictionSolver,
    std::array<std::shared_ptr<seissol::dr::friction_law::FrictionSolver>, seissol::multisim::NumSimulations>
        frictionSolverDevice,
    std::array<std::shared_ptr<dr::output::OutputManager>, seissol::multisim::NumSimulations>
        faultOutputManager,
    seissol::SeisSol& seissolInstance,
    LoopStatistics* loopStatistics,
    ActorStateStatistics* actorStateStatistics)
    : AbstractTimeCluster(maxTimeStepSize,
                          timeStepRate,
#ifdef ACL_DEVICE
                          i_clusterData->getNumberOfCells() >= deviceHostSwitch() ? Executor::Device
                                                                                  : Executor::Host

      seissolInstance.executionPlace(i_clusterData->getNumberOfCells())
#else
                          Executor::Host
#endif
                          ),
      // cluster ids
      usePlasticity(usePlasticity), seissolInstance(seissolInstance),
      m_globalDataOnHost(globalData.onHost), m_globalDataOnDevice(globalData.onDevice),
      m_clusterData(clusterData),
      // global data
      dynRupInteriorData(dynRupInteriorData), dynRupCopyData(dynRupCopyData), m_lts(lts),
      m_dynRup(dynRup), frictionSolver(frictionSolver),
      frictionSolverDevice(frictionSolverDevice), faultOutputManager(faultOutputManager),
      m_sourceCluster(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
      // cells
      m_loopStatistics(loopStatistics), actorStateStatistics(actorStateStatistics),
      m_receiverCluster(nullptr), layerType(layerType), printProgress(printProgress),
      m_clusterId(clusterId), m_globalClusterId(globalClusterId), m_profilingId(profilingId),
      dynamicRuptureScheduler(dynamicRuptureScheduler) {
  // assert all pointers are valid
  assert(m_clusterData != nullptr);
  assert(m_globalDataOnHost != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(m_globalDataOnDevice != nullptr);
  }

  // set timings to zero
  m_receiverTime                  = 0;

  m_timeKernel.setGlobalData(globalData);
  m_localKernel.setGlobalData(globalData);
  m_localKernel.setInitConds(&seissolInstance.getMemoryManager().getInitialConditions());
  m_localKernel.setGravitationalAcceleration(seissolInstance.getGravitationSetup().acceleration);
  m_neighborKernel.setGlobalData(globalData);
  m_dynamicRuptureKernel.setGlobalData(globalData);
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

#ifndef ACL_DEVICE
void seissol::time_stepping::TimeCluster::computeDynamicRupture(
    std::array<seissol::initializer::Layer*, seissol::multisim::NumSimulations>&
        layerData) {
  bool layerzero = true;

  for (unsigned int i = 0; i < seissol::multisim::NumSimulations; i++) {
    if (layerData[i]->getNumberOfCells() > 0) {
      layerzero = false;
    }
  }

  if (layerzero) {
    return;
  }

  SCOREP_USER_REGION_DEFINE(myRegionHandle)
  SCOREP_USER_REGION_BEGIN(myRegionHandle, "computeDynamicRuptureSpaceTimeInterpolation", SCOREP_USER_REGION_TYPE_COMMON )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  std::array<DRFaceInformation*, seissol::multisim::NumSimulations> faceInformation;
  std::array<DRGodunovData*, seissol::multisim::NumSimulations> godunovData;
  std::array<DREnergyOutput*, seissol::multisim::NumSimulations> drEnergyOutput;
  std::array<real**, seissol::multisim::NumSimulations> timeDerivativePlus;
  std::array<real**, seissol::multisim::NumSimulations> timeDerivativeMinus;
  std::array<real(*)[CONVERGENCE_ORDER][tensor::QInterpolated::size()], seissol::multisim::NumSimulations> qInterpolatedPlus;
  std::array<real(*)[CONVERGENCE_ORDER][tensor::QInterpolated::size()], seissol::multisim::NumSimulations> qInterpolatedMinus;  
  int dQ_DR_Size = 0;
  int dQ_Size = 0;
  for (unsigned int i = 0; i < CONVERGENCE_ORDER; i++) {
    dQ_DR_Size += tensor::dQ_DR::size(i);
    dQ_Size += tensor::dQ::size(i);
  }

  int numSimulationsIfAligned = dQ_Size/dQ_DR_Size;

  m_dynamicRuptureKernel.setTimeStepWidth(timeStepSize());

  for (unsigned int i=0; i < seissol::multisim::NumSimulations; i++){
    faceInformation[i] = layerData[i]->var(m_dynRup[i]->faceInformation);
    godunovData[i] = layerData[i]->var(m_dynRup[i]->godunovData);
    drEnergyOutput[i] = layerData[i]->var(m_dynRup[i]->drEnergyOutput);
    timeDerivativePlus[i] = layerData[i]->var(m_dynRup[i]->timeDerivativePlus); // These are interleaved -> hence need to pick the right values for the spacetime interpolation
    timeDerivativeMinus[i] = layerData[i]->var(m_dynRup[i]->timeDerivativeMinus);// These are interleaved -> hence need to pick the right values for the spacetime interpolation

    qInterpolatedPlus[i] = layerData[i]->var(m_dynRup[i]->qInterpolatedPlus); // This is normal
    qInterpolatedMinus[i] = layerData[i]->var(m_dynRup[i]->qInterpolatedMinus); // This is normal
    frictionSolver[i]->computeDeltaT(m_dynamicRuptureKernel.timePoints);
  }

#pragma omp parallel 
  {
  LIKWID_MARKER_START("computeDynamicRuptureSpaceTimeInterpolation");
  }
  for (unsigned int sim = 0; sim < seissol::multisim::NumSimulations; sim++) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned face = 0; face < layerData[sim]->getNumberOfCells(); ++face) {
      std::vector<real> timeDerivativePlusDR(dQ_DR_Size, 0.0);
      std::vector<real> timeDerivativeMinusDR(dQ_DR_Size, 0.0);
      for(unsigned int j = 0; j < dQ_DR_Size; j++){
        timeDerivativePlusDR[j] = timeDerivativePlus[sim][face][j*numSimulationsIfAligned + sim];
        timeDerivativeMinusDR[j] = timeDerivativeMinus[sim][face][j*numSimulationsIfAligned + sim];
      }
      
      unsigned prefetchFace = (face < layerData[sim]->getNumberOfCells() - 1) ? face + 1 : face;
      m_dynamicRuptureKernel.spaceTimeInterpolation(
          faceInformation[sim][face],
          m_globalDataOnHost,
          &godunovData[sim][face],
          &drEnergyOutput[sim][face],
          #ifdef MULTIPLE_SIMULATIONS
          timeDerivativePlusDR.data(),
          timeDerivativeMinusDR.data(),
          #else
          timeDerivativePlus[sim][face],
          timeDerivativeMinus[sim][face],
          #endif
          qInterpolatedPlus[sim][face], // DR part
          qInterpolatedMinus[sim][face],
          timeDerivativePlus[sim][prefetchFace],
          timeDerivativeMinus[sim][prefetchFace]);
    }
  }
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel 
  {
  LIKWID_MARKER_STOP("computeDynamicRuptureSpaceTimeInterpolation");
  LIKWID_MARKER_START("computeDynamicRuptureFrictionLaw");
  }

  SCOREP_USER_REGION_BEGIN(myRegionHandle, "computeDynamicRuptureFrictionLaw", SCOREP_USER_REGION_TYPE_COMMON )
  for (unsigned int i = 0; i < seissol::multisim::NumSimulations; i++) {
    frictionSolver[i]->evaluate(*layerData[i], m_dynRup[i].get(), ct.correctionTime, m_dynamicRuptureKernel.timeWeights, streamRuntime);
  }
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel 
  {
  LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }

  unsigned int numberofCells = 0;

  for(unsigned int i = 0 ; i < seissol::multisim::NumSimulations; i++){
    numberofCells += layerData[i]->getNumberOfCells();
  }

  m_loopStatistics->end(m_regionComputeDynamicRupture, numberofCells, m_profilingId);
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


void seissol::time_stepping::TimeCluster::computeDynamicRuptureFlops( 
  // seissol::initializer::Layer& layerData,
          std::array<seissol::initializer::Layer*, seissol::multisim::NumSimulations>& layerData,
                                                                      long long&                    nonZeroFlops,
                                                                      long long&                    hardwareFlops )
{
  nonZeroFlops = 0;
  hardwareFlops = 0;

  for(unsigned int i = 0; i < seissol::multisim::NumSimulations; i++)
{  
  DRFaceInformation* faceInformation = layerData[i]->var(m_dynRup[i]->faceInformation);

  for (unsigned face = 0; face < layerData[i]->getNumberOfCells(); ++face) {
    long long faceNonZeroFlops = 0;
    long long faceHardwareFlops = 0;
    m_dynamicRuptureKernel.flopsGodunovState(faceInformation[face], faceNonZeroFlops, faceHardwareFlops);
    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}
}

void seissol::time_stepping::TimeCluster::computeLocalIntegration(seissol::initializer::Layer& layerData, bool resetBuffers ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  real* bufferPointer = nullptr;
  kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);
  // local integration buffer
  alignas(Alignment) real integrationBuffer[tensor::I::size()];

  // pointer for the call of the ADER-function

  real** buffers = layerData.var(m_lts->buffers);
  real** derivatives = layerData.var(m_lts->derivatives);
  CellMaterialData* materialData = layerData.var(m_lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, layerData);


  auto timeStepLocal = timeStepSize();

#ifdef _OPENMP
  #pragma omp parallel for private(bufferPointer, integrationBuffer), firstprivate(tmp) schedule(static)
#endif
  for (unsigned int lCell = 0; lCell < layerData.getNumberOfCells(); lCell++) {

    auto data = loader.entry(lCell);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided = (data.cellInformation().ltsSetup >> 8) % 2 == 1; // buffers are provided
    const bool resetMyBuffers = buffersProvided && ( (data.cellInformation().ltsSetup >> 10) %2 == 0 || resetBuffers ); // they should be reset

    if (resetMyBuffers) {
      // assert presence of the buffer
      assert(buffers[lCell] != nullptr);

      bufferPointer = buffers[lCell];
    } else {
      // work on local buffer
      bufferPointer = integrationBuffer;
    }

    m_timeKernel.computeAder(timeStepLocal,
                             data,
                             tmp, // only for gravitational stuff
                             bufferPointer,
                             derivatives[lCell],
                             true);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping (*boundaryMapping)[4] = layerData.var(m_lts->boundaryMapping);
    m_localKernel.computeIntegral(bufferPointer,
                                  data,
                                  tmp,
                                  &materialData[lCell],
                                  &boundaryMapping[lCell],
                                  ct.correctionTime,
                                  timeStepLocal
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
        addVelocityKrnl.I = bufferPointer;
        addVelocityKrnl.execute(face);
      }
    }

    // TODO: Integrate this step into the kernel
    // We've used a temporary buffer -> need to accumulate update in
    // shared buffer.
    if (!resetMyBuffers && buffersProvided) {
      assert(buffers[lCell] != nullptr);

      for (unsigned int dof = 0; dof < tensor::I::size(); ++dof) {
        buffers[lCell][dof] += integrationBuffer[dof];
      }
    }
  }

  m_loopStatistics->end(m_regionComputeLocalIntegration, layerData.getNumberOfCells(), m_profilingId);
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
        displacementKrnl.integratedVelocities = const_cast<const real**>(
            entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
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
    unsigned cellNonZero = 0;
    unsigned cellHardware = 0;
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
    unsigned cellNonZero = 0;
    unsigned cellHardware = 0;
    long long cellDRNonZero = 0;
    long long cellDRHardware = 0;
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
  computeDynamicRuptureFlops(dynRupInteriorData,
                             m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)],
                             m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
  computeDynamicRuptureFlops(dynRupCopyData,
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
  if (m_clusterData->getNumberOfCells() == 0) {
    return;
  }

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
void TimeCluster::correct() { // Here, need to think what to do
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
  double subTimeStart = ct.correctionTime - lastSubTime; // This is significantly different, lastSubTime -> almost zero in fused, but a huge value in master

  // Note, if this is a copy layer actor, we need the FL_Copy and the FL_Int.
  // Otherwise, this is an interior layer actor, and we need only the FL_Int.
  // We need to avoid computing it twice.
  if (dynamicRuptureScheduler->hasDynamicRuptureFaces()) { // Only one scheduler because everything is concatenated
    if (dynamicRuptureScheduler->mayComputeInterior(ct.stepsSinceStart)) {
     computeDynamicRupture(dynRupInteriorData); // This method deals with fused simulations internally by looping for all
      seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(m_flops_nonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(m_flops_hardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      
      dynamicRuptureScheduler->setLastCorrectionStepsInterior(ct.stepsSinceStart);
    }
    if (layerType == Copy) {
      computeDynamicRupture(dynRupCopyData); // This method deals with fused simulations internally by looping for all
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
  auto timeStepLocal = timeStepSize();
  if (dynamicRuptureScheduler->isFirstClusterWithDynamicRuptureFaces() &&
      dynamicRuptureScheduler->mayComputeFaultOutput(ct.stepsSinceStart)) {
    for (unsigned int i = 0; i < seissol::multisim::NumSimulations; i++) {
      faultOutputManager[i]->writePickpointOutput(ct.correctionTime + timeStepLocal,
                                                  timeStepLocal);
    } /// \todo this method needs modification to incorporate fused simulations
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

std::pair<long, long>
    TimeCluster::computeNeighboringIntegration(seissol::initializer::Layer& layerData, double subTimeStart) {
  if (layerData.getNumberOfCells() == 0) {
    return {0, 0};
}
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  real*(*faceNeighbors)[4] = layerData.var(m_lts->faceNeighbors);
  CellDRMapping(*drMapping)[4] = layerData.var(m_lts->drMapping);
  CellLocalInformation* cellInformation = layerData.var(m_lts->cellInformation);
  auto plasticity = layerData.var(m_lts->plasticity);
  auto* pstrain = layerData.var(m_lts->pstrain);
  unsigned numberOTetsWithPlasticYielding = 0;

  kernels::NeighborData::Loader loader;
  loader.load(*m_lts, layerData);

  real* timeIntegrated[4];
  real* faceNeighborsPrefetch[4];
  auto timeStepLocal = timeStepSize();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) private(timeIntegrated,                   \
                                                                faceNeighborsPrefetch)            \
    shared(cellInformation,                                                                        \
               loader,                                                                             \
               faceNeighbors,                                                                      \
               pstrain,                                                                            \
               layerData,                                                                        \
               plasticity,                                                                         \
               drMapping,                                                                          \
               subTimeStart) reduction(+ : numberOTetsWithPlasticYielding)                         \
               firstprivate(timeStepLocal)
#endif

  for (unsigned int cell = 0; cell < layerData.getNumberOfCells(); cell++) {
    auto data = loader.entry(cell);
    seissol::kernels::TimeCommon::computeIntegrals(
        m_timeKernel,
        data.cellInformation().ltsSetup,
        data.cellInformation().faceTypes,
        subTimeStart,
        timeStepLocal,
        faceNeighbors[cell],
#ifdef _OPENMP
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(
            &(m_globalDataOnHost
                  ->integrationBufferLTS[omp_get_thread_num() * 4 * tensor::I::size()])),
#else
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(m_globalData->integrationBufferLTS),
#endif
        timeIntegrated);

    faceNeighborsPrefetch[0] =
        (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
            ? faceNeighbors[cell][1]
            : drMapping[cell][1].godunov[0]; // the prefetch does not change anything numerical,
                                               // just cache behaviour
    /// \todo think of a cleaner way to actually do the prefetch later
    faceNeighborsPrefetch[1] =
        (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
            ? faceNeighbors[cell][2]
            : drMapping[cell][2].godunov[0]; // the prefetch does not change anything numerical,
                                               // just cache behaviour
    /// \todo think of a cleaner way to actually do the prefetch later
    faceNeighborsPrefetch[2] =
        (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
            ? faceNeighbors[cell][3]
            : drMapping[cell][3].godunov[0]; // the prefetch does not change anything numerical,
                                               // just cache behaviour
    /// \todo think of a cleaner way to actually do the prefetch later

    // fourth face's prefetches
    if (cell < (layerData.getNumberOfCells() - 1)) {
      faceNeighborsPrefetch[3] =
          (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
              ? faceNeighbors[cell + 1][0]
              : drMapping[cell + 1][0].godunov[0]; // the prefetch does not change anything
                                                     // numerical, just cache behaviour
      /// \todo think of a cleaner way to actually do the prefetch later
    } else {
      faceNeighborsPrefetch[3] = faceNeighbors[cell][3];
    }

    m_neighborKernel.computeNeighborsIntegral(
        data,
        drMapping[cell],
        timeIntegrated,
        faceNeighborsPrefetch);

    if (usePlasticity) {
      updateRelaxTime();
      numberOTetsWithPlasticYielding +=
          seissol::kernels::Plasticity::computePlasticity(m_oneMinusIntegratingFactor,
                                                          timeStepLocal,
                                                          m_tv,
                                                          m_globalDataOnHost,
                                                          &plasticity[cell],
                                                          data.dofs(),
                                                          pstrain[cell]); // Plasticity needs to be managed now
    }
#ifdef INTEGRATE_QUANTITIES
    seissolInstance.postProcessor().integrateQuantities(
        m_timeStepWidth, layerData, cell, dofs[cell]);
#endif // INTEGRATE_QUANTITIES
  }

  const long long nonZeroFlopsPlasticity =
      layerData.getNumberOfCells() *
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)] +
      numberOTetsWithPlasticYielding *
          m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)];
  const long long hardwareFlopsPlasticity =
      layerData.getNumberOfCells() *
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)] +
      numberOTetsWithPlasticYielding *
          m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)];

  m_loopStatistics->end(
      m_regionComputeNeighboringIntegration, layerData.getNumberOfCells(), m_profilingId);

  return {nonZeroFlopsPlasticity, hardwareFlopsPlasticity};
}
#endif // ACL_DEVICE

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

