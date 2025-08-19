// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Allocator.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Kernels/Touch.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/GlobalData.h>
#include <Memory/MemoryAllocator.h>
#include <Memory/Tree/Layer.h>
#include <Memory/Tree/TimeCluster.h>
#include <cstddef>
#include <random>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef ACL_DEVICE
#include <Initializer/MemoryManager.h>
#endif

#ifdef USE_POROELASTIC
#include "Proxy/Constants.h"
#endif

namespace {
void fakeData(initializer::LTS& lts, initializer::Layer& layer, FaceType faceTp) {
  real(*dofs)[tensor::Q::size()] = layer.var(lts.dofs);
  real** buffers = layer.var(lts.buffers);
  real** derivatives = layer.var(lts.derivatives);
  real*(*faceNeighbors)[4] = layer.var(lts.faceNeighbors);
  auto* localIntegration = layer.var(lts.localIntegration);
  auto* neighboringIntegration = layer.var(lts.neighboringIntegration);
  auto* cellInformation = layer.var(lts.cellInformation);
  auto* secondaryInformation = layer.var(lts.secondaryInformation);
  real* bucket =
      static_cast<real*>(layer.var(lts.buffersDerivatives, initializer::AllocationPlace::Host));

  real** buffersDevice = layer.var(lts.buffersDevice);
  real** derivativesDevice = layer.var(lts.derivativesDevice);
  real*(*faceNeighborsDevice)[4] = layer.var(lts.faceNeighborsDevice);
  real* bucketDevice =
      static_cast<real*>(layer.var(lts.buffersDerivatives, initializer::AllocationPlace::Device));

  std::mt19937 rng(layer.size());
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<unsigned> orientationDist(0, 2);
  std::uniform_int_distribution<std::size_t> cellDist(0, layer.size() - 1);

  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    buffers[cell] = bucket + cell * tensor::I::size();
    derivatives[cell] = nullptr;
    buffersDevice[cell] = bucketDevice + cell * tensor::I::size();
    derivativesDevice[cell] = nullptr;

    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      cellInformation[cell].faceTypes[f] = faceTp;
      cellInformation[cell].faceRelations[f][0] = sideDist(rng);
      cellInformation[cell].faceRelations[f][1] = orientationDist(rng);
      secondaryInformation[cell].faceNeighborIds[f] = cellDist(rng);
    }
    cellInformation[cell].ltsSetup = 0;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      switch (faceTp) {
      case FaceType::FreeSurface:
        faceNeighbors[cell][f] = buffers[cell];
        faceNeighborsDevice[cell][f] = buffersDevice[cell];
        break;
      case FaceType::Periodic:
      case FaceType::Regular:
        faceNeighbors[cell][f] = buffers[secondaryInformation[cell].faceNeighborIds[f]];
        faceNeighborsDevice[cell][f] = buffersDevice[secondaryInformation[cell].faceNeighborIds[f]];
        break;
      default:
        faceNeighbors[cell][f] = nullptr;
        break;
      }
    }
  }

  kernels::fillWithStuff(reinterpret_cast<real*>(dofs), tensor::Q::size() * layer.size(), false);
  kernels::fillWithStuff(bucket, tensor::I::size() * layer.size(), false);
  kernels::fillWithStuff(reinterpret_cast<real*>(localIntegration),
                         sizeof(LocalIntegrationData) / sizeof(real) * layer.size(),
                         false);
  kernels::fillWithStuff(reinterpret_cast<real*>(neighboringIntegration),
                         sizeof(NeighboringIntegrationData) / sizeof(real) * layer.size(),
                         false);

#ifdef USE_POROELASTIC
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    localIntegration[cell].specific.typicalTimeStepWidth = seissol::proxy::Timestep;
  }
#endif

#ifdef ACL_DEVICE
  const auto& device = device::DeviceInstance::getInstance();
  layer.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                      device.api->getDefaultStream());
  device.api->syncDefaultStreamWithHost();
#endif
}
} // namespace

namespace seissol::proxy {

ProxyData::ProxyData(std::size_t cellCount, bool enableDR) : cellCount(cellCount) {
  initGlobalData();
  initDataStructures(enableDR);
  initDataStructuresOnDevice(enableDR);
}

void ProxyData::initGlobalData() {
  seissol::initializer::GlobalDataInitializerOnHost::init(
      globalDataOnHost, allocator, seissol::memory::Standard);

  CompoundGlobalData globalData{};
  globalData.onHost = &globalDataOnHost;
  globalData.onDevice = nullptr;
  if constexpr (seissol::isDeviceOn()) {
    seissol::initializer::GlobalDataInitializerOnDevice::init(
        globalDataOnDevice, allocator, seissol::memory::DeviceGlobalMemory);
    globalData.onDevice = &globalDataOnDevice;
  }
  spacetimeKernel.setGlobalData(globalData);
  timeKernel.setGlobalData(globalData);
  localKernel.setGlobalData(globalData);
  neighborKernel.setGlobalData(globalData);
  dynRupKernel.setGlobalData(globalData);
}

void ProxyData::initDataStructures(bool enableDR) {
  // init RNG
  lts.addTo(ltsTree, false); // proxy does not use plasticity
  ltsTree.setNumberOfTimeClusters(1);
  ltsTree.fixate();

  seissol::initializer::TimeCluster& cluster = ltsTree.child(0);
  cluster.child<Ghost>().setNumberOfCells(0);
  cluster.child<Copy>().setNumberOfCells(0);
  cluster.child<Interior>().setNumberOfCells(cellCount);

  seissol::initializer::Layer& layer = cluster.child<Interior>();
  layer.setEntrySize(lts.buffersDerivatives, sizeof(real) * tensor::I::size() * layer.size());

  ltsTree.allocateVariables();
  ltsTree.touchVariables();
  ltsTree.allocateBuckets();

  if (enableDR) {
    dynRup.addTo(dynRupTree);
    dynRupTree.setNumberOfTimeClusters(1);
    dynRupTree.fixate();

    seissol::initializer::TimeCluster& cluster = dynRupTree.child(0);
    cluster.child<Ghost>().setNumberOfCells(0);
    cluster.child<Copy>().setNumberOfCells(0);
    cluster.child<Interior>().setNumberOfCells(
        4 * cellCount); /// Every face is a potential dynamic rupture face

    dynRupTree.allocateVariables();
    dynRupTree.touchVariables();

    fakeDerivativesHost = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::Solver::DerivativesSize * sizeof(real),
        PagesizeHeap,
        seissol::memory::Standard));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
      const auto offset = omp_get_thread_num();
#else
      const auto offset = 0;
#endif
      std::mt19937 rng(cellCount + offset);
      std::uniform_real_distribution<real> urd;
      for (std::size_t cell = 0; cell < cellCount; ++cell) {
        for (std::size_t i = 0; i < seissol::kernels::Solver::DerivativesSize; i++) {
          fakeDerivativesHost[cell * seissol::kernels::Solver::DerivativesSize + i] = urd(rng);
        }
      }
    }

#ifdef ACL_DEVICE
    fakeDerivatives = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::Solver::DerivativesSize * sizeof(real),
        PagesizeHeap,
        seissol::memory::DeviceGlobalMemory));
    const auto& device = ::device::DeviceInstance::getInstance();
    device.api->copyTo(fakeDerivatives,
                       fakeDerivativesHost,
                       cellCount * seissol::kernels::Solver::DerivativesSize * sizeof(real));
#else
    fakeDerivatives = fakeDerivativesHost;
#endif
  }

  /* cell information and integration data*/
  fakeData(lts, layer, (enableDR) ? FaceType::DynamicRupture : FaceType::Regular);

  if (enableDR) {
    // From lts tree
    CellDRMapping(*drMapping)[4] =
        isDeviceOn() ? ltsTree.var(lts.drMappingDevice) : ltsTree.var(lts.drMapping);

    constexpr initializer::AllocationPlace Place =
        isDeviceOn() ? initializer::AllocationPlace::Device : initializer::AllocationPlace::Host;

    // From dynamic rupture tree
    seissol::initializer::Layer& interior = dynRupTree.child(0).child<Interior>();
    real(*imposedStatePlus)[seissol::tensor::QInterpolated::size()] =
        interior.var(dynRup.imposedStatePlus, Place);
    real(*fluxSolverPlus)[seissol::tensor::fluxSolver::size()] =
        interior.var(dynRup.fluxSolverPlus, Place);
    real** timeDerivativeHostPlus = interior.var(dynRup.timeDerivativePlus);
    real** timeDerivativeHostMinus = interior.var(dynRup.timeDerivativeMinus);
    real** timeDerivativePlus = isDeviceOn() ? interior.var(dynRup.timeDerivativePlusDevice)
                                             : interior.var(dynRup.timeDerivativePlus);
    real** timeDerivativeMinus = isDeviceOn() ? interior.var(dynRup.timeDerivativeMinusDevice)
                                              : interior.var(dynRup.timeDerivativeMinus);
    DRFaceInformation* faceInformation = interior.var(dynRup.faceInformation);

    std::mt19937 rng(cellCount);
    std::uniform_int_distribution<unsigned> sideDist(0, 3);
    std::uniform_int_distribution<unsigned> orientationDist(0, 2);
    std::uniform_int_distribution<std::size_t> drDist(0, interior.size() - 1);
    std::uniform_int_distribution<std::size_t> cellDist(0, cellCount - 1);

    /* init drMapping */
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        CellDRMapping& drm = drMapping[cell][face];
        const auto side = sideDist(rng);
        const auto orientation = orientationDist(rng);
        const auto drFace = drDist(rng);
        drm.side = side;
        drm.faceRelation = orientation;
        drm.godunov = imposedStatePlus[drFace];
        drm.fluxSolver = fluxSolverPlus[drFace];
      }
    }

    /* init dr godunov state */
    for (std::size_t face = 0; face < interior.size(); ++face) {
      const auto plusCell = cellDist(rng);
      const auto minusCell = cellDist(rng);
      timeDerivativeHostPlus[face] =
          &fakeDerivativesHost[plusCell * seissol::kernels::Solver::DerivativesSize];
      timeDerivativeHostMinus[face] =
          &fakeDerivativesHost[minusCell * seissol::kernels::Solver::DerivativesSize];
      timeDerivativePlus[face] =
          &fakeDerivatives[plusCell * seissol::kernels::Solver::DerivativesSize];
      timeDerivativeMinus[face] =
          &fakeDerivatives[minusCell * seissol::kernels::Solver::DerivativesSize];

      faceInformation[face].plusSide = sideDist(rng);
      faceInformation[face].minusSide = sideDist(rng);
      faceInformation[face].faceRelation = orientationDist(rng);
    }
  }
}

void ProxyData::initDataStructuresOnDevice(bool enableDR) {
#ifdef ACL_DEVICE
  const auto& device = ::device::DeviceInstance::getInstance();
  ltsTree.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                        device.api->getDefaultStream());
  device.api->syncDefaultStreamWithHost();

  seissol::initializer::TimeCluster& cluster = ltsTree.child(0);
  seissol::initializer::Layer& layer = cluster.child<Interior>();

  seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(false, ltsTree, lts);
  ltsTree.allocateScratchPads();

  seissol::initializer::recording::CompositeRecorder<seissol::initializer::LTS> recorder;
  recorder.addRecorder(new seissol::initializer::recording::LocalIntegrationRecorder);
  recorder.addRecorder(new seissol::initializer::recording::NeighIntegrationRecorder);

  recorder.addRecorder(new seissol::initializer::recording::PlasticityRecorder);
  recorder.record(lts, layer);
  if (enableDR) {
    dynRupTree.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                             device.api->getDefaultStream());
    device.api->syncDefaultStreamWithHost();
    seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForDr(dynRupTree, dynRup);
    dynRupTree.allocateScratchPads();

    CompositeRecorder<seissol::initializer::DynamicRupture> drRecorder;
    drRecorder.addRecorder(new DynamicRuptureRecorder);

    auto& drLayer = dynRupTree.child(0).child<Interior>();
    drRecorder.record(dynRup, drLayer);
  }
#endif // ACL_DEVICE
}

} // namespace seissol::proxy
