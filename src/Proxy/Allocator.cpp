// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Allocator.h"
#include <Common/Constants.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Kernels/Touch.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/GlobalData.h>
#include <Memory/MemoryAllocator.h>
#include <Memory/Tree/Layer.h>
#include <Memory/Tree/TimeCluster.h>
#include <Proxy/Constants.h>
#include <cstddef>
#include <stdlib.h>
#include <tensor.h>

#ifdef ACL_DEVICE
#include <Initializer/MemoryManager.h>
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
      static_cast<real*>(layer.bucket(lts.buffersDerivatives, initializer::AllocationPlace::Host));

  real** buffersDevice = layer.var(lts.buffersDevice);
  real** derivativesDevice = layer.var(lts.derivativesDevice);
  real*(*faceNeighborsDevice)[4] = layer.var(lts.faceNeighborsDevice);
  real* bucketDevice = static_cast<real*>(
      layer.bucket(lts.buffersDerivatives, initializer::AllocationPlace::Device));

  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    buffers[cell] = bucket + cell * tensor::I::size();
    derivatives[cell] = nullptr;
    buffersDevice[cell] = bucketDevice + cell * tensor::I::size();
    derivativesDevice[cell] = nullptr;

    for (unsigned f = 0; f < 4; ++f) {
      cellInformation[cell].faceTypes[f] = faceTp;
      cellInformation[cell].faceRelations[f][0] = ((unsigned int)lrand48() % 4);
      cellInformation[cell].faceRelations[f][1] = ((unsigned int)lrand48() % 3);
      secondaryInformation[cell].faceNeighborIds[f] =
          ((unsigned int)lrand48() % layer.getNumberOfCells());
    }
    cellInformation[cell].ltsSetup = 0;
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
    for (unsigned f = 0; f < 4; ++f) {
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

  kernels::fillWithStuff(
      reinterpret_cast<real*>(dofs), tensor::Q::size() * layer.getNumberOfCells(), false);
  kernels::fillWithStuff(bucket, tensor::I::size() * layer.getNumberOfCells(), false);
  kernels::fillWithStuff(reinterpret_cast<real*>(localIntegration),
                         sizeof(LocalIntegrationData) / sizeof(real) * layer.getNumberOfCells(),
                         false);
  kernels::fillWithStuff(reinterpret_cast<real*>(neighboringIntegration),
                         sizeof(NeighboringIntegrationData) / sizeof(real) *
                             layer.getNumberOfCells(),
                         false);

#ifdef USE_POROELASTIC
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned cell = 0; cell < layer.getNumberOfCells(); ++cell) {
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
  timeKernel.setGlobalData(globalData);
  localKernel.setGlobalData(globalData);
  neighborKernel.setGlobalData(globalData);
  dynRupKernel.setGlobalData(globalData);

  dynRupKernel.setTimeStepWidth(Timestep);
}

void ProxyData::initDataStructures(bool enableDR) {
  // init RNG
  srand48(cellCount);
  lts.addTo(ltsTree, false); // proxy does not use plasticity
  ltsTree.setNumberOfTimeClusters(1);
  ltsTree.fixate();

  seissol::initializer::TimeCluster& cluster = ltsTree.child(0);
  cluster.child<Ghost>().setNumberOfCells(0);
  cluster.child<Copy>().setNumberOfCells(0);
  cluster.child<Interior>().setNumberOfCells(cellCount);

  seissol::initializer::Layer& layer = cluster.child<Interior>();
  layer.setBucketSize(lts.buffersDerivatives,
                      sizeof(real) * tensor::I::size() * layer.getNumberOfCells());

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

    fakeDerivativesHost = reinterpret_cast<real*>(
        allocator.allocateMemory(cellCount * yateto::computeFamilySize<tensor::dQ>() * sizeof(real),
                                 PagesizeHeap,
                                 seissol::memory::Standard));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned cell = 0; cell < cellCount; ++cell) {
      for (unsigned i = 0; i < yateto::computeFamilySize<tensor::dQ>(); i++) {
        fakeDerivativesHost[cell * yateto::computeFamilySize<tensor::dQ>() + i] = (real)drand48();
      }
    }

#ifdef ACL_DEVICE
    fakeDerivatives = reinterpret_cast<real*>(
        allocator.allocateMemory(cellCount * yateto::computeFamilySize<tensor::dQ>() * sizeof(real),
                                 PagesizeHeap,
                                 seissol::memory::DeviceGlobalMemory));
    const auto& device = ::device::DeviceInstance::getInstance();
    device.api->copyTo(fakeDerivatives,
                       fakeDerivativesHost,
                       cellCount * yateto::computeFamilySize<tensor::dQ>() * sizeof(real));
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

    /* init drMapping */
    for (unsigned cell = 0; cell < cellCount; ++cell) {
      for (unsigned face = 0; face < 4; ++face) {
        CellDRMapping& drm = drMapping[cell][face];
        const unsigned side = (unsigned int)lrand48() % 4;
        const unsigned orientation = (unsigned int)lrand48() % 3;
        const unsigned drFace = (unsigned int)lrand48() % interior.getNumberOfCells();
        drm.side = side;
        drm.faceRelation = orientation;
        drm.godunov = imposedStatePlus[drFace];
        drm.fluxSolver = fluxSolverPlus[drFace];
      }
    }

    /* init dr godunov state */
    for (unsigned face = 0; face < interior.getNumberOfCells(); ++face) {
      const unsigned plusCell = (unsigned int)lrand48() % cellCount;
      const unsigned minusCell = (unsigned int)lrand48() % cellCount;
      timeDerivativeHostPlus[face] =
          &fakeDerivativesHost[plusCell * yateto::computeFamilySize<tensor::dQ>()];
      timeDerivativeHostMinus[face] =
          &fakeDerivativesHost[minusCell * yateto::computeFamilySize<tensor::dQ>()];
      timeDerivativePlus[face] =
          &fakeDerivatives[plusCell * yateto::computeFamilySize<tensor::dQ>()];
      timeDerivativeMinus[face] =
          &fakeDerivatives[minusCell * yateto::computeFamilySize<tensor::dQ>()];

      faceInformation[face].plusSide = (unsigned int)lrand48() % 4;
      faceInformation[face].minusSide = (unsigned int)lrand48() % 4;
      faceInformation[face].faceRelation = (unsigned int)lrand48() % 3;
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

  seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(ltsTree, lts);
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
