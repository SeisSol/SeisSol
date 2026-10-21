// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Allocator.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Kernels/Touch.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/GlobalData.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Colormap.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/OpenMP.h"
#include "Solver/Settings.h"

#include <algorithm>
#include <cstddef>
#include <random>
#include <stdlib.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/Recorders.h"
#include "Initializer/InitProcedure/Internal/Scratchpads.h"

#include <Device/device.h>
#endif

#ifdef SEISSOL_KERNELS_STP
#include "Proxy/Constants.h"
#endif

namespace seissol::proxy {

namespace {

void fakeData(LTS::Layer& layer, FaceType faceTp) {
  real(*dofs)[tensor::Q::size()] = layer.var<LTS::Dofs>();
  auto* faceTraces = layer.var<LTS::FaceTraces>();
  auto* faceBuffers = layer.var<LTS::FaceBuffers>();
  auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  auto* localIntegration = layer.var<LTS::LocalIntegration>();
  auto* neighboringIntegration = layer.var<LTS::NeighboringIntegration>();
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
  real* bucket = static_cast<real*>(layer.var<LTS::Buffers>(initializer::AllocationPlace::Host));

  auto* faceTracesDevice = layer.var<LTS::FaceTracesDevice>();
  auto* faceBuffersDevice = layer.var<LTS::FaceBuffersDevice>();
  auto* faceNeighborsDevice = layer.var<LTS::FaceNeighborsDevice>();
  real* bucketDevice =
      static_cast<real*>(layer.var<LTS::Buffers>(initializer::AllocationPlace::Device));

  std::mt19937 rng(layer.size());
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<unsigned> traceSideDist(0, 1);
  std::uniform_int_distribution<std::size_t> cellDist(0, layer.size() - 1);

  // every cell runs at the same time step and writes the trace of its step on every face, which
  // its neighbor reads from there
  constexpr auto TraceSize = kernels::Solver::FaceTraceSize;
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    cellInformation[cell].ltsSetup = LtsSetup();
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      faceTraces[cell][f] = bucket + (cell * Cell::NumFaces + f) * TraceSize;
      faceBuffers[cell][f] = faceTraces[cell][f];
      faceTracesDevice[cell][f] = bucketDevice + (cell * Cell::NumFaces + f) * TraceSize;
      faceBuffersDevice[cell][f] = faceTracesDevice[cell][f];

      cellInformation[cell].faceTypes[f] = faceTp;
      cellInformation[cell].faceRelations[f][0] = sideDist(rng);
      cellInformation[cell].faceRelations[f][1] = 0;
      cellInformation[cell].traceSides[f] = traceSideDist(rng);
      cellInformation[cell].ltsSetup.setNeighborBuffer(f, BufferType::StepIntegrals);
      cellInformation[cell].ltsSetup.setNeighborGTSRelation(f, true);
      cellInformation[cell].ltsSetup.setFaceBuffer(f, BufferType::StepIntegrals);

      const auto neighbor = cellDist(rng);
      secondaryInformation[cell].faceNeighbors[f].global = neighbor;
      secondaryInformation[cell].faceNeighbors[f].color = 0;
      secondaryInformation[cell].faceNeighbors[f].cell = neighbor;
    }
  }

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      const auto neighbor = secondaryInformation[cell].faceNeighbors[f].cell;
      const auto neighborFace = cellInformation[cell].faceRelations[f][0];
      switch (faceTp) {
      case FaceType::FreeSurface:
        faceNeighbors[cell][f] = faceBuffers[cell][f];
        faceNeighborsDevice[cell][f] = faceBuffersDevice[cell][f];
        break;
      case FaceType::Regular:
        faceNeighbors[cell][f] = faceBuffers[neighbor][neighborFace];
        faceNeighborsDevice[cell][f] = faceBuffersDevice[neighbor][neighborFace];
        break;
      default:
        faceNeighbors[cell][f] = nullptr;
        break;
      }
    }
  }

  kernels::fillWithStuff(reinterpret_cast<real*>(dofs), tensor::Q::size() * layer.size(), false);
  kernels::fillWithStuff(bucket, Cell::NumFaces * TraceSize * layer.size(), false);
  std::fill_n(layer.var<LTS::ZeroFaceTrace>(), TraceSize, real{0});
  kernels::fillWithStuff(reinterpret_cast<real*>(localIntegration),
                         sizeof(LocalIntegrationData) / sizeof(real) * layer.size(),
                         false);
  kernels::fillWithStuff(reinterpret_cast<real*>(neighboringIntegration),
                         sizeof(NeighboringIntegrationData) / sizeof(real) * layer.size(),
                         false);

#ifdef SEISSOL_KERNELS_STP

#pragma omp parallel for schedule(static)
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

ProxyData::ProxyData(std::size_t cellCount, bool enableDR) : cellCount(cellCount) {
  layerId =
      initializer::LayerIdentifier(HaloType::Interior, initializer::ConfigVariant{Config()}, 0);

  initGlobalData();
  initDataStructures(enableDR);
  initDataStructuresOnDevice(enableDR);
}

void ProxyData::initGlobalData() {
  seissol::initializer::GlobalDataInitializerOnHost::init(
      globalDataOnHost, allocator, seissol::memory::Memkind::Standard);

  CompoundGlobalData globalData{};
  globalData.onHost = &globalDataOnHost;
  globalData.onDevice = nullptr;
  if constexpr (seissol::isDeviceOn()) {
    seissol::initializer::GlobalDataInitializerOnDevice::init(
        globalDataOnDevice, allocator, seissol::memory::Memkind::DeviceGlobalMemory);
    globalData.onDevice = &globalDataOnDevice;
  }
  spacetimeKernel.setGlobalData(globalData);
  timeKernel.setGlobalData(globalData);
  localKernel.setGlobalData(globalData);
  neighborKernel.setGlobalData(globalData);
  dynRupKernel.setGlobalData(globalData);
}

void ProxyData::initDataStructures(bool enableDR) {
  const initializer::LTSColorMap map(
      initializer::EnumLayer<HaloType>({HaloType::Interior}),
      initializer::EnumLayer<std::size_t>({0}),
      initializer::TraitLayer<initializer::ConfigVariant>({initializer::ConfigVariant(Config())}));

  // init RNG
  const auto nullSettings = SimulationSettings(false, false, false, false);
  LTS::addTo(ltsStorage, nullSettings);
  ltsStorage.setLayerCount(map);
  ltsStorage.fixate();

  ltsStorage.layer(layerId).setNumberOfCells(cellCount);

  LTS::Layer& layer = ltsStorage.layer(layerId);
  layer.setEntrySize<LTS::Buffers>(sizeof(real) * Cell::NumFaces * kernels::Solver::FaceTraceSize *
                                   layer.size());
  layer.setEntrySize<LTS::ZeroFaceTrace>(sizeof(real) * kernels::Solver::FaceTraceSize);

  ltsStorage.allocateVariables();
  ltsStorage.touchVariables();
  ltsStorage.allocateBuckets();

  integrals = reinterpret_cast<real*>(
      allocator.allocateMemory(cellCount * seissol::kernels::Solver::IntegralsSize * sizeof(real),
                               PagesizeHeap,
                               seissol::memory::Memkind::Standard));
  kernels::fillWithStuff(integrals, seissol::kernels::Solver::IntegralsSize * cellCount, false);

  if (enableDR) {
    DynamicRupture dynRup;
    dynRup.addTo(drStorage);
    drStorage.setLayerCount(ltsStorage.getColorMap());
    drStorage.fixate();

    drStorage.layer(layerId).setNumberOfCells(4 * cellCount);

    drStorage.allocateVariables();
    drStorage.touchVariables();

    fakeDerivativesHost = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::Solver::FaceDerivativesSize * sizeof(real),
        PagesizeHeap,
        seissol::memory::Memkind::Standard));

#pragma omp parallel
    {
      const auto offset = OpenMP::threadId();
      std::mt19937 rng(cellCount + offset);
      std::uniform_real_distribution<real> urd;
      for (std::size_t cell = 0; cell < cellCount; ++cell) {
        for (std::size_t i = 0; i < seissol::kernels::Solver::FaceDerivativesSize; i++) {
          fakeDerivativesHost[cell * seissol::kernels::Solver::FaceDerivativesSize + i] = urd(rng);
        }
      }
    }

#ifdef ACL_DEVICE
    fakeDerivatives = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::Solver::FaceDerivativesSize * sizeof(real),
        PagesizeHeap,
        seissol::memory::Memkind::DeviceGlobalMemory));
    const auto& device = ::device::DeviceInstance::getInstance();
    device.api->copyTo(fakeDerivatives,
                       fakeDerivativesHost,
                       cellCount * seissol::kernels::Solver::FaceDerivativesSize * sizeof(real));
#else
    fakeDerivatives = fakeDerivativesHost;
#endif
  }

  /* cell information and integration data*/
  fakeData(layer, enableDR ? FaceType::DynamicRupture : FaceType::Regular);

  if (enableDR) {
    // From lts storage
    auto* drMapping =
        isDeviceOn() ? ltsStorage.var<LTS::DRMappingDevice>() : ltsStorage.var<LTS::DRMapping>();

    constexpr initializer::AllocationPlace Place =
        isDeviceOn() ? initializer::AllocationPlace::Device : initializer::AllocationPlace::Host;

    // From dynamic rupture storage
    auto& interior = drStorage.layer(layerId);
    real(*imposedStatePlus)[seissol::tensor::QInterpolated::size()] =
        interior.var<DynamicRupture::ImposedStatePlus>(Place);
    auto* fluxSolverPlus = interior.var<DynamicRupture::FluxSolverPlus>(Place);
    real** timeDerivativeHostPlus = interior.var<DynamicRupture::TimeDerivativePlus>();
    real** timeDerivativeHostMinus = interior.var<DynamicRupture::TimeDerivativeMinus>();
    real** timeDerivativePlus = isDeviceOn()
                                    ? interior.var<DynamicRupture::TimeDerivativePlusDevice>()
                                    : interior.var<DynamicRupture::TimeDerivativePlus>();
    real** timeDerivativeMinus = isDeviceOn()
                                     ? interior.var<DynamicRupture::TimeDerivativeMinusDevice>()
                                     : interior.var<DynamicRupture::TimeDerivativeMinus>();
    DRFaceInformation* faceInformation = interior.var<DynamicRupture::FaceInformation>();

    std::mt19937 rng(cellCount);
    std::uniform_int_distribution<unsigned> sideDist(0, 3);
    std::uniform_int_distribution<unsigned> orientationDist(0, 1);
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
          &fakeDerivativesHost[plusCell * seissol::kernels::Solver::FaceDerivativesSize];
      timeDerivativeHostMinus[face] =
          &fakeDerivativesHost[minusCell * seissol::kernels::Solver::FaceDerivativesSize];
      timeDerivativePlus[face] =
          &fakeDerivatives[plusCell * seissol::kernels::Solver::FaceDerivativesSize];
      timeDerivativeMinus[face] =
          &fakeDerivatives[minusCell * seissol::kernels::Solver::FaceDerivativesSize];

      faceInformation[face].plusSide = sideDist(rng);
      faceInformation[face].minusSide = sideDist(rng);
      // a fault face always addresses the minus side here
      faceInformation[face].faceRelation = 1;
    }
  }
}

void ProxyData::initDataStructuresOnDevice(bool enableDR) {
#ifdef ACL_DEVICE
  const auto& device = ::device::DeviceInstance::getInstance();
  ltsStorage.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                           device.api->getDefaultStream());
  device.api->syncDefaultStreamWithHost();

  auto& layer = ltsStorage.layer(layerId);

  seissol::initializer::internal::deriveRequiredScratchpadMemoryForWp(false, ltsStorage);
  ltsStorage.allocateScratchPads();

  seissol::recording::CompositeRecorder<LTS::LTSVarmap> recorder;
  recorder.addRecorder(new seissol::recording::LocalIntegrationRecorder(9.81));
  recorder.addRecorder(new seissol::recording::NeighIntegrationRecorder);

  recorder.addRecorder(new seissol::recording::PlasticityRecorder);
  recorder.record(layer);
  if (enableDR) {
    drStorage.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                            device.api->getDefaultStream());
    device.api->syncDefaultStreamWithHost();
    seissol::initializer::internal::deriveRequiredScratchpadMemoryForDr(drStorage);
    drStorage.allocateScratchPads();

    seissol::recording::CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
    drRecorder.addRecorder(new seissol::recording::DynamicRuptureRecorder);

    auto& drLayer = drStorage.layer(layerId);
    drRecorder.record(drLayer);
  }
#endif // ACL_DEVICE
}

} // namespace seissol::proxy
