// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Allocator.h"
#include "Parallel/OpenMP.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <Config.h>
#include <Equations/Datastructures.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Kernels/Touch.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/GlobalData.h>
#include <Memory/MemoryAllocator.h>
#include <Memory/Tree/Colormap.h>
#include <Memory/Tree/Layer.h>
#include <Model/CommonDatastructures.h>
#include <Proxy/Constants.h>
#include <cstddef>
#include <memory>
#include <random>
#include <stdlib.h>
#include <variant>

#ifdef ACL_DEVICE
#include <Initializer/MemoryManager.h>
#endif

namespace {
template <typename Cfg>
void fakeData(Cfg cfg, LTS::Layer& layer, FaceType faceTp) {
  using real = Real<Cfg>;

  real(*dofs)[tensor::Q<Cfg>::size()] = layer.var<LTS::Dofs>(cfg);
  real** buffers = layer.var<LTS::Buffers>(cfg);
  real** derivatives = layer.var<LTS::Derivatives>(cfg);
  void*(*faceNeighbors)[4] = layer.var<LTS::FaceNeighbors>();
  auto* localIntegration = layer.var<LTS::LocalIntegration>(cfg);
  auto* neighboringIntegration = layer.var<LTS::NeighboringIntegration>(cfg);
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
  real* bucket = static_cast<real*>(
      layer.var<LTS::BuffersDerivatives>(cfg, initializer::AllocationPlace::Host));

  real** buffersDevice = layer.var<LTS::BuffersDevice>(cfg);
  real** derivativesDevice = layer.var<LTS::DerivativesDevice>(cfg);
  void*(*faceNeighborsDevice)[4] = layer.var<LTS::FaceNeighborsDevice>();
  real* bucketDevice = static_cast<real*>(
      layer.var<LTS::BuffersDerivatives>(cfg, initializer::AllocationPlace::Device));

  std::mt19937 rng(layer.size());
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<unsigned> orientationDist(0, 2);
  std::uniform_int_distribution<std::size_t> cellDist(0, layer.size() - 1);

  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    buffers[cell] = bucket + cell * tensor::I<Cfg>::size();
    derivatives[cell] = nullptr;
    buffersDevice[cell] = bucketDevice + cell * tensor::I<Cfg>::size();
    derivativesDevice[cell] = nullptr;

    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      cellInformation[cell].faceTypes[f] = faceTp;
      cellInformation[cell].faceRelations[f][0] = sideDist(rng);
      cellInformation[cell].faceRelations[f][1] = orientationDist(rng);

      const auto neighbor = cellDist(rng);
      secondaryInformation[cell].faceNeighbors[f].global = neighbor;
      secondaryInformation[cell].faceNeighbors[f].color = 0;
      secondaryInformation[cell].faceNeighbors[f].cell = neighbor;
    }
    cellInformation[cell].ltsSetup = LtsSetup();
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
        faceNeighbors[cell][f] = buffers[secondaryInformation[cell].faceNeighbors[f].cell];
        faceNeighborsDevice[cell][f] =
            buffersDevice[secondaryInformation[cell].faceNeighbors[f].cell];
        break;
      default:
        faceNeighbors[cell][f] = nullptr;
        break;
      }
    }
  }

  kernels::fillWithStuff(
      reinterpret_cast<real*>(dofs), tensor::Q<Cfg>::size() * layer.size(), false);
  kernels::fillWithStuff(bucket, tensor::I<Cfg>::size() * layer.size(), false);
  kernels::fillWithStuff(reinterpret_cast<real*>(localIntegration),
                         sizeof(LocalIntegrationData<Cfg>) / sizeof(real) * layer.size(),
                         false);
  kernels::fillWithStuff(reinterpret_cast<real*>(neighboringIntegration),
                         sizeof(NeighboringIntegrationData<Cfg>) / sizeof(real) * layer.size(),
                         false);

  if constexpr (model::MaterialTT<Cfg>::Type == model::MaterialType::Poroelastic) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      localIntegration[cell].specific.typicalTimeStepWidth = seissol::proxy::Timestep;
    }
  }

#ifdef ACL_DEVICE
  const auto& device = device::DeviceInstance::getInstance();
  layer.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                      device.api->getDefaultStream());
  device.api->syncDefaultStreamWithHost();
#endif
}
} // namespace

namespace seissol::proxy {

template <typename Cfg>
ProxyDataImpl<Cfg>::ProxyDataImpl(std::size_t cellCount, bool enableDR) : cellCount(cellCount) {
  initGlobalData();
  initDataStructures(enableDR);
  initDataStructuresOnDevice(enableDR);
}

template <typename Cfg>
void ProxyDataImpl<Cfg>::initGlobalData() {
  globalData.init(ConfigVariant(Cfg()).index());

  spacetimeKernel.setGlobalData(globalData);
  timeKernel.setGlobalData(globalData);
  localKernel.setGlobalData(globalData);
  neighborKernel.setGlobalData(globalData);
  dynRupKernel.setGlobalData(globalData);
}

template <typename Cfg>
void ProxyDataImpl<Cfg>::initDataStructures(bool enableDR) {
  const initializer::LTSColorMap map(initializer::EnumLayer<HaloType>({HaloType::Interior}),
                                     initializer::EnumLayer<std::size_t>({0}),
                                     initializer::TraitLayer<ConfigVariant>({Cfg()}));

  Cfg cfg;

  // init RNG
  LTS::addTo(ltsStorage, false); // proxy does not use plasticity
  ltsStorage.setLayerCount(map);
  ltsStorage.fixate();

  ltsStorage.layer(layerId).setNumberOfCells(cellCount);

  LTS::Layer& layer = ltsStorage.layer(layerId);

  layer.setEntrySize<LTS::BuffersDerivatives>(sizeof(Real<Cfg>) * tensor::I<Cfg>::size() *
                                              layer.size());

  ltsStorage.allocateVariables();
  ltsStorage.touchVariables();
  ltsStorage.allocateBuckets();

  if (enableDR) {
    DynamicRupture dynRup;
    dynRup.addTo(drStorage);
    drStorage.setLayerCount(ltsStorage.getColorMap());
    drStorage.fixate();

    drStorage.layer(layerId).setNumberOfCells(4 * cellCount);

    drStorage.allocateVariables();
    drStorage.touchVariables();
  }

  // NOLINTNEXTLINE
  using real = Real<Cfg>;

  if (enableDR) {
    fakeDerivativesHost = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg> * sizeof(real),
        PagesizeHeap,
        seissol::memory::Standard));

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
      const auto offset = OpenMP::threadId();
      std::mt19937 rng(cellCount + offset);
      std::uniform_real_distribution<real> urd;
      for (std::size_t cell = 0; cell < cellCount; ++cell) {
        for (std::size_t i = 0; i < seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg>;
             i++) {
          fakeDerivativesHost[cell * seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg> +
                              i] = urd(rng);
        }
      }
    }

#ifdef ACL_DEVICE
    fakeDerivatives = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg> * sizeof(real),
        PagesizeHeap,
        seissol::memory::DeviceGlobalMemory));
    const auto& device = ::device::DeviceInstance::getInstance();
    device.api->copyTo(fakeDerivatives,
                       fakeDerivativesHost,
                       cellCount * seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg> *
                           sizeof(real));
#else
    fakeDerivatives = fakeDerivativesHost;
#endif
  }

  /* cell information and integration data*/
  fakeData(cfg, layer, (enableDR) ? FaceType::DynamicRupture : FaceType::Regular);

  if (enableDR) {
    // From lts storage
    CellDRMapping<Cfg>(*drMapping)[4] =
        isDeviceOn() ? layer.var<LTS::DRMappingDevice>(cfg) : layer.var<LTS::DRMapping>(cfg);

    constexpr initializer::AllocationPlace Place =
        isDeviceOn() ? initializer::AllocationPlace::Device : initializer::AllocationPlace::Host;

    // From dynamic rupture storage
    auto& interior = drStorage.layer(layerId);
    real(*imposedStatePlus)[seissol::tensor::QInterpolated<Cfg>::size()] =
        interior.var<DynamicRupture::ImposedStatePlus>(cfg, Place);
    real(*fluxSolverPlus)[seissol::tensor::fluxSolver<Cfg>::size()] =
        interior.var<DynamicRupture::FluxSolverPlus>(cfg, Place);
    real** timeDerivativeHostPlus = interior.var<DynamicRupture::TimeDerivativePlus>(cfg);
    real** timeDerivativeHostMinus = interior.var<DynamicRupture::TimeDerivativeMinus>(cfg);
    real** timeDerivativePlus = isDeviceOn()
                                    ? interior.var<DynamicRupture::TimeDerivativePlusDevice>(cfg)
                                    : interior.var<DynamicRupture::TimeDerivativePlus>(cfg);
    real** timeDerivativeMinus = isDeviceOn()
                                     ? interior.var<DynamicRupture::TimeDerivativeMinusDevice>(cfg)
                                     : interior.var<DynamicRupture::TimeDerivativeMinus>(cfg);
    DRFaceInformation* faceInformation = interior.var<DynamicRupture::FaceInformation>();

    std::mt19937 rng(cellCount);
    std::uniform_int_distribution<unsigned> sideDist(0, 3);
    std::uniform_int_distribution<unsigned> orientationDist(0, 2);
    std::uniform_int_distribution<std::size_t> drDist(0, interior.size() - 1);
    std::uniform_int_distribution<std::size_t> cellDist(0, cellCount - 1);

    /* init drMapping */
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        auto& drm = drMapping[cell][face];
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
          &fakeDerivativesHost[plusCell *
                               seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg>];
      timeDerivativeHostMinus[face] =
          &fakeDerivativesHost[minusCell *
                               seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg>];
      timeDerivativePlus[face] =
          &fakeDerivatives[plusCell * seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg>];
      timeDerivativeMinus[face] =
          &fakeDerivatives[minusCell *
                           seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg>];

      faceInformation[face].plusSide = sideDist(rng);
      faceInformation[face].minusSide = sideDist(rng);
      faceInformation[face].faceRelation = orientationDist(rng);
    }
  }
}

template <typename Cfg>
void ProxyDataImpl<Cfg>::initDataStructuresOnDevice(bool enableDR) {
#ifdef ACL_DEVICE
  const auto& device = ::device::DeviceInstance::getInstance();
  ltsStorage.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                           device.api->getDefaultStream());
  device.api->syncDefaultStreamWithHost();

  auto& layer = ltsStorage.layer(layerId);

  seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForWp(false, ltsStorage);
  ltsStorage.allocateScratchPads();

  seissol::initializer::recording::CompositeRecorder<LTS::LTSVarmap> recorder;
  recorder.addRecorder(new seissol::initializer::recording::LocalIntegrationRecorder);
  recorder.addRecorder(new seissol::initializer::recording::NeighIntegrationRecorder);

  recorder.addRecorder(new seissol::initializer::recording::PlasticityRecorder);
  recorder.record(layer);
  if (enableDR) {
    drStorage.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                            device.api->getDefaultStream());
    device.api->syncDefaultStreamWithHost();
    seissol::initializer::MemoryManager::deriveRequiredScratchpadMemoryForDr(drStorage);
    drStorage.allocateScratchPads();

    CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
    drRecorder.addRecorder(new DynamicRuptureRecorder);

    auto& drLayer = drStorage.layer(layerId);
    drRecorder.record(drLayer);
  }
#endif // ACL_DEVICE
}

std::shared_ptr<ProxyData>
    ProxyData::get(ConfigVariant variant, std::size_t cellCount, bool enableDR) {
  return std::visit(
      [&](auto cfg) -> std::shared_ptr<ProxyData> {
        using Cfg = decltype(cfg);
        return std::make_shared<ProxyDataImpl<Cfg>>(cellCount, enableDR);
      },
      variant);
}

} // namespace seissol::proxy
