// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Allocator.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Constants.h"
#include "Config.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Kernels/Touch.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/GlobalData.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Colormap.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/OpenMP.h"
#include "Proxy/Constants.h"
#include "Solver/Settings.h"

#include <cstddef>
#include <random>
#include <stdlib.h>

#ifdef ACL_DEVICE
#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Initializer/BatchRecorders/Recorders.h"
#include "Initializer/InitProcedure/Internal/Scratchpads.h"

#include <Device/device.h>
#include <memory>
#endif

namespace seissol::proxy {

namespace {

template <typename Cfg>
void fakeData(LTS::Layer& layer, FaceType faceTp) {
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  real(*dofs)[tensor::Q<Cfg>::size()] = layer.var<LTS::Dofs>(Cfg());
  real** buffers = layer.var<LTS::StepIntegrals>(Cfg());
  real** derivatives = layer.var<LTS::Derivatives>(Cfg());
  auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  auto* localIntegration = layer.var<LTS::LocalIntegration>(Cfg());
  auto* neighboringIntegration = layer.var<LTS::NeighboringIntegration>(Cfg());
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
  real* bucket =
      static_cast<real*>(layer.var<LTS::Buffers>(Cfg(), initializer::AllocationPlace::Host));

  real** buffersDevice = layer.var<LTS::StepIntegralsDevice>(Cfg());
  real** derivativesDevice = layer.var<LTS::DerivativesDevice>(Cfg());
  auto* faceNeighborsDevice = layer.var<LTS::FaceNeighborsDevice>();
  real* bucketDevice =
      static_cast<real*>(layer.var<LTS::Buffers>(Cfg(), initializer::AllocationPlace::Device));

  std::mt19937 rng(layer.size());
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<std::size_t> cellDist(0, layer.size() - 1);

  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    buffers[cell] = bucket + cell * kernels::SolverOf<Cfg>::IntegralsSize;
    derivatives[cell] = nullptr;
    buffersDevice[cell] = bucketDevice + cell * kernels::SolverOf<Cfg>::IntegralsSize;
    derivativesDevice[cell] = nullptr;

    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      cellInformation[cell].faceTypes[f] = faceTp;
      cellInformation[cell].faceRelations[f][0] = sideDist(rng);
      cellInformation[cell].faceRelations[f][1] = 0;
      cellInformation[cell].neighborConfigIds[f] = configIdOf<Cfg>();

      const auto neighbor = cellDist(rng);
      secondaryInformation[cell].faceNeighbors[f].global = neighbor;
      secondaryInformation[cell].faceNeighbors[f].color = 0;
      secondaryInformation[cell].faceNeighbors[f].cell = neighbor;
    }
    cellInformation[cell].ltsSetup = LtsSetup();
  }

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
      switch (faceTp) {
      case FaceType::FreeSurface:
        faceNeighbors[cell][f] = buffers[cell];
        faceNeighborsDevice[cell][f] = buffersDevice[cell];
        break;
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
  kernels::fillWithStuff(bucket, kernels::SolverOf<Cfg>::IntegralsSize * layer.size(), false);
  kernels::fillWithStuff(reinterpret_cast<real*>(localIntegration),
                         sizeof(LocalIntegrationData<Cfg>) / sizeof(real) * layer.size(),
                         false);
  kernels::fillWithStuff(reinterpret_cast<real*>(neighboringIntegration),
                         sizeof(NeighboringIntegrationData<Cfg>) / sizeof(real) * layer.size(),
                         false);

  if constexpr (Cfg::Solver == SolverType::STP) {
#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      localIntegration[cell].specific.typicalTimeStepWidth = seissol::proxy::Timestep;
    }
  }

#ifdef ACL_DEVICE
  const auto& device = device::DeviceInstance::instance();
  layer.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                      device.api().getDefaultStream());
  device.api().syncDefaultStreamWithHost();
#endif
}
} // namespace

ProxyData::ProxyData(std::size_t cellCount, ConfigId config)
    : cellCount(cellCount), config(config),
      layerId(initializer::LayerIdentifier(HaloType::Interior, config, 0)) {}

template <typename Cfg>
ProxyDataImpl<Cfg>::ProxyDataImpl(std::size_t cellCount, bool enableDR)
    : ProxyData(cellCount, configIdOf<Cfg>()) {
  initGlobalData();
  initDataStructures(enableDR);
  initDataStructuresOnDevice(enableDR);
}

template <typename Cfg>
void ProxyDataImpl<Cfg>::initGlobalData() {
  seissol::initializer::GlobalDataInitializerOnHost::init<Cfg>(
      globalDataOnHost, allocator, seissol::memory::Memkind::Standard);

  CompoundGlobalData<Cfg> globalData{};
  globalData.onHost = &globalDataOnHost;
  globalData.onDevice = nullptr;
  if constexpr (seissol::isDeviceOn()) {
    seissol::initializer::GlobalDataInitializerOnDevice::init<Cfg>(
        globalDataOnDevice, allocator, seissol::memory::Memkind::DeviceGlobalMemory);
    globalData.onDevice = &globalDataOnDevice;
  }
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
                                     initializer::EnumLayer<ConfigId>({configIdOf<Cfg>()}));

  // init RNG
  const auto nullSettings = SimulationSettings(false, false);
  LTS::addTo(ltsStorage, nullSettings);
  ltsStorage.setLayerCount(map);
  ltsStorage.fixate();

  ltsStorage.layer(layerId).setNumberOfCells(cellCount);

  LTS::Layer& layer = ltsStorage.layer(layerId);
  layer.setEntrySize<LTS::Buffers>(sizeof(real) * kernels::SolverOf<Cfg>::IntegralsSize *
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

    fakeDerivativesHost = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::SolverOf<Cfg>::DerivativesSize * sizeof(real),
        PagesizeHeap,
        seissol::memory::Memkind::Standard));

#pragma omp parallel
    {
      const auto offset = OpenMP::threadId();
      std::mt19937 rng(cellCount + offset);
      std::uniform_real_distribution<real> urd;
      for (std::size_t cell = 0; cell < cellCount; ++cell) {
        for (std::size_t i = 0; i < seissol::kernels::SolverOf<Cfg>::DerivativesSize; i++) {
          fakeDerivativesHost[cell * seissol::kernels::SolverOf<Cfg>::DerivativesSize + i] =
              urd(rng);
        }
      }
    }

#ifdef ACL_DEVICE
    fakeDerivatives = reinterpret_cast<real*>(allocator.allocateMemory(
        cellCount * seissol::kernels::SolverOf<Cfg>::DerivativesSize * sizeof(real),
        PagesizeHeap,
        seissol::memory::Memkind::DeviceGlobalMemory));
    const auto& device = ::device::DeviceInstance::instance();
    device.api().copyTo(fakeDerivatives,
                        fakeDerivativesHost,
                        cellCount * seissol::kernels::SolverOf<Cfg>::DerivativesSize *
                            sizeof(real));
#else
    fakeDerivatives = fakeDerivativesHost;
#endif
  }

  /* cell information and integration data*/
  fakeData<Cfg>(layer, enableDR ? FaceType::DynamicRupture : FaceType::Regular);

  if (enableDR) {
    // From lts storage
    auto* drMapping =
        isDeviceOn() ? layer.var<LTS::DRMappingDevice>(Cfg()) : layer.var<LTS::DRMapping>(Cfg());

    constexpr initializer::AllocationPlace Place =
        isDeviceOn() ? initializer::AllocationPlace::Device : initializer::AllocationPlace::Host;

    // From dynamic rupture storage
    DynamicRupture::Layer& interior = drStorage.layer(layerId);
    real(*imposedStatePlus)[seissol::tensor::QInterpolated<Cfg>::size()] =
        interior.var<DynamicRupture::ImposedStatePlus>(Cfg(), Place);
    real(*fluxSolverPlus)[seissol::tensor::fluxSolver<Cfg>::size()] =
        interior.var<DynamicRupture::FluxSolverPlus>(Cfg(), Place);
    real** timeDerivativeHostPlus = interior.var<DynamicRupture::TimeDerivativePlus>(Cfg());
    real** timeDerivativeHostMinus = interior.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
    real** timeDerivativePlus = isDeviceOn()
                                    ? interior.var<DynamicRupture::TimeDerivativePlusDevice>(Cfg())
                                    : interior.var<DynamicRupture::TimeDerivativePlus>(Cfg());
    real** timeDerivativeMinus =
        isDeviceOn() ? interior.var<DynamicRupture::TimeDerivativeMinusDevice>(Cfg())
                     : interior.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
    DRFaceInformation* faceInformation = interior.var<DynamicRupture::FaceInformation>();

    std::mt19937 rng(cellCount);
    std::uniform_int_distribution<unsigned> sideDist(0, 3);
    std::uniform_int_distribution<unsigned> orientationDist(0, 1);
    std::uniform_int_distribution<std::size_t> drDist(0, interior.size() - 1);
    std::uniform_int_distribution<std::size_t> cellDist(0, cellCount - 1);

    /* init drMapping */
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        auto& drm = drMapping[cell][face];
        const auto orientation = orientationDist(rng);
        const auto drFace = drDist(rng);
        // as in the solver, the mapping of a face addresses that face
        drm.side = face;
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
          &fakeDerivativesHost[plusCell * seissol::kernels::SolverOf<Cfg>::DerivativesSize];
      timeDerivativeHostMinus[face] =
          &fakeDerivativesHost[minusCell * seissol::kernels::SolverOf<Cfg>::DerivativesSize];
      timeDerivativePlus[face] =
          &fakeDerivatives[plusCell * seissol::kernels::SolverOf<Cfg>::DerivativesSize];
      timeDerivativeMinus[face] =
          &fakeDerivatives[minusCell * seissol::kernels::SolverOf<Cfg>::DerivativesSize];

      faceInformation[face].plusSide = sideDist(rng);
      faceInformation[face].minusSide = sideDist(rng);
      // a fault face always addresses the minus side here
      faceInformation[face].faceRelation = 1;
    }
  }
}

template <typename Cfg>
void ProxyDataImpl<Cfg>::initDataStructuresOnDevice(bool enableDR) {
#ifdef ACL_DEVICE
  const auto& device = ::device::DeviceInstance::instance();
  ltsStorage.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                           device.api().getDefaultStream());
  device.api().syncDefaultStreamWithHost();

  LTS::Layer& layer = ltsStorage.layer(layerId);

  seissol::initializer::internal::deriveRequiredScratchpadMemoryForWp(false, ltsStorage);
  ltsStorage.allocateScratchPads();

  seissol::recording::CompositeRecorder<LTS::LTSVarmap> recorder;
  recorder.addRecorder(new seissol::recording::LocalIntegrationRecorder(9.81));
  recorder.addRecorder(new seissol::recording::NeighIntegrationRecorder);

  recorder.addRecorder(new seissol::recording::PlasticityRecorder);
  recorder.record(layer);
  if (enableDR) {
    drStorage.synchronizeTo(seissol::initializer::AllocationPlace::Device,
                            device.api().getDefaultStream());
    device.api().syncDefaultStreamWithHost();
    seissol::initializer::internal::deriveRequiredScratchpadMemoryForDr(drStorage);
    drStorage.allocateScratchPads();

    seissol::recording::CompositeRecorder<DynamicRupture::DynrupVarmap> drRecorder;
    drRecorder.addRecorder(new seissol::recording::DynamicRuptureRecorder);

    auto& drLayer = drStorage.layer(layerId);
    drRecorder.record(drLayer);
  }
#endif // ACL_DEVICE
}

std::shared_ptr<ProxyData> makeProxyData(ConfigId config, std::size_t cellCount, bool enableDR) {
  return dispatchConfig(config, [&](auto cfg) -> std::shared_ptr<ProxyData> {
    using Cfg = decltype(cfg);
    return std::make_shared<ProxyDataImpl<Cfg>>(cellCount, enableDR);
  });
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg) template struct ProxyDataImpl<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

} // namespace seissol::proxy
