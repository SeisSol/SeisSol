// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SyntheticData.h"

#include "Alignment.h"
#include "Config.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Typedefs.h"
#include "Core/Memory/CellLocalInformation.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/GlobalData.h"
#include "Core/Memory/LtsSetup.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/SimulationSettings.h"
#include "Core/Memory/Tree/Colormap.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Parallel/OpenMP.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "Features/Plasticity/Solver/Storage.h"
#include "GeneratedCode/tensor.h"
#include "Simulation/Bench/LoopBench.h"
#include "Simulation/Bench/Proxy.h"
#include "Simulation/InitProcedure/Internal/Recording.h"
#include "WavePropagation/Solver/SolverSelector.h"
#include "WavePropagation/Solver/Storage.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <random>

#ifdef ACL_DEVICE
#include "Simulation/InitProcedure/Internal/Scratchpads.h"

#include <Device/device.h>
#endif

namespace seissol::proxy {

namespace {

/// Fills `buffer` with arbitrary finite values that are not zero.
template <typename RealT>
void fillWithStuff(RealT* buffer, std::size_t count) {
  // No real point for these numbers. Should be just something != 0 and != NaN and != Inf
  const auto stuff = [](std::size_t n) {
    return static_cast<RealT>((214013.0 * static_cast<double>(n) + 2531011.0) / 16777216.0);
  };
#pragma omp parallel for schedule(static)
  for (std::size_t n = 0; n < count; ++n) {
    buffer[n] = stuff(n);
  }
}

} // namespace

template <typename Cfg>
SyntheticData<Cfg>::SyntheticData(const BenchSetup& setup)
    : cellCount_(setup.cells), dynamicRupture_(setup.dynamicRupture), plasticity_(setup.plasticity),
      layerId_(initializer::LayerIdentifier(HaloType::Interior, configIdOf<Cfg>(), 0)) {
  initGlobalData();
  initCells(dynamicRupture_ ? FaceType::DynamicRupture : FaceType::Regular);
  if (plasticity_ != CellPlasticity::None) {
    initPlasticity();
  }
  if (dynamicRupture_) {
    initFaults();
  }
  initDevice();
}

template <typename Cfg>
void SyntheticData<Cfg>::initGlobalData() {
  initializer::GlobalDataInitializerOnHost::init<Cfg>(
      globalData_.template onHost<Cfg>(), allocator_, memory::Memkind::Standard);
  if constexpr (isDeviceOn()) {
    initializer::GlobalDataInitializerOnDevice::init<Cfg>(
        globalData_.template onDevice<Cfg>(), allocator_, memory::Memkind::DeviceGlobalMemory);
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initCells(FaceType faceType) {
  const initializer::LTSColorMap map(initializer::EnumLayer<HaloType>({HaloType::Interior}),
                                     initializer::EnumLayer<std::size_t>({0}),
                                     initializer::EnumLayer<ConfigId>({configIdOf<Cfg>()}));

  const bool plasticity = plasticity_ != CellPlasticity::None;
  const SimulationSettings settings(plasticity, false, Cfg::NumSimulations, Cfg::NumSimulations);
  LTS::addTo(ltsStorage_);
  LTS::addWavePropagation(ltsStorage_, settings);
  LTSPlasticity::addTo(ltsStorage_, plasticity);
  ltsStorage_.setLayerCount(map);
  ltsStorage_.fixate();

  LTS::Layer& layer = cells();
  layer.setNumberOfCells(cellCount_);
  layer.setEntrySize<LTS::Buffers>(sizeof(real) * kernels::SolverOf<Cfg>::IntegralsSize *
                                   cellCount_);

  ltsStorage_.allocateVariables();
  ltsStorage_.touchVariables();
  ltsStorage_.allocateBuckets();

  real(*dofs)[tensor::Q<Cfg>::size()] = layer.var<LTS::Dofs>(Cfg());
  real** buffers = layer.var<LTS::StepIntegrals>(Cfg());
  real** derivatives = layer.var<LTS::Derivatives>(Cfg());
  auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
  real* bucket =
      static_cast<real*>(layer.var<LTS::Buffers>(Cfg(), initializer::AllocationPlace::Host));

  // every cell keeps the time integral of its step in its buffer, which its neighbors read
  LtsSetup gts;
  gts.setHasBuffer(true, BufferType::StepIntegrals);
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    gts.setNeighborBuffer(face, BufferType::StepIntegrals).setNeighborGTSRelation(face, true);
  }

  std::mt19937 rng(layer.size());
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<std::size_t> cellDist(0, layer.size() - 1);

  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    buffers[cell] = bucket + cell * kernels::SolverOf<Cfg>::IntegralsSize;
    derivatives[cell] = nullptr;

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      cellInformation[cell].faceTypes[face] = faceType;
      cellInformation[cell].faceRelations[face][0] = sideDist(rng);
      cellInformation[cell].faceRelations[face][1] = 0;
      cellInformation[cell].neighborConfigIds[face] = configIdOf<Cfg>();

      const auto neighbor = cellDist(rng);
      secondaryInformation[cell].faceNeighbors[face].global = neighbor;
      secondaryInformation[cell].faceNeighbors[face].color = 0;
      secondaryInformation[cell].faceNeighbors[face].cell = neighbor;
    }
    cellInformation[cell].ltsSetup = gts;
  }

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (faceType == FaceType::Regular) {
        const auto neighbor = secondaryInformation[cell].faceNeighbors[face].cell;
        faceNeighbors[cell][face] = buffers[neighbor];
      } else {
        faceNeighbors[cell][face] = nullptr;
      }
    }
  }

  // the same addresses on the device (without one, these variables stay empty)
  if constexpr (isDeviceOn()) {
    real** buffersDevice = layer.var<LTS::StepIntegralsDevice>(Cfg());
    real** derivativesDevice = layer.var<LTS::DerivativesDevice>(Cfg());
    auto* faceNeighborsDevice = layer.var<LTS::FaceNeighborsDevice>();
    real* bucketDevice =
        static_cast<real*>(layer.var<LTS::Buffers>(Cfg(), initializer::AllocationPlace::Device));
#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      buffersDevice[cell] = bucketDevice + cell * kernels::SolverOf<Cfg>::IntegralsSize;
      derivativesDevice[cell] = nullptr;
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (faceType == FaceType::Regular) {
          const auto neighbor = secondaryInformation[cell].faceNeighbors[face].cell;
          faceNeighborsDevice[cell][face] =
              bucketDevice + neighbor * kernels::SolverOf<Cfg>::IntegralsSize;
        }
      }
    }
  }

  fillWithStuff(reinterpret_cast<real*>(dofs), tensor::Q<Cfg>::size() * layer.size());
  fillWithStuff(bucket, kernels::SolverOf<Cfg>::IntegralsSize * layer.size());
  // The local and the neighboring integration data (the star matrices, the flux solvers, ...)
  // stay zero: the kernels compute with them as with any other values, and the cells keep the
  // values of their unknowns from step to step. With arbitrary values there, these grew by about
  // 24 orders of magnitude in every step, and were not finite after a few (in single precision,
  // after the first); the plasticity then found nothing to yield.

  if constexpr (Cfg::Solver == SolverType::STP) {
    auto* localIntegration = layer.var<LTS::LocalIntegration>(Cfg());
#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      localIntegration[cell].specific.typicalTimeStepWidth = Timestep;
    }
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initPlasticity() {
  LTS::Layer& layer = cells();
  auto* cellInformation = layer.var<LTS::CellInformation>();
  auto* plasticity = layer.var<LTSPlasticity::Plasticity>(Cfg());

  // A cell without cohesion and friction yields wherever its stress has a deviator, which the
  // arbitrary values of the cells give it; with the slow relaxation (see RelaxationTime), it
  // stays beyond the yield surface. Otherwise, the cohesion lies far above any stress of the
  // cells.
  const auto cohesion =
      plasticity_ == CellPlasticity::Yield ? static_cast<real>(0) : static_cast<real>(1e20);
  // 1 / (2 mu), with the shear modulus of crustal rock (32 GPa)
  const auto muFactor = static_cast<real>(1.0 / (2.0 * 3.2e10));

#pragma omp parallel for schedule(static)
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    cellInformation[cell].plasticityEnabled = true;
    auto& data = plasticity[cell];
    std::fill(std::begin(data.initialLoading), std::end(data.initialLoading), 0);
    std::fill(std::begin(data.cohesionTimesCosAngularFriction),
              std::end(data.cohesionTimesCosAngularFriction),
              cohesion);
    std::fill(std::begin(data.sinAngularFriction), std::end(data.sinAngularFriction), 0);
    data.mufactor = muFactor;
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initFaults() {
  DynamicRupture dynamicRupture(&drParameters_);
  dynamicRupture.addTo(drStorage_);
  drStorage_.setLayerCount(ltsStorage_.getColorMap());
  drStorage_.fixate();

  DynamicRupture::Layer& faults = faces();
  faults.setNumberOfCells(Cell::NumFaces * cellCount_);

  drStorage_.allocateVariables();
  drStorage_.touchVariables();

  constexpr auto DerivativesSize = kernels::SolverOf<Cfg>::DerivativesSize;
  const auto derivativesBytes = cellCount_ * DerivativesSize * sizeof(real);
  faceDerivatives_ = static_cast<real*>(
      allocator_.allocateMemory(derivativesBytes, PagesizeHeap, memory::Memkind::Standard));

#pragma omp parallel
  {
    std::mt19937 rng(cellCount_ + OpenMP::threadId());
    std::uniform_real_distribution<real> urd;
#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < cellCount_; ++cell) {
      for (std::size_t i = 0; i < DerivativesSize; ++i) {
        faceDerivatives_[cell * DerivativesSize + i] = urd(rng);
      }
    }
  }

#ifdef ACL_DEVICE
  faceDerivativesOnDevice_ = static_cast<real*>(allocator_.allocateMemory(
      derivativesBytes, PagesizeHeap, memory::Memkind::DeviceGlobalMemory));
  device::DeviceInstance::instance().api().copyTo(
      faceDerivativesOnDevice_, faceDerivatives_, derivativesBytes);
#endif

  LTS::Layer& layer = cells();
  auto* drMapping = layer.var<LTS::DRMapping>(Cfg());

  real(*imposedStatePlus)[tensor::QInterpolated<Cfg>::size()] =
      faults.var<DynamicRupture::ImposedStatePlus>(Cfg());
  real(*fluxSolverPlus)[tensor::fluxSolver<Cfg>::size()] =
      faults.var<DynamicRupture::FluxSolverPlus>(Cfg());
  real** timeDerivativePlus = faults.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  real** timeDerivativeMinus = faults.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  auto* faceInformation = faults.var<DynamicRupture::FaceInformation>();

  std::mt19937 rng(cellCount_);
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<unsigned> orientationDist(0, 1);
  std::uniform_int_distribution<std::size_t> faceDist(0, faults.size() - 1);
  std::uniform_int_distribution<std::size_t> cellDist(0, cellCount_ - 1);

  // each face of a cell takes the flux of a random fault face; as in the solver, the mapping of a
  // face addresses that face
  for (std::size_t cell = 0; cell < cellCount_; ++cell) {
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      const auto orientation = orientationDist(rng);
      const auto faultFace = faceDist(rng);
      auto& mapping = drMapping[cell][face];
      mapping.side = static_cast<std::int8_t>(face);
      mapping.faceRelation = static_cast<std::int8_t>(orientation);
      mapping.godunov = imposedStatePlus[faultFace];
      mapping.fluxSolver = fluxSolverPlus[faultFace];
    }
  }

  for (std::size_t face = 0; face < faults.size(); ++face) {
    const auto plusCell = cellDist(rng);
    const auto minusCell = cellDist(rng);
    timeDerivativePlus[face] = faceDerivatives_ + plusCell * DerivativesSize;
    timeDerivativeMinus[face] = faceDerivatives_ + minusCell * DerivativesSize;

    faceInformation[face].plusSide = sideDist(rng);
    faceInformation[face].minusSide = sideDist(rng);
    // a fault face always addresses the minus side here
    faceInformation[face].faceRelation = 1;
  }

  // the same addresses on the device: at the same offsets into the data there
  if constexpr (isDeviceOn()) {
    auto* drMappingDevice = layer.var<LTS::DRMappingDevice>(Cfg());
    auto* imposedStatePlusDevice = reinterpret_cast<real*>(
        faults.var<DynamicRupture::ImposedStatePlus>(Cfg(), initializer::AllocationPlace::Device));
    auto* fluxSolverPlusDevice = reinterpret_cast<real*>(
        faults.var<DynamicRupture::FluxSolverPlus>(Cfg(), initializer::AllocationPlace::Device));
    real** timeDerivativePlusDevice = faults.var<DynamicRupture::TimeDerivativePlusDevice>(Cfg());
    real** timeDerivativeMinusDevice = faults.var<DynamicRupture::TimeDerivativeMinusDevice>(Cfg());
    for (std::size_t cell = 0; cell < cellCount_; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        const auto& mapping = drMapping[cell][face];
        auto& mappingDevice = drMappingDevice[cell][face];
        mappingDevice.side = mapping.side;
        mappingDevice.faceRelation = mapping.faceRelation;
        mappingDevice.godunov = imposedStatePlusDevice + (mapping.godunov - imposedStatePlus[0]);
        mappingDevice.fluxSolver = fluxSolverPlusDevice + (mapping.fluxSolver - fluxSolverPlus[0]);
      }
    }
    for (std::size_t face = 0; face < faults.size(); ++face) {
      timeDerivativePlusDevice[face] =
          faceDerivativesOnDevice_ + (timeDerivativePlus[face] - faceDerivatives_);
      timeDerivativeMinusDevice[face] =
          faceDerivativesOnDevice_ + (timeDerivativeMinus[face] - faceDerivatives_);
    }
  }
}

template <typename Cfg>
void SyntheticData<Cfg>::initDevice() {
  const bool plasticity = plasticity_ != CellPlasticity::None;
#ifdef ACL_DEVICE
  const auto& device = device::DeviceInstance::instance();
  ltsStorage_.synchronizeTo(initializer::AllocationPlace::Device, device.api().getDefaultStream());
  if (dynamicRupture_) {
    drStorage_.synchronizeTo(initializer::AllocationPlace::Device, device.api().getDefaultStream());
  }
  device.api().syncDefaultStreamWithHost();

  initializer::internal::deriveRequiredScratchpadMemoryForWp(plasticity, ltsStorage_);
  ltsStorage_.allocateScratchPads();
  if (dynamicRupture_) {
    initializer::internal::deriveRequiredScratchpadMemoryForDr(drStorage_);
    drStorage_.allocateScratchPads();
  }
#endif
  initializer::internal::setupRecorders(
      ltsStorage_, drStorage_, plasticity, GravitationalAcceleration);
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg) template class SyntheticData<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

} // namespace seissol::proxy
