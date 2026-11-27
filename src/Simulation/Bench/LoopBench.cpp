// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LoopBench.h"

#include "Alignment.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Real.h"
#include "Core/Common/Typedefs.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/CellLocalInformation.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/GlobalData.h"
#include "Core/Memory/LtsSetup.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/SimulationSettings.h"
#include "Core/Memory/Tree/Colormap.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/OpenMP.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Solver/DynamicRuptureCluster.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/CpuImpl/NoFault.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Features/DynamicRupture/Typedefs.h"
#include "GeneratedCode/tensor.h"
#include "Simulation/Bench/Proxy.h"
#include "Simulation/InitProcedure/Internal/Recording.h"
#include "WavePropagation/Solver/CellData.h"
#include "WavePropagation/Solver/SolverSelector.h"
#include "WavePropagation/Solver/Storage.h"
#include "WavePropagation/Solver/WaveCluster.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef LIKWID_PERFMON
#include "Core/Common/Instrumentation.h"
#endif

#ifdef ACL_DEVICE
#include "Features/DynamicRupture/Solver/FrictionLaws/GpuImpl/NoFault.h"
#include "Simulation/InitProcedure/Internal/Scratchpads.h"

#include <Device/device.h>
#endif

namespace seissol::proxy {

namespace {

/// the size of the time steps
constexpr double Timestep = 1.0;

constexpr double GravitationalAcceleration = 9.81;

/// a global time step from time 0
constexpr solver::StepParams Step{.time = 0,
                                  .timeStepSize = Timestep,
                                  .subTimeStart = 0,
                                  .neighborTimeStepSize = Timestep,
                                  .resetBuffers = true};

#ifdef ACL_DEVICE
template <typename Cfg>
using NoFaultOnDevice = dr::friction_law::gpu::NoFault<Cfg>;
#else
template <typename Cfg>
using NoFaultOnDevice = dr::friction_law::cpu::NoFault<Cfg>;
#endif

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

/// The loops of LoopBench in the configuration `Cfg`.
template <typename Cfg>
class Bench : public LoopBench {
  public:
  using real = Real<Cfg>;

  Bench(std::size_t cells, bool dynamicRupture, Executor executor);
  ~Bench() override;

  Bench(const Bench&) = delete;
  Bench(Bench&&) = delete;
  auto operator=(const Bench&) = delete;
  auto operator=(Bench&&) = delete;

  void run(Kernel kernel) override;
  void wait() override;
  [[nodiscard]] auto estimate(Kernel kernel) const -> PerformanceEstimate override;

  private:
  void initGlobalData();
  void initCells(FaceType faceType);
  void initFaults(const initializer::parameters::DRParameters& parameters);
  void initDevice(bool dynamicRupture);

  /// The cluster of the faults; throws without faults.
  [[nodiscard]] solver::DynamicRuptureCluster<Cfg>& faultCluster() const;

  /// `cluster`, once the work of the cluster that ran before, if another one, is done: each reads
  /// what the other one writes.
  template <typename ClusterT>
  ClusterT& use(ClusterT& cluster) {
    if (last_ != nullptr && last_ != &cluster) {
      waitFor(*last_);
    }
    last_ = &cluster;
    return cluster;
  }

  static void waitFor(solver::ComputeCluster& cluster);
  void waitForAll();

  /// Runs `work`, marked as the region `name` for LIKWID on the host.
  template <typename F>
  void marked(const char* name, F&& work);

  std::size_t cells_;
  Executor executor_;
  initializer::LayerIdentifier layerId_;

  memory::ManagedAllocator allocator_;
  GlobalDataStore globalData_;
  LTS::Storage ltsStorage_;
  DynamicRupture::Storage drStorage_;
  //! the derivatives of random cells, which the fault faces read
  real* faceDerivatives_{nullptr};
  real* faceDerivativesOnDevice_{nullptr};

  monitoring::FlopCounter flopCounter_;
  LoopStatistics loopStatistics_;

  std::unique_ptr<solver::WaveCluster<Cfg>> waveCluster_;
  std::unique_ptr<solver::DynamicRuptureCluster<Cfg>> faultCluster_;
  solver::ComputeCluster* last_{nullptr};
};

template <typename Cfg>
Bench<Cfg>::Bench(std::size_t cells, bool dynamicRupture, Executor executor)
    : cells_(cells), executor_(executor),
      layerId_(initializer::LayerIdentifier(HaloType::Interior, configIdOf<Cfg>(), 0)) {
  // the regions of the loops of the clusters, as the time manager adds them
  loopStatistics_.addRegion("computeLocalIntegration");
  loopStatistics_.addRegion("computeNeighboringIntegration");
  loopStatistics_.addRegion("computeDynamicRupture");

  // the parameters of the faults, for their storage and their friction law alike: no friction, and
  // the stress sources of a run by default
  const initializer::parameters::DRParameters drParameters{};

  initGlobalData();
  initCells(dynamicRupture ? FaceType::DynamicRupture : FaceType::Regular);
  if (dynamicRupture) {
    initFaults(drParameters);
  }
  initDevice(dynamicRupture);

  auto& layer = ltsStorage_.layer(layerId_);
  const typename solver::WaveCluster<Cfg>::Setup setup{.flopCounter = flopCounter_,
                                                       .globalData = globalData_,
                                                       .configs = {configIdOf<Cfg>()},
                                                       .initialConditions = nullptr,
                                                       .gravitationalAcceleration =
                                                           GravitationalAcceleration};
  waveCluster_ = std::make_unique<solver::WaveCluster<Cfg>>(/*clusterId=*/0,
                                                            /*globalClusterId=*/0,
                                                            layer.id(),
                                                            Timestep,
                                                            /*timeStepRate=*/1,
                                                            executor,
                                                            /*printProgress=*/false,
                                                            &layer,
                                                            setup,
                                                            &loopStatistics_,
                                                            nullptr,
                                                            nullptr);

  if (dynamicRupture) {
    const auto frictionLaw =
        dr::friction_law::makeFrictionSolverFactory<dr::friction_law::cpu::NoFault>(drParameters);
    const auto frictionLawDevice =
        dr::friction_law::makeFrictionSolverFactory<NoFaultOnDevice>(drParameters);
    auto& faces = drStorage_.layer(layerId_);
    faultCluster_ =
        std::make_unique<solver::DynamicRuptureCluster<Cfg>>(Timestep,
                                                             /*timeStepRate=*/1,
                                                             executor,
                                                             faces.id(),
                                                             std::vector<solver::PickpointPlan>{},
                                                             globalData_.template get<Cfg>(),
                                                             &faces,
                                                             frictionLaw,
                                                             frictionLawDevice,
                                                             flopCounter_,
                                                             &loopStatistics_,
                                                             nullptr,
                                                             nullptr);
  }
}

template <typename Cfg>
Bench<Cfg>::~Bench() {
  // as at the end of a run, the clusters give their streams and events back
  waitForAll();
  if (faultCluster_ != nullptr) {
    faultCluster_->finalize();
  }
  waveCluster_->finalize();
}

template <typename Cfg>
void Bench<Cfg>::initGlobalData() {
  initializer::GlobalDataInitializerOnHost::init<Cfg>(
      globalData_.template onHost<Cfg>(), allocator_, memory::Memkind::Standard);
  if constexpr (isDeviceOn()) {
    initializer::GlobalDataInitializerOnDevice::init<Cfg>(
        globalData_.template onDevice<Cfg>(), allocator_, memory::Memkind::DeviceGlobalMemory);
  }
}

template <typename Cfg>
void Bench<Cfg>::initCells(FaceType faceType) {
  const initializer::LTSColorMap map(initializer::EnumLayer<HaloType>({HaloType::Interior}),
                                     initializer::EnumLayer<std::size_t>({0}),
                                     initializer::EnumLayer<ConfigId>({configIdOf<Cfg>()}));

  const SimulationSettings settings(false, false, Cfg::NumSimulations, Cfg::NumSimulations);
  LTS::addTo(ltsStorage_);
  LTS::addWavePropagation(ltsStorage_, settings);
  ltsStorage_.setLayerCount(map);
  ltsStorage_.fixate();

  auto& layer = ltsStorage_.layer(layerId_);
  layer.setNumberOfCells(cells_);
  layer.setEntrySize<LTS::Buffers>(sizeof(real) * kernels::SolverOf<Cfg>::IntegralsSize * cells_);

  ltsStorage_.allocateVariables();
  ltsStorage_.touchVariables();
  ltsStorage_.allocateBuckets();

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
  fillWithStuff(reinterpret_cast<real*>(localIntegration),
                sizeof(LocalIntegrationData<Cfg>) / sizeof(real) * layer.size());
  fillWithStuff(reinterpret_cast<real*>(neighboringIntegration),
                sizeof(NeighboringIntegrationData<Cfg>) / sizeof(real) * layer.size());

  if constexpr (Cfg::Solver == SolverType::STP) {
#pragma omp parallel for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      localIntegration[cell].specific.typicalTimeStepWidth = Timestep;
    }
  }
}

template <typename Cfg>
void Bench<Cfg>::initFaults(const initializer::parameters::DRParameters& parameters) {
  DynamicRupture dynamicRupture(&parameters);
  dynamicRupture.addTo(drStorage_);
  drStorage_.setLayerCount(ltsStorage_.getColorMap());
  drStorage_.fixate();

  auto& faces = drStorage_.layer(layerId_);
  faces.setNumberOfCells(Cell::NumFaces * cells_);

  drStorage_.allocateVariables();
  drStorage_.touchVariables();

  constexpr auto DerivativesSize = kernels::SolverOf<Cfg>::DerivativesSize;
  const auto derivativesBytes = cells_ * DerivativesSize * sizeof(real);
  faceDerivatives_ = static_cast<real*>(
      allocator_.allocateMemory(derivativesBytes, PagesizeHeap, memory::Memkind::Standard));

#pragma omp parallel
  {
    std::mt19937 rng(cells_ + OpenMP::threadId());
    std::uniform_real_distribution<real> urd;
#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < cells_; ++cell) {
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

  auto& layer = ltsStorage_.layer(layerId_);
  auto* drMapping = layer.var<LTS::DRMapping>(Cfg());

  real(*imposedStatePlus)[tensor::QInterpolated<Cfg>::size()] =
      faces.var<DynamicRupture::ImposedStatePlus>(Cfg());
  real(*fluxSolverPlus)[tensor::fluxSolver<Cfg>::size()] =
      faces.var<DynamicRupture::FluxSolverPlus>(Cfg());
  real** timeDerivativePlus = faces.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  real** timeDerivativeMinus = faces.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  auto* faceInformation = faces.var<DynamicRupture::FaceInformation>();

  std::mt19937 rng(cells_);
  std::uniform_int_distribution<unsigned> sideDist(0, 3);
  std::uniform_int_distribution<unsigned> orientationDist(0, 1);
  std::uniform_int_distribution<std::size_t> faceDist(0, faces.size() - 1);
  std::uniform_int_distribution<std::size_t> cellDist(0, cells_ - 1);

  // each face of a cell takes the flux of a random fault face; as in the solver, the mapping of a
  // face addresses that face
  for (std::size_t cell = 0; cell < cells_; ++cell) {
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

  for (std::size_t face = 0; face < faces.size(); ++face) {
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
        faces.var<DynamicRupture::ImposedStatePlus>(Cfg(), initializer::AllocationPlace::Device));
    auto* fluxSolverPlusDevice = reinterpret_cast<real*>(
        faces.var<DynamicRupture::FluxSolverPlus>(Cfg(), initializer::AllocationPlace::Device));
    real** timeDerivativePlusDevice = faces.var<DynamicRupture::TimeDerivativePlusDevice>(Cfg());
    real** timeDerivativeMinusDevice = faces.var<DynamicRupture::TimeDerivativeMinusDevice>(Cfg());
    for (std::size_t cell = 0; cell < cells_; ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        const auto& mapping = drMapping[cell][face];
        auto& mappingDevice = drMappingDevice[cell][face];
        mappingDevice.side = mapping.side;
        mappingDevice.faceRelation = mapping.faceRelation;
        mappingDevice.godunov = imposedStatePlusDevice + (mapping.godunov - imposedStatePlus[0]);
        mappingDevice.fluxSolver = fluxSolverPlusDevice + (mapping.fluxSolver - fluxSolverPlus[0]);
      }
    }
    for (std::size_t face = 0; face < faces.size(); ++face) {
      timeDerivativePlusDevice[face] =
          faceDerivativesOnDevice_ + (timeDerivativePlus[face] - faceDerivatives_);
      timeDerivativeMinusDevice[face] =
          faceDerivativesOnDevice_ + (timeDerivativeMinus[face] - faceDerivatives_);
    }
  }
}

template <typename Cfg>
void Bench<Cfg>::initDevice([[maybe_unused]] bool dynamicRupture) {
#ifdef ACL_DEVICE
  const auto& device = device::DeviceInstance::instance();
  ltsStorage_.synchronizeTo(initializer::AllocationPlace::Device, device.api().getDefaultStream());
  if (dynamicRupture) {
    drStorage_.synchronizeTo(initializer::AllocationPlace::Device, device.api().getDefaultStream());
  }
  device.api().syncDefaultStreamWithHost();

  initializer::internal::deriveRequiredScratchpadMemoryForWp(false, ltsStorage_);
  ltsStorage_.allocateScratchPads();
  if (dynamicRupture) {
    initializer::internal::deriveRequiredScratchpadMemoryForDr(drStorage_);
    drStorage_.allocateScratchPads();
  }
#endif
  initializer::internal::setupRecorders(
      ltsStorage_, drStorage_, /*usePlasticity=*/false, GravitationalAcceleration);
}

template <typename Cfg>
solver::DynamicRuptureCluster<Cfg>& Bench<Cfg>::faultCluster() const {
  if (faultCluster_ == nullptr) {
    throw std::invalid_argument("The loops of the faces of dynamic rupture need the faces.");
  }
  return *faultCluster_;
}

template <typename Cfg>
void Bench<Cfg>::waitFor(solver::ComputeCluster& cluster) {
  for (auto* runtime : cluster.runtimes()) {
    runtime->wait();
  }
}

template <typename Cfg>
template <typename F>
void Bench<Cfg>::marked([[maybe_unused]] const char* name, F&& work) {
#ifdef LIKWID_PERFMON
  const bool host = executor_ == Executor::Host;
  if (host) {
#pragma omp parallel
    {
      LIKWID_MARKER_START(name);
    }
  }
#endif
  std::invoke(std::forward<F>(work));
#ifdef LIKWID_PERFMON
  if (host) {
#pragma omp parallel
    {
      LIKWID_MARKER_STOP(name);
    }
  }
#endif
}

template <typename Cfg>
void Bench<Cfg>::run(Kernel kernel) {
  switch (kernel) {
  case Kernel::Local:
    marked("local", [&] { use(*waveCluster_).predictOnce(Step); });
    break;
  case Kernel::Neighbor:
  case Kernel::NeighborDR:
    marked("neighboring", [&] { use(*waveCluster_).correctOnce(Step); });
    break;
  case Kernel::GodunovDR:
    use(faultCluster()).interactOnce(Step);
    break;
  case Kernel::All:
    run(Kernel::Local);
    run(Kernel::Neighbor);
    break;
  case Kernel::AllDR:
    run(Kernel::Local);
    run(Kernel::GodunovDR);
    run(Kernel::NeighborDR);
    break;
  }
}

template <typename Cfg>
void Bench<Cfg>::wait() {
  waitForAll();
}

template <typename Cfg>
void Bench<Cfg>::waitForAll() {
  waitFor(*waveCluster_);
  if (faultCluster_ != nullptr) {
    waitFor(*faultCluster_);
  }
}

template <typename Cfg>
auto Bench<Cfg>::estimate(Kernel kernel) const -> PerformanceEstimate {
  switch (kernel) {
  case Kernel::Local:
    return waveCluster_->predictionEstimate();
  case Kernel::Neighbor:
  case Kernel::NeighborDR:
    return waveCluster_->correctionEstimate();
  case Kernel::GodunovDR:
    return faultCluster().interactionEstimate();
  case Kernel::All:
    return estimate(Kernel::Local) + estimate(Kernel::Neighbor);
  case Kernel::AllDR:
    return estimate(Kernel::Local) + estimate(Kernel::GodunovDR) + estimate(Kernel::NeighborDR);
  }
  return {};
}

} // namespace

auto makeLoopBench(ConfigId config, std::size_t cells, bool dynamicRupture, Executor executor)
    -> std::unique_ptr<LoopBench> {
  return dispatchConfig(config, [&](auto cfg) -> std::unique_ptr<LoopBench> {
    using Cfg = decltype(cfg);
    return std::make_unique<Bench<Cfg>>(cells, dynamicRupture, executor);
  });
}

} // namespace seissol::proxy
