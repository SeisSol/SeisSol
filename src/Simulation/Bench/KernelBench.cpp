// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Alignment.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Instrumentation.h"
#include "Core/Common/Real.h"
#include "Core/Memory/CellLocalInformation.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/MemoryAllocator.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Numerical/Quadrature.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Features/DynamicRupture/Solver/DynamicRupture.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Features/Plasticity/Solver/Plasticity.h"
#include "Features/Plasticity/Solver/Storage.h"
#include "GeneratedCode/tensor.h"
#include "Simulation/Bench/LoopBench.h"
#include "Simulation/Bench/Proxy.h"
#include "Simulation/Bench/SyntheticData.h"
#include "WavePropagation/Solver/Interface.h"
#include "WavePropagation/Solver/Solver.h"
#include "WavePropagation/Solver/SolverSelector.h"
#include "WavePropagation/Solver/Storage.h"
#include "WavePropagation/Solver/TimeCommon.h"

#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>

#ifdef ACL_DEVICE
#include "Core/Memory/Batch/DeviceGraph.h"
#endif

namespace seissol::proxy {

namespace {

/**
 * The kernels of a time step in the configuration `Cfg`, each over all cells or faces of
 * SyntheticData on its own, as the proxy ran them before it ran the clusters: without what the
 * clusters do around them (the displacements, the receivers, the buffers of LTS, the friction law
 * after the space-time interpolation, ...).
 */
template <typename Cfg>
class KernelBench : public LoopBench {
  public:
  using real = Real<Cfg>;

  explicit KernelBench(const BenchSetup& setup);
  // (the device may still work on the data)
  ~KernelBench() override { runtime_.wait(); }

  KernelBench(const KernelBench&) = delete;
  KernelBench(KernelBench&&) = delete;
  auto operator=(const KernelBench&) = delete;
  auto operator=(KernelBench&&) = delete;

  void run(Kernel kernel) override;
  void wait() override { runtime_.wait(); }
  [[nodiscard]] auto estimate(Kernel kernel) const -> PerformanceEstimate override;

  private:
  void runAder();
  void runLocalWOAder();
  void runLocal();
  void runNeighbor();
  void runGodunovDR();
  void runPlasticity();

  void runAderDevice();
  void runLocalWOAderDevice();
  void runLocalDevice();
  void runNeighborDevice();
  void runGodunovDRDevice();
  void runPlasticityDevice();

  void computeEstimates();

  bool onDevice_;
  SyntheticData<Cfg> data_;

  kernels::Spacetime<Cfg> spacetimeKernel_;
  kernels::Time<Cfg> timeKernel_;
  kernels::Local<Cfg> localKernel_;
  kernels::Neighbor<Cfg> neighborKernel_;
  kernels::DynamicRupture<Cfg> dynRupKernel_;

  parallel::runtime::StreamRuntime runtime_;
  //! the cells that yielded (on the device)
  memory::MemkindArray<std::size_t> yielded_;

  PerformanceEstimate aderEstimate_{};
  PerformanceEstimate localWOAderEstimate_{};
  PerformanceEstimate neighborEstimate_{};
  PerformanceEstimate godunovEstimate_{};
  PerformanceEstimate plasticityEstimate_{};
};

template <typename Cfg>
KernelBench<Cfg>::KernelBench(const BenchSetup& setup)
    : onDevice_(setup.executor == Executor::Device), data_(setup),
      yielded_(1, isDeviceOn() ? memory::Memkind::DeviceGlobalMemory : memory::Memkind::Standard) {
  const auto globalData = data_.globalData().template get<Cfg>();
  spacetimeKernel_.setGlobalData(globalData);
  timeKernel_.setGlobalData(globalData);
  localKernel_.setGlobalData(globalData);
  neighborKernel_.setGlobalData(globalData);
  dynRupKernel_.setGlobalData(globalData);
  computeEstimates();
}

template <typename Cfg>
void KernelBench<Cfg>::computeEstimates() {
  const LTS::Layer& layer = data_.cells();
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  const auto* drMapping = layer.var<LTS::DRMapping>(Cfg());
  for (std::size_t cell = 0; cell < layer.size(); ++cell) {
    aderEstimate_ += spacetimeKernel_.metrics();
    localWOAderEstimate_ += localKernel_.metrics(cellInformation[cell].faceTypes);
    const auto [cellRegular, cellDR] = neighborKernel_.metrics(
        cellInformation[cell].faceTypes, cellInformation[cell].faceRelations, drMapping[cell]);
    neighborEstimate_ += cellRegular + cellDR;
  }

  if (data_.hasFaces()) {
    const DynamicRupture::Layer& faces = data_.faces();
    const auto* faceInformation = faces.var<DynamicRupture::FaceInformation>();
    for (std::size_t face = 0; face < faces.size(); ++face) {
      godunovEstimate_ += dynRupKernel_.metrics(faceInformation[face]);
    }
  }

  if (data_.plasticity() != CellPlasticity::None) {
    const auto [check, yield] = kernels::Plasticity<Cfg>::metrics();
    plasticityEstimate_ = check * layer.size();
    if (data_.plasticity() == CellPlasticity::Yield) {
      plasticityEstimate_ += yield * layer.size();
    }
  }
}

template <typename Cfg>
void KernelBench<Cfg>::run(Kernel kernel) {
  switch (kernel) {
  case Kernel::Ader:
    onDevice_ ? runAderDevice() : runAder();
    break;
  case Kernel::LocalWOAder:
    onDevice_ ? runLocalWOAderDevice() : runLocalWOAder();
    break;
  case Kernel::Local:
    onDevice_ ? runLocalDevice() : runLocal();
    break;
  case Kernel::Neighbor:
  case Kernel::NeighborDR:
    onDevice_ ? runNeighborDevice() : runNeighbor();
    break;
  case Kernel::GodunovDR:
    if (!data_.hasFaces()) {
      throw std::invalid_argument("The loops of the faces of dynamic rupture need the faces.");
    }
    onDevice_ ? runGodunovDRDevice() : runGodunovDR();
    break;
  case Kernel::Plasticity:
  case Kernel::PlasticityYield:
    onDevice_ ? runPlasticityDevice() : runPlasticity();
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
auto KernelBench<Cfg>::estimate(Kernel kernel) const -> PerformanceEstimate {
  switch (kernel) {
  case Kernel::Ader:
    return aderEstimate_;
  case Kernel::LocalWOAder:
    return localWOAderEstimate_;
  case Kernel::Local:
    return aderEstimate_ + localWOAderEstimate_;
  case Kernel::Neighbor:
  case Kernel::NeighborDR:
    return neighborEstimate_;
  case Kernel::GodunovDR:
    return godunovEstimate_;
  case Kernel::Plasticity:
  case Kernel::PlasticityYield:
    return plasticityEstimate_;
  case Kernel::All:
    return estimate(Kernel::Local) + estimate(Kernel::Neighbor);
  case Kernel::AllDR:
    return estimate(Kernel::Local) + estimate(Kernel::GodunovDR) + estimate(Kernel::NeighborDR);
  }
  return {};
}

template <typename Cfg>
void KernelBench<Cfg>::runAder() {
  LTS::Layer& layer = data_.cells();
  real* const* stepIntegrals = layer.var<LTS::StepIntegrals>(Cfg());
  real* const* derivatives = layer.var<LTS::Derivatives>(Cfg());

  const auto integrationCoeffs = kernels::timeBasis<Cfg>().integrate(0, Timestep, Timestep);

#pragma omp parallel
  {
    LIKWID_MARKER_START("ader");
    kernels::LocalTmp<Cfg> tmp(GravitationalAcceleration);

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); cell++) {
      auto local = layer.cellRef<Cfg>(cell);
      spacetimeKernel_.computeAder(integrationCoeffs.data(),
                                   Timestep,
                                   local,
                                   tmp,
                                   stepIntegrals[cell],
                                   derivatives[cell],
                                   /*updateDisplacement=*/false);
    }
    LIKWID_MARKER_STOP("ader");
  }
}

template <typename Cfg>
void KernelBench<Cfg>::runLocalWOAder() {
  LTS::Layer& layer = data_.cells();
  real* const* stepIntegrals = layer.var<LTS::StepIntegrals>(Cfg());

#pragma omp parallel
  {
    LIKWID_MARKER_START("localwoader");
    kernels::LocalTmp<Cfg> tmp(GravitationalAcceleration);

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); cell++) {
      auto local = layer.cellRef<Cfg>(cell);
      localKernel_.computeIntegral(stepIntegrals[cell], local, tmp, 0, 0);
    }
    LIKWID_MARKER_STOP("localwoader");
  }
}

template <typename Cfg>
void KernelBench<Cfg>::runLocal() {
  LTS::Layer& layer = data_.cells();
  real* const* stepIntegrals = layer.var<LTS::StepIntegrals>(Cfg());
  real* const* derivatives = layer.var<LTS::Derivatives>(Cfg());

  const auto integrationCoeffs = kernels::timeBasis<Cfg>().integrate(0, Timestep, Timestep);

#pragma omp parallel
  {
    LIKWID_MARKER_START("local");
    kernels::LocalTmp<Cfg> tmp(GravitationalAcceleration);

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); cell++) {
      auto local = layer.cellRef<Cfg>(cell);
      spacetimeKernel_.computeAder(integrationCoeffs.data(),
                                   Timestep,
                                   local,
                                   tmp,
                                   stepIntegrals[cell],
                                   derivatives[cell],
                                   /*updateDisplacement=*/false);
      localKernel_.computeIntegral(stepIntegrals[cell], local, tmp, 0, 0);
    }
    LIKWID_MARKER_STOP("local");
  }
}

template <typename Cfg>
void KernelBench<Cfg>::runNeighbor() {
  LTS::Layer& layer = data_.cells();
  const auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  const auto* drMapping = layer.var<LTS::DRMapping>(Cfg());
  const auto* cellInformation = layer.var<LTS::CellInformation>();

  const auto timeCoeffs = kernels::timeBasis<Cfg>().integrate(0, Timestep, Timestep);

  // the faces are regular, GTS, or all on the fault
#pragma omp parallel
  {
    LIKWID_MARKER_START("neighboring");
    std::array<real*, Cell::NumFaces> timeIntegrated{};
    std::array<real*, Cell::NumFaces> faceNeighborsPrefetch{};

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); cell++) {
      auto local = layer.cellRef<Cfg>(cell);

      // scratch for the neighbors integrated here, written before it is read
      alignas(Alignment)
          real integrationBuffer[Cell::NumFaces][kernels::SolverOf<Cfg>::IntegralsSize];
      std::array<real*, Cell::NumFaces> integrationBuffers{};
      for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
        integrationBuffers[i] = integrationBuffer[i];
      }

      kernels::TimeCommon<Cfg>::computeIntegrals(timeKernel_,
                                                 cellInformation[cell],
                                                 timeCoeffs.data(),
                                                 timeCoeffs.data(),
                                                 faceNeighbors[cell],
                                                 integrationBuffers,
                                                 timeIntegrated);

      for (std::size_t face = 0; face + 1 < Cell::NumFaces; ++face) {
        faceNeighborsPrefetch[face] =
            (cellInformation[cell].faceTypes[face + 1] != FaceType::DynamicRupture)
                ? static_cast<real*>(faceNeighbors[cell][face + 1])
                : drMapping[cell][face + 1].godunov;
      }
      // the fourth face prefetches the first face of the next cell
      if (cell + 1 < layer.size()) {
        faceNeighborsPrefetch[3] =
            (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
                ? static_cast<real*>(faceNeighbors[cell + 1][0])
                : drMapping[cell + 1][0].godunov;
      } else {
        faceNeighborsPrefetch[3] = static_cast<real*>(faceNeighbors[cell][3]);
      }

      neighborKernel_.computeNeighborsIntegral(local, timeIntegrated, faceNeighborsPrefetch);
    }
    LIKWID_MARKER_STOP("neighboring");
  }
}

template <typename Cfg>
void KernelBench<Cfg>::runGodunovDR() {
  DynamicRupture::Layer& faces = data_.faces();
  const auto* faceInformation = faces.var<DynamicRupture::FaceInformation>();
  const auto* godunovData = faces.var<DynamicRupture::GodunovData>(Cfg());
  real* const* timeDerivativePlus = faces.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  real* const* timeDerivativeMinus = faces.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  const auto [timePoints, timeWeights] =
      quadrature::shiftedGaussLegendre(Cfg::ConvergenceOrder, 0, Timestep);
  const auto coeffsCollocate = kernels::timeBasis<Cfg>().collocate(timePoints, Timestep);

#pragma omp parallel
  {
    // the interpolation of a face, which nothing reads here: one copy per thread
    alignas(Alignment)
        real qInterpolatedPlus[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()];
    alignas(Alignment)
        real qInterpolatedMinus[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()];

#pragma omp for schedule(static)
    for (std::size_t face = 0; face < faces.size(); ++face) {
      const std::size_t prefetchFace = (face + 1 < faces.size()) ? face + 1 : face;
      dynRupKernel_.spaceTimeInterpolation(faceInformation[face],
                                           &godunovData[face],
                                           timeDerivativePlus[face],
                                           timeDerivativeMinus[face],
                                           qInterpolatedPlus,
                                           qInterpolatedMinus,
                                           timeDerivativePlus[prefetchFace],
                                           timeDerivativeMinus[prefetchFace],
                                           coeffsCollocate.data());
    }
  }
}

template <typename Cfg>
void KernelBench<Cfg>::runPlasticity() {
  LTS::Layer& layer = data_.cells();
  const auto* plasticity = layer.var<LTSPlasticity::Plasticity>(Cfg());
  auto* pstrain = layer.var<LTSPlasticity::PStrain>(Cfg());
  real(*dofs)[tensor::Q<Cfg>::size()] = layer.var<LTS::Dofs>(Cfg());
  const auto* global = &data_.globalData().template onHost<Cfg>();
  const auto oneMinusIntegratingFactor =
      static_cast<real>(kernels::Plasticity<Cfg>::computeRelaxTime(RelaxationTime, Timestep));

#pragma omp parallel
  {
    LIKWID_MARKER_START("plasticity");
#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      kernels::Plasticity<Cfg>::computePlasticity(
          oneMinusIntegratingFactor, global, &plasticity[cell], dofs[cell], pstrain[cell]);
    }
    LIKWID_MARKER_STOP("plasticity");
  }
}

// On the device, every kernel records a graph of its own: the graphs of a layer are kept by their
// keys, which therefore differ from kernel to kernel here.

template <typename Cfg>
void KernelBench<Cfg>::runAderDevice() {
#ifdef ACL_DEVICE
  LTS::Layer& layer = data_.cells();
  auto& dataTable = layer.getConditionalTable();
  kernels::LocalTmp<Cfg> tmp(GravitationalAcceleration);
  const auto integrationCoeffs = kernels::timeBasis<Cfg>().integrate(0, Timestep, Timestep);

  const auto key = initializer::GraphKey(ComputeGraphType::AccumulatedVelocities, Timestep, false);
  runtime_.runGraph(key, layer, [&](auto& runtime) {
    spacetimeKernel_.computeBatchedAder(
        integrationCoeffs.data(), Timestep, layer, tmp, dataTable, false, runtime);
  });
#endif
}

template <typename Cfg>
void KernelBench<Cfg>::runLocalWOAderDevice() {
#ifdef ACL_DEVICE
  LTS::Layer& layer = data_.cells();
  auto& dataTable = layer.getConditionalTable();

  const auto key = initializer::GraphKey(ComputeGraphType::StreamedVelocities, Timestep, false);
  runtime_.runGraph(key, layer, [&](auto& runtime) {
    localKernel_.computeBatchedIntegral(dataTable, Timestep, runtime);
  });
#endif
}

template <typename Cfg>
void KernelBench<Cfg>::runLocalDevice() {
#ifdef ACL_DEVICE
  LTS::Layer& layer = data_.cells();
  auto& dataTable = layer.getConditionalTable();
  kernels::LocalTmp<Cfg> tmp(GravitationalAcceleration);
  const auto integrationCoeffs = kernels::timeBasis<Cfg>().integrate(0, Timestep, Timestep);

  const auto key = initializer::GraphKey(ComputeGraphType::AccumulatedVelocities, Timestep, true);
  runtime_.runGraph(key, layer, [&](auto& runtime) {
    spacetimeKernel_.computeBatchedAder(
        integrationCoeffs.data(), Timestep, layer, tmp, dataTable, false, runtime);
    localKernel_.computeBatchedIntegral(dataTable, Timestep, runtime);
  });
#endif
}

template <typename Cfg>
void KernelBench<Cfg>::runNeighborDevice() {
#ifdef ACL_DEVICE
  LTS::Layer& layer = data_.cells();
  auto& dataTable = layer.getConditionalTable();
  const auto timeCoeffs = kernels::timeBasis<Cfg>().integrate(0, Timestep, Timestep);

  kernels::TimeCommon<Cfg>::computeBatchedIntegrals(
      timeKernel_, timeCoeffs.data(), timeCoeffs.data(), dataTable, runtime_);

  const auto key = initializer::GraphKey(ComputeGraphType::NeighborIntegral);
  runtime_.runGraph(key, layer, [&](auto& runtime) {
    neighborKernel_.computeBatchedNeighborsIntegral(dataTable, runtime);
  });
#endif
}

template <typename Cfg>
void KernelBench<Cfg>::runGodunovDRDevice() {
#ifdef ACL_DEVICE
  DynamicRupture::Layer& faces = data_.faces();
  auto& dataTable = faces.getConditionalTable();
  const auto [timePoints, timeWeights] =
      quadrature::shiftedGaussLegendre(Cfg::ConvergenceOrder, 0, Timestep);
  const auto coeffsCollocate = kernels::timeBasis<Cfg>().collocate(timePoints, Timestep);

  const auto key = initializer::GraphKey(ComputeGraphType::DynamicRuptureInterface, 0.0);
  runtime_.runGraph(key, faces, [&](auto& runtime) {
    dynRupKernel_.batchedSpaceTimeInterpolation(dataTable, coeffsCollocate.data(), runtime);
  });
#endif
}

template <typename Cfg>
void KernelBench<Cfg>::runPlasticityDevice() {
#ifdef ACL_DEVICE
  LTS::Layer& layer = data_.cells();
  auto& dataTable = layer.getConditionalTable();
  auto* plasticity =
      layer.var<LTSPlasticity::Plasticity>(Cfg(), initializer::AllocationPlace::Device);
  auto* isAdjustableVector =
      layer.var<LTSPlasticity::FlagScratch>(initializer::AllocationPlace::Device);
  const auto* global = &data_.globalData().template onDevice<Cfg>();

  const auto key = initializer::GraphKey(ComputeGraphType::Plasticity, Timestep);
  runtime_.runGraph(key, layer, [&](auto& runtime) {
    kernels::Plasticity<Cfg>::computePlasticityBatched(Timestep,
                                                       RelaxationTime,
                                                       global,
                                                       dataTable,
                                                       plasticity,
                                                       yielded_.data(),
                                                       isAdjustableVector,
                                                       runtime);
  });
#endif
}

} // namespace

auto makeKernelBench(const BenchSetup& setup) -> std::unique_ptr<LoopBench> {
  return dispatchConfig(setup.config, [&](auto cfg) -> std::unique_ptr<LoopBench> {
    using Cfg = decltype(cfg);
    return std::make_unique<KernelBench<Cfg>>(setup);
  });
}

} // namespace seissol::proxy
