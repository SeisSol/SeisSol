// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "KernelHost.h"

#include "Alignment.h"
#include "Allocator.h"
#include "Common.h"
#include "Common/Constants.h"
#include "Constants.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/Typedefs.h"
#include "Kernel.h"
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/Instrumentation.h"
#include "Numerical/Quadrature.h"
#include "Parallel/OpenMP.h"
#include "Parallel/Runtime/Stream.h"

#include <cstddef>
#include <cstdint>
#include <memory>

namespace seissol::proxy {
void ProxyKernelHostAder::run(ProxyData& data,
                              seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  real* const* buffers = layer.var<LTS::Buffers>();
  real* const* derivatives = layer.var<LTS::Derivatives>();

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    LIKWID_MARKER_START("ader");
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);
      data.spacetimeKernel.computeAder(
          integrationCoeffs.data(), Timestep, local, tmp, buffers[cell], derivatives[cell]);
    }
    LIKWID_MARKER_STOP("ader");
  }
}
auto ProxyKernelHostAder::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0;
  ret.hardwareFlop = 0;

  // iterate over cells
  const auto nrOfCells = data.ltsStorage.layer(data.layerId).size();
  for (std::size_t cell = 0; cell < nrOfCells; ++cell) {
    std::uint64_t nonZeroFlops = 0;
    std::uint64_t hardwareFlops = 0;
    // get flops
    data.spacetimeKernel.flopsAder(nonZeroFlops, hardwareFlops);
    ret.nonzeroFlop += nonZeroFlops;
    ret.hardwareFlop += hardwareFlops;
  }

  ret.bytes = static_cast<std::size_t>(data.spacetimeKernel.bytesAder() * nrOfCells);

  return ret;
}
auto ProxyKernelHostAder::needsDR() const -> bool { return false; }

void ProxyKernelHostLocalWOAder::run(ProxyData& data,
                                     seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  real* const* buffers = layer.var<LTS::Buffers>();

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    LIKWID_MARKER_START("localwoader");
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);
      data.localKernel.computeIntegral(buffers[cell], local, tmp, nullptr, nullptr, 0, 0);
    }
    LIKWID_MARKER_STOP("localwoader");
  }
}
auto ProxyKernelHostLocalWOAder::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  const auto* cellInformation = layer.var<LTS::CellInformation>();
  for (std::size_t cell = 0; cell < nrOfCells; ++cell) {
    std::uint64_t nonZeroFlops = 0;
    std::uint64_t hardwareFlops = 0;
    data.localKernel.flopsIntegral(cellInformation[cell].faceTypes, nonZeroFlops, hardwareFlops);
    ret.nonzeroFlop += nonZeroFlops;
    ret.hardwareFlop += hardwareFlops;
  }

  const auto bytes = data.localKernel.bytesIntegral();

  ret.bytes = static_cast<std::size_t>(nrOfCells * bytes);

  return ret;
}
auto ProxyKernelHostLocalWOAder::needsDR() const -> bool { return false; }

void ProxyKernelHostLocal::run(ProxyData& data,
                               seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  real* const* buffers = layer.var<LTS::Buffers>();
  real* const* derivatives = layer.var<LTS::Derivatives>();

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    LIKWID_MARKER_START("local");
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);
      data.spacetimeKernel.computeAder(
          integrationCoeffs.data(), Timestep, local, tmp, buffers[cell], derivatives[cell]);
      data.localKernel.computeIntegral(buffers[cell], local, tmp, nullptr, nullptr, 0, 0);
    }
    LIKWID_MARKER_STOP("local");
  }
}

void ProxyKernelHostNeighbor::run(ProxyData& data,
                                  seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  real* const(*faceNeighbors)[4] = layer.var<LTS::FaceNeighbors>();
  const CellDRMapping(*drMapping)[4] = layer.var<LTS::DRMapping>();
  const CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();

  real* timeIntegrated[4];
  real* faceNeighborsPrefetch[4];

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, Timestep, Timestep);

  // note: we use GTS here, in all cases

#ifdef _OPENMP
#pragma omp parallel private(timeIntegrated, faceNeighborsPrefetch)
#endif
  {
    LIKWID_MARKER_START("neighboring");
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);
      seissol::kernels::TimeCommon::computeIntegrals(
          data.timeKernel,
          cellInformation[cell].ltsSetup,
          cellInformation[cell].faceTypes,
          timeCoeffs.data(),
          timeCoeffs.data(),
          faceNeighbors[cell],
          *reinterpret_cast<real(*)[4][tensor::I::size()]>(
              &data.globalDataOnHost
                   .integrationBufferLTS[OpenMP::threadId() *
                                         static_cast<size_t>(tensor::I::size()) * 4]),
          timeIntegrated);

      faceNeighborsPrefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
                                     ? faceNeighbors[cell][1]
                                     : drMapping[cell][1].godunov;
      faceNeighborsPrefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
                                     ? faceNeighbors[cell][2]
                                     : drMapping[cell][2].godunov;
      faceNeighborsPrefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
                                     ? faceNeighbors[cell][3]
                                     : drMapping[cell][3].godunov;

      // fourth face's prefetches
      if (cell + 1 < nrOfCells) {
        faceNeighborsPrefetch[3] =
            (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
                ? faceNeighbors[cell + 1][0]
                : drMapping[cell + 1][0].godunov;
      } else {
        faceNeighborsPrefetch[3] = faceNeighbors[cell][3];
      }

      data.neighborKernel.computeNeighborsIntegral(
          local, drMapping[cell], timeIntegrated, faceNeighborsPrefetch);
    }

    LIKWID_MARKER_STOP("neighboring");
  }
}
auto ProxyKernelHostNeighbor::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  // iterate over cells
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  const CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();
  const CellDRMapping(*drMapping)[4] = layer.var<LTS::DRMapping>();
  for (std::size_t cell = 0; cell < nrOfCells; cell++) {
    std::uint64_t nonZeroFlops = 0;
    std::uint64_t hardwareFlops = 0;
    std::uint64_t drNonZeroFlops = 0;
    std::uint64_t drHardwareFlops = 0;
    // get flops
    data.neighborKernel.flopsNeighborsIntegral(cellInformation[cell].faceTypes,
                                               cellInformation[cell].faceRelations,
                                               drMapping[cell],
                                               nonZeroFlops,
                                               hardwareFlops,
                                               drNonZeroFlops,
                                               drHardwareFlops);
    ret.nonzeroFlop += nonZeroFlops + drNonZeroFlops;
    ret.hardwareFlop += hardwareFlops + drHardwareFlops;
  }

  ret.bytes = static_cast<std::size_t>(data.neighborKernel.bytesNeighborsIntegral() * nrOfCells);

  return ret;
}
auto ProxyKernelHostNeighbor::needsDR() const -> bool { return false; }

auto ProxyKernelHostNeighborDR::needsDR() const -> bool { return true; }

void ProxyKernelHostGodunovDR::run(ProxyData& data,
                                   seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layerData = data.drStorage.layer(data.layerId);
  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();
  const DRGodunovData* godunovData = layerData.var<DynamicRupture::GodunovData>();
  DREnergyOutput* drEnergyOutput = layerData.var<DynamicRupture::DREnergyOutputVar>();
  real* const* timeDerivativePlus = layerData.var<DynamicRupture::TimeDerivativePlus>();
  real* const* timeDerivativeMinus = layerData.var<DynamicRupture::TimeDerivativeMinus>();
  alignas(Alignment) real qInterpolatedPlus[ConvergenceOrder][tensor::QInterpolated::size()];
  alignas(Alignment) real qInterpolatedMinus[ConvergenceOrder][tensor::QInterpolated::size()];
  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, Timestep);
  const auto coeffsCollocate = seissol::kernels::timeBasis().collocate(timePoints, Timestep);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(qInterpolatedPlus, qInterpolatedMinus)
#endif
  for (std::size_t face = 0; face < layerData.size(); ++face) {
    const std::size_t prefetchFace = (face + 1 < layerData.size()) ? face + 1 : face;
    data.dynRupKernel.spaceTimeInterpolation(faceInformation[face],
                                             &data.globalDataOnHost,
                                             &godunovData[face],
                                             &drEnergyOutput[face],
                                             timeDerivativePlus[face],
                                             timeDerivativeMinus[face],
                                             qInterpolatedPlus,
                                             qInterpolatedMinus,
                                             timeDerivativePlus[prefetchFace],
                                             timeDerivativeMinus[prefetchFace],
                                             coeffsCollocate.data());
  }
}
auto ProxyKernelHostGodunovDR::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  // iterate over cells
  auto& interior = data.drStorage.layer(data.layerId);
  const DRFaceInformation* faceInformation = interior.var<DynamicRupture::FaceInformation>();
  for (std::size_t face = 0; face < interior.size(); ++face) {
    std::uint64_t drNonZeroFlops = 0;
    std::uint64_t drHardwareFlops = 0;
    data.dynRupKernel.flopsGodunovState(faceInformation[face], drNonZeroFlops, drHardwareFlops);
    ret.nonzeroFlop += drNonZeroFlops;
    ret.hardwareFlop += drHardwareFlops;
  }

  return ret;
}
auto ProxyKernelHostGodunovDR::needsDR() const -> bool { return true; }

std::shared_ptr<ProxyKernel> getProxyKernelHost(Kernel kernel) {
  switch (kernel) {
  case Kernel::All:
    return std::make_shared<ProxyKernelHostAll>();
  case Kernel::AllDR:
    return std::make_shared<ProxyKernelHostAllDR>();
  case Kernel::Ader:
    return std::make_shared<ProxyKernelHostAder>();
  case Kernel::LocalWOAder:
    return std::make_shared<ProxyKernelHostLocalWOAder>();
  case Kernel::Local:
    return std::make_shared<ProxyKernelHostLocal>();
  case Kernel::Neighbor:
    return std::make_shared<ProxyKernelHostNeighbor>();
  case Kernel::NeighborDR:
    return std::make_shared<ProxyKernelHostNeighborDR>();
  case Kernel::GodunovDR:
    return std::make_shared<ProxyKernelHostGodunovDR>();
  }
  throw;
}

} // namespace seissol::proxy
