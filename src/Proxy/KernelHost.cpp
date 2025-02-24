// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "KernelHost.h"
#include "Allocator.h"
#include "Common.h"
#include "Constants.h"
#include "Kernel.h"

#include <Common/Constants.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/CellLocalInformation.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Kernels/TimeCommon.h>
#include <Memory/Tree/Layer.h>
#include <Monitoring/Instrumentation.h>
#include <Parallel/Runtime/Stream.h>
#include <memory>
#include <omp.h>
#include <tensor.h>

namespace seissol::proxy {
void ProxyKernelHostAder::run(ProxyData& data,
                              seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsTree.child(0).child<Interior>();
  const unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(data.lts.buffers);
  real** derivatives = layer.var(data.lts.derivatives);

  kernels::LocalData::Loader loader;
  loader.load(data.lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    LIKWID_MARKER_START("ader");
#endif
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned int cell = 0; cell < nrOfCells; cell++) {
      auto local = loader.entry(cell);
      data.timeKernel.computeAder(Timestep, local, tmp, buffers[cell], derivatives[cell]);
    }
#ifdef _OPENMP
    LIKWID_MARKER_STOP("ader");
  }
#endif
}
auto ProxyKernelHostAder::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0;
  ret.hardwareFlop = 0;

  // iterate over cells
  const unsigned nrOfCells = data.ltsTree.child(0).child<Interior>().getNumberOfCells();
  for (unsigned int cell = 0; cell < nrOfCells; ++cell) {
    unsigned int nonZeroFlops = 0;
    unsigned int hardwareFlops = 0;
    // get flops
    data.timeKernel.flopsAder(nonZeroFlops, hardwareFlops);
    ret.nonzeroFlop += nonZeroFlops;
    ret.hardwareFlop += hardwareFlops;
  }

  ret.bytes = data.timeKernel.bytesAder() * nrOfCells;

  return ret;
}
auto ProxyKernelHostAder::needsDR() const -> bool { return false; }

void ProxyKernelHostLocalWOAder::run(ProxyData& data,
                                     seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsTree.child(0).child<Interior>();
  const unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(data.lts.buffers);

  kernels::LocalData::Loader loader;
  loader.load(data.lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    LIKWID_MARKER_START("localwoader");
#endif
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned int cell = 0; cell < nrOfCells; cell++) {
      auto local = loader.entry(cell);
      data.localKernel.computeIntegral(buffers[cell], local, tmp, nullptr, nullptr, 0, 0);
    }
#ifdef _OPENMP
    LIKWID_MARKER_STOP("localwoader");
  }
#endif
}
auto ProxyKernelHostLocalWOAder::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  auto& layer = data.ltsTree.child(0).child<Interior>();
  const unsigned nrOfCells = layer.getNumberOfCells();
  CellLocalInformation* cellInformation = layer.var(data.lts.cellInformation);
  for (unsigned cell = 0; cell < nrOfCells; ++cell) {
    unsigned int nonZeroFlops = 0;
    unsigned int hardwareFlops = 0;
    data.localKernel.flopsIntegral(cellInformation[cell].faceTypes, nonZeroFlops, hardwareFlops);
    ret.nonzeroFlop += nonZeroFlops;
    ret.hardwareFlop += hardwareFlops;
  }

  const auto bytes = data.localKernel.bytesIntegral();

  ret.bytes = nrOfCells * bytes;

  return ret;
}
auto ProxyKernelHostLocalWOAder::needsDR() const -> bool { return false; }

void ProxyKernelHostLocal::run(ProxyData& data,
                               seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsTree.child(0).child<Interior>();
  const unsigned nrOfCells = layer.getNumberOfCells();
  real** buffers = layer.var(data.lts.buffers);
  real** derivatives = layer.var(data.lts.derivatives);

  kernels::LocalData::Loader loader;
  loader.load(data.lts, layer);

#ifdef _OPENMP
#pragma omp parallel
  {
    LIKWID_MARKER_START("local");
#endif
    kernels::LocalTmp tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (unsigned int cell = 0; cell < nrOfCells; cell++) {
      auto local = loader.entry(cell);
      data.timeKernel.computeAder(Timestep, local, tmp, buffers[cell], derivatives[cell]);
      data.localKernel.computeIntegral(buffers[cell], local, tmp, nullptr, nullptr, 0, 0);
    }
#ifdef _OPENMP
    LIKWID_MARKER_STOP("local");
  }
#endif
}

void ProxyKernelHostNeighbor::run(ProxyData& data,
                                  seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsTree.child(0).child<Interior>();
  const unsigned nrOfCells = layer.getNumberOfCells();
  real*(*faceNeighbors)[4] = layer.var(data.lts.faceNeighbors);
  CellDRMapping(*drMapping)[4] = layer.var(data.lts.drMapping);
  CellLocalInformation* cellInformation = layer.var(data.lts.cellInformation);

  kernels::NeighborData::Loader loader;
  loader.load(data.lts, layer);

  real* timeIntegrated[4];
  real* faceNeighborsPrefetch[4];

#ifdef _OPENMP
#pragma omp parallel private(timeIntegrated, faceNeighborsPrefetch)
  {
    LIKWID_MARKER_START("neighboring");
#pragma omp for schedule(static)
#endif
    for (unsigned cell = 0; cell < nrOfCells; cell++) {
      auto local = loader.entry(cell);
      seissol::kernels::TimeCommon::computeIntegrals(
          data.timeKernel,
          cellInformation[cell].ltsSetup,
          cellInformation[cell].faceTypes,
          0.0,
          Timestep,
          faceNeighbors[cell],
#ifdef _OPENMP
          *reinterpret_cast<real(*)[4][tensor::I::size()]>(
              &(data.globalDataOnHost
                    .integrationBufferLTS[omp_get_thread_num() * 4 * tensor::I::size()])),
#else
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(
            data.globalDataOnHost.integrationBufferLTS),
#endif
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
      if (cell < (nrOfCells - 1)) {
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

#ifdef _OPENMP
    LIKWID_MARKER_STOP("neighboring");
  }
#endif
}
auto ProxyKernelHostNeighbor::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  // iterate over cells
  auto& layer = data.ltsTree.child(0).child<Interior>();
  const unsigned nrOfCells = layer.getNumberOfCells();
  CellLocalInformation* cellInformation = layer.var(data.lts.cellInformation);
  CellDRMapping(*drMapping)[4] = layer.var(data.lts.drMapping);
  for (unsigned int cell = 0; cell < nrOfCells; cell++) {
    unsigned int nonZeroFlops = 0;
    unsigned int hardwareFlops = 0;
    long long drNonZeroFlops = 0;
    long long drHardwareFlops = 0;
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

  ret.bytes = data.neighborKernel.bytesNeighborsIntegral() * nrOfCells;

  return ret;
}
auto ProxyKernelHostNeighbor::needsDR() const -> bool { return false; }

auto ProxyKernelHostNeighborDR::needsDR() const -> bool { return true; }

void ProxyKernelHostGodunovDR::run(ProxyData& data,
                                   seissol::parallel::runtime::StreamRuntime& runtime) const {
  seissol::initializer::Layer& layerData = data.dynRupTree.child(0).child<Interior>();
  DRFaceInformation* faceInformation = layerData.var(data.dynRup.faceInformation);
  DRGodunovData* godunovData = layerData.var(data.dynRup.godunovData);
  DREnergyOutput* drEnergyOutput = layerData.var(data.dynRup.drEnergyOutput);
  real** timeDerivativePlus = layerData.var(data.dynRup.timeDerivativePlus);
  real** timeDerivativeMinus = layerData.var(data.dynRup.timeDerivativeMinus);
  alignas(Alignment) real qInterpolatedPlus[ConvergenceOrder][tensor::QInterpolated::size()];
  alignas(Alignment) real qInterpolatedMinus[ConvergenceOrder][tensor::QInterpolated::size()];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(qInterpolatedPlus, qInterpolatedMinus)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    const unsigned prefetchFace = (face < layerData.getNumberOfCells() - 1) ? face + 1 : face;
    data.dynRupKernel.spaceTimeInterpolation(faceInformation[face],
                                             &data.globalDataOnHost,
                                             &godunovData[face],
                                             &drEnergyOutput[face],
                                             timeDerivativePlus[face],
                                             timeDerivativeMinus[face],
                                             qInterpolatedPlus,
                                             qInterpolatedMinus,
                                             timeDerivativePlus[prefetchFace],
                                             timeDerivativeMinus[prefetchFace]);
  }
}
auto ProxyKernelHostGodunovDR::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  // iterate over cells
  seissol::initializer::Layer& interior = data.dynRupTree.child(0).child<Interior>();
  DRFaceInformation* faceInformation = interior.var(data.dynRup.faceInformation);
  for (unsigned face = 0; face < interior.getNumberOfCells(); ++face) {
    long long drNonZeroFlops = 0;
    long long drHardwareFlops = 0;
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
