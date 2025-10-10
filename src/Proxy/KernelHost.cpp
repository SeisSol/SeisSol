// SPDX-FileCopyrightText: 2013 SeisSol Group
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

#include "GeneratedCode/tensor.h"
#include "Parallel/OpenMP.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <Config.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/CellLocalInformation.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Kernels/TimeCommon.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <Monitoring/Instrumentation.h>
#include <Numerical/Quadrature.h>
#include <Parallel/Runtime/Stream.h>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <variant>

namespace seissol::proxy {
template <typename Cfg>
void ProxyKernelHostAder<Cfg>::run(ProxyData& predata,
                                   seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  Real<Cfg>* const* buffers = layer.template var<LTS::Buffers>(Cfg());
  Real<Cfg>* const* derivatives = layer.template var<LTS::Derivatives>(Cfg());

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    LIKWID_MARKER_START("ader");
    kernels::LocalTmp<Cfg> tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.template cellRef<Cfg>(cell);
      data.spacetimeKernel.computeAder(
          integrationCoeffs.data(), Timestep, local, tmp, buffers[cell], derivatives[cell]);
    }
    LIKWID_MARKER_STOP("ader");
  }
}

template <typename Cfg>
auto ProxyKernelHostAder<Cfg>::performanceEstimate(ProxyData& predata) const
    -> PerformanceEstimate {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
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

template <typename Cfg>
auto ProxyKernelHostAder<Cfg>::needsDR() const -> bool {
  return false;
}

template <typename Cfg>
void ProxyKernelHostLocalWOAder<Cfg>::run(
    ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  Real<Cfg>* const* buffers = layer.template var<LTS::Buffers>(Cfg());

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    LIKWID_MARKER_START("localwoader");
    kernels::LocalTmp<Cfg> tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.template cellRef<Cfg>(cell);
      data.localKernel.computeIntegral(buffers[cell], local, tmp, nullptr, nullptr, 0, 0);
    }
    LIKWID_MARKER_STOP("localwoader");
  }
}

template <typename Cfg>
auto ProxyKernelHostLocalWOAder<Cfg>::performanceEstimate(ProxyData& predata) const
    -> PerformanceEstimate {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  const auto* cellInformation = layer.template var<LTS::CellInformation>();
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

template <typename Cfg>
auto ProxyKernelHostLocalWOAder<Cfg>::needsDR() const -> bool {
  return false;
}

template <typename Cfg>
void ProxyKernelHostLocal<Cfg>::run(ProxyData& predata,
                                    seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  Real<Cfg>* const* buffers = layer.template var<LTS::Buffers>(Cfg());
  Real<Cfg>* const* derivatives = layer.template var<LTS::Derivatives>(Cfg());

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    LIKWID_MARKER_START("local");
    kernels::LocalTmp<Cfg> tmp(9.81);
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.template cellRef<Cfg>(cell);
      data.spacetimeKernel.computeAder(
          integrationCoeffs.data(), Timestep, local, tmp, buffers[cell], derivatives[cell]);
      data.localKernel.computeIntegral(buffers[cell], local, tmp, nullptr, nullptr, 0, 0);
    }
    LIKWID_MARKER_STOP("local");
  }
}

template <typename Cfg>
void ProxyKernelHostNeighbor<Cfg>::run(ProxyData& predata,
                                       seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  void* const(*faceNeighbors)[4] = layer.template var<LTS::FaceNeighbors>();
  const CellDRMapping<Cfg>(*drMapping)[4] = layer.template var<LTS::DRMapping>(Cfg());
  const CellLocalInformation* cellInformation = layer.template var<LTS::CellInformation>();

  Real<Cfg>* timeIntegrated[4];
  Real<Cfg>* faceNeighborsPrefetch[4];

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
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
      auto local = layer.template cellRef<Cfg>(cell);
      auto* const* faceNeighborsCell = reinterpret_cast<Real<Cfg>* const*>(faceNeighbors[cell]);
      auto* const* faceNeighborsCell1 = reinterpret_cast<Real<Cfg>* const*>(faceNeighbors[cell]);

      seissol::kernels::TimeCommon<Cfg>::computeIntegrals(
          data.timeKernel,
          cellInformation[cell],
          timeCoeffs.data(),
          timeCoeffs.data(),
          faceNeighbors[cell],
          *reinterpret_cast<Real<Cfg>(*)[4][tensor::I<Cfg>::size()]>(
              &data.globalData.template get<Cfg>()
                   .integrationBufferLTS[OpenMP::threadId() *
                                         static_cast<size_t>(tensor::I<Cfg>::size()) * 4]),
          timeIntegrated);

      faceNeighborsPrefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
                                     ? faceNeighborsCell[1]
                                     : drMapping[cell][1].godunov;
      faceNeighborsPrefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
                                     ? faceNeighborsCell[2]
                                     : drMapping[cell][2].godunov;
      faceNeighborsPrefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
                                     ? faceNeighborsCell[3]
                                     : drMapping[cell][3].godunov;

      // fourth face's prefetches
      if (cell + 1 < nrOfCells) {
        faceNeighborsPrefetch[3] =
            (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
                ? faceNeighborsCell1[0]
                : drMapping[cell + 1][0].godunov;
      } else {
        faceNeighborsPrefetch[3] = faceNeighborsCell[3];
      }

      data.neighborKernel.computeNeighborsIntegral(
          local, drMapping[cell], timeIntegrated, faceNeighborsPrefetch);
    }

    LIKWID_MARKER_STOP("neighboring");
  }
}

template <typename Cfg>
auto ProxyKernelHostNeighbor<Cfg>::performanceEstimate(ProxyData& predata) const
    -> PerformanceEstimate {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  // iterate over cells
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  const CellLocalInformation* cellInformation = layer.template var<LTS::CellInformation>();
  const CellDRMapping<Cfg>(*drMapping)[4] = layer.template var<LTS::DRMapping>(Cfg());
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

template <typename Cfg>
auto ProxyKernelHostNeighbor<Cfg>::needsDR() const -> bool {
  return false;
}

template <typename Cfg>
auto ProxyKernelHostNeighborDR<Cfg>::needsDR() const -> bool {
  return true;
}

template <typename Cfg>
void ProxyKernelHostGodunovDR<Cfg>::run(ProxyData& predata,
                                        seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layerData = data.drStorage.layer(data.layerId);
  const auto* faceInformation = layerData.template var<DynamicRupture::FaceInformation>();
  const auto* godunovData = layerData.template var<DynamicRupture::GodunovData>(Cfg());
  auto* drEnergyOutput = layerData.template var<DynamicRupture::DREnergyOutputVar>(Cfg());
  auto* const* timeDerivativePlus =
      layerData.template var<DynamicRupture::TimeDerivativePlus>(Cfg());
  auto* const* timeDerivativeMinus =
      layerData.template var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  alignas(Alignment) Real<Cfg> qInterpolatedPlus[dr::misc::TimeSteps<Cfg>]
                                                [tensor::QInterpolated<Cfg>::size()];
  alignas(Alignment) Real<Cfg> qInterpolatedMinus[dr::misc::TimeSteps<Cfg>]
                                                 [tensor::QInterpolated<Cfg>::size()];
  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, Timestep);
  const auto coeffsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, Timestep);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(qInterpolatedPlus, qInterpolatedMinus)
#endif
  for (std::size_t face = 0; face < layerData.size(); ++face) {
    const std::size_t prefetchFace = (face + 1 < layerData.size()) ? face + 1 : face;
    data.dynRupKernel.spaceTimeInterpolation(faceInformation[face],
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

template <typename Cfg>
auto ProxyKernelHostGodunovDR<Cfg>::performanceEstimate(ProxyData& predata) const
    -> PerformanceEstimate {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  PerformanceEstimate ret;
  ret.nonzeroFlop = 0.0;
  ret.hardwareFlop = 0.0;

  // iterate over cells
  auto& interior = data.drStorage.layer(data.layerId);
  const DRFaceInformation* faceInformation =
      interior.template var<DynamicRupture::FaceInformation>();
  for (std::size_t face = 0; face < interior.size(); ++face) {
    std::uint64_t drNonZeroFlops = 0;
    std::uint64_t drHardwareFlops = 0;
    data.dynRupKernel.flopsGodunovState(faceInformation[face], drNonZeroFlops, drHardwareFlops);
    ret.nonzeroFlop += drNonZeroFlops;
    ret.hardwareFlop += drHardwareFlops;
  }

  return ret;
}

template <typename Cfg>
auto ProxyKernelHostGodunovDR<Cfg>::needsDR() const -> bool {
  return true;
}

std::shared_ptr<ProxyKernel> getProxyKernelHost(Kernel kernel, ConfigVariant variant) {
  return std::visit(
      [&](auto cfg) -> std::shared_ptr<ProxyKernel> {
        using Cfg = decltype(cfg);
        switch (kernel) {
        case Kernel::All:
          return std::make_shared<ProxyKernelHostAll<Cfg>>();
        case Kernel::AllDR:
          return std::make_shared<ProxyKernelHostAllDR<Cfg>>();
        case Kernel::Ader:
          return std::make_shared<ProxyKernelHostAder<Cfg>>();
        case Kernel::LocalWOAder:
          return std::make_shared<ProxyKernelHostLocalWOAder<Cfg>>();
        case Kernel::Local:
          return std::make_shared<ProxyKernelHostLocal<Cfg>>();
        case Kernel::Neighbor:
          return std::make_shared<ProxyKernelHostNeighbor<Cfg>>();
        case Kernel::NeighborDR:
          return std::make_shared<ProxyKernelHostNeighborDR<Cfg>>();
        case Kernel::GodunovDR:
          return std::make_shared<ProxyKernelHostGodunovDR<Cfg>>();
        }
        throw;
      },
      variant);
}

} // namespace seissol::proxy
