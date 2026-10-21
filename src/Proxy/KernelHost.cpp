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
#include "Monitoring/Metric.h"
#include "Numerical/Quadrature.h"
#include "Parallel/Runtime/Stream.h"

#include <array>
#include <cstddef>
#include <memory>

namespace seissol::proxy {
void ProxyKernelHostAder::run(ProxyData& data,
                              seissol::parallel::runtime::StreamRuntime& /*runtime*/) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

#pragma omp parallel
  {
    LIKWID_MARKER_START("ader");
    kernels::LocalTmp tmp(9.81);

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);
      data.spacetimeKernel.computeAder(integrationCoeffs.data(),
                                       Timestep,
                                       local,
                                       tmp,
                                       data.integrals + cell * kernels::Solver::IntegralsSize);
    }
    LIKWID_MARKER_STOP("ader");
  }
}
auto ProxyKernelHostAder::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;

  const auto nrOfCells = data.ltsStorage.layer(data.layerId).size();
  for (std::size_t cell = 0; cell < nrOfCells; ++cell) {
    ret += data.spacetimeKernel.metrics();
  }

  return ret;
}
auto ProxyKernelHostAder::needsDR() const -> bool { return false; }

void ProxyKernelHostLocalWOAder::run(ProxyData& data,
                                     seissol::parallel::runtime::StreamRuntime& /*runtime*/) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();

#pragma omp parallel
  {
    LIKWID_MARKER_START("localwoader");
    kernels::LocalTmp tmp(9.81);

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);
      data.localKernel.computeIntegral(
          data.integrals + cell * kernels::Solver::IntegralsSize, local, tmp, 0, 0);
    }
    LIKWID_MARKER_STOP("localwoader");
  }
}
auto ProxyKernelHostLocalWOAder::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;

  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  const auto* cellInformation = layer.var<LTS::CellInformation>();

  for (std::size_t cell = 0; cell < nrOfCells; ++cell) {
    ret += data.localKernel.metrics(cellInformation[cell].faceTypes);
  }

  return ret;
}
auto ProxyKernelHostLocalWOAder::needsDR() const -> bool { return false; }

void ProxyKernelHostLocal::run(ProxyData& data,
                               seissol::parallel::runtime::StreamRuntime& /*runtime*/) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

#pragma omp parallel
  {
    LIKWID_MARKER_START("local");
    kernels::LocalTmp tmp(9.81);

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);
      real* integral = data.integrals + cell * kernels::Solver::IntegralsSize;
      data.spacetimeKernel.computeAder(integrationCoeffs.data(), Timestep, local, tmp, integral);
      data.localKernel.computeIntegral(integral, local, tmp, 0, 0);
    }
    LIKWID_MARKER_STOP("local");
  }
}

void ProxyKernelHostNeighbor::run(ProxyData& data,
                                  seissol::parallel::runtime::StreamRuntime& /*runtime*/) const {
  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  const auto* faceNeighbors = layer.var<LTS::FaceNeighbors>();
  const CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();
  const real* zeroTrace = layer.var<LTS::ZeroFaceTrace>();

  std::array<const real*, Cell::NumFaces> neighborTraces{};

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, Timestep, Timestep);

  // note: we use GTS here, in all cases

#pragma omp parallel private(neighborTraces)
  {
    LIKWID_MARKER_START("neighboring");

#pragma omp for schedule(static)
    for (std::size_t cell = 0; cell < nrOfCells; cell++) {
      auto local = layer.cellRef(cell);

      // See TimeCluster: scratch for the neighbours whose traces are evaluated here, written
      // before it is read, one copy per thread.
      alignas(Alignment) real traceBuffer[Cell::NumFaces][kernels::Solver::FaceTraceSize];
      std::array<real*, Cell::NumFaces> traceBuffers{};
      for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
        traceBuffers[i] = traceBuffer[i];
      }

      seissol::kernels::TimeCommon::computeNeighborTraces(data.timeKernel,
                                                          cellInformation[cell].ltsSetup,
                                                          cellInformation[cell].faceTypes,
                                                          timeCoeffs.data(),
                                                          faceNeighbors[cell],
                                                          zeroTrace,
                                                          traceBuffers,
                                                          neighborTraces);

      data.neighborKernel.computeNeighborsIntegral(local, neighborTraces);
    }

    LIKWID_MARKER_STOP("neighboring");
  }
}
auto ProxyKernelHostNeighbor::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  PerformanceEstimate ret;

  auto& layer = data.ltsStorage.layer(data.layerId);
  const auto nrOfCells = layer.size();
  const CellLocalInformation* cellInformation = layer.var<LTS::CellInformation>();
  const auto* drMapping = layer.var<LTS::DRMapping>();

  for (std::size_t cell = 0; cell < nrOfCells; cell++) {
    const auto [cellReg, cellDR] =
        data.neighborKernel.metrics(cellInformation[cell].faceTypes, drMapping[cell]);
    ret += cellReg + cellDR;
  }

  return ret;
}
auto ProxyKernelHostNeighbor::needsDR() const -> bool { return false; }

auto ProxyKernelHostNeighborDR::needsDR() const -> bool { return true; }

void ProxyKernelHostGodunovDR::run(ProxyData& data,
                                   seissol::parallel::runtime::StreamRuntime& /*runtime*/) const {
  auto& layerData = data.drStorage.layer(data.layerId);
  const DRGodunovData* godunovData = layerData.var<DynamicRupture::GodunovData>();
  real* const* timeDerivativePlus = layerData.var<DynamicRupture::TimeDerivativePlus>();
  real* const* timeDerivativeMinus = layerData.var<DynamicRupture::TimeDerivativeMinus>();
  alignas(Alignment) real qInterpolatedPlus[ConvergenceOrder][tensor::QInterpolated::size()];
  alignas(Alignment) real qInterpolatedMinus[ConvergenceOrder][tensor::QInterpolated::size()];
  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, Timestep);
  const auto coeffsCollocate = seissol::kernels::timeBasis().collocate(timePoints, Timestep);

#pragma omp parallel for schedule(static) private(qInterpolatedPlus, qInterpolatedMinus)
  for (std::size_t face = 0; face < layerData.size(); ++face) {
    data.dynRupKernel.spaceTimeInterpolation(&godunovData[face],
                                             timeDerivativePlus[face],
                                             timeDerivativeMinus[face],
                                             qInterpolatedPlus,
                                             qInterpolatedMinus,
                                             coeffsCollocate.data());
  }
}
auto ProxyKernelHostGodunovDR::performanceEstimate(ProxyData& data) const -> PerformanceEstimate {
  return data.dynRupKernel.metrics() * data.drStorage.layer(data.layerId).size();
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
