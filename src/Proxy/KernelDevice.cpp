// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "KernelDevice.h"
#include "Allocator.h"
#include "Common.h"
#include "Kernel.h"
#include <Parallel/Runtime/Stream.h>
#include <memory>

#ifdef ACL_DEVICE
#include "Constants.h"
#include <Kernels/TimeCommon.h>
#include <device.h>
#endif

namespace seissol::proxy {
#ifdef ACL_DEVICE
void ProxyKernelDeviceAder::run(ProxyData& data,
                                seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);

  kernels::<Cfg> tmp(9.81);

  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.getConditionalTable<inner_keys::Material>();

  const double timeStepWidth = static_cast<double>(Timestep);
  ComputeGraphType graphType{ComputeGraphType::AccumulatedVelocities};
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, false);

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.spacetimeKernel.computeBatchedAder(
        integrationCoeffs.data(), timeStepWidth, tmp, dataTable, materialTable, false, runtime);
  });
}

void ProxyKernelDeviceLocalWOAder::run(ProxyData& data,
                                       seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);

  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer.getConditionalTable<inner_keys::Indices>();

  const double timeStepWidth = 0.0;
  ComputeGraphType graphType{ComputeGraphType::AccumulatedVelocities};
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, false);

  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.localKernel.computeBatchedIntegral(
        dataTable, materialTable, indicesTable, timeStepWidth, runtime);
  });
}

void ProxyKernelDeviceLocal::run(ProxyData& data,
                                 seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);

  kernels::LocalTmp<Cfg> tmp(9.81);

  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer.getConditionalTable<inner_keys::Indices>();

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

  const double timeStepWidth = static_cast<double>(Timestep);
  ComputeGraphType graphType{ComputeGraphType::AccumulatedVelocities};
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, false);
  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.spacetimeKernel.computeBatchedAder(
        integrationCoeffs.data(), timeStepWidth, tmp, dataTable, materialTable, false, runtime);
    data.localKernel.computeBatchedIntegral(dataTable, materialTable, indicesTable, 0.0, runtime);
  });
}

void ProxyKernelDeviceNeighbor::run(ProxyData& data,
                                    seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.ltsStorage.layer(data.layerId);

  const double timeStepWidth = static_cast<double>(Timestep);
  auto& dataTable = layer.getConditionalTable<inner_keys::Wp>();

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, Timestep, Timestep);

  seissol::kernels::TimeCommon<Cfg>::computeBatchedIntegrals(
      data.timeKernel, timeCoeffs.data(), timeCoeffs.data(), dataTable, runtime);

  ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);
  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.neighborKernel.computeBatchedNeighborsIntegral(dataTable, runtime);
  });
}

void ProxyKernelDeviceGodunovDR::run(ProxyData& data,
                                     seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& layer = data.drStorage.layer(data.layerId);

  auto& dataTable = layer.getConditionalTable<inner_keys::Dr>();

  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, Timestep);
  const auto coeffsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, Timestep);

  ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
  auto computeGraphKey = initializer::GraphKey(graphType, 0.0);
  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.dynRupKernel.batchedSpaceTimeInterpolation(dataTable, coeffsCollocate.data(), runtime);
  });
}
#else
void ProxyKernelDeviceAder::run(ProxyData& data,
                                seissol::parallel::runtime::StreamRuntime& runtime) const {}

void ProxyKernelDeviceLocalWOAder::run(ProxyData& data,
                                       seissol::parallel::runtime::StreamRuntime& runtime) const {}

void ProxyKernelDeviceLocal::run(ProxyData& data,
                                 seissol::parallel::runtime::StreamRuntime& runtime) const {}

void ProxyKernelDeviceNeighbor::run(ProxyData& data,
                                    seissol::parallel::runtime::StreamRuntime& runtime) const {}

void ProxyKernelDeviceGodunovDR::run(ProxyData& data,
                                     seissol::parallel::runtime::StreamRuntime& runtime) const {}
#endif

auto ProxyKernelDeviceNeighborDR::needsDR() const -> bool { return true; }

std::shared_ptr<ProxyKernel> getProxyKernelDevice(Kernel kernel) {
  switch (kernel) {
  case Kernel::All:
    return std::make_shared<ProxyKernelDeviceAll>();
  case Kernel::AllDR:
    return std::make_shared<ProxyKernelDeviceAllDR>();
  case Kernel::Ader:
    return std::make_shared<ProxyKernelDeviceAder>();
  case Kernel::LocalWOAder:
    return std::make_shared<ProxyKernelDeviceLocalWOAder>();
  case Kernel::Local:
    return std::make_shared<ProxyKernelDeviceLocal>();
  case Kernel::Neighbor:
    return std::make_shared<ProxyKernelDeviceNeighbor>();
  case Kernel::NeighborDR:
    return std::make_shared<ProxyKernelDeviceNeighborDR>();
  case Kernel::GodunovDR:
    return std::make_shared<ProxyKernelDeviceGodunovDR>();
  }
  throw;
}
} // namespace seissol::proxy
