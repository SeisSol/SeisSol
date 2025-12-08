// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "KernelDevice.h"

#include "Allocator.h"
#include "Common.h"
#include "Config.h"
#include "Kernel.h"
#include "Parallel/Runtime/Stream.h"

#include <memory>
#include <variant>

#ifdef ACL_DEVICE
#include "Constants.h"
#include "Kernels/TimeCommon.h"

#include <Device/device.h>
#endif

using namespace seissol::recording;

namespace seissol::proxy {
#ifdef ACL_DEVICE
template <typename Cfg>
void ProxyKernelDeviceAder<Cfg>::run(ProxyData& predata,
                                     seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);

  kernels::LocalTmp<Cfg> tmp(9.81);

  auto& dataTable = layer.template getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.template getConditionalTable<inner_keys::Material>();

  ComputeGraphType graphType{ComputeGraphType::AccumulatedVelocities};
  auto computeGraphKey = initializer::GraphKey(graphType, Timestep, false);

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.spacetimeKernel.computeBatchedAder(
        integrationCoeffs.data(), Timestep, tmp, dataTable, materialTable, false, runtime);
  });
}

template <typename Cfg>
void ProxyKernelDeviceLocalWOAder<Cfg>::run(
    ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);

  auto& dataTable = layer.template getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.template getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer.template getConditionalTable<inner_keys::Indices>();

  ComputeGraphType graphType{ComputeGraphType::AccumulatedVelocities};
  auto computeGraphKey = initializer::GraphKey(graphType, Timestep, false);

  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.localKernel.computeBatchedIntegral(
        dataTable, materialTable, indicesTable, Timestep, runtime);
  });
}

template <typename Cfg>
void ProxyKernelDeviceLocal<Cfg>::run(ProxyData& predata,
                                      seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);

  kernels::LocalTmp<Cfg> tmp(9.81);

  auto& dataTable = layer.template getConditionalTable<inner_keys::Wp>();
  auto& materialTable = layer.template getConditionalTable<inner_keys::Material>();
  auto& indicesTable = layer.template getConditionalTable<inner_keys::Indices>();

  const auto integrationCoeffs = data.timeBasis.integrate(0, Timestep, Timestep);

  ComputeGraphType graphType{ComputeGraphType::AccumulatedVelocities};
  auto computeGraphKey = initializer::GraphKey(graphType, Timestep, false);
  runtime.runGraph(computeGraphKey, layer, [&](auto& runtime) {
    data.spacetimeKernel.computeBatchedAder(
        integrationCoeffs.data(), Timestep, tmp, dataTable, materialTable, false, runtime);
    data.localKernel.computeBatchedIntegral(dataTable, materialTable, indicesTable, 0.0, runtime);
  });
}

template <typename Cfg>
void ProxyKernelDeviceNeighbor<Cfg>::run(ProxyData& predata,
                                         seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.ltsStorage.layer(data.layerId);

  auto& dataTable = layer.template getConditionalTable<inner_keys::Wp>();

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

template <typename Cfg>
void ProxyKernelDeviceGodunovDR<Cfg>::run(
    ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const {
  auto& data = static_cast<ProxyDataImpl<Cfg>&>(predata);
  auto& layer = data.drStorage.layer(data.layerId);

  auto& dataTable = layer.template getConditionalTable<inner_keys::Dr>();

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
template <typename Cfg>
void ProxyKernelDeviceAder<Cfg>::run(ProxyData& predata,
                                     seissol::parallel::runtime::StreamRuntime& runtime) const {}

template <typename Cfg>
void ProxyKernelDeviceLocalWOAder<Cfg>::run(
    ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const {}

template <typename Cfg>
void ProxyKernelDeviceLocal<Cfg>::run(ProxyData& predata,
                                      seissol::parallel::runtime::StreamRuntime& runtime) const {}

template <typename Cfg>
void ProxyKernelDeviceNeighbor<Cfg>::run(ProxyData& predata,
                                         seissol::parallel::runtime::StreamRuntime& runtime) const {
}

template <typename Cfg>
void ProxyKernelDeviceGodunovDR<Cfg>::run(
    ProxyData& predata, seissol::parallel::runtime::StreamRuntime& runtime) const {}
#endif

template <typename Cfg>
auto ProxyKernelDeviceNeighborDR<Cfg>::needsDR() const -> bool {
  return true;
}

std::shared_ptr<ProxyKernel> getProxyKernelDevice(Kernel kernel, ConfigVariant variant) {
  return std::visit(
      [&](auto cfg) -> std::shared_ptr<ProxyKernel> {
        using Cfg = decltype(cfg);
        switch (kernel) {
        case Kernel::All:
          return std::make_shared<ProxyKernelDeviceAll<Cfg>>();
        case Kernel::AllDR:
          return std::make_shared<ProxyKernelDeviceAllDR<Cfg>>();
        case Kernel::Ader:
          return std::make_shared<ProxyKernelDeviceAder<Cfg>>();
        case Kernel::LocalWOAder:
          return std::make_shared<ProxyKernelDeviceLocalWOAder<Cfg>>();
        case Kernel::Local:
          return std::make_shared<ProxyKernelDeviceLocal<Cfg>>();
        case Kernel::Neighbor:
          return std::make_shared<ProxyKernelDeviceNeighbor<Cfg>>();
        case Kernel::NeighborDR:
          return std::make_shared<ProxyKernelDeviceNeighborDR<Cfg>>();
        case Kernel::GodunovDR:
          return std::make_shared<ProxyKernelDeviceGodunovDR<Cfg>>();
        }
        throw;
      },
      variant);
}
} // namespace seissol::proxy
