// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "Receiver.h"
#include "GeneratedCode/kernel.h"
#include "Monitoring/FlopCounter.h"
#include "Numerical/BasisFunction.h"
#include "SeisSol.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <Common/Executor.h>
#include <Equations/Datastructures.h>
#include <GeneratedCode/init.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Numerical/Transformation.h>
#include <Parallel/DataCollector.h>
#include <Parallel/Helper.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/MultipleSimulations.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yateto.h>

namespace seissol::kernels {

template <typename Cfg>
Receiver<Cfg>::Receiver(std::size_t pointId,
                        Eigen::Vector3d position,
                        const double* elementCoords[4],
                        LTS::Ref<Cfg> dataHost,
                        LTS::Ref<Cfg> dataDevice,
                        size_t reserved)
    : pointId(pointId), position(std::move(position)), dataHost(dataHost), dataDevice(dataDevice) {
  output.reserve(reserved);

  auto xiEtaZeta = seissol::transformations::tetrahedronGlobalToReference(
      elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], this->position);
  basisFunctions = basisFunction::SampledBasisFunctions<real>(
      Cfg::ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives = basisFunction::SampledBasisFunctionDerivatives<Cfg>(
      Cfg::ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives.transformToGlobalCoordinates(elementCoords);
}

template <typename Cfg>
ReceiverClusterImpl<Cfg>::ReceiverClusterImpl(seissol::SeisSol& seissolInstance)
    : m_samplingInterval(1.0e99), m_syncPointInterval(0.0), seissolInstance(seissolInstance) {}

template <typename Cfg>
ReceiverClusterImpl<Cfg>::ReceiverClusterImpl(
    const GlobalData& global,
    const std::vector<std::size_t>& quantities,
    double samplingInterval,
    double syncPointInterval,
    const std::vector<std::shared_ptr<DerivedReceiverQuantity<Cfg>>>& derivedQuantities,
    seissol::SeisSol& seissolInstance)
    : m_quantities(quantities), m_samplingInterval(samplingInterval),
      m_syncPointInterval(syncPointInterval), derivedQuantities(derivedQuantities),
      seissolInstance(seissolInstance) {
  timeKernel.setGlobalData(global);
  spacetimeKernel.setGlobalData(global);
  spacetimeKernel.flopsAder(m_nonZeroFlops, m_hardwareFlops);
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::addReceiver(std::size_t meshId,
                                           std::size_t pointId,
                                           const Eigen::Vector3d& point,
                                           const seissol::geometry::MeshReader& mesh,
                                           const LTS::Backmap& backmap) {
  const auto& elements = mesh.getElements();
  const auto& vertices = mesh.getVertices();

  const double* coords[Cell::NumVertices];
  for (std::size_t v = 0; v < Cell::NumVertices; ++v) {
    coords[v] = vertices[elements[meshId].vertices[v]].coords;
  }

  if (!extraRuntime.has_value()) {
    // use an extra stream if we have receivers
    extraRuntime.emplace(0);
  }

  // (time + number of quantities) * number of samples until sync point
  const size_t reserved = ncols() * (m_syncPointInterval / m_samplingInterval + 1);

  const auto position = backmap.get(meshId);
  auto& ltsStorage = seissolInstance.getMemoryManager().getLtsStorage();
  m_receivers.emplace_back(pointId,
                           point,
                           coords,
                           ltsStorage.lookupRef<Cfg>(position),
                           ltsStorage.lookupRef<Cfg>(position,
                                                     isDeviceOn()
                                                         ? initializer::AllocationPlace::Device
                                                         : initializer::AllocationPlace::Host),
                           reserved);
}

template <typename Cfg>
double ReceiverClusterImpl<Cfg>::calcReceivers(double time,
                                               double expansionPoint,
                                               double timeStepWidth,
                                               Executor executor,
                                               parallel::runtime::StreamRuntime& runtime) {

  double outReceiverTime = time;
  while (outReceiverTime < expansionPoint + timeStepWidth) {
    outReceiverTime += m_samplingInterval;
  }

  if (executor == Executor::Device) {
    // we need to sync with the new data copy (the rest can continue to run asynchronously)

    if (extraRuntime.has_value()) {
      runtime.eventSync(extraRuntime->eventRecord());
    }
    deviceCollector->gatherToHost(runtime.stream());
    if (extraRuntime.has_value()) {
      extraRuntime->eventSync(runtime.eventRecord());
    }
  }

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();

  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    const std::size_t recvCount = this->m_receivers.size();
    const auto receiverHandler =
        [this, timeBasis, timeStepWidth, time, expansionPoint, executor](std::size_t i) {
          alignas(Alignment) real timeEvaluated[tensor::Q<Cfg>::size()];
          alignas(Alignment) real timeEvaluatedAtPoint[tensor::QAtPoint<Cfg>::size()];
          alignas(Alignment)
              real timeEvaluatedDerivativesAtPoint[tensor::QDerivativeAtPoint<Cfg>::size()];
          alignas(PagesizeStack) real timeDerivatives[Solver<Cfg>::template DerivativesSize<Cfg>];

          kernels::LocalTmp<Cfg> tmp(seissolInstance.getGravitationSetup().acceleration);

          kernel::evaluateDOFSAtPoint<Cfg> krnl;
          krnl.QAtPoint = timeEvaluatedAtPoint;
          krnl.Q = timeEvaluated;
          kernel::evaluateDerivativeDOFSAtPoint<Cfg> derivativeKrnl;
          derivativeKrnl.QDerivativeAtPoint = timeEvaluatedDerivativesAtPoint;
          derivativeKrnl.Q = timeEvaluated;

          auto qAtPoint = init::QAtPoint<Cfg>::view::create(timeEvaluatedAtPoint);
          auto qDerivativeAtPoint =
              init::QDerivativeAtPoint<Cfg>::view::create(timeEvaluatedDerivativesAtPoint);

          auto& receiver = m_receivers[i];
          krnl.basisFunctionsAtPoint = receiver.basisFunctions.m_data.data();
          derivativeKrnl.basisFunctionDerivativesAtPoint =
              receiver.basisFunctionDerivatives.m_data.data();

          // Copy DOFs from device to host.
          auto tmpReceiverData{receiver.dataHost};

          if (executor == Executor::Device) {
            tmpReceiverData.template setPointer<LTS::Dofs>(
                reinterpret_cast<decltype(tmpReceiverData.template getPointer<LTS::Dofs>())>(
                    deviceCollector->get(deviceIndices[i])));
          }

          const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
          spacetimeKernel.computeAder(integrationCoeffs.data(),
                                      timeStepWidth,
                                      tmpReceiverData,
                                      tmp,
                                      timeEvaluated, // useless but the interface requires it
                                      timeDerivatives);

          seissolInstance.flopCounter().incrementNonZeroFlopsOther(m_nonZeroFlops);
          seissolInstance.flopCounter().incrementHardwareFlopsOther(m_hardwareFlops);

          double receiverTime = time;
          while (receiverTime < expansionPoint + timeStepWidth) {
            const auto coeffs = timeBasis.point(time - expansionPoint, timeStepWidth);

            timeKernel.evaluate(coeffs.data(), timeDerivatives, timeEvaluated);

            krnl.execute();
            derivativeKrnl.execute();

            // note: necessary receiver space is reserved in advance
            receiver.output.push_back(receiverTime);
            for (std::size_t sim = 0; sim < seissol::multisim::NumSimulations; ++sim) {
              for (auto quantity : m_quantities) {
                if (!std::isfinite(seissol::multisim::multisimWrap(qAtPoint, sim, quantity))) {
                  logError() << "Detected Inf/NaN in receiver output at" << receiver.position[0]
                             << "," << receiver.position[1] << "," << receiver.position[2]
                             << " in simulation" << sim << "."
                             << "Aborting.";
                }
                receiver.output.push_back(seissol::multisim::multisimWrap(qAtPoint, sim, quantity));
              }
              for (const auto& derived : derivedQuantities) {
                derived->compute(sim, receiver.output, qAtPoint, qDerivativeAtPoint);
              }
            }

            receiverTime += m_samplingInterval;
          }
        };

    auto& callRuntime = extraRuntime.has_value() ? extraRuntime.value() : runtime;
    callRuntime.enqueueLoop(recvCount, receiverHandler);
  }
  return outReceiverTime;
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::allocateData() {
  if constexpr (isDeviceOn()) {
    // collect all data pointers to transfer. If we have multiple receivers on the same cell, we
    // make sure to only transfer the related data once (hence, we use the `indexMap` here)
    deviceIndices.resize(m_receivers.size());
    std::vector<real*> dofs;
    std::unordered_map<real*, size_t> indexMap;
    for (size_t i = 0; i < m_receivers.size(); ++i) {
      real* currentDofs = m_receivers[i].dataDevice.template get<LTS::Dofs>(Cfg());
      if (indexMap.find(currentDofs) == indexMap.end()) {
        // point to the current array end
        indexMap[currentDofs] = dofs.size();
        dofs.push_back(currentDofs);
      }
      deviceIndices[i] = indexMap.at(currentDofs);
    }

    const bool hostAccessible = useUSM() && !extraRuntime.has_value();
    deviceCollector = std::make_unique<seissol::parallel::DataCollector<real>>(
        dofs, tensor::Q<Cfg>::size(), hostAccessible);
  }
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::freeData() {
  deviceCollector.reset(nullptr);
  extraRuntime.reset();
}

template <typename Cfg>
size_t ReceiverClusterImpl<Cfg>::ncols() const {
  size_t ncols = m_quantities.size();
  for (const auto& derived : derivedQuantities) {
    ncols += derived->quantities().size();
  }
  ncols *= seissol::multisim::NumSimulations;
  return 1 + ncols;
}

template <typename Cfg>
std::vector<std::string> ReceiverClusterImpl<Cfg>::header() const {
  std::vector<std::string> names;
  for (const auto& index : m_quantities) {
    names.emplace_back(model::MaterialTT<Cfg>::Quantities[index]);
  }
  for (const auto& derived : derivedQuantities) {
    auto derivedNames = derived->quantities();
    names.insert(names.end(), derivedNames.begin(), derivedNames.end());
  }

  if constexpr (seissol::multisim::MultisimEnabled) {
    std::vector<std::string> fusedNames;
    fusedNames.reserve(seissol::multisim::NumSimulations * names.size());
    for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
      for (const auto& name : names) {
        fusedNames.emplace_back(name + std::to_string(sim));
      }
    }
    return fusedNames;
  } else {
    return names;
  }
}

#define _H_(cfg) template class ReceiverClusterImpl<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::kernels
