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

Receiver::Receiver(unsigned pointId,
                   Eigen::Vector3d position,
                   const double* elementCoords[4],
                   LTS::Ref dataHost,
                   LTS::Ref dataDevice,
                   size_t reserved)
    : pointId(pointId), position(std::move(position)), dataHost(dataHost), dataDevice(dataDevice) {
  output.reserve(reserved);

  auto xiEtaZeta = seissol::transformations::tetrahedronGlobalToReference(
      elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], this->position);
  basisFunctions = basisFunction::SampledBasisFunctions<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives = basisFunction::SampledBasisFunctionDerivatives<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives.transformToGlobalCoordinates(elementCoords);
}

ReceiverCluster::ReceiverCluster(seissol::SeisSol& seissolInstance)
    : m_samplingInterval(1.0e99), m_syncPointInterval(0.0), seissolInstance(seissolInstance) {}

ReceiverCluster::ReceiverCluster(
    const CompoundGlobalData& global,
    const std::vector<unsigned>& quantities,
    double samplingInterval,
    double syncPointInterval,
    const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
    seissol::SeisSol& seissolInstance)
    : m_quantities(quantities), m_samplingInterval(samplingInterval),
      m_syncPointInterval(syncPointInterval), derivedQuantities(derivedQuantities),
      seissolInstance(seissolInstance) {
  timeKernel.setGlobalData(global);
  spacetimeKernel.setGlobalData(global);
  spacetimeKernel.flopsAder(m_nonZeroFlops, m_hardwareFlops);
}

void ReceiverCluster::addReceiver(unsigned meshId,
                                  unsigned pointId,
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
                           ltsStorage.lookupRef(position),
                           ltsStorage.lookupRef(position,
                                                isDeviceOn() ? initializer::AllocationPlace::Device
                                                             : initializer::AllocationPlace::Host),
                           reserved);
}

double ReceiverCluster::calcReceivers(double time,
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
    const std::size_t recvCount = m_receivers.size();
    const auto receiverHandler = [this, timeBasis, timeStepWidth, time, expansionPoint, executor](
                                     std::size_t i) {
      alignas(Alignment) real timeEvaluated[tensor::Q<Cfg>::size()];
      alignas(Alignment) real timeEvaluatedAtPoint[tensor::QAtPoint<Cfg>::size()];
      alignas(Alignment) real timeEvaluatedDerivativesAtPoint[tensor::QDerivativeAtPoint<Cfg>::size()];
      alignas(PagesizeStack) real timeDerivatives[Solver<Cfg>::DerivativesSize];

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
        tmpReceiverData.setPointer<LTS::Dofs>(
            reinterpret_cast<decltype(tmpReceiverData.getPointer<LTS::Dofs>())>(
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
        for (unsigned sim = seissol::multisim::MultisimStart; sim < seissol::multisim::MultisimEnd;
             ++sim) {
          for (auto quantity : m_quantities) {
            if (!std::isfinite(seissol::multisim::multisimWrap(qAtPoint, sim, quantity))) {
              logError() << "Detected Inf/NaN in receiver output at" << receiver.position[0] << ","
                         << receiver.position[1] << "," << receiver.position[2] << " in simulation"
                         << sim << "."
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

void ReceiverCluster::allocateData() {
  if constexpr (isDeviceOn()) {
    // collect all data pointers to transfer. If we have multiple receivers on the same cell, we
    // make sure to only transfer the related data once (hence, we use the `indexMap` here)
    deviceIndices.resize(m_receivers.size());
    std::vector<real*> dofs;
    std::unordered_map<real*, size_t> indexMap;
    for (size_t i = 0; i < m_receivers.size(); ++i) {
      real* currentDofs = m_receivers[i].dataDevice.get<LTS::Dofs>();
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
void ReceiverCluster::freeData() {
  deviceCollector.reset(nullptr);
  extraRuntime.reset();
}

size_t ReceiverCluster::ncols() const {
  size_t ncols = m_quantities.size();
  for (const auto& derived : derivedQuantities) {
    ncols += derived->quantities().size();
  }
  ncols *= seissol::multisim::MultisimEnd - seissol::multisim::MultisimStart;
  return 1 + ncols;
}

std::vector<std::string> ReceiverRotation::quantities() const { return {"rot1", "rot2", "rot3"}; }
void ReceiverRotation::compute(size_t sim,
                               std::vector<real>& output,
                               seissol::init::QAtPoint<Cfg>::view::type& qAtPoint,
                               seissol::init::QDerivativeAtPoint<Cfg>::view::type& qDerivativeAtPoint) {
  output.push_back(seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 8, 1) -
                   seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 7, 2));
  output.push_back(seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 6, 2) -
                   seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 8, 0));
  output.push_back(seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 7, 0) -
                   seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 6, 1));
}

std::vector<std::string> ReceiverStrain::quantities() const {
  return {"epsxx", "epsxy", "epsxz", "epsyy", "epsyz", "epszz"};
}
void ReceiverStrain::compute(size_t sim,
                             std::vector<real>& output,
                             seissol::init::QAtPoint<Cfg>::view::type& qAtPoint,
                             seissol::init::QDerivativeAtPoint<Cfg>::view::type& qDerivativeAtPoint) {
  // actually 9 quantities; 3 removed due to symmetry

  output.push_back(seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 6, 0));
  output.push_back((seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 6, 1) +
                    seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 7, 0)) /
                   2);
  output.push_back((seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 6, 2) +
                    seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 8, 0)) /
                   2);
  output.push_back(seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 7, 1));
  output.push_back((seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 7, 2) +
                    seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 8, 1)) /
                   2);
  output.push_back(seissol::multisim::multisimWrap(qDerivativeAtPoint, sim, 8, 2));
}

} // namespace seissol::kernels
