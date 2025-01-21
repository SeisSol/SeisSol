// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "Receiver.h"
#include "Monitoring/FlopCounter.h"
#include "Numerical/BasisFunction.h"
#include "SeisSol.h"
#include "generated_code/kernel.h"
#include <Common/Constants.h>
#include <Common/Executor.h>
#include <Kernels/Common.h>
#include <Kernels/Interface.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Memory/Tree/Lut.h>
#include <Numerical/Transformation.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <init.h>
#include <memory>
#include <omp.h>
#include <string>
#include <tensor.h>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "device.h"
#include <Parallel/Helper.h>
#include <unordered_map>
#endif

namespace {
#ifdef MULTIPLE_SIMULATIONS
template <typename T, typename F, typename... Args>
T multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), sim, std::forward<Args>(args)...);
}
constexpr size_t MultisimStart = init::QAtPoint::Start[0];
constexpr size_t MultisimEnd = init::QAtPoint::Stop[0];
#else
template <typename F, typename... Args>
auto multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), std::forward<Args>(args)...);
}
constexpr size_t MultisimStart = 0;
constexpr size_t MultisimEnd = 1;
#endif
} // namespace

namespace seissol::kernels {

Receiver::Receiver(unsigned pointId,
                   Eigen::Vector3d position,
                   const double* elementCoords[4],
                   kernels::LocalData dataHost,
                   kernels::LocalData dataDevice,
                   size_t reserved)
    : pointId(pointId), position(std::move(position)), dataHost(dataHost), dataDevice(dataDevice) {
  output.reserve(reserved);

  auto xiEtaZeta = seissol::transformations::tetrahedronGlobalToReference(
      elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], position);
  basisFunctions = basisFunction::SampledBasisFunctions<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives = basisFunction::SampledBasisFunctionDerivatives<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives.transformToGlobalCoordinates(elementCoords);
}

ReceiverCluster::ReceiverCluster(seissol::SeisSol& seissolInstance)
    : m_samplingInterval(1.0e99), m_syncPointInterval(0.0), seissolInstance(seissolInstance) {}

ReceiverCluster::ReceiverCluster(
    const GlobalData* global,
    const std::vector<unsigned>& quantities,
    double samplingInterval,
    double syncPointInterval,
    const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
    seissol::SeisSol& seissolInstance)
    : m_quantities(quantities), m_samplingInterval(samplingInterval),
      m_syncPointInterval(syncPointInterval), derivedQuantities(derivedQuantities),
      seissolInstance(seissolInstance) {
  m_timeKernel.setHostGlobalData(global);
  m_timeKernel.flopsAder(m_nonZeroFlops, m_hardwareFlops);
}

void ReceiverCluster::addReceiver(unsigned meshId,
                                  unsigned pointId,
                                  const Eigen::Vector3d& point,
                                  const seissol::geometry::MeshReader& mesh,
                                  const seissol::initializer::Lut& ltsLut,
                                  seissol::initializer::LTS const& lts) {
  const auto& elements = mesh.getElements();
  const auto& vertices = mesh.getVertices();

  const double* coords[4];
  for (unsigned v = 0; v < 4; ++v) {
    coords[v] = vertices[elements[meshId].vertices[v]].coords;
  }

  // (time + number of quantities) * number of samples until sync point
  const size_t reserved = ncols() * (m_syncPointInterval / m_samplingInterval + 1);
  m_receivers.emplace_back(
      pointId,
      point,
      coords,
      kernels::LocalData::lookup(lts, ltsLut, meshId, initializer::AllocationPlace::Host),
      kernels::LocalData::lookup(lts,
                                 ltsLut,
                                 meshId,
                                 isDeviceOn() ? initializer::AllocationPlace::Device
                                              : initializer::AllocationPlace::Host),
      reserved);
}

double ReceiverCluster::calcReceivers(
    double time, double expansionPoint, double timeStepWidth, Executor executor, void* stream) {

  double outReceiverTime = time;
  while (outReceiverTime < expansionPoint + timeStepWidth) {
    outReceiverTime += m_samplingInterval;
  }

#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    deviceCollector->gatherToHost(device::DeviceInstance::getInstance().api->getDefaultStream());
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif

  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    // heuristic; to avoid the overhead from the parallel region
    const std::size_t threshold = std::max(1000, omp_get_num_threads() * 100);
    const std::size_t recvCount = m_receivers.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (recvCount >= threshold)
#endif
    for (size_t i = 0; i < recvCount; ++i) {
      alignas(Alignment) real timeEvaluated[tensor::Q::size()];
      alignas(Alignment) real timeEvaluatedAtPoint[tensor::QAtPoint::size()];
      alignas(Alignment) real timeEvaluatedDerivativesAtPoint[tensor::QDerivativeAtPoint::size()];
#ifdef USE_STP
      alignas(PagesizeStack) real stp[tensor::spaceTimePredictor::size()];
      kernel::evaluateDOFSAtPointSTP krnl;
      krnl.QAtPoint = timeEvaluatedAtPoint;
      krnl.spaceTimePredictor = stp;
      kernel::evaluateDerivativeDOFSAtPointSTP derivativeKrnl;
      derivativeKrnl.QDerivativeAtPoint = timeEvaluatedDerivativesAtPoint;
      derivativeKrnl.spaceTimePredictor = stp;
#else
      alignas(Alignment) real timeDerivatives[yateto::computeFamilySize<tensor::dQ>()];
      kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);

      kernel::evaluateDOFSAtPoint krnl;
      krnl.QAtPoint = timeEvaluatedAtPoint;
      krnl.Q = timeEvaluated;
      kernel::evaluateDerivativeDOFSAtPoint derivativeKrnl;
      derivativeKrnl.QDerivativeAtPoint = timeEvaluatedDerivativesAtPoint;
      derivativeKrnl.Q = timeEvaluated;
#endif

      auto qAtPoint = init::QAtPoint::view::create(timeEvaluatedAtPoint);
      auto qDerivativeAtPoint =
          init::QDerivativeAtPoint::view::create(timeEvaluatedDerivativesAtPoint);

      auto& receiver = m_receivers[i];
      krnl.basisFunctionsAtPoint = receiver.basisFunctions.m_data.data();
      derivativeKrnl.basisFunctionDerivativesAtPoint =
          receiver.basisFunctionDerivatives.m_data.data();

      // Copy DOFs from device to host.
      LocalData tmpReceiverData{receiver.dataHost};
#ifdef ACL_DEVICE
      if (executor == Executor::Device) {
        tmpReceiverData.dofs_ptr = reinterpret_cast<decltype(tmpReceiverData.dofs_ptr)>(
            deviceCollector->get(deviceIndices[i]));
      }
#endif

#ifdef USE_STP
      m_timeKernel.executeSTP(timeStepWidth, tmpReceiverData, timeEvaluated, stp);
#else
      m_timeKernel.computeAder(timeStepWidth,
                               tmpReceiverData,
                               tmp,
                               timeEvaluated, // useless but the interface requires it
                               timeDerivatives);
#endif
      seissolInstance.flopCounter().incrementNonZeroFlopsOther(m_nonZeroFlops);
      seissolInstance.flopCounter().incrementHardwareFlopsOther(m_hardwareFlops);

      double receiverTime = time;
      while (receiverTime < expansionPoint + timeStepWidth) {
#ifdef USE_STP
        // eval time basis
        const double tau = (time - expansionPoint) / timeStepWidth;
        seissol::basisFunction::SampledTimeBasisFunctions<real> timeBasisFunctions(ConvergenceOrder,
                                                                                   tau);
        krnl.timeBasisFunctionsAtPoint = timeBasisFunctions.m_data.data();
        derivativeKrnl.timeBasisFunctionsAtPoint = timeBasisFunctions.m_data.data();
#else
        m_timeKernel.computeTaylorExpansion(
            receiverTime, expansionPoint, timeDerivatives, timeEvaluated);
#endif

        krnl.execute();
        derivativeKrnl.execute();

        // note: necessary receiver space is reserved in advance
        receiver.output.push_back(receiverTime);
        for (unsigned sim = MultisimStart; sim < MultisimEnd; ++sim) {
          for (auto quantity : m_quantities) {
            if (!std::isfinite(multisimWrap(qAtPoint, sim, quantity))) {
              logError() << "Detected Inf/NaN in receiver output at" << receiver.position[0] << ","
                         << receiver.position[1] << "," << receiver.position[2] << " in simulation"
                         << sim << "."
                         << "Aborting.";
            }
            receiver.output.push_back(multisimWrap(qAtPoint, sim, quantity));
          }
          for (const auto& derived : derivedQuantities) {
            derived->compute(sim, receiver.output, qAtPoint, qDerivativeAtPoint);
          }
        }

        receiverTime += m_samplingInterval;
      }
    }
  }
  return outReceiverTime;
}

void ReceiverCluster::allocateData() {
#ifdef ACL_DEVICE
  // collect all data pointers to transfer. If we have multiple receivers on the same cell, we make
  // sure to only transfer the related data once (hence, we use the `indexMap` here)
  deviceIndices.resize(m_receivers.size());
  std::vector<real*> dofs;
  std::unordered_map<real*, size_t> indexMap;
  for (size_t i = 0; i < m_receivers.size(); ++i) {
    real* currentDofs = m_receivers[i].dataDevice.dofs();
    if (indexMap.find(currentDofs) == indexMap.end()) {
      // point to the current array end
      indexMap[currentDofs] = dofs.size();
      dofs.push_back(currentDofs);
    }
    deviceIndices[i] = indexMap.at(currentDofs);
  }
  deviceCollector =
      std::make_unique<seissol::parallel::DataCollector>(dofs, tensor::Q::size(), useUSM());
#endif
}
void ReceiverCluster::freeData() {
#ifdef ACL_DEVICE
  deviceCollector.reset(nullptr);
#endif
}

size_t ReceiverCluster::ncols() const {
  size_t ncols = m_quantities.size();
  for (const auto& derived : derivedQuantities) {
    ncols += derived->quantities().size();
  }
  ncols *= MultisimEnd - MultisimStart;
  return 1 + ncols;
}

std::vector<std::string> ReceiverRotation::quantities() const { return {"rot1", "rot2", "rot3"}; }
void ReceiverRotation::compute(size_t sim,
                               std::vector<real>& output,
                               seissol::init::QAtPoint::view::type& qAtPoint,
                               seissol::init::QDerivativeAtPoint::view::type& qDerivativeAtPoint) {
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 8, 1) -
                   multisimWrap(qDerivativeAtPoint, sim, 7, 2));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6, 2) -
                   multisimWrap(qDerivativeAtPoint, sim, 8, 0));
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 7, 0) -
                   multisimWrap(qDerivativeAtPoint, sim, 6, 1));
}

std::vector<std::string> ReceiverStrain::quantities() const {
  return {"epsxx", "epsxy", "epsxz", "epsyy", "epsyz", "epszz"};
}
void ReceiverStrain::compute(size_t sim,
                             std::vector<real>& output,
                             seissol::init::QAtPoint::view::type& qAtPoint,
                             seissol::init::QDerivativeAtPoint::view::type& qDerivativeAtPoint) {
  // actually 9 quantities; 3 removed due to symmetry

  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 6, 0));
  output.push_back(
      (multisimWrap(qDerivativeAtPoint, sim, 6, 1) + multisimWrap(qDerivativeAtPoint, sim, 7, 0)) /
      2);
  output.push_back(
      (multisimWrap(qDerivativeAtPoint, sim, 6, 2) + multisimWrap(qDerivativeAtPoint, sim, 8, 0)) /
      2);
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 7, 1));
  output.push_back(
      (multisimWrap(qDerivativeAtPoint, sim, 7, 2) + multisimWrap(qDerivativeAtPoint, sim, 8, 1)) /
      2);
  output.push_back(multisimWrap(qDerivativeAtPoint, sim, 8, 2));
}

} // namespace seissol::kernels
