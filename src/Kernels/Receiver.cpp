// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "Receiver.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Executor.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Metric.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Transformation.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Helper.h"
#include "Parallel/Runtime/Stream.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yateto.h>

namespace seissol::kernels {

Receiver::Receiver(std::size_t pointId,
                   Eigen::Vector3d position,
                   const double* elementCoords[4],
                   size_t reserved)
    : pointId(pointId), position(std::move(position)) {
  output.reserve(reserved);

  auto xiEtaZeta = seissol::transformations::tetrahedronGlobalToReference(
      elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], this->position);
  basisFunctions = basisFunction::SampledBasisFunctions<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives = basisFunction::SampledBasisFunctionDerivatives<real>(
      ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives.transformToGlobalCoordinates(elementCoords);
}

ReceiverCell::ReceiverCell(std::size_t meshId, LTS::Ref dataHost, LTS::Ref dataDevice)
    : meshId(meshId), dataHost(dataHost), dataDevice(dataDevice) {}

ReceiverCluster::ReceiverCluster(seissol::SeisSol& seissolInstance)
    : samplingInterval_(1.0e99), syncPointInterval_(0.0), seissolInstance_(seissolInstance) {}

ReceiverCluster::ReceiverCluster(
    const CompoundGlobalData& global,
    const std::vector<std::size_t>& quantities,
    double samplingInterval,
    double syncPointInterval,
    const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
    seissol::SeisSol& seissolInstance)
    : quantities_(quantities), samplingInterval_(samplingInterval),
      syncPointInterval_(syncPointInterval), derivedQuantities_(derivedQuantities),
      seissolInstance_(seissolInstance) {
  timeKernel_.setGlobalData(global);
  spacetimeKernel_.setGlobalData(global);

  estimatePerCell_ = spacetimeKernel_.metrics();
  estimatePerCellStep_ = timeKernel_.metrics();
  estimatePerPoint_ = PerformanceEstimate::fromKernel<kernel::evaluateDOFSAtPoint>() +
                      PerformanceEstimate::fromKernel<kernel::evaluateDerivativeDOFSAtPoint>();

  perfHandle_ = seissolInstance_.flopCounter().addMetric("receiver", "WP");
}

void ReceiverCluster::addReceiver(std::size_t meshId,
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

  if (!extraRuntime_.has_value()) {
    // use an extra stream if we have receivers
    extraRuntime_.emplace(0);
  }

  // (time + number of quantities) * number of samples until sync point
  const size_t reserved = ncols() * (syncPointInterval_ / samplingInterval_ + 1);

  if (meshToReceiverCell_.find(meshId) == meshToReceiverCell_.end()) {
    const auto position = backmap.get(meshId);
    auto& ltsStorage = seissolInstance_.memoryManager().ltsStorage();

    meshToReceiverCell_[meshId] = receiverCells_.size();

    auto& cell = receiverCells_.emplace_back(
        meshId,
        ltsStorage.lookupRef(position),
        ltsStorage.lookupRef(position,
                             isDeviceOn() ? initializer::AllocationPlace::Device
                                          : initializer::AllocationPlace::Host));
    cell.ltsPosition = position.global;
  }

  receiverCells_[meshToReceiverCell_.at(meshId)].receiverIds.emplace_back(receivers_.size());

  receivers_.emplace_back(pointId, point, coords, reserved);
}

double ReceiverCluster::calcReceivers(double time,
                                      double expansionPoint,
                                      double timeStepWidth,
                                      Executor executor,
                                      parallel::runtime::StreamRuntime& runtime) {

  double outReceiverTime = time;
  std::size_t samplingSteps = 0;
  while (outReceiverTime < expansionPoint + timeStepWidth) {
    outReceiverTime += samplingInterval_;
    ++samplingSteps;
  }

  // copy dofs from the device to the host.
  if (executor == Executor::Device) {
    // we need to sync with the new data copy (the rest can continue to run asynchronously)

    if (extraRuntime_.has_value()) {
      runtime.eventSync(extraRuntime_->eventRecord());
    }
    deviceCollector_->gatherToHost(runtime.stream());
    if (extraRuntime_.has_value()) {
      extraRuntime_->eventSync(runtime.eventRecord());
    }
  }

  const auto timeBasis = seissol::kernels::timeBasis();

  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    const std::size_t cellCount = receiverCells_.size();
    const auto receiverHandler = [this, timeBasis, timeStepWidth, time, expansionPoint, executor](
                                     std::size_t i) {
      alignas(Alignment) real timeEvaluated[tensor::Q::size()]{};
      alignas(Alignment) real timeEvaluatedAtPoint[tensor::QAtPoint::size()]{};
      alignas(Alignment) real timeEvaluatedDerivativesAtPoint[tensor::QDerivativeAtPoint::size()]{};
      alignas(PagesizeStack) real timeDerivatives[Solver::DerivativesSize]{};

      kernels::LocalTmp tmp(seissolInstance_.gravitationSetup().acceleration);

      kernel::evaluateDOFSAtPoint krnl;
      krnl.QAtPoint = timeEvaluatedAtPoint;
      krnl.Q = timeEvaluated;
      kernel::evaluateDerivativeDOFSAtPoint derivativeKrnl;
      derivativeKrnl.QDerivativeAtPoint = timeEvaluatedDerivativesAtPoint;
      derivativeKrnl.Q = timeEvaluated;

      auto qAtPoint = init::QAtPoint::view::create(timeEvaluatedAtPoint);
      auto qDerivativeAtPoint =
          init::QDerivativeAtPoint::view::create(timeEvaluatedDerivativesAtPoint);

      auto& receiverCell = receiverCells_[i];

      // Use device pointers where required.
      auto tmpReceiverData{receiverCell.dataHost};

      if (executor == Executor::Device) {
        tmpReceiverData.setPointer<LTS::Dofs>(
            reinterpret_cast<decltype(tmpReceiverData.getPointer<LTS::Dofs>())>(
                deviceCollector_->get(i)));
      }

      const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
      spacetimeKernel_.computeAder(integrationCoeffs.data(),
                                   timeStepWidth,
                                   tmpReceiverData,
                                   tmp,
                                   timeEvaluated, // useless but the interface requires it
                                   timeDerivatives);

      double receiverTime = time;
      while (receiverTime < expansionPoint + timeStepWidth) {
        const auto coeffs = timeBasis.point(receiverTime - expansionPoint, timeStepWidth);

        timeKernel_.evaluate(coeffs.data(), timeDerivatives, timeEvaluated);

        for (const auto& receiverId : receiverCell.receiverIds) {

          auto& receiver = receivers_[receiverId];

          krnl.basisFunctionsAtPoint = receiver.basisFunctions.data().data();
          derivativeKrnl.basisFunctionDerivativesAtPoint =
              receiver.basisFunctionDerivatives.data().data();

          krnl.execute();
          derivativeKrnl.execute();

          // note: necessary receiver space is reserved in advance
          receiver.output.push_back(receiverTime);
          for (auto sim = seissol::multisim::MultisimStart; sim < seissol::multisim::MultisimEnd;
               ++sim) {
            for (auto quantity : quantities_) {
              if (!std::isfinite(seissol::multisim::multisimWrap(qAtPoint, sim, quantity))) {
                logError() << "Detected Inf/NaN in receiver output at" << receiver.position[0]
                           << "," << receiver.position[1] << "," << receiver.position[2]
                           << " in simulation" << sim << "."
                           << "Aborting.";
              }
              receiver.output.push_back(seissol::multisim::multisimWrap(qAtPoint, sim, quantity));
            }
            for (const auto& derived : derivedQuantities_) {
              derived->compute(sim, receiver.output, qAtPoint, qDerivativeAtPoint);
            }
          }
        }

        receiverTime += samplingInterval_;
      }
    };

    auto& callRuntime = extraRuntime_.has_value() ? extraRuntime_.value() : runtime;
    callRuntime.enqueueLoop(cellCount, receiverHandler);

    const auto recvCount = receivers_.size();

    seissolInstance_.flopCounter().incrementMetric(
        perfHandle_,
        estimatePerCell_ * cellCount + estimatePerCellStep_ * cellCount * samplingSteps +
            estimatePerPoint_ * recvCount * samplingSteps);
  }
  return outReceiverTime;
}

void ReceiverCluster::allocateData() {
  // Visit the cells in storage order, so that both the host loop and the device gather read the
  // DOFs sequentially instead of in the order the receivers happened to appear in the parameter
  // file. This is only safe because a receiver does not refer back to its cell.
  std::vector<std::size_t> order(receiverCells_.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [this](std::size_t a, std::size_t b) {
    return receiverCells_[a].ltsPosition < receiverCells_[b].ltsPosition;
  });

  std::vector<ReceiverCell> sorted;
  sorted.reserve(receiverCells_.size());
  for (const auto cellId : order) {
    sorted.push_back(receiverCells_[cellId]);
  }
  receiverCells_ = std::move(sorted);

  if constexpr (isDeviceOn()) {
    // one entry per cell; the gather index equals the cell index
    std::vector<real*> dofs;
    dofs.reserve(receiverCells_.size());
    for (const auto& receiverCell : receiverCells_) {
      dofs.push_back(receiverCell.dataDevice.get<LTS::Dofs>());
    }

    const bool hostAccessible = useUSM() && !extraRuntime_.has_value();
    deviceCollector_ = std::make_unique<seissol::parallel::DataCollector<real>>(
        dofs, tensor::Q::size(), hostAccessible);
  }

  meshToReceiverCell_ = {};
}
void ReceiverCluster::freeData() {
  deviceCollector_.reset(nullptr);
  extraRuntime_.reset();
}

size_t ReceiverCluster::ncols() const {
  size_t ncols = quantities_.size();
  for (const auto& derived : derivedQuantities_) {
    ncols += derived->quantities().size();
  }
  ncols *= seissol::multisim::MultisimEnd - seissol::multisim::MultisimStart;
  return 1 + ncols;
}

std::vector<std::string> ReceiverRotation::quantities() const { return {"rot1", "rot2", "rot3"}; }
void ReceiverRotation::compute(size_t sim,
                               std::vector<real>& output,
                               seissol::init::QAtPoint::view::type& /*qAtPoint*/,
                               seissol::init::QDerivativeAtPoint::view::type& qDerivativeAtPoint) {
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
                             seissol::init::QAtPoint::view::type& /*qAtPoint*/,
                             seissol::init::QDerivativeAtPoint::view::type& qDerivativeAtPoint) {
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
