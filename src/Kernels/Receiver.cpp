// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "Receiver.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Common/Executor.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Metric.h"
#include "Monitoring/Trace.h"
#include "Numerical/BasisFunction.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Helper.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "Parallel/Runtime/Stream.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
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

#ifdef ACL_DEVICE
#include "Kernels/DeviceAux/HostCopyAux.h"

#include <Device/device.h>
#endif

namespace seissol::kernels {

Receiver::Receiver(std::size_t pointId, Eigen::Vector3d position, size_t reserved)
    : pointId(pointId), position(std::move(position)) {
  output.reserve(reserved);
}

template <typename Cfg>
ReceiverBasis<Cfg>::ReceiverBasis(const Eigen::Vector3d& position,
                                  const seissol::geometry::CellTransform& transform) {
  const auto xiEtaZeta = transform.spaceToRef(position);
  basisFunctions = basisFunction::SampledBasisFunctions<Real<Cfg>>(
      Cfg::ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives = basisFunction::SampledBasisFunctionDerivatives<Real<Cfg>>(
      Cfg::ConvergenceOrder, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
  basisFunctionDerivatives.transformToGlobalCoordinates(
      transform, xiEtaZeta[0], xiEtaZeta[1], xiEtaZeta[2]);
}

template <typename Cfg>
ReceiverCell<Cfg>::ReceiverCell(std::size_t meshId,
                                LTS::Ref<Cfg> dataHost,
                                LTS::Ref<Cfg> dataDevice)
    : meshId(meshId), dataHost(dataHost), dataDevice(dataDevice) {}

template <typename Cfg>
ReceiverClusterImpl<Cfg>::ReceiverClusterImpl(seissol::SeisSol& seissolInstance)
    : samplingInterval_(1.0e99), syncPointInterval_(0.0), seissolInstance_(seissolInstance) {}

template <typename Cfg>
ReceiverClusterImpl<Cfg>::ReceiverClusterImpl(
    const CompoundGlobalData<Cfg>& global,
    const std::vector<std::size_t>& quantities,
    double samplingInterval,
    double syncPointInterval,
    const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
    seissol::SeisSol& seissolInstance)
    : quantities_(quantities), estimatePerCell_(spacetimeKernel_.metrics()),
      estimatePerCellStep_(timeKernel_.metrics()),
      estimatePerPoint_(
          PerformanceEstimate::fromKernel<kernel::evaluateDOFSAtPoint<Cfg>>() +
          PerformanceEstimate::fromKernel<kernel::evaluateDerivativeDOFSAtPoint<Cfg>>()),
      samplingInterval_(samplingInterval), syncPointInterval_(syncPointInterval),
      derivedQuantities_(derivedQuantities), seissolInstance_(seissolInstance) {
  timeKernel_.setGlobalData(global);
  spacetimeKernel_.setGlobalData(global);

  perfHandle_ = seissolInstance_.flopCounter().addMetric("receiver", "WP");
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::addReceiver(std::size_t meshId,
                                           std::size_t pointId,
                                           const Eigen::Vector3d& point,
                                           const seissol::geometry::MeshReader& mesh,
                                           const LTS::Backmap& backmap) {
  const auto transform = seissol::geometry::AffineTransform::fromMeshCell(meshId, mesh);

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
        ltsStorage.template lookupRef<Cfg>(position),
        ltsStorage.template lookupRef<Cfg>(position,
                                           isDeviceOn() ? initializer::AllocationPlace::Device
                                                        : initializer::AllocationPlace::Host));
    cell.ltsPosition = position.global;
  }

  receiverCells_[meshToReceiverCell_.at(meshId)].receiverIds.emplace_back(receivers_.size());

  receivers_.emplace_back(pointId, point, reserved);
  receiverBases_.emplace_back(point, transform);
}

template <typename Cfg>
ReceiverCluster::Sampling ReceiverClusterImpl<Cfg>::planSampling(double time,
                                                                 double expansionPoint,
                                                                 double timeStepWidth) const {
  Sampling sampling;
  sampling.time = time;
  sampling.nextTime = time;
  while (sampling.nextTime < expansionPoint + timeStepWidth) {
    sampling.nextTime += samplingInterval_;
    ++sampling.steps;
  }
  sampling.due =
      !receivers_.empty() && time >= expansionPoint && time < expansionPoint + timeStepWidth;
  return sampling;
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::sampleReceiver(
    std::size_t cell, double time, double expansionPoint, double timeStepWidth, Executor executor) {
  using Multisim = seissol::multisim::MultisimHelperWrapper<Cfg>;

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  alignas(Alignment) real timeEvaluated[tensor::Q<Cfg>::size()]{};
  alignas(Alignment) real timeEvaluatedAtPoint[tensor::QAtPoint<Cfg>::size()]{};
  alignas(Alignment)
      real timeEvaluatedDerivativesAtPoint[tensor::QDerivativeAtPoint<Cfg>::size()]{};
  alignas(PagesizeStack) real timeDerivatives[SolverOf<Cfg>::DerivativesSize]{};

  kernels::LocalTmp<Cfg> tmp(seissolInstance_.gravitationSetup().acceleration);

  constexpr auto Variant = configIdOf<Cfg>();
  runtime::kernel::evaluateDOFSAtPoint krnl;
  krnl.QAtPoint = runtime::init::QAtPoint::view(Variant, timeEvaluatedAtPoint);
  krnl.Q = runtime::init::Q::view(Variant, timeEvaluated);
  runtime::kernel::evaluateDerivativeDOFSAtPoint derivativeKrnl;
  derivativeKrnl.QDerivativeAtPoint =
      runtime::init::QDerivativeAtPoint::view(Variant, timeEvaluatedDerivativesAtPoint);
  derivativeKrnl.Q = runtime::init::Q::view(Variant, timeEvaluated);

  auto qAtPoint = init::QAtPoint<Cfg>::view::create(timeEvaluatedAtPoint);
  auto qDerivativeAtPoint =
      init::QDerivativeAtPoint<Cfg>::view::create(timeEvaluatedDerivativesAtPoint);

  auto& receiverCell = receiverCells_[cell];

  // Use device pointers where required.
  auto tmpReceiverData{receiverCell.dataHost};

  if (executor == Executor::Device) {
    tmpReceiverData.template setPointer<LTS::Dofs>(
        reinterpret_cast<decltype(tmpReceiverData.template getPointer<LTS::Dofs>())>(
            deviceCollector_->get(cell)));
    if constexpr (kernels::size<tensor::Qane<Cfg>>() > 0) {
      tmpReceiverData.template setPointer<LTS::DofsAne>(
          reinterpret_cast<decltype(tmpReceiverData.template getPointer<LTS::DofsAne>())>(
              deviceCollectorAne_->get(cell)));
    }
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
      const auto& basis = receiverBases_[receiverId];

      krnl.basisFunctionsAtPoint =
          runtime::init::basisFunctionsAtPoint::view(Variant, basis.basisFunctions.data().data());
      derivativeKrnl.basisFunctionDerivativesAtPoint =
          runtime::init::basisFunctionDerivativesAtPoint::view(
              Variant, basis.basisFunctionDerivatives.data().data());

      krnl.execute(Variant);
      derivativeKrnl.execute(Variant);

      // note: necessary receiver space is reserved in advance
      receiver.output.push_back(receiverTime);
      for (auto sim = Multisim::MultisimStart; sim < Multisim::MultisimEnd; ++sim) {
        for (auto quantity : quantities_) {
          if (!std::isfinite(seissol::multisim::multisimWrap<Cfg>(qAtPoint, sim, quantity))) {
            logError() << "Detected Inf/NaN in receiver output at" << receiver.position[0] << ","
                       << receiver.position[1] << "," << receiver.position[2] << " in simulation"
                       << sim << "."
                       << "Aborting.";
          }
          receiver.output.push_back(seissol::multisim::multisimWrap<Cfg>(qAtPoint, sim, quantity));
        }
        if (!derivedQuantities_.empty()) {
          const auto gradient = velocityGradient<Cfg>(qDerivativeAtPoint, sim);
          for (const auto& derived : derivedQuantities_) {
            derived->compute(receiver.output, gradient);
          }
        }
      }
    }

    receiverTime += samplingInterval_;
  }
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::countSamples(std::size_t samplingSteps) {
  const auto cellCount = receiverCells_.size();
  const auto recvCount = receivers_.size();
  seissolInstance_.flopCounter().incrementMetric(
      perfHandle_,
      estimatePerCell_ * cellCount + estimatePerCellStep_ * cellCount * samplingSteps +
          estimatePerPoint_ * recvCount * samplingSteps);
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::sample(const Sampling& sampling,
                                      double expansionPoint,
                                      double timeStepWidth,
                                      Executor executor,
                                      parallel::runtime::StreamRuntime& runtime) {
  if (!sampling.due) {
    return;
  }
  const double time = sampling.time;

  // copy dofs from the device to the host.
  if (executor == Executor::Device) {
    // we need to sync with the new data copy (the rest can continue to run asynchronously)
    if (extraRuntime_.has_value()) {
      runtime.eventSync(extraRuntime_->eventRecord());
    }
    gatherDofs(runtime);
  }

  const std::size_t cellCount = receiverCells_.size();
  const auto receiverHandler =
      [this, timeStepWidth, time, expansionPoint, executor](std::size_t cell) {
        sampleReceiver(cell, time, expansionPoint, timeStepWidth, executor);
      };

  if (executor == Executor::Host) {
    // A cluster that runs on the host goes on to integrate right after this and overwrites the
    // DOFs the sampling reads, so it samples right here. (On CUDA and SYCL, enqueueLoop would
    // leave the sampling to a host function on a stream.)
    static const auto TraceName = monitoring::Trace::name("receivers");
    const parallel::runtime::HostRuntime::LoopName loopName(TraceName);
    parallel::runtime::HostRuntime::forEach(cellCount, receiverHandler);
  } else {
    auto& callRuntime = extraRuntime_.has_value() ? extraRuntime_.value() : runtime;
    static const auto TraceName = monitoring::Trace::name("receivers");
    const parallel::runtime::HostRuntime::LoopName loopName(TraceName);
    callRuntime.enqueueLoop(cellCount, receiverHandler);
  }

  countSamples(sampling.steps);
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::sampleAtRunTime([[maybe_unused]] const double* deviceClock,
                                               double hostTime,
                                               double timeStepWidth,
                                               Executor executor,
                                               parallel::runtime::StreamRuntime& runtime) {
  if (receivers_.empty()) {
    return;
  }

  const auto sampleStep = [this, timeStepWidth, executor]() {
    const double expansionPoint = *stepStart_;
    const auto sampling = planSampling(nextSampleTime_, expansionPoint, timeStepWidth);
    nextSampleTime_ = sampling.nextTime;
    if (sampling.due) {
      static const auto TraceName = monitoring::Trace::name("receivers");
      const parallel::runtime::HostRuntime::LoopName loopName(TraceName);
      parallel::runtime::HostRuntime::forEach(receiverCells_.size(), [&](std::size_t cell) {
        sampleReceiver(cell, sampling.time, expansionPoint, timeStepWidth, executor);
      });
      countSamples(sampling.steps);
    }
  };

  if (executor == Executor::Host) {
    // A cluster that runs on the host goes on to integrate right after this and overwrites the
    // DOFs the sampling reads, so it samples right here, like sample() does.
    *stepStart_ = hostTime;
    sampleStep();
    return;
  }

#ifdef ACL_DEVICE
  // wait until the samples of the last step have been taken
  if (extraRuntime_.has_value()) {
    runtime.eventSync(extraRuntime_->eventRecord());
  }
  // by a kernel, as the copy engines may be busy with a large copy of another stream
  kernels::device::aux::hostcopy::copyToHost(
      stepStartOnDevice_, deviceClock, sizeof(double), runtime.stream());
  gatherDofs(runtime);
#endif

  auto& callRuntime = extraRuntime_.has_value() ? extraRuntime_.value() : runtime;
  callRuntime.enqueueHost(sampleStep);
}

void* ReceiverCluster::sampleStream() {
  if (extraRuntime_.has_value()) {
    return extraRuntime_->stream();
  }
  return nullptr;
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::gatherDofs(parallel::runtime::StreamRuntime& runtime) {
  // in the order of the stream of the cluster, which goes on with the DOFs; the samples get taken
  // from the copies, on a stream of their own
  deviceCollector_->gatherToHost(runtime.stream());
  if constexpr (kernels::size<tensor::Qane<Cfg>>() > 0) {
    deviceCollectorAne_->gatherToHost(runtime.stream());
  }
  if (extraRuntime_.has_value()) {
    extraRuntime_->eventSync(runtime.eventRecord());
  }
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::allocateData() {
#ifdef ACL_DEVICE
  if (stepStart_ == &stepStartHost_) {
    stepStart_ = static_cast<double*>(
        ::device::DeviceInstance::instance().api().allocPinnedMem(sizeof(double)));
    stepStartOnDevice_ =
        static_cast<double*>(::device::DeviceInstance::instance().api().devicePointer(stepStart_));
  }
#endif

  // Visit the cells in storage order, so that both the host loop and the device gather read the
  // DOFs sequentially instead of in the order the receivers happened to appear in the parameter
  // file. This is only safe because a receiver does not refer back to its cell.
  std::vector<std::size_t> order(receiverCells_.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [this](std::size_t a, std::size_t b) {
    return receiverCells_[a].ltsPosition < receiverCells_[b].ltsPosition;
  });

  std::vector<ReceiverCell<Cfg>> sorted;
  sorted.reserve(receiverCells_.size());
  for (const auto cellId : order) {
    sorted.push_back(receiverCells_[cellId]);
  }
  receiverCells_ = std::move(sorted);

  if constexpr (isDeviceOn()) {
    // one entry per cell; the gather index equals the cell index
    std::vector<real*> dofs;
    dofs.reserve(receiverCells_.size());
    for (auto& receiverCell : receiverCells_) {
      dofs.push_back(receiverCell.dataDevice.template get<LTS::Dofs>());
    }

    const bool hostAccessible = useUSM() && !extraRuntime_.has_value();
    deviceCollector_ = std::make_unique<seissol::parallel::DataCollector<real>>(
        dofs, tensor::Q<Cfg>::size(), hostAccessible);

    if constexpr (kernels::size<tensor::Qane<Cfg>>() > 0) {
      std::vector<real*> dofsAne;
      dofsAne.reserve(receiverCells_.size());
      for (auto& receiverCell : receiverCells_) {
        dofsAne.push_back(receiverCell.dataDevice.template get<LTS::DofsAne>());
      }
      deviceCollectorAne_ = std::make_unique<seissol::parallel::DataCollector<real>>(
          dofsAne, kernels::size<tensor::Qane<Cfg>>(), hostAccessible);
    }
  }

  meshToReceiverCell_ = {};
}
template <typename Cfg>
void ReceiverClusterImpl<Cfg>::freeData() {
  // a handler still running would read the collector, the step start, and write the outputs
  waitForSamples();
  deviceCollector_.reset(nullptr);
  deviceCollectorAne_.reset(nullptr);
  extraRuntime_.reset();
#ifdef ACL_DEVICE
  if (stepStart_ != &stepStartHost_) {
    ::device::DeviceInstance::instance().api().freePinnedMem(stepStart_);
    stepStart_ = &stepStartHost_;
    stepStartOnDevice_ = nullptr;
  }
#endif
}

void ReceiverCluster::waitForSamples() {
  // On CUDA and SYCL, sample() and sampleAtRunTime() leave the sampling of a device cluster to a
  // host function on a stream, which appends to the output of the receivers once the gathered DOFs
  // are there. Nothing in the time stepping waits for the one enqueued last before a
  // synchronization point, so whoever reads the output has to.
  if (extraRuntime_.has_value()) {
    extraRuntime_->wait();
  }
}

template <typename Cfg>
size_t ReceiverClusterImpl<Cfg>::ncols() const {
  size_t ncols = quantities_.size();
  for (const auto& derived : derivedQuantities_) {
    ncols += derived->quantities().size();
  }
  ncols *= seissol::multisim::MultisimHelperWrapper<Cfg>::MultisimEnd -
           seissol::multisim::MultisimHelperWrapper<Cfg>::MultisimStart;
  return 1 + ncols;
}

template <typename Cfg>
std::size_t ReceiverClusterImpl<Cfg>::simulationCount() const {
  using Multisim = seissol::multisim::MultisimHelperWrapper<Cfg>;
  return Multisim::MultisimEnd - Multisim::MultisimStart;
}

template <typename Cfg>
std::vector<std::string> ReceiverClusterImpl<Cfg>::quantityNames() const {
  std::vector<std::string> names;
  names.reserve(quantities_.size());
  for (const auto quantity : quantities_) {
    names.emplace_back(seissol::model::MaterialOf<Cfg>::Quantities[quantity]);
  }
  for (const auto& derived : derivedQuantities_) {
    auto derivedNames = derived->quantities();
    names.insert(names.end(), derivedNames.begin(), derivedNames.end());
  }
  return names;
}

// The derived quantities differentiate the particle velocity. Where it sits among the quantities
// depends on the material: behind the six stresses for the solids (the solid velocity, for
// poroelastic ones), right behind the pressure for the (visco)acoustic ones.
template <typename Cfg>
VelocityGradient velocityGradient(
    const typename seissol::init::QDerivativeAtPoint<Cfg>::view::type& qDerivativeAtPoint,
    std::size_t sim) {
  static_assert(seissol::model::MaterialOf<Cfg>::VelocityOffset + Cell::Dim <=
                    tensor::QDerivativeAtPoint<Cfg>::Shape[seissol::multisim::BasisDim<Cfg>],
                "The velocity has to lie within the point derivatives.");
  VelocityGradient gradient{};
  for (std::size_t i = 0; i < Cell::Dim; ++i) {
    for (std::size_t j = 0; j < Cell::Dim; ++j) {
      gradient[i][j] = seissol::multisim::multisimWrap<Cfg>(
          qDerivativeAtPoint, sim, seissol::model::MaterialOf<Cfg>::VelocityOffset + i, j);
    }
  }
  return gradient;
}

std::vector<std::string> ReceiverRotation::quantities() const { return {"rot1", "rot2", "rot3"}; }
void ReceiverRotation::compute(std::vector<double>& output,
                               const VelocityGradient& gradient) const {
  output.push_back(gradient[2][1] - gradient[1][2]);
  output.push_back(gradient[0][2] - gradient[2][0]);
  output.push_back(gradient[1][0] - gradient[0][1]);
}

std::vector<std::string> ReceiverStrain::quantities() const {
  return {"epsxx", "epsxy", "epsxz", "epsyy", "epsyz", "epszz"};
}
void ReceiverStrain::compute(std::vector<double>& output, const VelocityGradient& gradient) const {
  // actually 9 quantities; 3 removed due to symmetry
  output.push_back(gradient[0][0]);
  output.push_back((gradient[0][1] + gradient[1][0]) / 2);
  output.push_back((gradient[0][2] + gradient[2][0]) / 2);
  output.push_back(gradient[1][1]);
  output.push_back((gradient[1][2] + gradient[2][1]) / 2);
  output.push_back(gradient[2][2]);
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template struct ReceiverBasis<Cfg>;                                                              \
  template struct ReceiverCell<Cfg>;                                                               \
  template class ReceiverClusterImpl<Cfg>;                                                         \
  template VelocityGradient velocityGradient<Cfg>(                                                 \
      const typename seissol::init::QDerivativeAtPoint<Cfg>::view::type& qDerivativeAtPoint,       \
      std::size_t sim);
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
