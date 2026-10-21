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
#include "Initializer/MemoryManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Metric.h"
#include "Numerical/BasisFunction.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Helper.h"
#include "Parallel/Runtime/Stream.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"
#include "Solver/Settings.h"

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
                                std::vector<LTS::Ref<Cfg>> dataHost,
                                std::vector<LTS::Ref<Cfg>> dataDevice)
    : meshId(meshId), dataHost(std::move(dataHost)), dataDevice(std::move(dataDevice)) {}

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

  const auto& settings = seissolInstance_.memoryManager().simulationSettings();
  simulationCount_ = settings.simulationCount;
  copyCount_ = settings.copyCount();

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
    auto& ltsStorage = seissolInstance_.memoryManager().ltsStorage();

    meshToReceiverCell_[meshId] = receiverCells_.size();

    std::vector<LTS::Ref<Cfg>> dataHost;
    std::vector<LTS::Ref<Cfg>> dataDevice;
    for (std::size_t copy = 0; copy < copyCount_; ++copy) {
      const auto position = backmap.get(meshId, copy);
      dataHost.emplace_back(ltsStorage.template lookupRef<Cfg>(position));
      dataDevice.emplace_back(
          ltsStorage.template lookupRef<Cfg>(position,
                                             isDeviceOn() ? initializer::AllocationPlace::Device
                                                          : initializer::AllocationPlace::Host));
    }

    auto& cell = receiverCells_.emplace_back(meshId, std::move(dataHost), std::move(dataDevice));
    cell.ltsPosition = backmap.get(meshId, 0).global;
  }

  receiverCells_[meshToReceiverCell_.at(meshId)].receiverIds.emplace_back(receivers_.size());

  receivers_.emplace_back(pointId, point, reserved);
  receiverBases_.emplace_back(point, transform);
}

template <typename Cfg>
double ReceiverClusterImpl<Cfg>::calcReceivers(double time,
                                               double expansionPoint,
                                               double timeStepWidth,
                                               Executor executor,
                                               parallel::runtime::StreamRuntime& runtime) {
  using Multisim = seissol::multisim::MultisimHelperWrapper<Cfg>;

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
    if constexpr (kernels::size<tensor::Qane<Cfg>>() > 0) {
      deviceCollectorAne_->gatherToHost(runtime.stream());
    }
    if (extraRuntime_.has_value()) {
      extraRuntime_->eventSync(runtime.eventRecord());
    }
  }

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();

  if (time >= expansionPoint && time < expansionPoint + timeStepWidth) {
    const std::size_t cellCount = receiverCells_.size();
    const auto receiverHandler =
        [this, timeBasis, timeStepWidth, time, expansionPoint, executor, samplingSteps](
            std::size_t i) {
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

          auto& receiverCell = receiverCells_[i];

          // A sample takes the time, then the quantities of one simulation after the other. The
          // copies of the cell fill in the simulations they hold, so the samples of this time step
          // are laid out first. (only the handler of the cell touches its receivers)
          const auto columns = ncols();
          const auto columnsPerSimulation = (columns - 1) / simulationCount_;
          std::vector<std::size_t> firstSample(receiverCell.receiverIds.size());
          for (std::size_t j = 0; j < receiverCell.receiverIds.size(); ++j) {
            auto& output = receivers_[receiverCell.receiverIds[j]].output;
            firstSample[j] = output.size();
            output.resize(firstSample[j] + samplingSteps * columns);
          }

          std::vector<double> derivedValues;

          for (std::size_t copy = 0; copy < receiverCell.dataHost.size(); ++copy) {
            // Use device pointers where required.
            auto tmpReceiverData{receiverCell.dataHost[copy]};

            if (executor == Executor::Device) {
              const auto gatherIndex = i * receiverCell.dataHost.size() + copy;
              tmpReceiverData.template setPointer<LTS::Dofs>(
                  reinterpret_cast<decltype(tmpReceiverData.template getPointer<LTS::Dofs>())>(
                      deviceCollector_->get(gatherIndex)));
              if constexpr (kernels::size<tensor::Qane<Cfg>>() > 0) {
                tmpReceiverData.template setPointer<LTS::DofsAne>(
                    reinterpret_cast<decltype(tmpReceiverData.template getPointer<LTS::DofsAne>())>(
                        deviceCollectorAne_->get(gatherIndex)));
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
            for (std::size_t step = 0; step < samplingSteps; ++step) {
              const auto coeffs = timeBasis.point(receiverTime - expansionPoint, timeStepWidth);

              timeKernel_.evaluate(coeffs.data(), timeDerivatives, timeEvaluated);

              for (std::size_t j = 0; j < receiverCell.receiverIds.size(); ++j) {
                const auto receiverId = receiverCell.receiverIds[j];
                auto& receiver = receivers_[receiverId];
                const auto& basis = receiverBases_[receiverId];

                krnl.basisFunctionsAtPoint = runtime::init::basisFunctionsAtPoint::view(
                    Variant, basis.basisFunctions.data().data());
                derivativeKrnl.basisFunctionDerivativesAtPoint =
                    runtime::init::basisFunctionDerivativesAtPoint::view(
                        Variant, basis.basisFunctionDerivatives.data().data());

                krnl.execute(Variant);
                derivativeKrnl.execute(Variant);

                auto* sample = receiver.output.data() + firstSample[j] + step * columns;
                sample[0] = receiverTime;
                for (auto fused = Multisim::MultisimStart; fused < Multisim::MultisimEnd; ++fused) {
                  const auto sim = copy * Cfg::NumSimulations + (fused - Multisim::MultisimStart);
                  if (sim >= simulationCount_) {
                    // the fused simulations past the last one only fill up the last copy
                    break;
                  }
                  auto* target = sample + 1 + sim * columnsPerSimulation;
                  for (auto quantity : quantities_) {
                    const auto value =
                        seissol::multisim::multisimWrap<Cfg>(qAtPoint, fused, quantity);
                    if (!std::isfinite(value)) {
                      logError() << "Detected Inf/NaN in receiver output at" << receiver.position[0]
                                 << "," << receiver.position[1] << "," << receiver.position[2]
                                 << " in simulation" << sim << "."
                                 << "Aborting.";
                    }
                    *target++ = value;
                  }
                  if (!derivedQuantities_.empty()) {
                    const auto gradient = velocityGradient<Cfg>(qDerivativeAtPoint, fused);
                    derivedValues.clear();
                    for (const auto& derived : derivedQuantities_) {
                      derived->compute(derivedValues, gradient);
                    }
                    std::copy(derivedValues.begin(), derivedValues.end(), target);
                  }
                }
              }

              receiverTime += samplingInterval_;
            }
          }
        };

    if (executor == Executor::Host) {
      // A cluster that runs on the host goes on to integrate right after this and overwrites the
      // DOFs the sampling reads, so it samples right here. (On CUDA and SYCL, enqueueLoop would
      // leave the sampling to a host function on a stream.)
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < cellCount; ++i) {
        receiverHandler(i);
      }
    } else {
      auto& callRuntime = extraRuntime_.has_value() ? extraRuntime_.value() : runtime;
      callRuntime.enqueueLoop(cellCount, receiverHandler);
    }

    // every copy of a cell is evaluated
    const auto copyCellCount = cellCount * copyCount_;
    const auto recvCount = receivers_.size() * copyCount_;

    seissolInstance_.flopCounter().incrementMetric(
        perfHandle_,
        estimatePerCell_ * copyCellCount + estimatePerCellStep_ * copyCellCount * samplingSteps +
            estimatePerPoint_ * recvCount * samplingSteps);
  }
  return outReceiverTime;
}

template <typename Cfg>
void ReceiverClusterImpl<Cfg>::allocateData() {
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
    // one entry per cell and copy; the gather index is the cell index times the copy count plus
    // the copy
    std::vector<real*> dofs;
    dofs.reserve(receiverCells_.size() * copyCount_);
    for (auto& receiverCell : receiverCells_) {
      for (auto& data : receiverCell.dataDevice) {
        dofs.push_back(data.template get<LTS::Dofs>());
      }
    }

    const bool hostAccessible = useUSM() && !extraRuntime_.has_value();
    deviceCollector_ = std::make_unique<seissol::parallel::DataCollector<real>>(
        dofs, tensor::Q<Cfg>::size(), hostAccessible);

    if constexpr (kernels::size<tensor::Qane<Cfg>>() > 0) {
      std::vector<real*> dofsAne;
      dofsAne.reserve(receiverCells_.size() * copyCount_);
      for (auto& receiverCell : receiverCells_) {
        for (auto& data : receiverCell.dataDevice) {
          dofsAne.push_back(data.template get<LTS::DofsAne>());
        }
      }
      deviceCollectorAne_ = std::make_unique<seissol::parallel::DataCollector<real>>(
          dofsAne, kernels::size<tensor::Qane<Cfg>>(), hostAccessible);
    }
  }

  meshToReceiverCell_ = {};
}
template <typename Cfg>
void ReceiverClusterImpl<Cfg>::freeData() {
  // a handler still running would read the collector and write the outputs
  waitForSamples();
  deviceCollector_.reset(nullptr);
  deviceCollectorAne_.reset(nullptr);
  extraRuntime_.reset();
}

void ReceiverCluster::waitForSamples() {
  // On CUDA and SYCL, calcReceivers leaves the sampling of a device cluster to a host function on a
  // stream, which appends to the output of the receivers once the gathered DOFs are there. Nothing
  // in the time stepping waits for the one enqueued last before a synchronization point, so whoever
  // reads the output has to.
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
  ncols *= simulationCount_;
  return 1 + ncols;
}

template <typename Cfg>
std::size_t ReceiverClusterImpl<Cfg>::simulationCount() const {
  return simulationCount_;
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
