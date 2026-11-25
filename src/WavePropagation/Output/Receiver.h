// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_RECEIVER_H_
#define SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_RECEIVER_H_

#include "Core/Common/Executor.h"
#include "Core/Common/Real.h"
#include "Core/Geometry/CellTransform.h"
#include "Core/Geometry/MeshDefinition.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Geometry/PointMapper.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/DataCollector.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/SimulationSettings.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Numerical/BasisFunction.h"
#include "Core/Numerical/Transformation.h"
#include "Core/Parallel/Runtime/HostActor.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Run/Context.h"
#include "GeneratedCode/init.h"
#include "WavePropagation/Solver/Interface.h"
#include "WavePropagation/Solver/Solver.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::kernels {

/// The cells the receivers of a run lie in, and how they compute: their storage and the global
/// data of their configurations, the simulations, and the gravitational acceleration at the free
/// surface.
struct ReceiverCells {
  LTS::Storage& storage;
  GlobalDataStore& globalData;
  SimulationSettings settings;
  double gravitationalAcceleration{};
};

/// What a receiver records, whatever the configuration of the cell it sits in.
struct Receiver {
  Receiver(std::size_t pointId, Eigen::Vector3d position, size_t reserved);
  std::size_t pointId;
  Eigen::Vector3d position;
  std::vector<double> output;
};

/// The basis functions of the cell of a receiver, and their derivatives, at the receiver; for the
/// configuration `Cfg` of the cell.
template <typename Cfg>
struct ReceiverBasis {
  ReceiverBasis(const Eigen::Vector3d& position, const seissol::geometry::CellTransform& transform);
  basis_function::SampledBasisFunctions<Real<Cfg>> basisFunctions;
  basis_function::SampledBasisFunctionDerivatives<Real<Cfg>> basisFunctionDerivatives;
};

/**
  A cell carrying at least one receiver. The time evaluation runs once per cell and copy, and only
  the point evaluation is repeated for each receiver in it.
 */
template <typename Cfg>
struct ReceiverCell {
  ReceiverCell(std::size_t meshId,
               std::vector<LTS::Ref<Cfg>> dataHost,
               std::vector<LTS::Ref<Cfg>> dataDevice);
  std::size_t meshId{};
  std::size_t ltsPosition{};
  // the data of the copies of the cell, which hold the simulations of the run
  std::vector<LTS::Ref<Cfg>> dataHost;
  std::vector<LTS::Ref<Cfg>> dataDevice;
  std::vector<std::size_t> receiverIds;
};

/// The velocity gradient of one simulation at a point, gradient[i][j] = d_j v_i.
using VelocityGradient = std::array<std::array<double, Cell::Dim>, Cell::Dim>;

/// The velocity gradient of the simulation `sim` at a point, from the derivatives of the quantities
/// there in the configuration `Cfg`; in double, as the receivers record what is derived from it.
template <typename Cfg>
VelocityGradient velocityGradient(
    const typename seissol::init::QDerivativeAtPoint<Cfg>::view::type& qDerivativeAtPoint,
    std::size_t sim);

/// A quantity that receivers record in addition to the quantities of the material, derived from
/// the velocity gradient.
struct DerivedReceiverQuantity {
  virtual ~DerivedReceiverQuantity() = default;
  [[nodiscard]] virtual std::vector<std::string> quantities() const = 0;
  /// Appends the values of the quantities for one simulation to `output`.
  virtual void compute(std::vector<double>& output, const VelocityGradient& gradient) const = 0;
};

struct ReceiverRotation : public DerivedReceiverQuantity {
  ~ReceiverRotation() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override;
  void compute(std::vector<double>& output, const VelocityGradient& gradient) const override;
};

struct ReceiverStrain : public DerivedReceiverQuantity {
  ~ReceiverStrain() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override;
  void compute(std::vector<double>& output, const VelocityGradient& gradient) const override;
};

/**
  The receivers of a time cluster, whatever the configuration of its cells: what the time cluster
  and the receiver writer see of them.
 */
class ReceiverCluster {
  public:
  virtual ~ReceiverCluster() = default;

  virtual void addReceiver(std::size_t meshId,
                           std::size_t pointId,
                           const Eigen::Vector3d& point,
                           const seissol::geometry::MeshReader& mesh,
                           const LTS::Backmap& backmap) = 0;

  /**
   * The receiver samples of one time step.
   */
  struct Sampling {
    /// whether samples fall into the step
    bool due{false};
    /// the first sample time
    double time{0};
    /// the next sample time after the step
    double nextTime{0};
    /// the number of sample times up to the end of the step
    std::size_t steps{0};
  };

  /**
   * Determines which samples fall into the step that starts at `expansionPoint`, given the next
   * sample time `time`; only on the host, without touching any data.
   */
  [[nodiscard]] virtual Sampling
      planSampling(double time, double expansionPoint, double timeStepWidth) const = 0;

  /**
   * Takes the samples of a step as planned.
   */
  virtual void sample(const Sampling& sampling,
                      double expansionPoint,
                      double timeStepWidth,
                      Executor executor,
                      parallel::runtime::StreamRuntime& runtime) = 0;

  /**
   * Takes the samples of the step that starts at the time on the clock, and decides about them
   * only when the enqueued work runs; the next sample time is then kept here. Replaying the
   * enqueued work (e.g. from a graph) thus takes the samples of the step it runs in. On the device,
   * the receiver data gets copied to the host in every step; on the host, the step starts at
   * `hostTime`.
   */
  virtual void sampleAtRunTime(const double* deviceClock,
                               double hostTime,
                               double timeStepWidth,
                               Executor executor,
                               parallel::runtime::StreamRuntime& runtime) = 0;

  /**
   * Sets the next sample time for sampleAtRunTime().
   */
  void setNextSampleTime(double time) { nextSampleTime_ = time; }

  std::vector<Receiver>::iterator begin() { return receivers_.begin(); }

  std::vector<Receiver>::iterator end() { return receivers_.end(); }

  [[nodiscard]] virtual size_t ncols() const = 0;

  //! @brief The number of simulations a sample holds: all simulations of the run.
  [[nodiscard]] virtual std::size_t simulationCount() const = 0;

  //! @brief What one simulation records, without the time. A sample takes the time, then this for
  //! one simulation after the other.
  [[nodiscard]] virtual std::vector<std::string> quantityNames() const = 0;

  virtual void allocateData() = 0;
  virtual void freeData() = 0;

  //! @brief Waits for the samples taken so far to be in the output of the receivers.
  void waitForSamples();

  /**
   * The stream runtime the samples of device clusters are taken on, if it is not the one of the
   * cluster; null otherwise.
   */
  [[nodiscard]] parallel::runtime::StreamRuntime* sampleRuntime();

  protected:
  std::optional<parallel::runtime::StreamRuntime> extraRuntime_;
  // with the host actor (instead of the extra runtime): the samples of device clusters, taken once
  // the device has gathered what they read
  std::unique_ptr<parallel::runtime::HostActor::Queue> hostQueue_;
  std::vector<Receiver> receivers_;

  // for sampleAtRunTime(): the next sample time
  double nextSampleTime_{0};
};

/// The receivers of a time cluster whose cells compute in the configuration `Cfg`.
template <typename Cfg>
class ReceiverClusterImpl : public ReceiverCluster {
  public:
  using real = Real<Cfg>;

  ReceiverClusterImpl(
      const CompoundGlobalData<Cfg>& global,
      const std::vector<std::size_t>& quantities,
      double samplingInterval,
      double syncPointInterval,
      const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
      RunContext& run,
      const ReceiverCells& cells);

  void addReceiver(std::size_t meshId,
                   std::size_t pointId,
                   const Eigen::Vector3d& point,
                   const seissol::geometry::MeshReader& mesh,
                   const LTS::Backmap& backmap) override;

  [[nodiscard]] Sampling
      planSampling(double time, double expansionPoint, double timeStepWidth) const override;

  void sample(const Sampling& sampling,
              double expansionPoint,
              double timeStepWidth,
              Executor executor,
              parallel::runtime::StreamRuntime& runtime) override;

  void sampleAtRunTime(const double* deviceClock,
                       double hostTime,
                       double timeStepWidth,
                       Executor executor,
                       parallel::runtime::StreamRuntime& runtime) override;

  [[nodiscard]] size_t ncols() const override;

  [[nodiscard]] std::size_t simulationCount() const override;

  [[nodiscard]] std::vector<std::string> quantityNames() const override;

  void allocateData() override;
  void freeData() override;

  private:
  std::unique_ptr<seissol::parallel::DataCollector<real>> deviceCollector_{nullptr};
  // anelastic DOFs (LinearCKAnelastic only); their host copy is stale between sync points or,
  // with USM, written by the device concurrently
  std::unique_ptr<seissol::parallel::DataCollector<real>> deviceCollectorAne_{nullptr};
  std::vector<ReceiverBasis<Cfg>> receiverBases_;
  std::vector<ReceiverCell<Cfg>> receiverCells_;
  std::unordered_map<std::size_t, std::size_t> meshToReceiverCell_;
  seissol::kernels::Spacetime<Cfg> spacetimeKernel_;
  seissol::kernels::Time<Cfg> timeKernel_;
  std::vector<std::size_t> quantities_;
  PerformanceEstimate estimatePerCell_{};
  PerformanceEstimate estimatePerCellStep_{};
  PerformanceEstimate estimatePerPoint_{};
  std::size_t perfHandle_{};
  double samplingInterval_;

  /**
   * Samples the receivers of one receiver cell in the step that starts at `expansionPoint`, from
   * the sample time `time` on. A device cluster reads the DOFs from `gathered` (a buffer that
   * gatherForHost() has filled) if given, and from the buffers of its collectors otherwise.
   */
  void sampleReceiver(std::size_t cell,
                      double time,
                      double expansionPoint,
                      double timeStepWidth,
                      Executor executor,
                      const void* gathered = nullptr);

  /// samples all receiver cells as sampleReceiver() does, and counts it
  void sampleAll(const Sampling& sampling,
                 double expansionPoint,
                 double timeStepWidth,
                 Executor executor,
                 const void* gathered);

  /**
   * For the host actor: gathers the DOFs of the receiver cells, and the start of the step from the
   * device clock if given, into a buffer of their own (see HostActor::Queue::take()) on the stream
   * of the cluster; the host reads them from there once the device is done.
   */
  parallel::runtime::HostActor::Staging gatherForHost(parallel::runtime::StreamRuntime& runtime,
                                                      const double* deviceClock);

  /// counts the operations of sampling all receiver cells in `samplingSteps` sample times
  void countSamples(std::size_t samplingSteps);

  /// copies the DOFs of the receiver cells to the host on the stream of the cluster, for the stream
  /// of the samples
  void gatherDofs(parallel::runtime::StreamRuntime& runtime);

  // with the host actor: the buffers gatherForHost() fills hold the start of the step at their
  // beginning, then the DOFs at these offsets
  std::size_t dofsOffset_{0};
  std::size_t dofsAneOffset_{0};

  // for sampleAtRunTime(): the start of the step the samples are taken in, and its address on the
  // device, which writes it
  double stepStartHost_{0};
  double* stepStart_{&stepStartHost_};
  double* stepStartOnDevice_{nullptr};
  double syncPointInterval_;
  // the simulations of the run, and the copies of the cells which hold them
  std::size_t simulationCount_{Cfg::NumSimulations};
  std::size_t copyCount_{1};
  std::vector<std::shared_ptr<DerivedReceiverQuantity>> derivedQuantities_;
  RunContext& run_;
  LTS::Storage& storage_;
  double gravitationalAcceleration_;
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_RECEIVER_H_
