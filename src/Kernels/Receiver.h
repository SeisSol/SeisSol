// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_RECEIVER_H_
#define SEISSOL_SRC_KERNELS_RECEIVER_H_

#include "Common/Executor.h"
#include "GeneratedCode/init.h"
#include "Geometry/CellTransform.h"
#include "Geometry/MeshReader.h"
#include "Initializer/PointMapper.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Interface.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Backmap.h"
#include "Monitoring/Metric.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Transformation.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Runtime/Stream.h"

#include <Eigen/Dense>
#include <optional>
#include <unordered_map>
#include <vector>

namespace seissol {
class SeisSol;

namespace kernels {
struct Receiver {
  Receiver(std::size_t pointId,
           Eigen::Vector3d position,
           const seissol::geometry::CellTransform& transform,
           size_t reserved);
  std::size_t pointId;
  Eigen::Vector3d position;
  basisFunction::SampledBasisFunctions<real> basisFunctions;
  basisFunction::SampledBasisFunctionDerivatives<real> basisFunctionDerivatives;
  std::vector<real> output;
};

/**
  A cell carrying at least one receiver. The time evaluation runs once per cell, and only the
  point evaluation is repeated for each receiver in it.
 */
struct ReceiverCell {
  ReceiverCell(std::size_t meshId, LTS::Ref dataHost, LTS::Ref dataDevice);
  std::size_t meshId{};
  std::size_t ltsPosition{};
  LTS::Ref dataHost;
  LTS::Ref dataDevice;
  std::vector<std::size_t> receiverIds;
};

struct DerivedReceiverQuantity {
  virtual ~DerivedReceiverQuantity() = default;
  [[nodiscard]] virtual std::vector<std::string> quantities() const = 0;
  virtual void compute(size_t sim,
                       std::vector<real>&,
                       seissol::init::QAtPoint::view::type&,
                       seissol::init::QDerivativeAtPoint::view::type&) = 0;
};

struct ReceiverRotation : public DerivedReceiverQuantity {
  ~ReceiverRotation() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override;
  void compute(size_t sim,
               std::vector<real>& /*output*/,
               seissol::init::QAtPoint::view::type& /*qAtPoint*/,
               seissol::init::QDerivativeAtPoint::view::type& /*qDerivativeAtPoint*/) override;
};

struct ReceiverStrain : public DerivedReceiverQuantity {
  ~ReceiverStrain() override = default;
  [[nodiscard]] std::vector<std::string> quantities() const override;
  void compute(size_t sim,
               std::vector<real>& /*output*/,
               seissol::init::QAtPoint::view::type& /*qAtPoint*/,
               seissol::init::QDerivativeAtPoint::view::type& /*qDerivativeAtPoint*/) override;
};

class ReceiverCluster {
  public:
  explicit ReceiverCluster(seissol::SeisSol& seissolInstance);

  ReceiverCluster(const CompoundGlobalData& global,
                  const std::vector<std::size_t>& quantities,
                  double samplingInterval,
                  double syncPointInterval,
                  const std::vector<std::shared_ptr<DerivedReceiverQuantity>>& derivedQuantities,
                  seissol::SeisSol& seissolInstance);

  void addReceiver(std::size_t meshId,
                   std::size_t pointId,
                   const Eigen::Vector3d& point,
                   const seissol::geometry::MeshReader& mesh,
                   const LTS::Backmap& backmap);

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
  [[nodiscard]] Sampling
      planSampling(double time, double expansionPoint, double timeStepWidth) const;

  /**
   * Takes the samples of a step as planned.
   */
  void sample(const Sampling& sampling,
              double expansionPoint,
              double timeStepWidth,
              Executor executor,
              parallel::runtime::StreamRuntime& runtime);

  /**
   * Takes the samples of the step that starts at the time on the clock, and decides about them
   * only when the enqueued work runs; the next sample time is then kept here. Replaying the
   * enqueued work (e.g. from a graph) thus takes the samples of the step it runs in. On the device,
   * the receiver data gets copied to the host in every step; on the host, the step starts at
   * `hostTime`.
   */
  void sampleAtRunTime(const double* deviceClock,
                       double hostTime,
                       double timeStepWidth,
                       Executor executor,
                       parallel::runtime::StreamRuntime& runtime);

  /**
   * Sets the next sample time for sampleAtRunTime().
   */
  void setNextSampleTime(double time) { nextSampleTime_ = time; }

  std::vector<Receiver>::iterator begin() { return receivers_.begin(); }

  std::vector<Receiver>::iterator end() { return receivers_.end(); }

  [[nodiscard]] size_t ncols() const;

  void allocateData();
  void freeData();

  //! @brief Waits for the samples taken so far to be in the output of the receivers.
  void waitForSamples();

  private:
  std::optional<parallel::runtime::StreamRuntime> extraRuntime_;
  std::unique_ptr<seissol::parallel::DataCollector<real>> deviceCollector_{nullptr};
  // anelastic DOFs (LinearCKAnelastic only); their host copy is stale between sync points or,
  // with USM, written by the device concurrently
  std::unique_ptr<seissol::parallel::DataCollector<real>> deviceCollectorAne_{nullptr};
  std::vector<Receiver> receivers_;
  std::vector<ReceiverCell> receiverCells_;
  std::unordered_map<std::size_t, std::size_t> meshToReceiverCell_;
  seissol::kernels::Spacetime spacetimeKernel_;
  seissol::kernels::Time timeKernel_;
  std::vector<std::size_t> quantities_;
  PerformanceEstimate estimatePerCell_{};
  PerformanceEstimate estimatePerCellStep_{};
  PerformanceEstimate estimatePerPoint_{};
  std::size_t perfHandle_{};
  double samplingInterval_;

  /**
   * Samples the receivers of one receiver cell in the step that starts at `expansionPoint`, from
   * the sample time `time` on.
   */
  void sampleReceiver(std::size_t cell,
                      double time,
                      double expansionPoint,
                      double timeStepWidth,
                      Executor executor);

  /// counts the operations of sampling all receiver cells in `samplingSteps` sample times
  void countSamples(std::size_t samplingSteps);

  // for sampleAtRunTime(): the next sample time, and the start of the step the samples are taken in
  double nextSampleTime_{0};
  double stepStartHost_{0};
  double* stepStart_{&stepStartHost_};
  double syncPointInterval_;
  std::vector<std::shared_ptr<DerivedReceiverQuantity>> derivedQuantities_;
  seissol::SeisSol& seissolInstance_;
};
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_RECEIVER_H_
