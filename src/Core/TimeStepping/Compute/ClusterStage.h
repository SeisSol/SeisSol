// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CLUSTERSTAGE_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CLUSTERSTAGE_H_

#include "Core/Common/Executor.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/StepParams.h"

#include <cstddef>
#include <string>

namespace seissol::solver {

/// When a stage runs within the step of a cell cluster.
enum class StagePhase {
  //! after the prediction of the cluster, i.e. after its local integration
  AfterPrediction,
};

/**
 * Work that a feature adds to the steps of a cell cluster over items of its own, such as point
 * sources: it runs in its phase, on the host or on the device as the cluster does, on the stream of
 * the cluster, and counts in the loop statistics under its name.
 *
 * (Work over the cells of the cluster themselves, fused into the loops of the cluster, will be a
 * kind of stage of its own.)
 */
class ClusterStage {
  public:
  ClusterStage() = default;
  virtual ~ClusterStage() = default;

  ClusterStage(const ClusterStage&) = delete;
  ClusterStage(ClusterStage&&) = delete;
  auto operator=(const ClusterStage&) = delete;
  auto operator=(ClusterStage&&) = delete;

  [[nodiscard]] virtual StagePhase phase() const = 0;

  /// the region of the loop statistics the stage counts in
  [[nodiscard]] virtual std::string name() const = 0;

  /// whether the stage has work on `executor`
  [[nodiscard]] virtual bool hasWork(Executor executor) const = 0;

  /// the number of items the stage works on with `executor`, for the loop statistics
  [[nodiscard]] virtual std::size_t size(Executor executor) const = 0;

  /// lets the work on the device read the start of each step from the clock of the cluster
  virtual void setClock(const double* /*deviceClock*/) {}

  /// the work of the step `params` with `executor`, on the stream of the cluster
  virtual void run(const StepParams& params,
                   Executor executor,
                   parallel::runtime::StreamRuntime& runtime) = 0;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CLUSTERSTAGE_H_
