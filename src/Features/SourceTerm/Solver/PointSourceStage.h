// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_SOURCETERM_SOLVER_POINTSOURCESTAGE_H_
#define SEISSOL_SRC_FEATURES_SOURCETERM_SOLVER_POINTSOURCESTAGE_H_

#include "Core/Common/Executor.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/ClusterStage.h"
#include "Features/SourceTerm/Solver/PointSourceCluster.h"

#include <cstddef>
#include <string>
#include <utility>

namespace seissol::kernels {

/// The point sources of a cell cluster, added to its cells after each prediction.
class PointSourceStage : public seissol::solver::ClusterStage {
  public:
  explicit PointSourceStage(PointSourceClusterPair sources) : sources_(std::move(sources)) {}

  [[nodiscard]] seissol::solver::StagePhase phase() const override {
    return seissol::solver::StagePhase::AfterPrediction;
  }

  [[nodiscard]] std::string name() const override { return "computePointSources"; }

  [[nodiscard]] bool hasWork(Executor executor) const override {
    return sourcesOn(executor) != nullptr;
  }

  [[nodiscard]] std::size_t size(Executor executor) const override {
    const auto* sources = sourcesOn(executor);
    return sources == nullptr ? 0 : sources->size();
  }

  void setClock(const double* deviceClock) override {
    if (sources_.device != nullptr) {
      sources_.device->setClock(deviceClock);
    }
  }

  void run(const seissol::solver::StepParams& params,
           Executor executor,
           parallel::runtime::StreamRuntime& runtime) override {
    if (auto* sources = sourcesOn(executor); sources != nullptr) {
      sources->addTimeIntegratedPointSources(params.time, params.timeStepSize, runtime);
    }
  }

  private:
  [[nodiscard]] PointSourceCluster* sourcesOn(Executor executor) const {
    return executor == Executor::Device ? sources_.device.get() : sources_.host.get();
  }

  PointSourceClusterPair sources_;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_FEATURES_SOURCETERM_SOLVER_POINTSOURCESTAGE_H_
