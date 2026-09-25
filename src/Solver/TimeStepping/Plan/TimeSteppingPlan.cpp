// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "TimeSteppingPlan.h"

#include "Solver/TimeStepping/Actor/ActorState.h"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

namespace seissol::time_stepping {

namespace {

enum class Phase { Correction = 0, Prediction = 1, FaceWork = 2 };

struct Key {
  long time;
  Phase phase;
  int priority;
  std::size_t cluster;
  long step;

  bool operator<(const Key& other) const {
    return std::tie(time, phase, priority, cluster, step) <
           std::tie(other.time, other.phase, other.priority, other.cluster, other.step);
  }
};

} // namespace

std::vector<PlannedAction> planTimeSteps(const std::vector<PlannedCluster>& clusters) {
  std::vector<std::pair<Key, PlannedAction>> actions;
  for (std::size_t cluster = 0; cluster < clusters.size(); ++cluster) {
    const auto& info = clusters[cluster];
    const auto priority = info.priority == ActorPriority::High ? 0 : 1;
    const bool face = info.dataReadiness == DataReadiness::AfterCorrection;
    const auto steps = (info.stepsUntilSync + info.timeStepRate - 1) / info.timeStepRate;
    for (long step = 0; step < steps; ++step) {
      const auto start = step * info.timeStepRate;
      const auto end = std::min(start + info.timeStepRate, info.stepsUntilSync);
      actions.emplace_back(Key{start, Phase::Prediction, priority, cluster, step},
                           PlannedAction{cluster, ActorAction::Predict, step});
      const auto correction = face ? Key{start, Phase::FaceWork, priority, cluster, step}
                                   : Key{end, Phase::Correction, priority, cluster, step};
      actions.emplace_back(correction, PlannedAction{cluster, ActorAction::Correct, step});
    }
  }
  std::sort(actions.begin(), actions.end(), [](const auto& a, const auto& b) {
    return a.first < b.first;
  });

  std::vector<PlannedAction> plan;
  plan.reserve(actions.size());
  for (const auto& [key, action] : actions) {
    plan.push_back(action);
  }
  return plan;
}

} // namespace seissol::time_stepping
