// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_FEATURE_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_FEATURE_H_

#include "Core/IO/Instance/Checkpoint/CheckpointManager.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/TimeStepping/Compute/FaceClusterFactory.h"

namespace seissol::solver {

class TimeManager;

/**
 * A feature of a simulation, such as the plasticity or the dynamic rupture: what it adds to the
 * storage of the cells and to the clusters of the time stepping. The composition of a run makes
 * its features from the parameters, and calls each hook on all of them in their order.
 */
class Feature {
  public:
  Feature() = default;
  virtual ~Feature() = default;

  Feature(const Feature&) = delete;
  Feature(Feature&&) = delete;
  auto operator=(const Feature&) = delete;
  auto operator=(Feature&&) = delete;

  /// the name of the feature, for the log
  [[nodiscard]] virtual const char* name() const = 0;

  /// adds the variables of the feature to the storage of the cells, before it gets allocated
  virtual void addStorage(LTS::Storage& /*cells*/) {}

  /// the clusters of faces the feature makes, if any; asked for when the clusters get added
  [[nodiscard]] virtual FaceClusterFactory* faceClusters() { return nullptr; }

  /// adds the work of the feature to the clusters of the cells of `cells`, once they exist
  virtual void addStages(TimeManager& /*timeManager*/, LTS::Storage& /*cells*/) {}

  /// registers what the checkpoints hold of the feature in `checkpoint`, after the variables of
  /// the cells of `cells` (the tree "lts"), e.g. its own variables of the cells or a tree of its
  /// own
  virtual void registerCheckpoint(io::instance::checkpoint::CheckpointManager& /*checkpoint*/,
                                  LTS::Storage& /*cells*/) {}
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_FEATURE_H_
