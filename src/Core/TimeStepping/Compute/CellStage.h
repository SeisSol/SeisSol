// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLSTAGE_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLSTAGE_H_

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/StepParams.h"

#include <cstddef>

namespace seissol::solver {

/// What a cell cluster sees of a cell stage, whatever the configuration of its cells.
class CellStageBase {
  public:
  CellStageBase() = default;
  virtual ~CellStageBase() = default;

  CellStageBase(const CellStageBase&) = delete;
  CellStageBase(CellStageBase&&) = delete;
  auto operator=(const CellStageBase&) = delete;
  auto operator=(CellStageBase&&) = delete;

  /// the configuration of the cells the stage works on
  [[nodiscard]] virtual ConfigId config() const = 0;

  /// the correction of the step `params` begins, on the host, before its cells
  virtual void beginCorrection(const StepParams& /*params*/) {}

  /// the correction ends, on the host, with the sum of the counts its cells returned
  virtual void endCorrection(std::size_t /*count*/) {}

  /// the work on all cells after the correction of the step `params`, on the device; `recurring`
  /// tells whether the width of the step recurs, so that a recording of its work pays off
  virtual void correctDevice(const StepParams& /*params*/,
                             bool /*recurring*/,
                             parallel::runtime::StreamRuntime& /*runtime*/) {}

  /// a correction gets planned, on the host and also for the device, e.g. to count its FLOPs
  virtual void countCorrection() {}

  /// the end of a phase of the simulation, e.g. to collect what the device counted
  virtual void finishPhase() {}
};

/**
 * Work that a feature adds to the corrections of a cell cluster over its cells: on the host right
 * after the neighbor integral of each cell, while the data of the cell are in the cache, and on the
 * device after the neighbor integral of all cells. E.g. the plasticity.
 */
template <typename Cfg>
class CellStage : public CellStageBase {
  public:
  [[nodiscard]] ConfigId config() const final { return configIdOf<Cfg>(); }

  /// the work on the cell `cell`, with its data `data`, after its correction on the host; returns
  /// a count that the cluster sums up over the cells for endCorrection()
  virtual std::size_t correctCell(std::size_t cell, LTS::Ref<Cfg>& data) = 0;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLSTAGE_H_
