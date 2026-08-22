// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_GRIDUPDATEMODULE_H_
#define SEISSOL_SRC_READER_DATAFIELD_GRIDUPDATEMODULE_H_

#include "Modules/Module.h"
#include "Reader/Datafield/Grid.h"

#include <optional>

namespace seissol::reader::datafield {

/// Advances the resident window of every time-dependent grid so it still covers
/// the current simulated time.
///
/// The module carries no cadence of its own. GridStore::update() is idempotent
/// and cheap when nothing needs reloading, and the interval comes from
/// suggestedSyncInterval(), which derives it from the file's own time spacing
/// and the memory budget. There is deliberately no user-facing sampling
/// interval, for the reason Grid.h gives: a second, independent way to get the
/// cadence wrong.
class GridUpdateModule : public Module {
  public:
  explicit GridUpdateModule(GridStore& store) : store_(&store) {}

  /// Bring the windows up to date before the first step, including after a
  /// restart, where the first time the solver asks for data is not t = 0.
  void simulationStart(std::optional<double> checkpointTime) override;

  void syncPoint(double currentTime) override;

  /// Register for the two hooks and adopt `interval`. A member because
  /// Module::setSyncInterval is protected -- the framework's way of saying that
  /// a module's cadence is its own business and not something a registrar sets
  /// from outside.
  void enable(double interval);

  private:
  GridStore* store_;
};

/// Register a GridUpdateModule for the shared store, if and only if some grid
/// in it is time-dependent.
///
/// ORDERING: call this AFTER every model has been loaded. The store is filled by
/// makeKernel interning a program's grids, which happens when a CompiledReader
/// prepares -- so before the models are read, the store is empty, every grid is
/// trivially static, and this would decline to register a module that turns out
/// to be needed. Registering unconditionally instead would be worse: a module
/// with no time-dependent grid still forces a synchronisation point on every
/// interval, which costs a global barrier for nothing.
///
/// Idempotent: calling it twice registers one module.
void registerGridUpdateModule();

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_GRIDUPDATEMODULE_H_
