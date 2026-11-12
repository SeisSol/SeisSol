// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_DATAFIELD_GRIDFACTORY_H_
#define SEISSOL_SRC_READER_DATAFIELD_GRIDFACTORY_H_

// Backend construction, kept out of Grid.cpp so that GridStore does not drag
// HDF5 and MPI into every translation unit that schedules a sampling update —
// and so that the store can be unit-tested against a synthetic backend.

#include "Grid.h"

#include <cstddef>
#include <memory>
#include <string>

namespace seissol::reader::datafield {

/// Opens the backend named by `desc` and reads its metadata. Collective over the
/// rank's communicator. A static grid is fully loaded here; a time-dependent one
/// is not, because how much of it to load is a question only answerable once the
/// geometry is known — see resizeWindow().
std::unique_ptr<Grid> makeGrid(const GridDesc& desc);

/// The SHA-256 of the file at `path` (expr::sha256OfFile()), as rank 0 of the rank's communicator
/// reads it, given to every rank. Collective where MPI runs. Empty where the file cannot be read,
/// or where SEISSOL_EXPR_GRID_HASH is off.
[[nodiscard]] std::string fileHashOf(const std::string& path);

/// Allocates and fills the resident time window. Collective. Called by
/// GridStore::load() after the budget has been turned into a slice count. A no-op
/// for a grid whose backend keeps no window.
void resizeWindow(Grid& grid, std::size_t residentSlices);

/// Moves a time-dependent grid's resident window so that it covers `time`.
/// Idempotent and cheap when the window already does — which is the normal case,
/// since it is called at every synchronisation point.
void advanceWindow(Grid& grid, double time);

/// Closes the file of the grid while HDF5 is still there -- at the end of the run
/// (GridStore::finish()), since the store may live until the process ends, past
/// HDF5. The grid keeps its array and samples on: a script evaluated after the end
/// of the simulation (the analysis of the initial condition, say) reads it. Its
/// window does not move any more.
void closeGridFile(Grid& grid);

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_GRIDFACTORY_H_
