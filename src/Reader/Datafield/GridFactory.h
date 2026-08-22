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

namespace seissol::reader::datafield {

/// Opens the backend named by `desc` and reads its metadata. Collective over the
/// rank's communicator. A static grid is fully loaded here; a time-dependent one
/// is not, because how much of it to load is a question only answerable once the
/// geometry is known — see resizeWindow().
std::unique_ptr<Grid> makeGrid(const GridDesc& desc);

/// Allocates and fills the resident time window. Collective. Called by
/// GridStore::load() after the budget has been turned into a slice count.
void resizeWindow(Grid& grid, std::size_t residentSlices);

/// Moves a time-dependent grid's resident window so that it covers `time`.
/// Idempotent and cheap when the window already does — which is the normal case,
/// since it is called at every synchronisation point.
void advanceWindow(Grid& grid, double time);

} // namespace seissol::reader::datafield

#endif // SEISSOL_SRC_READER_DATAFIELD_GRIDFACTORY_H_
