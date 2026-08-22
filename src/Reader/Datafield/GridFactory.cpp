// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "GridFactory.h"

#include "AsagiLite.h"
#include "Grid.h"
#include "Parallel/MPI.h"
#include "utils/logger.h"

#include <cstddef>
#include <memory>

namespace seissol::reader::datafield {

namespace {

const char* describe(AsagiLiteGrid::Error error) {
  switch (error) {
  case AsagiLiteGrid::Error::Success:
    return "success";
  case AsagiLiteGrid::Error::IoError:
    return "the file or variable could not be read";
  case AsagiLiteGrid::Error::BadRank:
    return "the variable has an unsupported shape";
  case AsagiLiteGrid::Error::BadElementType:
    return "the variable is not 32- or 64-bit floating point";
  case AsagiLiteGrid::Error::AllocationFailed:
    return "memory for the resident window could not be allocated";
  }
  return "an unknown error occurred";
}

/// Backends are constructed here and nowhere else, so the downcast is safe by
/// construction. It exists because window management is inherently backend
/// specific — a distributed or texture-backed grid would slide differently —
/// while the sampling interface is not, and putting resizeWindow()/advanceTo()
/// on Grid would force every future backend to answer a question that may not
/// apply to it.
AsagiLiteGrid& asAsagiLite(Grid& grid) { return dynamic_cast<AsagiLiteGrid&>(grid); }

} // namespace

std::unique_ptr<Grid> makeGrid(const GridDesc& desc) {
  switch (desc.kind) {
  case GridKind::AsagiLite: {
    auto grid = std::make_unique<AsagiLiteGrid>();
    grid->setComm(seissol::Mpi::mpi.comm());
    const auto error =
        grid->open(desc.path, desc.variable, desc.interpolation, desc.boundary, desc.timeAxis);
    if (error != AsagiLiteGrid::Error::Success) {
      logError() << "datafield: could not open" << desc.variable << "in" << desc.path << "--"
                 << describe(error) << ".";
    }
    return grid;
  }
  case GridKind::Scec:
    logError() << "datafield: the SCEC backend is not implemented; requested for" << desc.path
               << ".";
    return nullptr;
  case GridKind::Distributed:
    logError() << "datafield: the distributed backend is not implemented; requested for"
               << desc.path << ".";
    return nullptr;
  }
  logError() << "datafield: unknown grid backend requested for" << desc.path << ".";
  return nullptr;
}

void resizeWindow(Grid& grid, std::size_t residentSlices) {
  auto& backend = asAsagiLite(grid);
  if (backend.resizeWindow(residentSlices) != AsagiLiteGrid::Error::Success) {
    logError() << "datafield: could not allocate a" << residentSlices
               << "slice window for a time-dependent grid.";
  }
}

void advanceWindow(Grid& grid, double time) {
  auto& backend = asAsagiLite(grid);
  if (backend.advanceTo(time) != AsagiLiteGrid::Error::Success) {
    logError() << "datafield: could not slide the time window to t =" << time << ".";
  }
}

} // namespace seissol::reader::datafield
