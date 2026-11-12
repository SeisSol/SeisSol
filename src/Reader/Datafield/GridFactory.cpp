// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "GridFactory.h"

#include "AsagiLite.h"
#include "Expr/Sha256.h"
#include "Grid.h"
#include "InlineGrid.h"
#include "Parallel/MPI.h"
#include "ScecGrid.h"
#include "utils/logger.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

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
  case AsagiLiteGrid::Error::BadKnots:
    return "its axes are not as declared (see the warning before)";
  }
  return "an unknown error occurred";
}

/// The SCEC file `path`, as rank 0 of the rank's communicator reads it, given to every rank (where
/// MPI runs); nothing, with an error, where it cannot be read.
std::optional<ScecData> scecDataOf(const std::string& path) {
  int running = 0;
  MPI_Initialized(&running);
  const bool reads = running == 0 || seissol::Mpi::mpi.rank() == 0;
  ScecData data;
  std::string problem;
  if (reads) {
    try {
      data = readScecFile(path);
    } catch (const std::exception& error) {
      problem = error.what();
    }
  }
  if (running != 0) {
    const MPI_Comm comm = seissol::Mpi::mpi.comm();
    int length = static_cast<int>(problem.size());
    MPI_Bcast(&length, 1, MPI_INT, 0, comm);
    problem.resize(static_cast<std::size_t>(length));
    MPI_Bcast(problem.data(), length, MPI_CHAR, 0, comm);
    if (problem.empty()) {
      MPI_Bcast(data.nodes.data(), 2, MPI_UNSIGNED, 0, comm);
      MPI_Bcast(data.lengths.data(), 2, MPI_DOUBLE, 0, comm);
      data.values.resize(std::size_t{data.nodes[0]} * data.nodes[1] * ScecData::Components);
      // in pieces an int counts
      constexpr std::size_t Piece = std::size_t{1} << 28U;
      for (std::size_t at = 0; at < data.values.size(); at += Piece) {
        const auto count = static_cast<int>(std::min(Piece, data.values.size() - at));
        MPI_Bcast(data.values.data() + at, count, MPI_DOUBLE, 0, comm);
      }
    }
  }
  if (!problem.empty()) {
    logError() << "datafield:" << problem;
    return std::nullopt;
  }
  return data;
}

/// The backend with a time window, or null. Window management is inherently
/// backend specific -- a distributed or texture-backed grid would slide
/// differently -- while the sampling interface is not, and putting
/// resizeWindow()/advanceTo() on Grid would force every future backend to answer
/// a question that may not apply to it. A grid without a window (such as one
/// injected by a test) has nothing to resize or slide.
AsagiLiteGrid* asAsagiLite(Grid& grid) { return dynamic_cast<AsagiLiteGrid*>(&grid); }

} // namespace

std::unique_ptr<Grid> makeGrid(const GridDesc& desc) {
  switch (desc.kind) {
  case GridKind::AsagiLite: {
    auto grid = std::make_unique<AsagiLiteGrid>();
    grid->setComm(seissol::Mpi::mpi.comm());
    grid->setLookupsOnly(!desc.interpolations.empty() || !desc.boundaries.empty());
    // the window moves by the stencil of the time axis
    const Interpolation interpolation =
        desc.timeAxis.has_value() ? desc.interpolationOf(static_cast<std::size_t>(*desc.timeAxis))
                                  : desc.interpolation;
    const auto error = grid->open(
        desc.path, desc.variable, interpolation, desc.boundary, desc.timeAxis, desc.axes);
    if (error != AsagiLiteGrid::Error::Success) {
      logError() << "datafield: could not open" << desc.variable << "in" << desc.path << "--"
                 << describe(error) << ".";
    }
    // an interpolation or a boundary per axis: one per axis of the file
    const auto perAxis = [&](std::size_t given, const char* what) {
      if (given > 1 && given != grid->dimensions()) {
        logError() << "datafield: the grid of" << desc.variable << "in" << desc.path << "gives"
                   << given << what << "for" << grid->dimensions() << "axes.";
      }
    };
    perAxis(desc.interpolations.size(), "interpolations");
    perAxis(desc.boundaries.size(), "boundaries");
    return grid;
  }
  case GridKind::Inline: {
    if (desc.data == nullptr || desc.data->knots.empty()) {
      logError() << "datafield: an inline grid without data.";
      return nullptr;
    }
    std::size_t nodes = 1;
    for (const auto& knots : desc.data->knots) {
      nodes *= std::max<std::size_t>(1, knots.size());
    }
    try {
      return std::make_unique<InlineGrid>(desc.data, desc.data->values.size() / nodes);
    } catch (const std::invalid_argument& error) {
      logError() << error.what();
    }
    return nullptr;
  }
  case GridKind::Scec: {
    // rank 0 reads the file, and gives it to the others
    auto data = scecDataOf(desc.path);
    if (!data.has_value()) {
      return nullptr;
    }
    return std::make_unique<ScecGrid>(
        std::move(*data), desc.path, desc.interpolation, desc.boundary);
  }
  case GridKind::Distributed:
    logError() << "datafield: the distributed backend is not implemented; requested for"
               << desc.path << ".";
    return nullptr;
  }
  logError() << "datafield: unknown grid backend requested for" << desc.path << ".";
  return nullptr;
}

std::string fileHashOf(const std::string& path) {
  if (const char* env = std::getenv("SEISSOL_EXPR_GRID_HASH"); env != nullptr) {
    const std::string value(env);
    if (value == "off" || value == "0" || value == "none") {
      return {};
    }
  }
  int running = 0;
  MPI_Initialized(&running);
  if (running == 0) {
    return expr::sha256OfFile(path);
  }
  // 64 digits and the end, or nothing
  std::array<char, 65> hash{};
  if (seissol::Mpi::mpi.rank() == 0) {
    const std::string hex = expr::sha256OfFile(path);
    std::copy_n(hex.c_str(), std::min<std::size_t>(hex.size(), 64), hash.data());
  }
  MPI_Bcast(hash.data(), static_cast<int>(hash.size()), MPI_CHAR, 0, seissol::Mpi::mpi.comm());
  return hash.data();
}

void resizeWindow(Grid& grid, std::size_t residentSlices) {
  auto* backend = asAsagiLite(grid);
  if (backend == nullptr) {
    return;
  }
  if (backend->resizeWindow(residentSlices) != AsagiLiteGrid::Error::Success) {
    logError() << "datafield: could not allocate a" << residentSlices
               << "slice window for a time-dependent grid.";
  }
}

void closeGridFile(Grid& grid) {
  auto* backend = asAsagiLite(grid);
  if (backend != nullptr) {
    backend->closeFile();
  }
}

void advanceWindow(Grid& grid, double time) {
  auto* backend = asAsagiLite(grid);
  if (backend == nullptr) {
    return;
  }
  if (backend->advanceTo(time) != AsagiLiteGrid::Error::Success) {
    logError() << "datafield: could not slide the time window to t =" << time << ".";
  }
}

} // namespace seissol::reader::datafield
