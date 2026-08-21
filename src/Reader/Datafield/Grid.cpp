// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Grid.h"

#include "GridFactory.h"
#include "utils/logger.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace seissol::reader::datafield {

bool GridDesc::operator==(const GridDesc& other) const {
  return path == other.path && variable == other.variable && kind == other.kind &&
         interpolation == other.interpolation && boundary == other.boundary &&
         timeAxis == other.timeAxis;
}

void Grid::bounds(double* min, double* max) const {
  GridGeometry geom;
  geometry(geom);
  for (unsigned d = 0; d < geom.dimensions; ++d) {
    min[d] = geom.min[d];
    max[d] = geom.min[d] + static_cast<double>(geom.num[d] - 1) * geom.delta[d];
  }
}

std::size_t GridStore::intern(const GridDesc& desc) {
  for (std::size_t i = 0; i < descs_.size(); ++i) {
    if (descs_[i] == desc) {
      return i;
    }
  }
  descs_.push_back(desc);
  grids_.emplace_back();
  return descs_.size() - 1;
}

const Grid& GridStore::get(std::size_t id) const {
  if (id >= grids_.size() || grids_[id] == nullptr) {
    logError() << "datafield: grid" << id << "requested before GridStore::load().";
  }
  return *grids_[id];
}

void GridStore::setResidentSlices(std::optional<std::size_t> slices) {
  if (slices.has_value() && *slices <= MaxStencilWidth) {
    logError() << "datafield: residentSlices must exceed the widest stencil (" << MaxStencilWidth
               << "), got" << *slices << ".";
  }
  residentSlicesOverride_ = slices;
}

void GridStore::setWindowMemoryBudget(std::size_t bytes) { windowMemoryBudget_ = bytes; }

std::size_t GridStore::residentSlicesFor(std::size_t id) const {
  const GridDesc& desc = descs_.at(id);
  const std::size_t width = kernelOf(desc.interpolation).width;
  // One more than the stencil width is the narrowest useful window: it leaves
  // exactly one cell of forward headroom, so the window is reloaded on every
  // step but is never invalid.
  const std::size_t floorSlices = width + 1;

  if (!desc.timeDependent()) {
    return 1;
  }

  const Grid& grid = get(id);
  GridGeometry geom;
  grid.geometry(geom);
  const int axis = *desc.timeAxis;
  const std::size_t numSlices =
      (axis >= 0 && static_cast<unsigned>(axis) < geom.dimensions) ? geom.num[axis] : 1;

  if (residentSlicesOverride_.has_value()) {
    return std::min(*residentSlicesOverride_, std::max(numSlices, floorSlices));
  }

  std::size_t timeDependent = 0;
  for (const auto& d : descs_) {
    timeDependent += d.timeDependent() ? 1 : 0;
  }
  const std::size_t sliceBytes = grid.bytesPerTimeSlice();
  if (timeDependent == 0 || sliceBytes == 0) {
    // Backend cannot price a slice: take the narrowest correct window rather
    // than guessing at a number that would be wrong in whichever direction the
    // guess happens to fall.
    return std::min(std::max(floorSlices, std::size_t{1}), std::max(numSlices, floorSlices));
  }

  const std::size_t share = windowMemoryBudget_ / timeDependent;
  std::size_t slices = share / sliceBytes;
  if (slices < floorSlices) {
    logWarning() << "datafield:" << desc.path << "needs" << (floorSlices * sliceBytes)
                 << "bytes for the narrowest window its interpolation scheme allows, but the"
                 << "budget share is only" << share
                 << "bytes. Using the narrowest window and exceeding the budget; raise"
                 << "windowmemory or use a cheaper scheme.";
    slices = floorSlices;
  }
  return std::min(slices, std::max(numSlices, floorSlices));
}

void GridStore::setDefaultTimeSpacing(double dt) {
  if (!(dt > 0.0)) {
    logError() << "datafield: default time spacing must be positive, got" << dt << ".";
  }
  defaultTimeSpacing_ = dt;
}

std::optional<double> GridStore::timeSpacing(std::size_t id) const {
  const GridDesc& desc = descs_.at(id);
  if (!desc.timeDependent()) {
    return std::nullopt;
  }

  const int axis = *desc.timeAxis;
  GridGeometry geom;
  get(id).geometry(geom);
  if (axis >= 0 && static_cast<unsigned>(axis) < geom.dimensions) {
    const double dt = geom.delta[axis];
    if (std::isfinite(dt) && dt > 0.0) {
      return dt;
    }
  }

  // The file could not say. Fall back to the configured global value, but name
  // the file: a global constant standing in for a per-file property is exactly
  // the kind of thing that is fine until the one run where it is not.
  if (defaultTimeSpacing_ > 0.0) {
    logWarning() << "datafield:" << desc.path
                 << "has no usable spacing on its time axis; falling back to the configured"
                 << "default of" << defaultTimeSpacing_ << ". Sampling cadence for this grid is"
                 << "a guess.";
    return defaultTimeSpacing_;
  }
  logError() << "datafield:" << desc.path
             << "has no usable spacing on its time axis and no default time spacing is"
             << "configured.";
  return std::nullopt;
}

void GridStore::load() {
  // Two phases on purpose. Window sizing needs the geometry — how many slices
  // there are and what one costs — and the geometry is only known once the file
  // is open. So makeGrid() opens and reads metadata, and only then can the
  // budget be turned into a slice count and the window actually allocated.
  for (std::size_t i = 0; i < descs_.size(); ++i) {
    if (grids_[i] == nullptr) {
      grids_[i] = makeGrid(descs_[i]);
    }
  }
  for (std::size_t i = 0; i < descs_.size(); ++i) {
    if (descs_[i].timeDependent()) {
      resizeWindow(*grids_[i], residentSlicesFor(i));
    }
  }
}

void GridStore::injectForTesting(std::size_t id, std::unique_ptr<Grid> grid) {
  grids_.at(id) = std::move(grid);
}

void GridStore::update(double time) {
  for (std::size_t i = 0; i < grids_.size(); ++i) {
    if (descs_[i].timeDependent() && grids_[i] != nullptr) {
      advanceWindow(*grids_[i], time);
    }
  }
}

std::optional<double> GridStore::suggestedSyncInterval() const {
  std::optional<double> interval;
  for (std::size_t i = 0; i < descs_.size(); ++i) {
    if (!descs_[i].timeDependent()) {
      continue;
    }
    const auto dtFile = timeSpacing(i);
    if (!dtFile.has_value()) {
      continue;
    }
    const std::size_t width = kernelOf(descs_[i].interpolation).width;
    const std::size_t resident = residentSlicesFor(i);
    if (resident <= width) {
      throw std::invalid_argument("datafield: " + descs_[i].path + " needs more than " +
                                  std::to_string(width) +
                                  " resident slices for its interpolation scheme, but only " +
                                  std::to_string(resident) + " are available");
    }
    const double candidate = static_cast<double>(resident - width) * *dtFile;
    interval = interval.has_value() ? std::min(*interval, candidate) : candidate;
  }
  return interval;
}

} // namespace seissol::reader::datafield
