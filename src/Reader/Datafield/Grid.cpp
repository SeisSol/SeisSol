// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Grid.h"

#include "GridFactory.h"
#include "InlineGrid.h"
#include "utils/logger.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace seissol::reader::datafield {

bool GridData::operator==(const GridData& other) const {
  const auto same = [](const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() &&
           (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  };
  if (knots.size() != other.knots.size() || !same(values, other.values)) {
    return false;
  }
  for (std::size_t d = 0; d < knots.size(); ++d) {
    if (!same(knots[d], other.knots[d])) {
      return false;
    }
  }
  return true;
}

bool GridDesc::operator==(const GridDesc& other) const {
  const bool sameData =
      data == other.data || (data != nullptr && other.data != nullptr && *data == *other.data);
  return path == other.path && variable == other.variable && kind == other.kind &&
         interpolation == other.interpolation && boundary == other.boundary &&
         interpolations == other.interpolations && boundaries == other.boundaries &&
         axes == other.axes && sameData && timeAxis == other.timeAxis;
}

bool GridDesc::failsAnywhere() const {
  return boundary == BoundaryMode::Fail ||
         std::find(boundaries.begin(), boundaries.end(), BoundaryMode::Fail) != boundaries.end();
}

namespace {

const char* token(GridKind kind) {
  switch (kind) {
  case GridKind::AsagiLite:
    return "asagilite";
  case GridKind::Scec:
    return "scec";
  case GridKind::Distributed:
    return "distributed";
  case GridKind::Inline:
    return "inline";
  }
  return "unknown";
}

const char* token(Interpolation interp) { return spelling(interp); }
const char* token(BoundaryMode boundary) { return spelling(boundary); }

/// FNV-1a over the bits of `values`.
void hashInto(std::uint64_t& hash, const std::vector<double>& values) {
  for (const double value : values) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    for (int byte = 0; byte < 8; ++byte) {
      hash ^= (bits >> (8 * byte)) & 0xffU;
      hash *= 0x100000001b3ULL;
    }
  }
}

} // namespace

const char* spelling(GridKind kind) { return kind == GridKind::AsagiLite ? "asagi" : token(kind); }

const char* spelling(Interpolation interpolation) {
  switch (interpolation) {
  case Interpolation::Nearest:
    return "nearest";
  case Interpolation::Linear:
    return "linear";
  case Interpolation::Cubic:
    return "cubic";
  case Interpolation::Previous:
    return "previous";
  case Interpolation::Next:
    return "next";
  }
  return "unknown";
}

const char* spelling(BoundaryMode boundary) {
  switch (boundary) {
  case BoundaryMode::Clamp:
    return "clamp";
  case BoundaryMode::Fail:
    return "fail";
  }
  return "unknown";
}

const char* spelling(AxisKind axis) { return axis == AxisKind::Knots ? "knots" : "uniform"; }

std::optional<GridKind> parseGridKind(const std::string& text) {
  if (text == "asagi" || text == "asagilite") {
    return GridKind::AsagiLite;
  }
  if (text == "scec") {
    return GridKind::Scec;
  }
  if (text == "inline") {
    return GridKind::Inline;
  }
  return std::nullopt;
}

std::optional<Interpolation> parseInterpolation(const std::string& text) {
  for (const auto interpolation : {Interpolation::Nearest,
                                   Interpolation::Linear,
                                   Interpolation::Cubic,
                                   Interpolation::Previous,
                                   Interpolation::Next}) {
    if (text == spelling(interpolation)) {
      return interpolation;
    }
  }
  return std::nullopt;
}

std::optional<BoundaryMode> parseBoundaryMode(const std::string& text) {
  if (text == "clamp") {
    return BoundaryMode::Clamp;
  }
  if (text == "fail") {
    return BoundaryMode::Fail;
  }
  return std::nullopt;
}

std::optional<AxisKind> parseAxisKind(const std::string& text) {
  if (text == "uniform") {
    return AxisKind::Uniform;
  }
  if (text == "knots") {
    return AxisKind::Knots;
  }
  return std::nullopt;
}

GridDesc gridDescOf(const GridSpec& spec) {
  const auto fail = [](const std::string& what) { throw std::invalid_argument(what); };
  GridDesc desc;
  const auto kind = parseGridKind(spec.kind);
  if (!kind.has_value()) {
    fail("the kind of a grid is \"asagi\", \"scec\" or \"inline\", not \"" + spec.kind + "\"");
  }
  desc.kind = *kind;
  if (spec.components.empty()) {
    fail("a grid names its components");
  }
  desc.componentNames = spec.components;
  const bool inlined = desc.kind == GridKind::Inline;
  // an inline grid has its data, a file its path, and none both
  if (inlined) {
    if (!spec.path.empty() || !spec.variable.empty()) {
      fail("an inline grid has its knots and values, and no path or variable");
    }
    if (!spec.axes.empty()) {
      fail("every axis of an inline grid is knots; it gives no axes");
    }
    if (spec.timeAxis.has_value()) {
      fail("an inline grid has no time axis");
    }
    auto data = std::make_shared<GridData>();
    data->knots = spec.knots;
    data->values = spec.values;
    if (const auto problem = InlineGrid::problemOf(*data, spec.components.size());
        !problem.empty()) {
      fail("the data of the inline grid: " + problem);
    }
    desc.data = std::move(data);
  } else {
    if (!spec.knots.empty() || !spec.values.empty()) {
      fail("a grid of a file reads its knots and values from the file");
    }
    desc.path = spec.path;
    desc.variable = spec.variable;
    desc.timeAxis = spec.timeAxis;
    for (const auto& axis : spec.axes) {
      const auto parsed = parseAxisKind(axis);
      if (!parsed.has_value()) {
        fail("an axis is \"uniform\" or \"knots\", not \"" + axis + "\"");
      }
      desc.axes.push_back(*parsed);
    }
    if (spec.timeAxis.has_value() && static_cast<std::size_t>(*spec.timeAxis) < desc.axes.size() &&
        desc.axes[static_cast<std::size_t>(*spec.timeAxis)] == AxisKind::Knots) {
      fail("the time axis of a grid is uniform");
    }
    // a SCEC file has a uniform grid of the fault, and no time
    if (desc.kind == GridKind::Scec) {
      if (spec.timeAxis.has_value()) {
        fail("a SCEC file has no time axis");
      }
      for (const auto axis : desc.axes) {
        if (axis == AxisKind::Knots) {
          fail("the axes of a SCEC file are uniform");
        }
      }
    }
  }
  // one for all axes, or one per axis, of those it gives (a file says how many it has)
  const std::size_t axes = inlined ? spec.knots.size() : spec.axes.size();
  const auto perAxis = [&](const std::vector<std::string>& given, const char* what) {
    if (given.empty()) {
      fail(std::string("a grid gives its ") + what);
    }
    if (given.size() > 1 && axes > 0 && given.size() != axes) {
      fail(std::string("a grid of ") + std::to_string(axes) + " axes gives " +
           std::to_string(given.size()) + " " + what + "s");
    }
  };
  perAxis(spec.interpolation, "interpolation");
  perAxis(spec.boundary, "boundary");
  for (const auto& text : spec.interpolation) {
    const auto parsed = parseInterpolation(text);
    if (!parsed.has_value()) {
      fail("the interpolation of a grid is \"nearest\", \"linear\", \"cubic\", \"previous\" or "
           "\"next\", not \"" +
           text + "\"");
    }
    desc.interpolations.push_back(*parsed);
  }
  for (const auto& text : spec.boundary) {
    const auto parsed = parseBoundaryMode(text);
    if (!parsed.has_value()) {
      fail("the boundary of a grid is \"clamp\" or \"fail\", not \"" + text + "\"");
    }
    desc.boundaries.push_back(*parsed);
  }
  desc.interpolation = desc.interpolations.front();
  desc.boundary = desc.boundaries.front();
  if (desc.interpolations.size() == 1) {
    desc.interpolations.clear();
  }
  if (desc.boundaries.size() == 1) {
    desc.boundaries.clear();
  }
  for (std::size_t d = 0; d < axes || d < desc.axes.size(); ++d) {
    if (desc.axisOf(d) == AxisKind::Knots && desc.interpolationOf(d) == Interpolation::Cubic) {
      fail("an axis of knots interpolates nearest, linear, previous or next, not cubic");
    }
  }
  return desc;
}

std::string GridDesc::canonicalKey() const {
  // The time axis is part of the key because it changes the sampling arity the
  // kernel emits, not merely which slices are resident.
  std::string key =
      std::string(token(kind)) + '\x1f' + token(interpolation) + '\x1f' + token(boundary) + '\x1f';
  key += timeAxis.has_value() ? ('t' + std::to_string(*timeAxis)) : std::string("static");
  // per axis, where it is said per axis
  for (const auto each : interpolations) {
    key += '\x1f' + std::string("i:") + token(each);
  }
  for (const auto each : boundaries) {
    key += '\x1f' + std::string("b:") + token(each);
  }
  for (const auto each : axes) {
    key += '\x1f' + std::string("a:") + spelling(each);
  }
  if (data != nullptr) {
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    for (const auto& knots : data->knots) {
      hashInto(hash, knots);
      hash ^= 0xff;
      hash *= 0x100000001b3ULL;
    }
    hashInto(hash, data->values);
    key += '\x1f' + std::string("data:") + std::to_string(data->knots.size()) + ':' +
           std::to_string(data->values.size()) + ':' + std::to_string(hash);
  }
  return key;
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
  hashes_.emplace_back();
  windowSlices_.push_back(0);
  mirrors_.emplace_back();
  return descs_.size() - 1;
}

const Grid& GridStore::get(std::size_t id) const {
  if (id >= grids_.size() || grids_[id] == nullptr) {
    logError() << "datafield: grid" << id << "requested before GridStore::load().";
  }
  return *grids_[id];
}

void GridStore::setResidentSlices(std::optional<std::size_t> slices) {
  // Whether a window is wide enough depends on the interpolation scheme of each grid, which
  // suggestedSyncInterval() checks; only a window that no scheme can sample is refused here.
  if (slices.has_value() && *slices < 2) {
    throw std::invalid_argument(
        "datafield: a time window needs at least two resident slices, got " +
        std::to_string(*slices));
  }
  residentSlicesOverride_ = slices;
}

void GridStore::setWindowMemoryBudget(std::size_t bytes) { windowMemoryBudget_ = bytes; }

std::size_t GridStore::residentSlicesFor(std::size_t id) const {
  const GridDesc& desc = descs_.at(id);
  // The narrowest useful window leaves exactly one cell of forward headroom, so
  // the window is reloaded on every step but is never invalid.
  const std::size_t floorSlices = narrowestWindow(
      kernelOf(desc.interpolationOf(static_cast<std::size_t>(desc.timeAxis.value_or(0)))));

  if (!desc.timeDependent()) {
    return 1;
  }

  const Grid& grid = get(id);
  // never more than the file: a window that holds all of it never moves (suggestedSyncInterval())
  const std::size_t numSlices = timeSlicesOf(id);

  if (residentSlicesOverride_.has_value()) {
    return std::min(*residentSlicesOverride_, numSlices);
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
    return std::min(floorSlices, numSlices);
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
  return std::min(slices, numSlices);
}

std::size_t GridStore::timeSlicesOf(std::size_t id) const {
  const GridDesc& desc = descs_.at(id);
  if (!desc.timeDependent()) {
    return 1;
  }
  // the file, not the window that geometry() describes
  const std::size_t slices = get(id).timeSlices();
  if (slices != 0) {
    return slices;
  }
  GridGeometry geom;
  get(id).geometry(geom);
  const int axis = *desc.timeAxis;
  return (axis >= 0 && static_cast<unsigned>(axis) < geom.dimensions) ? geom.num[axis] : 1;
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

GridStore& sharedGridStore() {
  // Function-local static: constructed on first use, after MPI is up, which a
  // namespace-scope object could not promise.
  static GridStore store;
  return store;
}

void GridStore::load() {
  // Two phases on purpose. Window sizing needs the geometry — how many slices
  // there are and what one costs — and the geometry is only known once the file
  // is open. So makeGrid() opens and reads metadata, and only then can the
  // budget be turned into a slice count and the window actually allocated.
  //
  // Both phases skip what is done already: opening and sizing are collective, and
  // a call that finds nothing to do must not reach the communicator, or a kernel
  // made on some ranks only would wait for the others forever (see Grid.h). A
  // window is sized again where its share of the budget changed, i.e. where a
  // time-dependent grid was added -- on all ranks alike, as the scripts are.
  for (std::size_t i = 0; i < descs_.size(); ++i) {
    if (grids_[i] == nullptr) {
      grids_[i] = makeGrid(descs_[i]);
      // what the evaluation records name the file by; an inline grid is in its program
      if (descs_[i].kind != GridKind::Inline) {
        hashes_[i] = fileHashOf(descs_[i].path);
      }
    }
  }
  for (std::size_t i = 0; i < descs_.size(); ++i) {
    if (descs_[i].timeDependent()) {
      const std::size_t slices = residentSlicesFor(i);
      if (windowSlices_[i] != slices) {
        resizeWindow(*grids_[i], slices);
        windowSlices_[i] = slices;
      }
    }
  }
  refreshMirrors();
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
  refreshMirrors();
}

GridStore::Mirror* GridStore::mirror(std::size_t id, int key) const {
  for (const auto& [kind, mirror] : mirrors_.at(id)) {
    if (kind == key) {
      return mirror.get();
    }
  }
  return nullptr;
}

GridStore::Mirror&
    GridStore::attachMirror(std::size_t id, int key, std::unique_ptr<Mirror> mirror) {
  auto& mirrors = mirrors_.at(id);
  for (auto& [kind, attached] : mirrors) {
    if (kind == key) {
      attached = std::move(mirror);
      return *attached;
    }
  }
  mirrors.emplace_back(key, std::move(mirror));
  return *mirrors.back().second;
}

void GridStore::refreshMirrors() {
  for (std::size_t i = 0; i < grids_.size(); ++i) {
    if (grids_[i] == nullptr) {
      continue;
    }
    for (auto& [kind, mirror] : mirrors_[i]) {
      mirror->refresh(*grids_[i]);
    }
  }
}

void GridStore::finish() {
  for (std::size_t i = 0; i < grids_.size(); ++i) {
    if (grids_[i] == nullptr) {
      continue;
    }
    for (auto& [kind, mirror] : mirrors_[i]) {
      mirror->finish(*grids_[i]);
    }
    // the file of the grid now, while HDF5 is there: the store may live until the process ends;
    // the array stays, for what samples it after the end of the simulation
    closeGridFile(*grids_[i]);
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
    const std::size_t resident = residentSlicesFor(i);
    if (resident >= timeSlicesOf(i)) {
      // the window holds the whole file and never moves: any interval is safe
      continue;
    }
    // the stencil and the slices to spare on either side of it
    const std::size_t span =
        kernelOf(descs_[i].interpolationOf(static_cast<std::size_t>(*descs_[i].timeAxis))).width +
        2 * TimeWindowSlack;
    if (resident <= span) {
      throw std::invalid_argument("datafield: " + descs_[i].path + " needs more than " +
                                  std::to_string(span) +
                                  " resident slices for its interpolation scheme, but only " +
                                  std::to_string(resident) + " are available");
    }
    const double candidate = static_cast<double>(resident - span) * *dtFile;
    interval = interval.has_value() ? std::min(*interval, candidate) : candidate;
  }
  return interval;
}

} // namespace seissol::reader::datafield
