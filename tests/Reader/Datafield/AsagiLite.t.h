// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Reader/Datafield/AsagiLite.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <netcdf.h>
#include <optional>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace seissol::unit_test {

namespace {

namespace df = reader::datafield;

/// A NetCDF file of z(t, x) = t + 10 x at t = 0, 1, ... and x = 0, 1: linear in both, which linear
/// interpolation reproduces. Time is the slowest dimension of the file, as AsagiLiteGrid wants it,
/// so the last axis of the grid. Removed at the end of the scope.
class LinearInTime {
  public:
  explicit LinearInTime(std::size_t slices)
      : path_("seissol-test-asagilite-" + std::to_string(::getpid()) + ".nc") {
    int file = 0;
    REQUIRE(nc_create(path_.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
    int dims[2] = {};
    REQUIRE(nc_def_dim(file, "t", slices, &dims[0]) == NC_NOERR);
    REQUIRE(nc_def_dim(file, "x", 2, &dims[1]) == NC_NOERR);
    int times = 0;
    int points = 0;
    int values = 0;
    REQUIRE(nc_def_var(file, "t", NC_DOUBLE, 1, &dims[0], &times) == NC_NOERR);
    REQUIRE(nc_def_var(file, "x", NC_DOUBLE, 1, &dims[1], &points) == NC_NOERR);
    REQUIRE(nc_def_var(file, "z", NC_DOUBLE, 2, dims, &values) == NC_NOERR);
    REQUIRE(nc_enddef(file) == NC_NOERR);
    std::vector<double> t(slices);
    std::vector<double> z(slices * 2);
    for (std::size_t i = 0; i < slices; ++i) {
      t[i] = static_cast<double>(i);
      z[i * 2] = static_cast<double>(i);
      z[i * 2 + 1] = static_cast<double>(i) + 10.0;
    }
    const std::vector<double> x = {0.0, 1.0};
    REQUIRE(nc_put_var_double(file, times, t.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(file, points, x.data()) == NC_NOERR);
    REQUIRE(nc_put_var_double(file, values, z.data()) == NC_NOERR);
    REQUIRE(nc_close(file) == NC_NOERR);
  }
  ~LinearInTime() { std::remove(path_.c_str()); }
  LinearInTime(const LinearInTime&) = delete;
  LinearInTime& operator=(const LinearInTime&) = delete;
  LinearInTime(LinearInTime&&) = delete;
  LinearInTime& operator=(LinearInTime&&) = delete;

  [[nodiscard]] const std::string& path() const { return path_; }

  private:
  std::string path_;
};

} // namespace

TEST_CASE("Datafield AsagiLite: a sliding window serves the times until the next update") {
  // The update module advances the window at most (W - w - 2) cells apart
  // (GridStore::suggestedSyncInterval()); every time up to the next update samples at full
  // accuracy, and never a slice clamped onto the end of the window.
  constexpr std::size_t Slices = 40;
  constexpr std::size_t Resident = 6;
  const LinearInTime file(Slices);
  df::AsagiLiteGrid grid;
  REQUIRE(grid.open(file.path(), "z", df::Interpolation::Linear, df::BoundaryMode::Clamp, 1) ==
          df::AsagiLiteGrid::Error::Success);
  REQUIRE(grid.resizeWindow(Resident) == df::AsagiLiteGrid::Error::Success);
  REQUIRE(grid.arrayView() != nullptr);

  const double interval = static_cast<double>(Resident - 2 - 2);
  std::uint64_t revision = grid.revision();
  for (double update = 0.0; update + interval <= static_cast<double>(Slices - 1);
       update += interval) {
    CAPTURE(update);
    REQUIRE(grid.advanceTo(update) == df::AsagiLiteGrid::Error::Success);
    if (update > 0.0) {
      // a window that moved is a new array for a copy of it
      CHECK(grid.revision() != revision);
    }
    revision = grid.revision();
    // the view counts the slices of the file around the window (at t = 0, 1, ...), against which a
    // written-out lookup checks that the window serves a time; the other axes are all resident
    const df::ArrayView* view = grid.arrayView();
    REQUIRE(view != nullptr);
    CHECK(static_cast<double>(view->before[1]) == view->min[1]);
    CHECK(view->before[1] + view->num[1] + view->after[1] == Slices);
    CHECK(view->before[0] == 0);
    CHECK(view->after[0] == 0);
    for (double time = update; time <= update + interval; time += 0.25) {
      CAPTURE(time);
      const double point[2] = {0.5, time};
      double value = 0;
      grid.sample(point, &value);
      CHECK(value == time + 5.0);
    }
  }
  // an update that stays in its cell reads nothing, and the copies stay current
  REQUIRE(grid.advanceTo(2.0) == df::AsagiLiteGrid::Error::Success);
  revision = grid.revision();
  REQUIRE(grid.advanceTo(2.5) == df::AsagiLiteGrid::Error::Success);
  CHECK(grid.revision() == revision);
}

TEST_CASE("Datafield AsagiLite: a store moves the window of a file it holds in part") {
  // The geometry of the grid is that of its window; the store sizes the window, and asks for
  // synchronization points, by the slices of the file (Grid::timeSlices()).
  const LinearInTime file(12);
  df::GridDesc desc;
  desc.path = file.path();
  desc.variable = "z";
  desc.interpolation = df::Interpolation::Linear;
  desc.timeAxis = 1;
  for (const std::size_t resident : {6, 12, 20}) {
    CAPTURE(resident);
    df::GridStore store;
    store.setResidentSlices(resident);
    const std::size_t id = store.intern(desc);
    store.load();
    CHECK(store.timeSlicesOf(id) == 12);
    CHECK(store.residentSlicesFor(id) == std::min<std::size_t>(resident, 12));
    const auto interval = store.suggestedSyncInterval();
    if (resident < 12) {
      // (W - w - 2) cells of the file's time step of 1
      REQUIRE(interval.has_value());
      CHECK(*interval == doctest::Approx(static_cast<double>(resident - 2 - 2)));
    } else {
      // the whole file is resident, and the window never moves
      CHECK_FALSE(interval.has_value());
    }
    store.finish();
  }
}

TEST_CASE("Datafield AsagiLite: the grids of a store sample on after the end of the run") {
  // The end of the simulation closes the files (GridStore::finish()); what samples the grids after
  // it -- the analysis of the initial condition evaluates its script -- still finds them, and
  // loading again finds nothing to do.
  const LinearInTime file(12);
  df::GridStore store;
  store.setResidentSlices(4);
  df::GridDesc desc;
  desc.path = file.path();
  desc.variable = "z";
  desc.interpolation = df::Interpolation::Linear;
  desc.timeAxis = 1;
  const std::size_t id = store.intern(desc);
  store.load();
  store.update(2.0);
  store.finish();
  CHECK_NOTHROW(store.load());
  const double point[2] = {1.0, 2.5};
  double value = 0;
  store.get(id).sample(point, &value);
  CHECK(value == 2.5 + 10.0);
}

TEST_CASE("Datafield AsagiLite: the dimensions of the file and the attributes of their axes") {
  const std::string path = "seissol-test-asagilite-axes-" + std::to_string(::getpid()) + ".nc";
  int file = 0;
  REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
  int dims[2] = {};
  REQUIRE(nc_def_dim(file, "depth", 3, &dims[0]) == NC_NOERR);
  REQUIRE(nc_def_dim(file, "x", 2, &dims[1]) == NC_NOERR);
  int depth = 0;
  int x = 0;
  int values = 0;
  REQUIRE(nc_def_var(file, "depth", NC_DOUBLE, 1, &dims[0], &depth) == NC_NOERR);
  REQUIRE(nc_def_var(file, "x", NC_DOUBLE, 1, &dims[1], &x) == NC_NOERR);
  REQUIRE(nc_def_var(file, "rho", NC_DOUBLE, 2, dims, &values) == NC_NOERR);
  // text of a fixed length, and a string of a variable length
  REQUIRE(nc_put_att_text(file, x, "units", 1, "m") == NC_NOERR);
  REQUIRE(nc_put_att_text(file, x, "standard_name", 23, "projection_x_coordinate") == NC_NOERR);
  REQUIRE(nc_put_att_text(file, depth, "units", 2, "km") == NC_NOERR);
  const char* positive = "down";
  REQUIRE(nc_put_att_string(file, depth, "positive", 1, &positive) == NC_NOERR);
  REQUIRE(nc_enddef(file) == NC_NOERR);
  const std::vector<double> depths = {0.0, 1.0, 2.0};
  const std::vector<double> xs = {0.0, 1.0};
  const std::vector<double> rho(6, 2.0);
  REQUIRE(nc_put_var_double(file, depth, depths.data()) == NC_NOERR);
  REQUIRE(nc_put_var_double(file, x, xs.data()) == NC_NOERR);
  REQUIRE(nc_put_var_double(file, values, rho.data()) == NC_NOERR);
  REQUIRE(nc_close(file) == NC_NOERR);

  df::AsagiLiteGrid grid;
  REQUIRE(
      grid.open(path, "rho", df::Interpolation::Linear, df::BoundaryMode::Clamp, std::nullopt) ==
      df::AsagiLiteGrid::Error::Success);
  // in the order of the axes of the grid, the opposite of the dimensions of the file
  const auto axes = grid.axisInfo();
  std::remove(path.c_str());
  REQUIRE(axes.size() == 2);
  CHECK(axes[0].dimension == "x");
  CHECK(axes[0].units == "m");
  CHECK(axes[0].standardName == "projection_x_coordinate");
  CHECK(axes[0].positive.empty());
  CHECK(axes[1].dimension == "depth");
  CHECK(axes[1].units == "km");
  CHECK(axes[1].positive == "down");
  CHECK(axes[1].standardName.empty());
}

TEST_CASE("Datafield AsagiLite: an axis declared knots reads its coordinates whole") {
  // rho(depth, x) over the depths 0, 100 and 1000, which are not uniform, and x = 0 and 1; and
  // up(down, x) over an axis that decreases
  const std::string path = "seissol-test-asagilite-knots-" + std::to_string(::getpid()) + ".nc";
  int file = 0;
  REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
  int dims[3] = {};
  REQUIRE(nc_def_dim(file, "depth", 3, &dims[0]) == NC_NOERR);
  REQUIRE(nc_def_dim(file, "x", 2, &dims[1]) == NC_NOERR);
  REQUIRE(nc_def_dim(file, "down", 3, &dims[2]) == NC_NOERR);
  int depth = 0;
  int x = 0;
  int down = 0;
  int rho = 0;
  int up = 0;
  REQUIRE(nc_def_var(file, "depth", NC_DOUBLE, 1, &dims[0], &depth) == NC_NOERR);
  REQUIRE(nc_def_var(file, "x", NC_DOUBLE, 1, &dims[1], &x) == NC_NOERR);
  REQUIRE(nc_def_var(file, "down", NC_DOUBLE, 1, &dims[2], &down) == NC_NOERR);
  REQUIRE(nc_def_var(file, "rho", NC_DOUBLE, 2, dims, &rho) == NC_NOERR);
  const int upDims[2] = {dims[2], dims[1]};
  REQUIRE(nc_def_var(file, "up", NC_DOUBLE, 2, upDims, &up) == NC_NOERR);
  REQUIRE(nc_enddef(file) == NC_NOERR);
  const std::vector<double> depths = {0.0, 100.0, 1000.0};
  const std::vector<double> xs = {0.0, 1.0};
  const std::vector<double> downs = {3.0, 2.0, 1.0};
  const std::vector<double> values = {2000.0, 2100.0, 2500.0, 2650.0, 2900.0, 3300.0};
  REQUIRE(nc_put_var_double(file, depth, depths.data()) == NC_NOERR);
  REQUIRE(nc_put_var_double(file, x, xs.data()) == NC_NOERR);
  REQUIRE(nc_put_var_double(file, down, downs.data()) == NC_NOERR);
  REQUIRE(nc_put_var_double(file, rho, values.data()) == NC_NOERR);
  REQUIRE(nc_put_var_double(file, up, values.data()) == NC_NOERR);
  REQUIRE(nc_close(file) == NC_NOERR);

  using df::AxisKind;
  const auto open =
      [&](df::AsagiLiteGrid& grid, const std::string& variable, const std::vector<AxisKind>& axes) {
        return grid.open(
            path, variable, df::Interpolation::Linear, df::BoundaryMode::Clamp, std::nullopt, axes);
      };
  {
    // the grid axes in the opposite order of the dimensions of the file: x, then depth
    df::AsagiLiteGrid grid;
    REQUIRE(open(grid, "rho", {AxisKind::Uniform, AxisKind::Knots}) ==
            df::AsagiLiteGrid::Error::Success);
    const df::ArrayView* view = grid.arrayView();
    REQUIRE(view != nullptr);
    CHECK(view->knots[0] == nullptr);
    REQUIRE(view->knots[1] != nullptr);
    REQUIRE(view->inverseSpacings[1] != nullptr);
    CHECK(view->num[1] == 3);
    CHECK(view->knots[1][0] == 0.0);
    CHECK(view->knots[1][1] == 100.0);
    CHECK(view->knots[1][2] == 1000.0);
    CHECK(view->inverseSpacings[1][0] == 1.0 / 100.0);
    CHECK(view->inverseSpacings[1][1] == 1.0 / 900.0);
    CHECK(view->min[1] == 0.0);
    // knots have no one spacing: the mean
    df::GridGeometry geometry;
    grid.geometry(geometry);
    CHECK(geometry.delta[1] == 500.0);
    CHECK(geometry.delta[0] == 1.0);
    // only the written-out lookup searches knots
    const double point[2] = {0.5, 50.0};
    double out = 0.0;
    CHECK_THROWS_AS(grid.sample(point, &out), std::invalid_argument);
  }
  {
    // undeclared, the axis is uniform, by the mean spacing between its ends (as ASAGI has it)
    df::AsagiLiteGrid grid;
    REQUIRE(open(grid, "rho", {}) == df::AsagiLiteGrid::Error::Success);
    CHECK(grid.arrayView()->knots[1] == nullptr);
    CHECK(grid.arrayView()->deltaInv[1] == 1.0 / 500.0);
  }
  {
    // declarations of another number of axes than the file has, and knots that decrease
    df::AsagiLiteGrid grid;
    CHECK(open(grid, "rho", {AxisKind::Knots}) == df::AsagiLiteGrid::Error::BadKnots);
    CHECK(open(grid, "up", {AxisKind::Uniform, AxisKind::Knots}) ==
          df::AsagiLiteGrid::Error::BadKnots);
    CHECK(open(grid, "up", {AxisKind::Knots, AxisKind::Uniform}) ==
          df::AsagiLiteGrid::Error::Success);
  }
  std::remove(path.c_str());
}

} // namespace seissol::unit_test
