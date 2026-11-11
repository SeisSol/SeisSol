// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Reader/Datafield/AsagiLite.h"
#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <netcdf.h>
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
  // The update module advances the window at most (W - w) cells apart
  // (GridStore::suggestedSyncInterval()); every time up to the next update samples at full
  // accuracy, and never a slice clamped onto the end of the window.
  constexpr std::size_t Slices = 40;
  constexpr std::size_t Resident = 4;
  const LinearInTime file(Slices);
  df::AsagiLiteGrid grid;
  REQUIRE(grid.open(file.path(), "z", df::Interpolation::Linear, df::BoundaryMode::Clamp, 1) ==
          df::AsagiLiteGrid::Error::Success);
  REQUIRE(grid.resizeWindow(Resident) == df::AsagiLiteGrid::Error::Success);
  REQUIRE(grid.arrayView() != nullptr);

  const double interval = static_cast<double>(Resident - 2);
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

} // namespace seissol::unit_test
