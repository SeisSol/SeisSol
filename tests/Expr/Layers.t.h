// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Core.h"
#include "Expr/Interp.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/Record.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/InlineGrid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/EasiReader.h"

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <netcdf.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace seissol::expr::test {

namespace {

/// Three layers of rho and mu at the depths 1000, 2000 and 3000, as easi's layered model has them;
/// one grid per interpolation.
const char* const LayeredModel =
    "#!sderiv 1\nin z\n"
    "grid m (kind = \"inline\", components = (\"rho\", \"mu\"), knots = ((1000, 2000, 3000)), "
    "values = (2600, 10, 2700, 20, 3000, 40), interpolation = \"linear\", boundary = \"clamp\")\n"
    "grid p (kind = \"inline\", components = (\"rho\", \"mu\"), knots = ((1000, 2000, 3000)), "
    "values = (2600, 10, 2700, 20, 3000, 40), interpolation = \"previous\", boundary = "
    "\"clamp\")\n"
    "grid n (kind = \"inline\", components = (\"rho\", \"mu\"), knots = ((1000, 2000, 3000)), "
    "values = (2600, 10, 2700, 20, 3000, 40), interpolation = \"next\", boundary = \"clamp\")\n"
    "grid q (kind = \"inline\", components = (\"rho\", \"mu\"), knots = ((1000, 2000, 3000)), "
    "values = (2600, 10, 2700, 20, 3000, 40), interpolation = \"nearest\", boundary = "
    "\"clamp\")\n"
    "def t0 = m.rho(z)\ndef t1 = p.rho(z)\ndef t2 = n.rho(z)\ndef t3 = q.rho(z)\ndef t4 = m.mu(z)\n"
    "out linear = t0\nout lower = t1\nout upper = t2\nout nearest = t3\nout mu = t4\n";

/// The outputs of `program` at the points of `columns`, one per input of the program in its order,
/// [output][point], on `backend`; the backend that ran in `used`. `grid`, where given, is the grid
/// of the program's first grid.
std::vector<double> evaluateAt(const Program& program,
                               const std::vector<std::vector<double>>& columns,
                               BackendKind backend = BackendKind::Interpreter,
                               BackendKind* used = nullptr,
                               std::unique_ptr<reader::datafield::Grid> grid = nullptr) {
  REQUIRE(columns.size() == program.inputs().size());
  const std::size_t count = columns.front().size();
  reader::scripting::DataTable table(count);
  for (std::size_t i = 0; i < columns.size(); ++i) {
    table.bindViewConst<double>(
        program.inputs()[i].name, reader::scripting::Direction::In, columns[i].data());
  }
  std::vector<double> out(program.outputs().size() * count, -1.0);
  for (std::size_t o = 0; o < program.outputs().size(); ++o) {
    table.bindView<double>(
        program.outputs()[o].name, reader::scripting::Direction::Out, out.data() + o * count);
  }
  auto binding = Binding::bind(program, table);
  reader::datafield::GridStore store;
  if (grid != nullptr) {
    store.injectForTesting(store.intern(program.grids()[0]), std::move(grid));
  }
  BackendOptions options;
  options.quiet = true;
  options.preferred = backend;
  const auto kernel = makeKernel(program, binding, store, options);
  kernel->precompute(table);
  kernel->run(table);
  if (used != nullptr) {
    *used = kernel->kind();
  }
  return out;
}

/// The outputs of `program`, of the one input z, at the points of `z`.
std::vector<double> evaluate(const Program& program,
                             const std::vector<double>& z,
                             BackendKind backend = BackendKind::Interpreter,
                             BackendKind* used = nullptr) {
  return evaluateAt(program, {z}, backend, used);
}

bool bitwiseEqual(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

} // namespace

TEST_SUITE("expr layers") {
  TEST_CASE("layers interpolate linearly, and step as easi's lower and upper") {
    const Program program = loadCore(LayeredModel, "");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> z = {500, 1000, 1500, 2000, 2500, 3000, 3500, nan, 1499, 1501};
    const auto out = evaluate(program, z);
    const std::size_t n = z.size();
    const auto at = [&](std::size_t output, std::size_t point) { return out[output * n + point]; };
    // linear: the value of a knot on it, constant beyond the ends, NaN at the first knot
    const std::vector<double> linear = {
        2600, 2600, 2650, 2700, 2850, 3000, 3000, 2600, 2649.9, 2650.1};
    // previous (easi's lower): a point on a knot belongs to the layer above it
    const std::vector<double> previous = {
        2600, 2600, 2600, 2600, 2700, 2700, 3000, 2600, 2600, 2600};
    // next (easi's upper)
    const std::vector<double> next = {2600, 2600, 2700, 2700, 3000, 3000, 3000, 2600, 2700, 2700};
    // nearest: half up
    const std::vector<double> nearest = {
        2600, 2600, 2700, 2700, 3000, 3000, 3000, 2600, 2600, 2700};
    for (std::size_t k = 0; k < n; ++k) {
      CAPTURE(z[k]);
      CHECK(at(0, k) == doctest::Approx(linear[k]).epsilon(1e-15));
      CHECK(at(1, k) == previous[k]);
      CHECK(at(2, k) == next[k]);
      CHECK(at(3, k) == nearest[k]);
    }
    // on a knot exactly the value of the knot
    CHECK(at(0, 1) == 2600.0);
    CHECK(at(0, 3) == 2700.0);
    CHECK(at(0, 5) == 3000.0);
    CHECK(at(4, 3) == 20.0);
  }

  TEST_CASE("the compiled kernel samples layers bit for bit as the interpreter does") {
    const Program program = loadCore(LayeredModel, "");
    std::vector<double> z;
    for (int i = -50; i < 450; ++i) {
      z.push_back(10.0 * i);
    }
    z.push_back(std::numeric_limits<double>::quiet_NaN());
    z.push_back(std::numeric_limits<double>::infinity());
    z.push_back(-std::numeric_limits<double>::infinity());
    z.push_back(std::nextafter(2000.0, 0.0));
    z.push_back(std::nextafter(2000.0, 3000.0));
    const auto interpreted = evaluate(program, z);
    BackendKind used = BackendKind::Interpreter;
    const auto compiled = evaluate(program, z, BackendKind::RtcCpu, &used);
    if (used != BackendKind::RtcCpu) {
      WARN_MESSAGE(false, "no usable C++ compiler; the compiled CPU kernel was not executed");
      return;
    }
    CHECK(bitwiseEqual(interpreted, compiled));
  }

  TEST_CASE("a grid of two axes of knots, each with an interpolation of its own") {
    const Program program = loadCore(
        "#!sderiv 1\nin x, z\n"
        "grid g (kind = \"inline\", components = (\"v\"), knots = ((0, 1, 3), (10, 20)), values = "
        "(0, 1, 3, 100, 101, 103), interpolation = (\"linear\", \"previous\"), boundary = "
        "\"clamp\")\n"
        "def t0 = g.v(x, z)\nout v = t0\n",
        "");
    const auto valueAt = [&](double x, double z) {
      reader::scripting::DataTable table(1);
      table.bindViewConst<double>("x", reader::scripting::Direction::In, &x);
      table.bindViewConst<double>("z", reader::scripting::Direction::In, &z);
      double out = -1.0;
      table.bindView<double>("v", reader::scripting::Direction::Out, &out);
      auto binding = Binding::bind(program, table);
      reader::datafield::GridStore store;
      BackendOptions options;
      options.quiet = true;
      const auto kernel = makeKernel(program, binding, store, options);
      kernel->precompute(table);
      kernel->run(table);
      return out;
    };
    // linear along x between the knots 1 and 3; z steps: up to 20 the layer of 10
    CHECK(valueAt(2.0, 15.0) == 2.0);
    CHECK(valueAt(2.0, 20.0) == 2.0);
    CHECK(valueAt(2.0, 25.0) == 102.0);
    CHECK(valueAt(0.5, 25.0) == 100.5);
    CHECK(valueAt(-1.0, 5.0) == 0.0);
  }

  TEST_CASE("an inline grid reads back from the normal form, of the same fingerprint") {
    const Program program = loadCore(LayeredModel, "");
    const std::string printed = printCore(program);
    CAPTURE(printed);
    CHECK(printed.find("kind = \"inline\", components = (\"rho\", \"mu\"), knots = ((1000, 2000, "
                       "3000)), values = (2600, 10, 2700, 20, 3000, 40), interpolation = "
                       "\"previous\"") != std::string::npos);
    const Program again = loadCore(printed, "");
    CHECK(again.fingerprint() == program.fingerprint());
    CHECK(printCore(again) == printed);
    // other values are another program
    std::string other = LayeredModel;
    other.replace(other.find("2700, 20"), 8, "2701, 20");
    CHECK(loadCore(other, "").fingerprint() != program.fingerprint());
  }

  TEST_CASE("a grid that does not form one is an error where it is declared") {
    const auto failure = [](const std::string& grid) {
      try {
        static_cast<void>(
            loadCore("#!sderiv 1\nin z\n" + grid + "\ndef t0 = m.a(z)\nout a = t0\n", ""));
      } catch (const CoreError& error) {
        return std::string(error.what());
      }
      return std::string();
    };
    CHECK(failure("grid m (kind = \"inline\", components = (\"a\"), knots = ((2, 1)), values = "
                  "(0, 1), interpolation = \"linear\", boundary = \"clamp\")")
              .find("do not increase strictly") != std::string::npos);
    CHECK(failure("grid m (kind = \"inline\", components = (\"a\"), knots = ((1, 2)), values = "
                  "(0, 1, 2), interpolation = \"linear\", boundary = \"clamp\")")
              .find("3 values for 2 nodes") != std::string::npos);
    CHECK(failure("grid m (kind = \"inline\", components = (\"a\"), knots = ((1, 2)), values = "
                  "(0, 1), interpolation = \"cubic\", boundary = \"clamp\")")
              .find("not cubic") != std::string::npos);
    CHECK(failure("grid m (kind = \"inline\", components = (\"a\"), knots = ((1, 2)), values = "
                  "(0, 1), interpolation = (\"linear\", \"next\"), boundary = \"clamp\")")
              .find("2 interpolations") != std::string::npos);
    CHECK(failure("grid m (kind = \"inline\", path = \"a.nc\", components = (\"a\"), knots = "
                  "((1, 2)), values = (0, 1), interpolation = \"linear\", boundary = \"clamp\")")
              .find("no path or variable") != std::string::npos);
  }

  TEST_CASE("layers step at the knots as easi's layered model does") {
    // the same layers through easi, lower and upper, at the knots and between them
    const std::string path = "seissol-test-layers-" + std::to_string(::getpid()) + ".yaml";
    const std::vector<double> z = {500, 1000, 1500, 2000, 2500, 3000, 3500, 1999.75, 2000.25};
    const Program program = loadCore(LayeredModel, "");
    const auto ours = evaluate(program, z);
    for (const auto& [interpolation, column] : std::vector<std::pair<std::string, std::size_t>>{
             {"lower", 1}, {"upper", 2}, {"linear", 0}}) {
      CAPTURE(interpolation);
      {
        std::ofstream file(path);
        file << "!LayeredModel\nmap: !AffineMap\n  matrix:\n    depth: [0.0, 0.0, 1.0]\n"
                "  translation:\n    depth: 0\ninterpolation: "
             << interpolation
             << "\nparameters: [rho, mu]\nnodes:\n  1000: [2600, 10]\n  2000: [2700, 20]\n"
                "  3000: [3000, 40]\n";
      }
      reader::scripting::EasiReader reader(path, {"x", "y", "z"});
      std::vector<double> zeros(z.size(), 0.0);
      std::vector<double> rho(z.size(), -1.0);
      std::vector<double> mu(z.size(), -1.0);
      reader::scripting::DataTable table(z.size());
      table.bindViewConst<double>("x", reader::scripting::Direction::In, zeros.data());
      table.bindViewConst<double>("y", reader::scripting::Direction::In, zeros.data());
      table.bindViewConst<double>("z", reader::scripting::Direction::In, z.data());
      table.bindView<double>("rho", reader::scripting::Direction::Out, rho.data());
      table.bindView<double>("mu", reader::scripting::Direction::Out, mu.data());
      reader.prepare(table);
      reader.call(table);
      for (std::size_t k = 0; k < z.size(); ++k) {
        CAPTURE(z[k]);
        const double value = ours[column * z.size() + k];
        if (interpolation == "linear") {
          // easi as slope and intercept per layer: a few ulp apart
          CHECK(value == doctest::Approx(rho[k]).epsilon(1e-14));
        } else {
          CHECK(value == rho[k]);
        }
      }
    }
    std::remove(path.c_str());
  }

  TEST_CASE("a uniform axis steps as an axis of knots does") {
    // the ramp i + 10 j + 100 k over 8 nodes per axis, from 0 in steps of 1: previous along x,
    // next along y, linear along z
    const Program program =
        loadCore("#!sderiv 1\nin x, y, z\n"
                 "grid m (kind = \"asagi\", path = \"ramp.nc\", variable = \"data\", components = "
                 "(\"a\", \"b\"), "
                 "interpolation = (\"previous\", \"next\", \"linear\"), boundary = \"clamp\")\n"
                 "def t0 = m.a(x, y, z)\nout a = t0\n",
                 "");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> x = {1.5, 2.0, 0.0, -3.0, 7.0, 6.5, nan, 9.0};
    const std::vector<double> y = {1.0, 1.5, 0.0, 9.0, 7.0, 0.5, nan, 0.0};
    const std::vector<double> z = {0.25, 3.0, 0.0, 7.5, 1.0, 2.0, 0.0, 0.0};
    // x: the last node below it; y: the first at it or above it
    const std::vector<double> expected = {36, 321, 0, 770, 176, 216, 0, 7};
    const auto out = evaluateAt(
        program, {x, y, z}, BackendKind::Interpreter, nullptr, std::make_unique<LinearRampGrid>());
    for (std::size_t k = 0; k < x.size(); ++k) {
      CAPTURE(k);
      CHECK(out[k] == expected[k]);
    }
    BackendKind used = BackendKind::Interpreter;
    const auto compiled = evaluateAt(
        program, {x, y, z}, BackendKind::RtcCpu, &used, std::make_unique<LinearRampGrid>());
    if (used == BackendKind::RtcCpu) {
      CHECK(bitwiseEqual(out, compiled));
    }
  }

  TEST_CASE("a file whose axis is declared knots samples it as an inline grid of its data") {
    // rho(depth, x) over the depths 0, 100 and 1000 and x = 0 and 1: the depth axis is not
    // uniform; declared knots, the file samples bit for bit as the same data inline does
    const std::string path = "seissol-test-layers-" + std::to_string(::getpid()) + ".nc";
    const std::vector<double> depths = {0.0, 100.0, 1000.0};
    const std::vector<double> xs = {0.0, 1.0};
    // depth-major in the file: rho(depth, x)
    const std::vector<double> rho = {2000.0, 2100.0, 2500.0, 2650.0, 2900.0, 3300.0};
    {
      int file = 0;
      REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
      int dims[2] = {};
      REQUIRE(nc_def_dim(file, "depth", depths.size(), &dims[0]) == NC_NOERR);
      REQUIRE(nc_def_dim(file, "x", xs.size(), &dims[1]) == NC_NOERR);
      int depth = 0;
      int x = 0;
      int values = 0;
      REQUIRE(nc_def_var(file, "depth", NC_DOUBLE, 1, &dims[0], &depth) == NC_NOERR);
      REQUIRE(nc_def_var(file, "x", NC_DOUBLE, 1, &dims[1], &x) == NC_NOERR);
      REQUIRE(nc_def_var(file, "rho", NC_DOUBLE, 2, dims, &values) == NC_NOERR);
      REQUIRE(nc_enddef(file) == NC_NOERR);
      REQUIRE(nc_put_var_double(file, depth, depths.data()) == NC_NOERR);
      REQUIRE(nc_put_var_double(file, x, xs.data()) == NC_NOERR);
      REQUIRE(nc_put_var_double(file, values, rho.data()) == NC_NOERR);
      REQUIRE(nc_close(file) == NC_NOERR);
    }
    const std::string schemes =
        "interpolation = (\"linear\", \"previous\"), boundary = \"clamp\")\n";
    const Program fromFile =
        loadCore("#!sderiv 1\nin x, z\ngrid m (kind = \"asagi\", path = \"" + path +
                     "\", variable = \"rho\", components = (\"rho\"), axes = "
                     "(\"uniform\", \"knots\"), " +
                     schemes + "def t0 = m.rho(x, z)\nout rho = t0\n",
                 "");
    // the grid axes in the opposite order of the dimensions of the file: x, then depth
    const Program inlined = loadCore(
        "#!sderiv 1\nin x, z\ngrid m (kind = \"inline\", components = (\"rho\"), knots = ((0, 1), "
        "(0, 100, 1000)), values = (2000, 2100, 2500, 2650, 2900, 3300), " +
            schemes + "def t0 = m.rho(x, z)\nout rho = t0\n",
        "");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> x = {0.5, 0.25, 0.0, 1.0, 0.75, -1.0, 2.0, nan, 0.5};
    const std::vector<double> z = {50.0, 100.0, 550.0, 1000.0, 2000.0, -5.0, 100.5, 0.0, nan};
    const auto expected = evaluateAt(inlined, {x, z});
    // previous along the depth: the layer above, a point on a knot in the layer above it
    CHECK(expected[0] == 2050.0);
    CHECK(expected[1] == 2025.0);
    CHECK(expected[2] == 2500.0);
    CHECK(expected[3] == 2650.0);
    const auto interpreted = evaluateAt(fromFile, {x, z});
    CHECK(bitwiseEqual(interpreted, expected));
    BackendKind used = BackendKind::Interpreter;
    const auto compiled = evaluateAt(fromFile, {x, z}, BackendKind::RtcCpu, &used);
    if (used == BackendKind::RtcCpu) {
      CHECK(bitwiseEqual(compiled, expected));
    }
    std::remove(path.c_str());
  }

  TEST_CASE("a record names an inline grid, and its axes of knots, and replays it") {
    const Program program = loadCore(
        "#!sderiv 1\nin z\n"
        "grid m (kind = \"inline\", components = (\"rho\", \"mu\"), knots = ((1000, 2000, "
        "3000)), values = (2600, 10, 2700, 20, 3000, 40), interpolation = \"linear\", boundary = "
        "\"clamp\")\n"
        "def t0 = m.rho(z)\nout rho = t0\n",
        "");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> z = {500, 1000, 1500, 2000, 2500, 3000, 3500, nan};
    std::vector<double> rho(z.size(), 0.0);
    reader::scripting::DataTable table(z.size());
    table.bindViewConst<double>("z", reader::scripting::Direction::In, z.data());
    table.bindView<double>("rho", reader::scripting::Direction::Out, rho.data());
    auto binding = Binding::bind(program, table);
    reader::datafield::GridStore store;
    BackendOptions options;
    options.quiet = true;
    const auto kernel = makeKernel(program, binding, store, options);
    kernel->precompute(table);
    kernel->run(table);
    const auto comparison = compareWithInterpreter(program, binding, store, &table, nullptr);
    REQUIRE(comparison.has_value());
    const std::string text = printRecord(recordOf(program, "interpreter", *comparison, {}, 0));
    CAPTURE(text);
    CHECK(text.find("grid m = inline, components \"rho\" \"mu\"\n") != std::string::npos);
    CHECK(text.find("grid m axis 0 = knots, min 1000, max 3000, nodes 3\n") != std::string::npos);
    // below, above, NaN
    CHECK(text.find("clamped m = 1 1 1 of 8 points\n") != std::string::npos);
    const Record read = parseRecord(text);
    REQUIRE(read.grids.size() == 1);
    CHECK(read.grids[0] == comparison->grids.at(0));
    CHECK(read.grids[0].inlined);
    CHECK(read.grids[0].axes.at(0).knots);
    CHECK(printRecord(read) == text);
    // the grid is in the program of the record
    reader::datafield::GridStore again;
    const Replay replay = replayRecord(read, again);
    CHECK(replay.points > 0);
    CHECK(replay.differing == 0);
  }

  TEST_CASE("an axis of knots that fails reports the points beyond its knots") {
    // two axes of knots: x clamps, z fails; a point beyond the knots of z, or NaN there, is
    // reported, one beyond those of x is not
    const Program program = loadCore(
        "#!sderiv 1\nin x, z\n"
        "grid g (kind = \"inline\", components = (\"v\"), knots = ((0, 1, 3), (10, 20)), values = "
        "(0, 1, 3, 100, 101, 103), interpolation = (\"linear\", \"next\"), boundary = "
        "(\"clamp\", \"fail\"))\n"
        "def t0 = g.v(x, z)\nout v = t0\n",
        "");
    REQUIRE(program.grids().size() == 1);
    const auto& desc = program.grids()[0];
    REQUIRE(desc.data != nullptr);

    class Reporter final : public GridArrays {
  public:
      explicit Reporter(const reader::datafield::GridDesc& desc) : grid_(desc.data, 1) {}
      [[nodiscard]] const reader::datafield::ArrayView& array(GridId /*grid*/) const override {
        return *grid_.arrayView();
      }
      void report(GridId /*grid*/, CheckKind kind) const override { reports.push_back(kind); }
      mutable std::vector<CheckKind> reports;

  private:
      reader::datafield::InlineGrid grid_;
    };

    const LoweredProgram lowered = lower(program);
    std::vector<CheckKind> reports;
    const auto sample = [&](double x, double z) {
      const Reporter grids(desc);
      const std::vector<double> inputs = {x, z};
      std::vector<double> outputs(1, 0.0);
      ArrayTileIo<double> io(inputs, outputs, 1, 2, 1);
      TileInterpreter<double> interp(program, lowered, &grids, InterpreterOptions{});
      std::vector<double> persistent(
          std::max<std::size_t>(static_cast<std::size_t>(lowered.persistentSlotCount()), 1), 0.0);
      interp.run(io, 1, persistent.data());
      reports = grids.reports;
      return outputs[0];
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // within the knots, on them, and beyond those of x
    CHECK(sample(2.0, 15.0) == 102.0);
    CHECK(reports.empty());
    CHECK(sample(2.0, 10.0) == 2.0);
    CHECK(reports.empty());
    CHECK(sample(2.0, 20.0) == 102.0);
    CHECK(reports.empty());
    CHECK(sample(-5.0, 15.0) == 100.0);
    CHECK(reports.empty());
    CHECK(sample(nan, 15.0) == 100.0);
    CHECK(reports.empty());
    // beyond those of z, or NaN there
    for (const double z : {9.999, 20.001, nan, -std::numeric_limits<double>::infinity()}) {
      CAPTURE(z);
      sample(2.0, z);
      const bool outside = reports == std::vector<CheckKind>{CheckKind::Outside};
      CHECK(outside);
    }
  }
}

} // namespace seissol::expr::test
