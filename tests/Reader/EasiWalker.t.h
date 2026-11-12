// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Interp.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Datafield/AsagiModule.h"
#include "Reader/Datafield/AsagiReader.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/ScecGrid.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/EasiReader.h"
#include "Reader/Scripting/EasiWalker.h"
#include "Reader/Scripting/ImpalaToLua.h"
#include "Reader/Scripting/ReaderBuilder.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <easi/Component.h>
#include <easi/YAMLParser.h>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <netcdf.h>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <utils/env.h>
#include <vector>

#if __has_include(<easi/util/ImpalaHandler.h>)
#include <easi/util/ImpalaHandler.h>
#define SEISSOL_TEST_EASI_IMPALA_TO_LUA
#endif

namespace seissol::unit_test {

namespace easiwalkertest {

namespace sc = reader::scripting;

/// Whether easi samples a grid as the walker does: easi 1.7 and later (InterpolationKernel.h, which
/// Reader/Datafield/Interpolation.h ports), with cubic interpolation. An earlier easi interpolates
/// linearly within a few ulp of it.
#if defined(EASI_VERSION_MAJOR) && defined(EASI_VERSION_MINOR) &&                                  \
    (EASI_VERSION_MAJOR > 1 || (EASI_VERSION_MAJOR == 1 && EASI_VERSION_MINOR >= 7))
constexpr bool EasiKernel = true;
#else
constexpr bool EasiKernel = false;
#endif
/// The ulp the linear interpolation of a grid is from easi's.
constexpr double GridUlps = EasiKernel ? 0.0 : 4.0;

/// A directory of easi files for one test, removed at the end of the scope.
class Files {
  public:
  Files() {
    directory_ =
        std::filesystem::temp_directory_path() /
        ("seissol-easiwalker-" + std::to_string(::getpid()) + "-" + std::to_string(counter()++));
    std::filesystem::create_directories(directory_);
  }
  ~Files() {
    std::error_code error;
    std::filesystem::remove_all(directory_, error);
  }
  Files(const Files&) = delete;
  Files& operator=(const Files&) = delete;
  Files(Files&&) = delete;
  Files& operator=(Files&&) = delete;

  std::string write(const std::string& name, const std::string& text) const {
    const auto path = directory_ / name;
    std::ofstream(path) << text;
    return path.string();
  }
  [[nodiscard]] std::string path(const std::string& name) const {
    return (directory_ / name).string();
  }

  private:
  static int& counter() {
    static int count = 0;
    return count;
  }
  std::filesystem::path directory_;
};

/// The points: x, y and z over [-lo, hi] each, the corners and the planes of the filters among
/// them, and the groups 1 to 3.
struct Points {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> group;
  [[nodiscard]] std::size_t size() const { return x.size(); }
};

inline Points pointsIn(double lo, double hi, const std::vector<double>& planes, std::size_t count) {
  Points points;
  std::uint64_t state = 0x9e3779b97f4a7c15ULL;
  const auto next = [&state] {
    state ^= state << 13U;
    state ^= state >> 7U;
    state ^= state << 17U;
    return static_cast<double>(state >> 11U) * 0x1.0p-53;
  };
  for (std::size_t k = 0; k < count; ++k) {
    std::array<double, 3> c{};
    for (auto& value : c) {
      value = lo + (hi - lo) * next();
      // every fifth on a plane of a filter or a node, where a point belongs to the one before
      if (!planes.empty() && k % 5 == 0) {
        value = planes[static_cast<std::size_t>(next() * static_cast<double>(planes.size())) %
                       planes.size()];
      }
    }
    points.x.push_back(c[0]);
    points.y.push_back(c[1]);
    points.z.push_back(c[2]);
    points.group.push_back(static_cast<double>(1 + k % 3));
  }
  return points;
}

/// The values of `parameters` at `points`, [parameter][point], by `reader`.
inline std::vector<std::vector<double>> valuesOf(sc::DataReader& reader,
                                                 const Points& points,
                                                 const std::vector<std::string>& parameters) {
  std::vector<std::vector<double>> values(parameters.size(),
                                          std::vector<double>(points.size(), -1.0));
  sc::DataTable table(points.size());
  table.bindViewConst<double>("x", sc::Direction::In, points.x.data());
  table.bindViewConst<double>("y", sc::Direction::In, points.y.data());
  table.bindViewConst<double>("z", sc::Direction::In, points.z.data());
  table.bindViewConst<double>("group", sc::Direction::In, points.group.data());
  for (std::size_t p = 0; p < parameters.size(); ++p) {
    table.bindView<double>(parameters[p], sc::Direction::Out, values[p].data());
  }
  reader.prepare(table);
  reader.call(table);
  return values;
}

inline bool sameBits(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

/// Units in the last place between `a` and `b`.
inline double ulpsBetween(double a, double b) {
  if (sameBits(a, b)) {
    return 0.0;
  }
  const double ulp = std::max(std::fabs(std::nextafter(a, INFINITY) - a),
                              std::fabs(std::nextafter(b, INFINITY) - b));
  return std::fabs(a - b) / ulp;
}

/// That the walker gives easi's values of `parameters` at `points`: the same bits, or at most
/// `ulps` apart.
inline void sameAsEasi(const std::string& path,
                       const std::vector<std::string>& parameters,
                       const Points& points,
                       double ulps = 0.0) {
  sc::EasiReader easi(path, {"x", "y", "z"});
  sc::EasiProgramReader walker(path, {"x", "y", "z"}, nullptr);
  const auto expected = valuesOf(easi, points, parameters);
  const auto actual = valuesOf(walker, points, parameters);
  for (std::size_t p = 0; p < parameters.size(); ++p) {
    CAPTURE(parameters[p]);
    double worst = 0.0;
    std::size_t at = 0;
    for (std::size_t k = 0; k < points.size(); ++k) {
      const double off = ulpsBetween(actual[p][k], expected[p][k]);
      if (!(off <= worst)) {
        worst = off;
        at = k;
      }
    }
    INFO("at (" << points.x[at] << ", " << points.y[at] << ", " << points.z[at] << ") in group "
                << points.group[at] << ": " << actual[p][at] << " for easi's " << expected[p][at]);
    CHECK(worst <= ulps);
  }
}

/// What easi's ASAGI reads its settings from, in a process without SeisSol's modules.
inline void asagiModule() {
  static utils::Env env("SEISSOL_");
  static const bool Made = [] {
    asagi::AsagiModule::initInstance(env);
    return true;
  }();
  static_cast<void>(Made);
}

/// A NetCDF file of an ASAGI grid over x and y (0 to 0.8 by 0.1 in binary32, 0 to 3 by 1) of the
/// components
/// rho and mu, as ASAGI has them: a compound of their values per node; values that no
/// interpolation reproduces exactly.
inline void writeAsagi(const std::string& path) {
  struct Material {
    float rho;
    float mu;
  };
  int file = 0;
  REQUIRE(nc_create(path.c_str(), NC_NETCDF4 | NC_CLOBBER, &file) == NC_NOERR);
  constexpr std::size_t Nx = 9;
  constexpr std::size_t Ny = 4;
  int material = 0;
  REQUIRE(nc_def_compound(file, sizeof(Material), "material", &material) == NC_NOERR);
  REQUIRE(nc_insert_compound(file, material, "rho", offsetof(Material, rho), NC_FLOAT) == NC_NOERR);
  REQUIRE(nc_insert_compound(file, material, "mu", offsetof(Material, mu), NC_FLOAT) == NC_NOERR);
  int dims[2] = {};
  REQUIRE(nc_def_dim(file, "y", Ny, &dims[0]) == NC_NOERR);
  REQUIRE(nc_def_dim(file, "x", Nx, &dims[1]) == NC_NOERR);
  int y = 0;
  int x = 0;
  int data = 0;
  REQUIRE(nc_def_var(file, "y", NC_FLOAT, 1, &dims[0], &y) == NC_NOERR);
  REQUIRE(nc_def_var(file, "x", NC_FLOAT, 1, &dims[1], &x) == NC_NOERR);
  REQUIRE(nc_def_var(file, "data", material, 2, dims, &data) == NC_NOERR);
  REQUIRE(nc_enddef(file) == NC_NOERR);
  std::vector<float> ys(Ny);
  std::vector<float> xs(Nx);
  for (std::size_t j = 0; j < Ny; ++j) {
    ys[j] = static_cast<float>(j);
  }
  // not exact in binary: the spacing is the mean between the ends, as ASAGI derives it
  for (std::size_t i = 0; i < Nx; ++i) {
    xs[i] = 0.1F * static_cast<float>(i);
  }
  std::vector<Material> values(Ny * Nx);
  for (std::size_t j = 0; j < Ny; ++j) {
    for (std::size_t i = 0; i < Nx; ++i) {
      const auto a = static_cast<double>(i);
      const auto b = static_cast<double>(j);
      values[j * Nx + i] = Material{static_cast<float>(2600.0 + 37.3 * std::sin(0.7 * a + 1.3 * b)),
                                    static_cast<float>(3.1e10 + 1.7e9 * std::cos(0.4 * a * b))};
    }
  }
  REQUIRE(nc_put_var_float(file, y, ys.data()) == NC_NOERR);
  REQUIRE(nc_put_var_float(file, x, xs.data()) == NC_NOERR);
  REQUIRE(nc_put_var(file, data, values.data()) == NC_NOERR);
  REQUIRE(nc_close(file) == NC_NOERR);
}

/// Evaluates `program` with the interpreter at `points`, and gives the points that fail a check,
/// by the check.
inline std::vector<std::pair<std::int32_t, std::size_t>> failuresOf(const expr::Program& program,
                                                                    const Points& points) {
  class Io final : public expr::TileIo<double> {
public:
    Io(const expr::Program& program, const Points& points) : program_(program), points_(points) {}
    void gather(std::size_t first, std::size_t count, double* dst) const override {
      for (std::size_t i = 0; i < program_.inputs().size(); ++i) {
        const auto& name = program_.inputs()[i].name;
        const auto& column = name == "x"   ? points_.x
                             : name == "y" ? points_.y
                             : name == "z" ? points_.z
                                           : points_.group;
        for (std::size_t l = 0; l < count; ++l) {
          dst[i * count + l] = column[first + l];
        }
      }
    }
    void scatter(std::size_t /*first*/, std::size_t /*count*/, const double* /*src*/) override {}

private:
    const expr::Program& program_;
    const Points& points_;
  };
  class Reporter final : public expr::GridArrays {
public:
    [[nodiscard]] const reader::datafield::ArrayView& array(expr::GridId /*grid*/) const override {
      throw std::logic_error("no grid here");
    }
    void report(expr::GridId /*grid*/, expr::CheckKind /*kind*/) const override {}
    void failed(std::int32_t check,
                std::size_t point,
                const std::vector<double>& /*inputs*/) const override {
      failures.emplace_back(check, point);
    }
    mutable std::vector<std::pair<std::int32_t, std::size_t>> failures;
  };
  const expr::LoweredProgram lowered = expr::lower(program);
  Io io(program, points);
  const Reporter reporter;
  expr::InterpreterOptions options;
  options.tileSize = 1;
  expr::TileInterpreter<double> interpreter(program, lowered, &reporter, options);
  std::vector<double> persistent(
      std::max<std::size_t>(static_cast<std::size_t>(lowered.persistentSlotCount()), 1), 0.0);
  interpreter.run(io, points.size(), persistent.data());
  return reporter.failures;
}

/// Whether the easi of this build reads the file `path` (an easi of another version may not: one
/// without cubic interpolation, say); EasiReader stops the run where it does not.
inline bool easiReads(const std::string& path) {
  asagi::AsagiReader asagi;
  easi::YAMLParser parser(std::set<std::string>{"x", "y", "z"}, &asagi);
  try {
    const std::unique_ptr<easi::Component> root(parser.parse(path));
    return root != nullptr;
  } catch (const std::exception& /*error*/) {
    return false;
  }
}

/// The points among `points` that meet every check of `program` (the points easi evaluates), with
/// the grids of `store`.
inline std::vector<std::size_t> acceptedBy(const expr::Program& program,
                                           reader::datafield::GridStore& store,
                                           const Points& points) {
  class Io final : public expr::TileIo<double> {
public:
    Io(const expr::Program& program, const Points& points) : program_(program), points_(points) {}
    void gather(std::size_t first, std::size_t count, double* dst) const override {
      for (std::size_t i = 0; i < program_.inputs().size(); ++i) {
        const auto& name = program_.inputs()[i].name;
        const auto& column = name == "x"   ? points_.x
                             : name == "y" ? points_.y
                             : name == "z" ? points_.z
                                           : points_.group;
        for (std::size_t l = 0; l < count; ++l) {
          dst[i * count + l] = column[first + l];
        }
      }
    }
    void scatter(std::size_t /*first*/, std::size_t /*count*/, const double* /*src*/) override {}

private:
    const expr::Program& program_;
    const Points& points_;
  };
  class Arrays final : public expr::GridArrays {
public:
    Arrays(reader::datafield::GridStore& store, std::vector<std::size_t> ids)
        : store_(store), ids_(std::move(ids)) {}
    [[nodiscard]] const reader::datafield::ArrayView& array(expr::GridId grid) const override {
      return *store_.get(ids_.at(static_cast<std::size_t>(grid))).arrayView();
    }
    void report(expr::GridId /*grid*/, expr::CheckKind /*kind*/) const override {}
    void failed(std::int32_t /*check*/,
                std::size_t point,
                const std::vector<double>& /*inputs*/) const override {
      failing.push_back(point);
    }
    mutable std::vector<std::size_t> failing;

private:
    reader::datafield::GridStore& store_;
    std::vector<std::size_t> ids_;
  };
  const expr::LoweredProgram lowered = expr::lower(program);
  Io io(program, points);
  const Arrays arrays(store, expr::internGrids(program, store));
  expr::InterpreterOptions options;
  options.tileSize = 1;
  expr::TileInterpreter<double> interpreter(program, lowered, &arrays, options);
  std::vector<double> persistent(
      std::max<std::size_t>(static_cast<std::size_t>(lowered.persistentSlotCount()), 1), 0.0);
  interpreter.run(io, points.size(), persistent.data());
  std::vector<std::size_t> accepted;
  for (std::size_t k = 0; k < points.size(); ++k) {
    if (std::find(arrays.failing.begin(), arrays.failing.end(), k) == arrays.failing.end()) {
      accepted.push_back(k);
    }
  }
  return accepted;
}

} // namespace easiwalkertest

using namespace easiwalkertest;

TEST_CASE("easi walker: the files of SEISSOL_EASI_EXAMPLES agree with easi where it translates "
          "them" *
          doctest::test_suite("reader")) {
  // opt-in: a directory of easi files, such as SeisSol/Examples and easi's examples, whose data
  // files lie next to them. A file is translated, and compared with easi at points of every scale
  // that easi evaluates (those that meet every check of the program); a file the walker does not
  // translate is counted.
  const char* directory = std::getenv("SEISSOL_EASI_EXAMPLES");
  if (directory == nullptr) {
    return;
  }
  asagiModule();
  std::vector<std::string> paths;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
    if (entry.path().extension() == ".yaml") {
      paths.push_back(entry.path().string());
    }
  }
  std::sort(paths.begin(), paths.end());
  std::size_t translated = 0;
  std::map<std::string, std::size_t> refused;
  for (const auto& path : paths) {
    CAPTURE(path);
    std::vector<std::string> parameters;
    expr::Program program;
    try {
      parameters = sc::easiParameters(path, {"x", "y", "z"});
      sc::EasiQuery query;
      query.dimensions = {"x", "y", "z"};
      query.parameters = parameters;
      program = sc::easiProgram(path, query);
    } catch (const std::exception& error) {
      const std::string what = error.what();
      const auto tag = what.find(" translate ");
      ++refused[tag == std::string::npos
                    ? what.substr(what.find(": ") + 2, 60)
                    : what.substr(tag + 11, what.find(' ', tag + 11) - tag - 11)];
      continue;
    }
    if (parameters.empty()) {
      continue;
    }
    ++translated;
    // points of every scale, the groups 1 to 3, where easi evaluates them
    Points all;
    for (const auto& [lo, hi] : std::vector<std::pair<double, double>>{{-1.0, 1.0},
                                                                       {0.0, 4.0},
                                                                       {-4.0, 0.0},
                                                                       {-2.0, 4.0},
                                                                       {-20.0, 20.0},
                                                                       {-200.0, 200.0},
                                                                       {-2.0e4, 2.0e4},
                                                                       {-2.0e5, 2.0e5}}) {
      const auto some = pointsIn(lo, hi, {}, 300);
      all.x.insert(all.x.end(), some.x.begin(), some.x.end());
      all.y.insert(all.y.end(), some.y.begin(), some.y.end());
      all.z.insert(all.z.end(), some.z.begin(), some.z.end());
      all.group.insert(all.group.end(), some.group.begin(), some.group.end());
    }
    reader::datafield::GridStore store;
    // the parameters easi evaluates somewhere: of those easi says the file supplies, a !Switch
    // gives those its branches name only, where a branch supplies more (an !IdentityMap: all it
    // is given), and fails for the others everywhere
    std::vector<std::string> evaluated;
    for (const auto& parameter : parameters) {
      sc::EasiQuery query;
      query.dimensions = {"x", "y", "z"};
      query.parameters = {parameter};
      if (!acceptedBy(sc::easiProgram(path, query), store, all).empty()) {
        evaluated.push_back(parameter);
      }
    }
    if (evaluated.empty()) {
      MESSAGE(path << ": easi evaluates no parameter at these points");
      continue;
    }
    if (evaluated.size() != parameters.size()) {
      parameters = evaluated;
      sc::EasiQuery query;
      query.dimensions = {"x", "y", "z"};
      query.parameters = parameters;
      program = sc::easiProgram(path, query);
    }
    Points points;
    for (const auto k : acceptedBy(program, store, all)) {
      points.x.push_back(all.x[k]);
      points.y.push_back(all.y[k]);
      points.z.push_back(all.z[k]);
      points.group.push_back(all.group[k]);
    }
    MESSAGE(path << ": " << points.size() << " of " << all.size() << " points evaluated, of "
                 << parameters.size() << " parameters");
    if (!easiReads(path)) {
      MESSAGE(path << ": the easi of this build does not read it; not compared");
      continue;
    }
    if (points.size() > 0) {
      sameAsEasi(path, parameters, points, GridUlps);
    }
  }
  std::string summary;
  for (const auto& [why, count] : refused) {
    summary += "\n  " + std::to_string(count) + "x " + why;
  }
  MESSAGE(translated << " of " << paths.size() << " files translated; the others:" << summary);
}

TEST_CASE("easi walker: maps, filters and switches give easi's bits" *
          doctest::test_suite("reader")) {
  const Files files;
  // an affine map and a polynomial on it; a switch; the first filter of a list that accepts a
  // point, with the boundaries of the boxes among the points; a group filter; a sphere; and the
  // map that takes the rest
  const std::string path = files.write("model.yaml",
                                       "!Switch\n"
                                       "[rho, mu]: !AffineMap\n"
                                       "  matrix:\n"
                                       "    xi: [0.3, -1.25, 0.0]\n"
                                       "    zeta: [0.0, 0.1, 1.0]\n"
                                       "  translation:\n"
                                       "    xi: 0.7\n"
                                       "    zeta: -1.5\n"
                                       "  components:\n"
                                       "    - !AxisAlignedCuboidalDomainFilter\n"
                                       "      limits:\n"
                                       "        xi: [-.inf, 0.0]\n"
                                       "        zeta: [-1.0, 1.0]\n"
                                       "      components: !AffineMap\n"
                                       "        matrix:\n"
                                       "          s: [1.0, 0.5]\n"
                                       "        translation:\n"
                                       "          s: 0.25\n"
                                       "        components: !PolynomialMap\n"
                                       "          map:\n"
                                       "            rho: [1.5, -0.25, 2.0, 2600.0]\n"
                                       "            mu: [0.0, 0.0, 3.0e9, 3.0e10]\n"
                                       "    - !SphericalDomainFilter\n"
                                       "      radius: 1.5\n"
                                       "      centre:\n"
                                       "        xi: 1.0\n"
                                       "        zeta: 0.5\n"
                                       "      components: !ConstantMap\n"
                                       "        map:\n"
                                       "          rho: 2500.0\n"
                                       "          mu: 2.0e10\n"
                                       "    - !ConstantMap\n"
                                       "      map:\n"
                                       "        rho: 2400.0\n"
                                       "        mu: 1.0e10\n"
                                       "[lambda]:\n"
                                       "  !Any\n"
                                       "  components:\n"
                                       "    - !GroupFilter\n"
                                       "      groups: [1, 3]\n"
                                       "      components: !ConstantMap\n"
                                       "        map:\n"
                                       "          lambda: 4.0e10\n"
                                       "    - !GroupFilter\n"
                                       "      groups: 2\n"
                                       "      components: !AffineMap\n"
                                       "        matrix:\n"
                                       "          lambda: [1.0e9, 2.0e9, 3.0e9]\n"
                                       "        translation:\n"
                                       "          lambda: 5.0e10\n");
  const auto points = pointsIn(-3.0, 3.0, {0.0, 1.0, -1.0, 0.5}, 2000);
  sameAsEasi(path, {"rho", "mu", "lambda"}, points);
  sameAsEasi(path, {"lambda"}, points);
  CHECK(sc::easiParameters(path, {"x", "y", "z"}) ==
        std::vector<std::string>{"lambda", "mu", "rho"});
}

TEST_CASE("easi walker: a layered model steps and interpolates as easi" *
          doctest::test_suite("reader")) {
  const Files files;
  for (const std::string interpolation : {"lower", "upper", "linear"}) {
    CAPTURE(interpolation);
    const std::string path = files.write("layers-" + interpolation + ".yaml",
                                         "!LayeredModel\n"
                                         "map: !AffineMap\n"
                                         "  matrix:\n"
                                         "    depth: [0.0, 0.0, -1.0]\n"
                                         "  translation:\n"
                                         "    depth: 0.0\n"
                                         "interpolation: " +
                                             interpolation +
                                             "\n"
                                             "parameters: [rho, mu, lambda]\n"
                                             "nodes:\n"
                                             "  -1.0: [2500.0, 1.0e10, 1.5e10]\n"
                                             "  0.5: [2600.0, 2.0e10, 2.5e10]\n"
                                             "  1.75: [2700.0, 3.0e10, 3.5e10]\n"
                                             "  2.0: [3000.0, 3.3e10, 3.6e10]\n");
    const auto points = pointsIn(-3.0, 3.0, {1.0, -0.5, -1.75, -2.0, 0.0, -0.0}, 2000);
    sameAsEasi(path, {"rho", "mu", "lambda"}, points);
  }
}

TEST_CASE("easi walker: an evaluated model and an included file" * doctest::test_suite("reader")) {
  const Files files;
  files.write("inner.yaml",
              "!ConstantMap\n"
              "map:\n"
              "  mu: 3.2e10\n");
  const std::string path = files.write("outer.yaml",
                                       "!Switch\n"
                                       "[rho]: !EvalModel\n"
                                       "  parameters: [elevation, z]\n"
                                       "  model: !Switch\n"
                                       "    [elevation]: !ConstantMap\n"
                                       "      map:\n"
                                       "        elevation: 100.0\n"
                                       "    [z]: !AffineMap\n"
                                       "      matrix:\n"
                                       "        z: [0.0, 0.0, 1.0]\n"
                                       "      translation:\n"
                                       "        z: 0.0\n"
                                       "  components: !AffineMap\n"
                                       "    matrix:\n"
                                       "      d: [-1.0, 1.0]\n"
                                       "    translation:\n"
                                       "      d: 0.0\n"
                                       "    components: !PolynomialMap\n"
                                       "      map:\n"
                                       "        rho: [-0.5, 2700.0]\n"
                                       "[mu]: !Include inner.yaml\n");
  sameAsEasi(path, {"rho", "mu"}, pointsIn(-200.0, 200.0, {}, 500));
}

TEST_CASE("easi walker: a grid of ASAGI samples as easi's" * doctest::test_suite("reader")) {
  asagiModule();
  const Files files;
  writeAsagi(files.path("grid.nc"));
  for (const std::string interpolation : {"nearest", "linear", "cubic"}) {
    if (interpolation == "cubic" && !EasiKernel) {
      continue;
    }
    CAPTURE(interpolation);
    const std::string path = files.write("asagi-" + interpolation + ".yaml",
                                         "!AffineMap\n"
                                         "matrix:\n"
                                         "  x: [1.0, 0.0, 0.0]\n"
                                         "  y: [0.0, 1.0, 0.0]\n"
                                         "translation:\n"
                                         "  x: 0.0\n"
                                         "  y: 0.0\n"
                                         "components:\n"
                                         "  - !ASAGI\n"
                                         "    file: grid.nc\n"
                                         "    parameters: [rho, mu]\n"
                                         "    var: data\n"
                                         "    interpolation: " +
                                             interpolation +
                                             "\n"
                                             "  - !ConstantMap\n"
                                             "    map:\n"
                                             "      rho: 1.0\n"
                                             "      mu: 2.0\n");
    // within the grid, beyond it (where the next component takes the point), and on its nodes
    // and its bounds
    sameAsEasi(path,
               {"rho", "mu"},
               pointsIn(-0.5, 3.5, {0.0, 0.1, 0.3, 0.8, 1.0, 3.0, static_cast<double>(0.8F)}, 2000),
               interpolation == "nearest" ? 0.0 : GridUlps);
  }
}

TEST_CASE("easi walker: Impala as easi reads it, the grouping of its chains included" *
          doctest::test_suite("reader")) {
  const Files files;
  // chains that easi groups to the right; branches and early returns; variables; the functions of
  // Lua's math with lmathx; integers of Lua, whose zero has no sign; NaN in comparisons
  const std::string path =
      files.write("impala.yaml",
                  "!FunctionMap\n"
                  "map:\n"
                  "  a: return x - y + z - 1.5 * x / y * 2;\n"
                  "  b: |\n"
                  "    if (x > 0 && y <= 1 || z == 2) {\n"
                  "      r = sqrt(x * x + 1) - log(y * y + 2, 10);\n"
                  "      if (z != 0) { return r * pow(2.0, x); }\n"
                  "      r = r + exp(0.5 * z);\n"
                  "    } else {\n"
                  "      r = atan(x, y) + min(x, y) - max(y, z, x) + fmin(x, z);\n"
                  "    }\n"
                  "    return r - floor(z) + ceil(0 - 0.5) + abs(0 - 3) * 0;\n"
                  "  c: return 0 * (0 - 1) + deg(x) + rad(y) + trunc(z) + log(x * x + 1, 2) + "
                  "log(y * y + 1, 3);\n"
                  "  d: |\n"
                  "    if (x < y) { return 1; } else { if (x == y) { return 2; } }\n"
                  "    return 3 + tanh(z) + erf(y) + asinh(x) + cosh(0.1 * z) + round(x);\n");
  auto points = pointsIn(-3.0, 3.0, {0.0, 1.0, 2.0, -1.0}, 2000);
  // and NaN, which no comparison holds
  points.x.push_back(std::numeric_limits<double>::quiet_NaN());
  points.y.push_back(1.0);
  points.z.push_back(2.0);
  points.group.push_back(1.0);
  sameAsEasi(path, {"a", "b", "c", "d"}, points);

  // `a - k` with a small integer constant k, which Lua computes as `a + -k`: at -0.0, 0.0
  const std::string zero = files.write("zero.yaml",
                                       "!FunctionMap\n"
                                       "map:\n"
                                       "  a: return x - 0;\n"
                                       "  b: return y - 0.0;\n"
                                       "  c: return z - -3;\n"
                                       "  d: return x - 1 - 1;\n");
  Points zeros;
  zeros.x = {-0.0, 0.0, 1.0};
  zeros.y = {-0.0, 0.0, -1.0};
  zeros.z = {-0.0, -3.0, 2.0};
  zeros.group = {1.0, 1.0, 1.0};
  sameAsEasi(zero, {"a", "b", "c", "d"}, zeros);

  // a variable read where it may have no value, which easi takes from an earlier point
  const std::string earlier = files.write("earlier.yaml",
                                          "!FunctionMap\n"
                                          "map:\n"
                                          "  a: |\n"
                                          "    if (x > 0) { v = 1; }\n"
                                          "    return v;\n");
  sc::EasiQuery query;
  query.dimensions = {"x", "y", "z"};
  query.parameters = {"a"};
  CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(earlier, query)), sc::EasiUnsupported);
}

#ifdef SEISSOL_TEST_EASI_IMPALA_TO_LUA
TEST_CASE("easi walker: the Lua made of Impala is easi's, token for token" *
          doctest::test_suite("reader")) {
  // easi's conversion (easi >= 1.5) and the walker's, of the same entries in the same order: the
  // same text but for the space, or both fail
  const std::vector<std::map<std::string, std::string>> maps = {
      {{"a", "return x - y + z - 1.5 * x / y * 2;"},
       {"b",
        "if (x > 0 && y <= 1 || z == 2) {\n  r = sqrt(x * x + 1) - log(y * y + 2, 10);\n"
        "  if (z != 0) { return r * pow(2.0, x); }\n} else {\n  r = atan(x, y);\n}\n"
        "return r - floor(z) + ceil(0 - 0.5) + abs(0 - 3) * 0;"},
       {"c", "return 0 * (0 - 1) + deg(x) + rad(y) + trunc(z) + log(x * x + 1, 2);"},
       {"d", "if (x < y) { return 1; } else { if (x == y) { return 2; } }\nreturn 3 + tanh(z);"}},
      {{"rho", "return 2600.;"}, {"mu", "return 007 + x;"}, {"lambda", "return -x - -y;"}},
      {{"s", "v = 1; v = v + x; if ((x > 0) || (y < 0 && z >= 1)) { v = v / 2; } return v;"}},
      // a variable read before any assignment, a unary +, an expression for a statement, `else if`,
      // a comparison without an operator, an exponent, a character easi does not read, a brace
      // that ends the code early
      {{"a", "return q;"}},
      {{"a", "return +x;"}},
      {{"a", "x + 1; return x;"}},
      {{"a", "if (x > 0) { return 1; } else if (x < 0) { return 2; } return 3;"}},
      {{"a", "if (x) { return 1; } return 2;"}},
      {{"a", "return 1e3 * x;"}},
      {{"a", "return x % 2;"}},
      {{"a", "return x; } return y;"}},
      {}};
  const std::vector<std::string> in{"x", "y", "z"};
  const auto withoutSpace = [](std::string text) {
    text.erase(std::remove_if(
                   text.begin(), text.end(), [](unsigned char c) { return std::isspace(c) != 0; }),
               text.end());
    return text;
  };
  for (const auto& precode : maps) {
    const std::unordered_map<std::string, std::string> code(precode.begin(), precode.end());
    std::vector<sc::ImpalaEntry> entries;
    for (const auto& [name, text] : code) {
      entries.push_back(sc::ImpalaEntry{name, text, 1, true});
    }
    std::string theirs;
    bool theyFail = false;
    try {
      theirs = easi::convertImpalaToLua(code, in);
    } catch (const std::exception&) {
      theyFail = true;
    }
    std::string ours;
    bool weFail = false;
    try {
      ours = sc::impalaToLua(entries, {in.begin(), in.end()}).code;
    } catch (const sc::ImpalaError&) {
      weFail = true;
    }
    CAPTURE(theirs);
    CAPTURE(ours);
    CHECK(theyFail == weFail);
    CHECK(withoutSpace(ours) == withoutSpace(theirs));
  }
}
#endif

TEST_CASE("easi walker: a message about Impala names the line of the file" *
          doctest::test_suite("reader")) {
  const Files files;
  const std::string path = files.write("lines.yaml",
                                       "!FunctionMap\n"
                                       "map:\n"
                                       "  a: return x;\n"
                                       "  b: |\n"
                                       "    r = x;\n"
                                       "    return r + q;\n");
  const std::string calls = files.write("calls.yaml",
                                        "!FunctionMap\n"
                                        "map:\n"
                                        "  a: |\n"
                                        "    r = x;\n"
                                        "\n"
                                        "    return nofunction(r);\n");
  const std::string earlier = files.write("earlier.yaml",
                                          "!FunctionMap\n"
                                          "map:\n"
                                          "  a: return x;\n"
                                          "  b: |\n"
                                          "    if (x > 0) { v = 1; }\n"
                                          "    return v;\n");
  sc::EasiQuery query;
  query.dimensions = {"x", "y", "z"};
  query.parameters = {"a", "b"};
  const auto messageOf = [&](const std::string& file) {
    try {
      static_cast<void>(sc::easiProgram(file, query));
    } catch (const std::exception& error) {
      return std::string(error.what());
    }
    return std::string();
  };
  // a variable easi's conversion does not know; a function Lua's math does not have, which fails
  // at every point; a variable of an earlier point
  CHECK(messageOf(path).find("lines.yaml@6: Unknown variable: q") != std::string::npos);
  CHECK(messageOf(calls).find("calls.yaml@6: ") != std::string::npos);
  CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(calls, query)), sc::EasiModelError);
  CHECK(messageOf(earlier).find("earlier.yaml@6: ") != std::string::npos);
  CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(earlier, query)), sc::EasiUnsupported);
}

/// A SCEC file of 4 by 3 cells over 300 by 100, its nodes in reverse order: every value another.
inline std::string scecFile() {
  std::ostringstream text;
  text.precision(17);
  text << "200 3\n4 3 3.000000E+02 1.000000E+02\n2 1 1.5E+02 5.0E+01 1.0 1.0\n";
  for (int j = 3; j >= 0; --j) {
    for (int i = 4; i >= 0; --i) {
      text << i << " " << j << " " << 75.0 * i << " " << 100.0 / 3.0 * j;
      for (int c = 0; c < 10; ++c) {
        text << " " << (c + 1) * 1000.0 + 13.7 * i - 7.3 * j + 0.123 * i * j;
      }
      text << "\n";
    }
  }
  return text.str();
}

TEST_CASE("easi walker: a SCEC file samples as easi's" * doctest::test_suite("reader")) {
  const Files files;
  files.write("fault.txt", scecFile());
  const std::vector<std::string> parameters{
      "cohesion", "d_c", "forced_rupture_time", "mu_d", "mu_s", "s_dip", "s_normal", "s_strike"};
  for (const std::string interpolation : {"nearest", "linear", "cubic"}) {
    if (interpolation == "cubic" && !EasiKernel) {
      continue;
    }
    CAPTURE(interpolation);
    const std::string path = files.write("scec-" + interpolation + ".yaml",
                                         "!AffineMap\n"
                                         "matrix:\n"
                                         "  x: [1.0, 0.0, 0.0]\n"
                                         "  z: [0.0, 0.3, 0.0]\n"
                                         "translation:\n"
                                         "  x: 0.0\n"
                                         "  z: 0.0\n"
                                         "components: !SCECFile\n"
                                         "  file: fault.txt\n"
                                         "  interpolation: " +
                                             interpolation + "\n");
    // within the fault, beyond it, and on its nodes
    sameAsEasi(path,
               parameters,
               pointsIn(-50.0, 350.0, {0.0, 75.0, 150.0, 300.0, 100.0, 1000.0 / 9.0}, 2000),
               interpolation == "nearest" ? 0.0 : GridUlps);
  }
  // as easi's example of TPV16 has it: stresses and friction from the file, rotated into the
  // fault, by a switch of a constant, a function and the identity
  const std::string tpv16 =
      files.write("tpv16.yaml",
                  "!Switch\n"
                  "[s_xx, s_zz, s_xz, s_yy, s_xy, s_yz, cohesion, mu_s, mu_d, "
                  "d_c, forced_rupture_time]: !AffineMap\n"
                  "  matrix:\n"
                  "    x: [1.0, 0.0,  0.0]\n"
                  "    z: [0.0, 0.0, -1.0]\n"
                  "  translation:\n"
                  "    x: 0.0\n"
                  "    z: 0.0\n"
                  "  components: !SCECFile\n"
                  "    file: fault.txt\n"
                  "    components: !Switch\n"
                  "      [s_xx, s_zz, s_xz]: !ConstantMap\n"
                  "        map:\n"
                  "          s_xx: 0.0\n"
                  "          s_zz: 0.0\n"
                  "          s_xz: 0.0\n"
                  "      [s_yy, s_xy, s_yz, cohesion]: !FunctionMap\n"
                  "        map:\n"
                  "          s_yy: return -s_normal;\n"
                  "          cohesion: return -cohesion;\n"
                  "          s_xy: return s_strike;\n"
                  "          s_yz: return s_dip;\n"
                  "      [mu_s, mu_d, d_c, forced_rupture_time]: !IdentityMap {}\n");
  sameAsEasi(tpv16,
             {"s_xx",
              "s_zz",
              "s_xz",
              "s_yy",
              "s_xy",
              "s_yz",
              "cohesion",
              "mu_s",
              "mu_d",
              "d_c",
              "forced_rupture_time"},
             pointsIn(-150.0, 350.0, {0.0, 75.0, 300.0, -100.0}, 2000),
             GridUlps);

  // easi's errors: three dimensions, no file, an interpolation easi does not know
  sc::EasiQuery query;
  query.dimensions = {"x", "y", "z"};
  query.parameters = parameters;
  const auto fails = [&](const std::string& name, const std::string& text) {
    CAPTURE(name);
    CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(files.write(name, text), query)),
                    sc::EasiModelError);
  };
  fails("three.yaml", "!SCECFile\nfile: fault.txt\n");
  fails("missing.yaml",
        "!AffineMap\nmatrix:\n  x: [1.0, 0.0, 0.0]\n  z: [0.0, 0.0, 1.0]\ntranslation:\n  x: 0.0\n"
        "  z: 0.0\ncomponents: !SCECFile\n  file: nofault.txt\n");
  fails("next.yaml",
        "!AffineMap\nmatrix:\n  x: [1.0, 0.0, 0.0]\n  z: [0.0, 0.0, 1.0]\ntranslation:\n  x: 0.0\n"
        "  z: 0.0\ncomponents: !SCECFile\n  file: fault.txt\n  interpolation: next\n");

  // the file: its nodes by their indices, and every node once
  const auto data = reader::datafield::readScecFile(files.path("fault.txt"));
  CHECK(data.nodes == std::array<unsigned, 2>{5, 4});
  CHECK(data.spacing(0) == 75.0);
  CHECK(data.spacing(1) == 100.0 / 3.0);
  CHECK(data.values[(2 + 3 * 5) * 10 + 4] == 5 * 1000.0 + 13.7 * 2 - 7.3 * 3 + 0.123 * 2 * 3);
  const auto problem = [&](const std::string& text) {
    try {
      static_cast<void>(reader::datafield::readScecFile(files.write("bad.txt", text)));
    } catch (const std::runtime_error& error) {
      return std::string(error.what());
    }
    return std::string();
  };
  const std::string head = "1 2\n1 1 10.0 10.0\n0 0 0.0 0.0 0.0 0.0\n";
  const std::string values = " 0 0 1 2 3 4 5 6 7 8 9 10\n";
  CHECK(problem(head + "0 0" + values + "1 0" + values + "0 1" + values + "0 1" + values)
            .find("bad.txt:7: the node (0, 1) is given twice") != std::string::npos);
  CHECK(problem(head + "0 0" + values + "2 0" + values).find("bad.txt:5: the index along strike") !=
        std::string::npos);
  CHECK(problem(head + "0 0" + values + "1 0" + values + "0 1" + values).find("the file ends") !=
        std::string::npos);
  CHECK(problem("1 2\n1 1 -10.0 10.0\n0 0 0.0 0.0 0.0 0.0\n")
            .find("Degenerate grid in dimension 0") != std::string::npos);
}

TEST_CASE("easi walker: the special components give easi's bits" * doctest::test_suite("reader")) {
  const Files files;
  const std::vector<std::string> stress{"b_xx", "b_yy", "b_zz", "b_xy", "b_yz", "b_xz"};
  // the stress from strike and dip of Aochi and Madariaga, of strike slip, of dip slip, and of a
  // kind that depends on the point (below epsilon strike slip)
  for (const std::string faulting : {"0.0", "1.0", "y"}) {
    CAPTURE(faulting);
    const std::string path =
        files.write("am-" + std::to_string(faulting.size()) + faulting.substr(0, 1) + ".yaml",
                    "!FunctionMap\n"
                    "map:\n"
                    "  strike: return 300.0 + 10.0 * x;\n"
                    "  DipSlipFaulting: return " +
                        faulting +
                        ";\n"
                        "components: !STRESS_STR_DIP_SLIP_AM\n"
                        "  constants:\n"
                        "    mu_d:     0.25\n"
                        "    mu_s:     0.3\n"
                        "    dip:      8\n"
                        "    s_zz:     555562000\n"
                        "    cohesion: 0.4e6\n"
                        "    R:        0.7\n"
                        "    s2ratio:  0.5\n");
    sameAsEasi(path, stress, pointsIn(-3.0, 3.0, {0.0, 1.0e-17, -1.0e-17, 1.0}, 2000));
  }
  // the optimal stress of Ulrich et al., strike, dip and rake of the point
  const std::string optimal = files.write("optimal.yaml",
                                          "!AffineMap\n"
                                          "matrix:\n"
                                          "  strike: [10.0, 0.0, 0.0]\n"
                                          "  dip: [0.0, 5.0, 0.0]\n"
                                          "  rake: [0.0, 0.0, 20.0]\n"
                                          "translation:\n"
                                          "  strike: 30.0\n"
                                          "  dip: 60.0\n"
                                          "  rake: -90.0\n"
                                          "components: !OptimalStress\n"
                                          "  constants:\n"
                                          "    mu_d: 0.1\n"
                                          "    mu_s: 0.6\n"
                                          "    effectiveConfiningStress: -1.0e8\n"
                                          "    cohesion: 0.0\n"
                                          "    R: 0.7\n"
                                          "    s2ratio: 0.5\n");
  sameAsEasi(optimal, stress, pointsIn(-3.0, 3.0, {0.0}, 2000));
  // the Andersonian stress, with the vertical principal stress of the point: 1, 2 or 3
  const std::string andersonian = files.write("andersonian.yaml",
                                              "!FunctionMap\n"
                                              "map:\n"
                                              "  SH_max: return 10.0 * y + 30.0;\n"
                                              "  S_v: |\n"
                                              "    if (x < -1) { return 1; }\n"
                                              "    if (x < 1) { return 2; }\n"
                                              "    return 3;\n"
                                              "  sig_zz: return -100000000.0 + 1000000.0 * z;\n"
                                              "components: !AndersonianStress\n"
                                              "  constants:\n"
                                              "    mu_d: 0.2\n"
                                              "    mu_s: 0.6\n"
                                              "    cohesion: 1.0e6\n"
                                              "    S: 0.5\n"
                                              "    s2ratio: 0.3\n");
  sameAsEasi(andersonian, stress, pointsIn(-3.0, 3.0, {-1.0, 1.0}, 2000));
  // every input a constant: the stress is a number
  const std::string constant = files.write("constant.yaml",
                                           "!OptimalStress\n"
                                           "constants:\n"
                                           "  mu_d: 0.1\n"
                                           "  mu_s: 0.6\n"
                                           "  strike: 30.0\n"
                                           "  dip: 60.0\n"
                                           "  rake: -90.0\n"
                                           "  effectiveConfiningStress: -1.0e8\n"
                                           "  cohesion: 0.0\n"
                                           "  R: 0.7\n"
                                           "  s2ratio: 0.5\n");
  // easi requires the inputs of a component to be the dimensions: none of x, y and z
  sc::EasiQuery query;
  query.dimensions = {"x", "y", "z"};
  query.parameters = stress;
  CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(constant, query)), sc::EasiModelError);
  const std::string unknown = files.write("unknown.yaml",
                                          "!AffineMap\n"
                                          "matrix:\n"
                                          "  strike: [10.0, 0.0, 0.0]\n"
                                          "translation:\n"
                                          "  strike: 30.0\n"
                                          "components: !OptimalStress\n"
                                          "  constants:\n"
                                          "    friction: 0.1\n");
  CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(unknown, query)), sc::EasiModelError);
}

TEST_CASE("easi walker: a point no component accepts fails a check, as easi fails there" *
          doctest::test_suite("reader")) {
  const Files files;
  // x in [0, 1] or in group 2 -- the other points easi does not evaluate; a component the walker
  // does not reach fails nothing, whatever it lacks
  const std::string path = files.write("partial.yaml",
                                       "!Any\n"
                                       "components:\n"
                                       "  - !AxisAlignedCuboidalDomainFilter\n"
                                       "    limits:\n"
                                       "      x: [0.0, 1.0]\n"
                                       "      y: [-.inf, .inf]\n"
                                       "      z: [-.inf, .inf]\n"
                                       "    components: !ConstantMap\n"
                                       "      map:\n"
                                       "        rho: 2600.0\n"
                                       "  - !GroupFilter\n"
                                       "    groups: [2]\n"
                                       "    components: !ConstantMap\n"
                                       "      map:\n"
                                       "        rho: 2700.0\n"
                                       "  - !GroupFilter\n"
                                       "    groups: [7]\n"
                                       "    components: !ConstantMap\n"
                                       "      map:\n"
                                       "        mu: 1.0\n");
  sc::EasiQuery query;
  query.dimensions = {"x", "y", "z"};
  query.parameters = {"rho"};
  const auto program = sc::easiProgram(path, query);
  Points points;
  points.x = {0.5, 2.0, 2.0, std::numeric_limits<double>::quiet_NaN(), 1.0};
  points.y = {0.0, 0.0, 0.0, 0.0, 0.0};
  points.z = {0.0, 0.0, 0.0, 0.0, 0.0};
  points.group = {1.0, 1.0, 2.0, 1.0, 3.0};
  const auto failures = failuresOf(program, points);
  std::vector<std::size_t> failing;
  for (const auto& [check, point] : failures) {
    failing.push_back(point);
  }
  std::sort(failing.begin(), failing.end());
  failing.erase(std::unique(failing.begin(), failing.end()), failing.end());
  // x = 2 outside group 2, and NaN, which no box contains
  CHECK(failing == std::vector<std::size_t>{1, 3});
}

TEST_CASE("easi walker: a !LuaMap gives easi's bits" * doctest::test_suite("reader")) {
  const Files files;
  // a layer search with a break, branches on the point, `c and x or y`, Lua's max with an integer,
  // lmathx's floor and fmax, `%`, `^`; x[i], x.name and x["name"]
  const std::string path =
      files.write("luamap.yaml",
                  "!LuaMap\n"
                  "returns: [rho, mu, lambda, layer]\n"
                  "function: |\n"
                  "  local tops = {1.0, 0.0, -1.0, -2.0}\n"
                  "  function f(x)\n"
                  "    local z = x[\"z\"]\n"
                  "    local layer = #tops + 1\n"
                  "    for i, top in ipairs(tops) do\n"
                  "      if z > top then\n"
                  "        layer = i\n"
                  "        break\n"
                  "      end\n"
                  "    end\n"
                  "    local rho = 2600 + 10 * layer\n"
                  "    if x.x > 0 and x[2] < 1 then rho = rho + math.floor(x.x * 3) end\n"
                  "    local vs = x.y > 0 and 3000 or math.max(x.y, -1) * 100\n"
                  "    return {rho = rho, mu = rho * vs ^ 2, lambda = math.fmax(x.x, z) % 2,\n"
                  "            layer = layer}\n"
                  "  end\n");
  sameAsEasi(
      path, {"rho", "mu", "lambda", "layer"}, pointsIn(-3.0, 3.0, {0.0, 1.0, -1.0, -2.0}, 2000));

  // after a map, whose outputs are its dimensions, in easi's order: x[1] is eta
  const std::string mapped = files.write("mapped.yaml",
                                         "!AffineMap\n"
                                         "matrix:\n"
                                         "  xi: [1.0, 0.0, 0.0]\n"
                                         "  eta: [0.0, 1.0, 1.0]\n"
                                         "translation:\n"
                                         "  xi: 0.5\n"
                                         "  eta: -1.0\n"
                                         "components: !LuaMap\n"
                                         "  returns: [rho, mu]\n"
                                         "  function: |\n"
                                         "    function f(x)\n"
                                         "      return {rho = x.xi * 2 + x[2], mu = x[1] - x.eta}\n"
                                         "    end\n");
  sameAsEasi(mapped, {"rho", "mu"}, pointsIn(-3.0, 3.0, {}, 500));

  // where Lua fails at some points, a check fails there, with the line in the file
  const std::string failing = files.write("failing.yaml",
                                          "!LuaMap\n"
                                          "returns: [rho]\n"
                                          "function: |\n"
                                          "  function f(x)\n"
                                          "    if x.z < -2.5 then\n"
                                          "      error('too deep')\n"
                                          "    end\n"
                                          "    return {rho = 2600}\n"
                                          "  end\n");
  sc::EasiQuery query;
  query.dimensions = {"x", "y", "z"};
  query.parameters = {"rho"};
  const auto program = sc::easiProgram(failing, query);
  REQUIRE(program.checks().size() == 1);
  CHECK(program.checks().front().message.find("failing.yaml@6: too deep") != std::string::npos);
  Points points;
  points.x = {0.0, 0.0, 0.0};
  points.y = {0.0, 0.0, 0.0};
  points.z = {-3.0, -2.5, 0.0};
  points.group = {1.0, 1.0, 1.0};
  const auto failures = failuresOf(program, points);
  REQUIRE(failures.size() == 1);
  CHECK(failures.front().second == 0);

  // f keeps a value from one point to the next, which the walker refuses
  const std::string counting = files.write("counting.yaml",
                                           "!LuaMap\n"
                                           "returns: [rho]\n"
                                           "function: |\n"
                                           "  function f(x)\n"
                                           "    n = (n or 0) + 1\n"
                                           "    return {rho = n}\n"
                                           "  end\n");
  CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(counting, query)), sc::EasiUnsupported);
}

TEST_CASE("easi walker: what it does not translate, easi evaluates" *
          doctest::test_suite("reader")) {
  const Files files;
  const std::string lua = files.write("lua.yaml",
                                      "!LuaMap\n"
                                      "returns: [rho]\n"
                                      "function: |\n"
                                      "  function f(x)\n"
                                      "    return {rho = 2600.0 + math.random()}\n"
                                      "  end\n");
  sc::EasiQuery query;
  query.dimensions = {"x", "y", "z"};
  query.parameters = {"rho"};
  CHECK_THROWS_AS(static_cast<void>(sc::easiProgram(lua, query)), sc::EasiUnsupported);
  const auto reader = sc::buildReader(lua, {"x", "y", "z"});
  CHECK(dynamic_cast<sc::EasiReader*>(reader.get()) != nullptr);
  const std::string constant = files.write("constant.yaml", "!ConstantMap\nmap:\n  rho: 2.0\n");
  const auto walked = sc::buildReader(constant, {"x", "y", "z"});
  CHECK(dynamic_cast<sc::EasiProgramReader*>(walked.get()) != nullptr);
}

} // namespace seissol::unit_test
