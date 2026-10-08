// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/elastic/Model/Datastructures.h"
#include "Expr/SderivFrontend.h"
#include "Initializer/Typedefs.h"
#include "Physics/NonlinearDirichlet.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::unit_test::nonlinear_dirichlet {

const std::vector<std::string> Quantities{
    "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

/// A script written to a file of its own, removed again when done.
class ScriptFile {
  public:
  ScriptFile(const std::string& name, const std::string& content)
      : path_(
            std::filesystem::temp_directory_path() /
            ("seissol-nonlineardirichlet-" + std::to_string(std::random_device{}()) + "-" + name)) {
    std::ofstream(path_) << content;
  }
  ~ScriptFile() { std::filesystem::remove(path_); }
  ScriptFile(const ScriptFile&) = delete;
  ScriptFile& operator=(const ScriptFile&) = delete;
  ScriptFile(ScriptFile&&) = delete;
  ScriptFile& operator=(ScriptFile&&) = delete;

  [[nodiscard]] std::string path() const { return path_.string(); }

  private:
  std::filesystem::path path_;
};

/// The ghost state of `condition` at two points, at the time 2 in the simulation 1, from an inner
/// state whose quantity j is j + 1 at the first point and -(j + 1) at the second.
inline std::vector<double> ghostOf(const physics::NonlinearDirichlet& condition) {
  model::ElasticMaterial material({2.5, 1.5, 3.0});
  CellMaterialData materialData;
  materialData.local = &material;
  const std::vector<std::array<double, 3>> points = {{1.0, 2.0, 3.0}, {-1.0, 0.5, -2.0}};
  std::vector<double> inner(Quantities.size() * points.size());
  for (std::size_t j = 0; j < Quantities.size(); ++j) {
    inner[j * points.size() + 0] = static_cast<double>(j + 1);
    inner[j * points.size() + 1] = -static_cast<double>(j + 1);
  }
  std::vector<double> ghost(inner.size(), -100.0);
  condition.evaluate(
      2.0, 1, points.data(), points.size(), materialData, inner.data(), ghost.data());
  return ghost;
}

TEST_CASE("NonlinearDirichlet: the ghost state is a function of the inner state" *
          doctest::test_suite("physics")) {
  // v1 is 7 (-7) inside, s_xy 4 (-4); the script reads the inner v1 next to its own definition
  const auto check = [](const physics::NonlinearDirichlet& condition) {
    CHECK(!condition.faceAligned());
    const auto ghost = ghostOf(condition);
    CHECK(ghost[6 * 2 + 0] == -7.0);
    CHECK(ghost[6 * 2 + 1] == 7.0);
    CHECK(ghost[3 * 2 + 0] == doctest::Approx(4.0 + 0.5 * 49.0 + 1.0));
    CHECK(ghost[3 * 2 + 1] == doctest::Approx(-4.0 + 0.5 * 49.0 - 1.0));
    // the time, the simulation and the material
    CHECK(ghost[7 * 2 + 0] == doctest::Approx(2.0 * 2.5 + 1.0 + 1.5 + 3.0));
    // what the script does not give is the inner state
    for (const std::size_t j : {0, 1, 2, 4, 5, 8}) {
      CHECK(ghost[j * 2 + 0] == static_cast<double>(j + 1));
      CHECK(ghost[j * 2 + 1] == -static_cast<double>(j + 1));
    }
  };

  const ScriptFile sderiv("ghost.sderiv",
                          "out def v1 = -v1\n"
                          "out def s_xy = s_xy + 0.5 * v1 * v1 + x\n"
                          "out def v2 = t * rho + sim + mu + lambda\n");
  check(physics::NonlinearDirichlet("sderiv:" + sderiv.path(), Quantities));

  const ScriptFile lua("ghost.lua",
                       "local M = {}\n"
                       "function M.evaluate(fields, x, t, sim, rho, mu, lambda, v1, s_xy)\n"
                       "  return {v1 = -v1, s_xy = s_xy + 0.5 * v1 * v1 + x,\n"
                       "          v2 = t * rho + sim + mu + lambda}\n"
                       "end\n"
                       "return M\n");
  check(physics::NonlinearDirichlet("lua:" + lua.path(), Quantities));
}

TEST_CASE("NonlinearDirichlet: the frame states the condition in the face-aligned basis" *
          doctest::test_suite("physics")) {
  const ScriptFile script("wall.sderiv", "out def frame = 1.0\nout def v1 = -v1\n");
  const physics::NonlinearDirichlet condition("sderiv:" + script.path(), Quantities);
  CHECK(condition.faceAligned());
  CHECK(ghostOf(condition)[6 * 2 + 0] == -7.0);
}

TEST_CASE("NonlinearDirichlet: a script that cannot be a boundary condition says why" *
          doctest::test_suite("physics")) {
  const auto rejects = [](const std::string& source) {
    expr::SderivOptions options;
    options.inputs.insert(Quantities.begin(), Quantities.end());
    CHECK_THROWS_AS(
        physics::NonlinearDirichlet(expr::compileSderivModule(source, options), "test", Quantities),
        std::invalid_argument);
  };
  // a name that is no quantity
  rejects("out def v4 = v1\n");
  // an input the boundary does not have
  rejects("out def v1 = pressure\n");
  // a frame that varies, or is neither 0 nor 1
  rejects("out def frame = select(lt(x, 0.0), 1.0, 0.0)\nout def v1 = -v1\n");
  rejects("out def frame = 2.0\nout def v1 = -v1\n");
  // state
  rejects("state s = 0.0\nout def s = s + v1\nout def v1 = s\n");
}

} // namespace seissol::unit_test::nonlinear_dirichlet
