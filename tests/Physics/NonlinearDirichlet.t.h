// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Expr/SderivFrontend.h"
#include "Numerical/Quadrature.h"
#include "Physics/GhostState.h"
#include "Physics/InitialField.h"
#include "Physics/NonlinearDirichlet.h"
#include "Physics/ScriptField.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::unit_test::nonlinear_dirichlet {

const std::vector<std::string> Quantities{
    "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

/// The parameters of an elastic material, as scripts read them.
const std::vector<std::string> Parameters{"rho", "mu", "lambda"};

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

/// The ghost state of `condition`, integrated over a step with one time point, at two points: the
/// step [1.5, 2.5] has its time point at 2 and the weight 1. In the simulation 1, beyond a cell
/// whose material has the parameters `parameters` (rho = 2.5, mu = 1.5, lambda = 3), from an inner
/// state whose quantity j is j + 1 at the first point and -(j + 1) at the second.
inline std::vector<double> ghostOf(const physics::NonlinearDirichlet& condition,
                                   const std::vector<double>& parameters = std::vector<double>{
                                       2.5, 1.5, 3.0}) {
  const std::vector<std::array<double, 3>> points = {{1.0, 2.0, 3.0}, {-1.0, 0.5, -2.0}};
  const auto& integral = condition.integral();
  REQUIRE(integral.timePointCount() == 1);
  REQUIRE(integral.innerCount() == Quantities.size());
  REQUIRE(integral.componentCount() == Quantities.size());
  std::vector<double> inner(Quantities.size() * points.size());
  for (std::size_t j = 0; j < Quantities.size(); ++j) {
    inner[j * points.size() + 0] = static_cast<double>(j + 1);
    inner[j * points.size() + 1] = -static_cast<double>(j + 1);
  }
  std::vector<double> ghost(inner.size(), -100.0);
  integral.evaluate(
      1.5, 1.0, 1, points.data(), points.size(), parameters.data(), inner.data(), ghost.data());
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
  check(physics::NonlinearDirichlet("sderiv:" + sderiv.path(), Quantities, Parameters, 1));

  const ScriptFile lua("ghost.lua",
                       "local M = {}\n"
                       "function M.evaluate(fields, x, t, sim, rho, mu, lambda, v1, s_xy)\n"
                       "  return {v1 = -v1, s_xy = s_xy + 0.5 * v1 * v1 + x,\n"
                       "          v2 = t * rho + sim + mu + lambda}\n"
                       "end\n"
                       "return M\n");
  check(physics::NonlinearDirichlet("lua:" + lua.path(), Quantities, Parameters, 1));
}

TEST_CASE("NonlinearDirichlet: the frame states the condition in the face-aligned basis" *
          doctest::test_suite("physics")) {
  const ScriptFile script("wall.sderiv", "out def frame = 1.0\nout def v1 = -v1\n");
  const physics::NonlinearDirichlet condition("sderiv:" + script.path(), Quantities, Parameters, 1);
  CHECK(condition.faceAligned());
  CHECK(ghostOf(condition)[6 * 2 + 0] == -7.0);
}

TEST_CASE("NonlinearDirichlet: the script reads every parameter of the material" *
          doctest::test_suite("physics")) {
  // the parameters of a poroelastic material, by name; the values in that order
  const std::vector<std::string> parameters{"rho", "porosity", "viscosity"};
  const ScriptFile script("poro.sderiv", "out def v1 = porosity * v1 + viscosity\n");
  const physics::NonlinearDirichlet condition("sderiv:" + script.path(), Quantities, parameters, 1);
  const auto ghost = ghostOf(condition, {1000.0, 0.25, 2.0});
  CHECK(ghost[6 * 2 + 0] == doctest::Approx(0.25 * 7.0 + 2.0));
  CHECK(ghost[6 * 2 + 1] == doctest::Approx(-0.25 * 7.0 + 2.0));
}

TEST_CASE("NonlinearDirichlet: a script that cannot be a boundary condition says why" *
          doctest::test_suite("physics")) {
  const auto rejects = [](const std::string& source) {
    expr::SderivOptions options;
    options.inputs.insert(Quantities.begin(), Quantities.end());
    CHECK_THROWS_AS(
        physics::NonlinearDirichlet(
            expr::compileSderivModule(source, options), "test", Quantities, Parameters, 1),
        std::invalid_argument);
  };
  // a name that is no quantity
  rejects("out def v4 = v1\n");
  // an input the boundary does not have, also a parameter of another material
  rejects("out def v1 = pressure\n");
  rejects("out def v1 = porosity * v1\n");
  // a frame that varies, or is neither 0 nor 1
  rejects("out def frame = select(lt(x, 0.0), 1.0, 0.0)\nout def v1 = -v1\n");
  rejects("out def frame = 2.0\nout def v1 = -v1\n");
  // state
  rejects("state s = 0.0\nout def s = s + v1\nout def v1 = s\n");
}

TEST_CASE("GhostStateIntegral: the script integrated over the time points of a step" *
          doctest::test_suite("physics")) {
  expr::SderivOptions options;
  options.inputs.insert("v2");
  const auto program = expr::compileSderivModule(
      "out def v1 = t * t * x\nout def v2 = -v2\nout def s_xx = sim\n", options);
  // s_xx is not asked for; v3 is no inner component, and the script does not give it
  const physics::GhostStateIntegral integral(program,
                                             "test",
                                             {"v1", "v2", "v3"},
                                             {"v2"},
                                             Parameters,
                                             physics::GhostStateIntegral::Missing::Inner,
                                             3);
  REQUIRE(integral.timePointCount() == 3);
  REQUIRE(integral.innerCount() == 1);
  REQUIRE(integral.componentCount() == 3);

  std::array<double, 3> quadPoints{};
  std::array<double, 3> quadWeights{};
  quadrature::GaussLegendre(quadPoints.data(), quadWeights.data(), 3);
  std::array<double, 3> times{};
  integral.timePoints(1.0, 0.5, times.data());

  const std::array<double, 3> parameters{1.0, 1.0, 1.0};
  const std::vector<std::array<double, 3>> points = {{2.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
  // v2 inside is 10 + i at time point i at the first point, and -(10 + i) at the second
  std::vector<double> inner(3 * points.size());
  double expected = 0;
  for (std::size_t i = 0; i < 3; ++i) {
    CHECK(times[i] == doctest::Approx(1.25 + 0.25 * quadPoints[i]));
    inner[i * points.size() + 0] = 10.0 + static_cast<double>(i);
    inner[i * points.size() + 1] = -10.0 - static_cast<double>(i);
    expected -= 0.25 * quadWeights[i] * (10.0 + static_cast<double>(i));
  }
  std::vector<double> ghost(3 * points.size(), -100.0);
  integral.evaluate(
      1.0, 0.5, 0, points.data(), points.size(), parameters.data(), inner.data(), ghost.data());
  // the quadrature is exact for t^2: the integral over [1, 1.5] is (1.5^3 - 1) / 3
  CHECK(ghost[0] == doctest::Approx(2.0 * (3.375 - 1.0) / 3.0));
  CHECK(ghost[1] == doctest::Approx(-1.0 * (3.375 - 1.0) / 3.0));
  CHECK(ghost[2] == doctest::Approx(expected));
  CHECK(ghost[3] == doctest::Approx(-expected));
  CHECK(ghost[4] == 0.0);
  CHECK(ghost[5] == 0.0);
}

TEST_CASE("GhostStateIntegral: a quantity the script does not give is the integral of the state "
          "inside, as a device takes it" *
          doctest::test_suite("physics")) {
  // a device integrates what the script gives (Missing::Zero) and takes the rest from the
  // integrated state of the cell; the quadrature is exact for the Taylor series of the state, so
  // that is what the host computes (Missing::Inner)
  constexpr std::size_t TimePoints = 3;
  expr::SderivOptions options;
  options.inputs.insert(Quantities.begin(), Quantities.end());
  const auto program = expr::compileSderivModule("out def v1 = -v1 * x + s_xx\n", options);
  const physics::GhostStateIntegral host(program,
                                         "test",
                                         Quantities,
                                         Quantities,
                                         Parameters,
                                         physics::GhostStateIntegral::Missing::Inner,
                                         TimePoints);
  const physics::GhostStateIntegral device(program,
                                           "test",
                                           Quantities,
                                           Quantities,
                                           Parameters,
                                           physics::GhostStateIntegral::Missing::Zero,
                                           TimePoints);
  REQUIRE(device.program().outputs().size() == 1);
  // the device program reads the state inside only where the script does
  for (std::size_t i = 0; i < device.program().inputs().size(); ++i) {
    if (device.input(i) == physics::GhostStateIntegral::Input::Inner) {
      const auto component = device.inputIndex(i) % device.innerCount();
      CHECK((Quantities[component] == "v1" || Quantities[component] == "s_xx"));
    }
  }

  const double start = 0.75;
  const double width = 0.5;
  std::array<double, TimePoints> times{};
  host.timePoints(start, width, times.data());
  const std::vector<std::array<double, 3>> points = {{1.0, 2.0, 3.0}, {-0.5, 0.25, 4.0}};
  const std::size_t count = points.size();
  // the state inside is a polynomial of degree 2 in the time into the step, per quantity and point
  const auto coefficient = [](std::size_t k, std::size_t p, std::size_t power) {
    return 0.1 * static_cast<double>(k + 1) - 0.2 * static_cast<double>(p) +
           0.3 * static_cast<double>(power);
  };
  std::vector<double> inner(TimePoints * Quantities.size() * count);
  for (std::size_t i = 0; i < TimePoints; ++i) {
    const double tau = times[i] - start;
    for (std::size_t k = 0; k < Quantities.size(); ++k) {
      for (std::size_t p = 0; p < count; ++p) {
        inner[(i * Quantities.size() + k) * count + p] =
            coefficient(k, p, 0) + coefficient(k, p, 1) * tau + coefficient(k, p, 2) * tau * tau;
      }
    }
  }
  const std::array<double, 3> parameters{1.0, 2.0, 3.0};
  std::vector<double> onHost(Quantities.size() * count, -100.0);
  std::vector<double> onDevice(Quantities.size() * count, -100.0);
  host.evaluate(
      start, width, 0, points.data(), count, parameters.data(), inner.data(), onHost.data());
  device.evaluate(
      start, width, 0, points.data(), count, parameters.data(), inner.data(), onDevice.data());
  for (std::size_t k = 0; k < Quantities.size(); ++k) {
    for (std::size_t p = 0; p < count; ++p) {
      CAPTURE(k);
      CAPTURE(p);
      if (Quantities[k] == "v1") {
        CHECK(onDevice[k * count + p] == onHost[k * count + p]);
      } else {
        // the integral of the state inside, which the integrated state of the cell gives
        const double integral = coefficient(k, p, 0) * width +
                                coefficient(k, p, 1) * width * width / 2.0 +
                                coefficient(k, p, 2) * width * width * width / 3.0;
        CHECK(onDevice[k * count + p] == 0.0);
        CHECK(onHost[k * count + p] == doctest::Approx(integral).epsilon(1e-13));
      }
    }
  }
}

TEST_CASE("AnalyticalGhostState: the fields of one script share one program" *
          doctest::test_suite("physics")) {
  const auto program = expr::compileSderivModule("out def v1 = sim + t\n");
  std::vector<std::unique_ptr<physics::InitialField>> fields;
  for (std::size_t sim = 0; sim < 2; ++sim) {
    fields.push_back(
        std::make_unique<physics::ScriptField>(program, "test", Quantities, Parameters, sim));
  }
  const auto state =
      physics::AnalyticalGhostState::fromFields(fields, Quantities, Parameters, 2, 2);
  REQUIRE(state != nullptr);
  CHECK(&state->integral(0) == &state->integral(1));
  CHECK(state->scriptSimulation(0) == 0);
  CHECK(state->scriptSimulation(1) == 1);
  CHECK(state->integral(1).innerCount() == 0);

  // the integral of 1 + t over [0, 2] in the simulation 1
  const std::array<double, 3> parameters{1.0, 1.0, 1.0};
  const std::array<double, 3> point{0.0, 0.0, 0.0};
  std::vector<double> ghost(Quantities.size(), -100.0);
  state->integral(1).evaluate(0.0, 2.0, 1, &point, 1, parameters.data(), nullptr, ghost.data());
  CHECK(ghost[6] == doctest::Approx(4.0));
  CHECK(ghost[0] == 0.0);

  // a field that is no script leaves the boundary to evaluate the fields one time after the other
  std::vector<std::unique_ptr<physics::InitialField>> other;
  other.push_back(std::make_unique<physics::ZeroField>());
  CHECK(physics::AnalyticalGhostState::fromFields(other, Quantities, Parameters, 1, 2) == nullptr);
}

} // namespace seissol::unit_test::nonlinear_dirichlet
