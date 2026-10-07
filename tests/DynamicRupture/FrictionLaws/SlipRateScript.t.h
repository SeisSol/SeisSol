// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_SLIPRATESCRIPT_T_H_
#define SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_SLIPRATESCRIPT_T_H_

#include <doctest.h>

#include "DynamicRupture/FrictionLaws/SlipRateScript.h"
#include "Expr/SderivFrontend.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "Reader/Scripting/DataTable.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace sliprate {

using seissol::dr::friction_law::SlipRateEvaluator;
using seissol::dr::friction_law::SlipRateScript;
using Row = SlipRateScript::Row;

/// The rows of a layer of `points` points, filled per row by `value(row, point)`, and the slip
/// rates of `steps` sub-steps.
struct Layer {
  std::size_t points{0};
  std::vector<double> rows;
  std::vector<double> slipRates;

  Layer(const SlipRateScript& script,
        std::size_t points,
        std::size_t steps,
        const std::function<double(Row, const std::string&, std::size_t)>& value)
      : points(points), rows(script.rows() * points), slipRates(2 * steps * points, -1.0) {
    for (std::size_t row = 0; row < script.rows(); ++row) {
      for (std::size_t p = 0; p < points; ++p) {
        rows[row * points + p] = value(script.rowKinds()[row], script.rowNames()[row], p);
      }
    }
  }

  [[nodiscard]] double rate(std::size_t step, std::size_t direction, std::size_t point) const {
    return slipRates[(2 * step + direction) * points + point];
  }
};

/// Evaluates `script` on `layer` for a step starting at `start` with sub-steps `deltaT`.
inline void evaluate(const std::shared_ptr<const SlipRateScript>& script,
              Layer& layer,
              double start,
              const std::vector<double>& deltaT) {
  SlipRateEvaluator evaluator(script, reader::scripting::DataType::F64, deltaT.size());
  SlipRateEvaluator::LayerData data;
  data.parameters = layer.rows.data();
  data.slipRates = layer.slipRates.data();
  evaluator.setLayer(layer.points, data, data);
  CHECK(!evaluator.evaluate(start, deltaT, nullptr, false));
}

} // namespace sliprate

TEST_CASE("The slip of a script imposes its increment over each sub-step, rotated into the face") {
  using namespace sliprate;
  const auto script = std::make_shared<const SlipRateScript>(
      expr::compileSderivModule("out def slip_strike = strike_slip * t * t\n"
                                "out def slip_dip = dip_slip * t\n"),
      "test");
  CHECK(script->cumulative());
  // the fault parameters and the rotation; the times of a sub-step are no rows
  REQUIRE(script->rows() == 4);
  CHECK(std::count(script->rowKinds().begin(), script->rowKinds().end(), Row::Parameter) == 2);
  CHECK(std::count(script->rowKinds().begin(), script->rowKinds().end(), Row::Cos) == 1);
  CHECK(std::count(script->rowKinds().begin(), script->rowKinds().end(), Row::Sin) == 1);

  constexpr std::size_t Points = 2500; // several host calls
  const std::vector<double> deltaT = {0.1, 0.25, 0.05};
  const double start = 1.5;
  const auto value = [](Row kind, const std::string& name, std::size_t p) {
    const double angle = 0.001 * static_cast<double>(p);
    switch (kind) {
    case Row::Cos:
      return std::cos(angle);
    case Row::Sin:
      return std::sin(angle);
    default:
      return name == "strike_slip" ? 1.0 + 0.01 * static_cast<double>(p)
                                   : -2.0 + 0.003 * static_cast<double>(p);
    }
  };
  Layer layer(*script, Points, deltaT.size(), value);
  evaluate(script, layer, start, deltaT);

  double begin = start;
  for (std::size_t step = 0; step < deltaT.size(); ++step) {
    const double end = begin + deltaT[step];
    for (std::size_t p = 0; p < Points; p += 7) {
      const double strike =
          value(Row::Parameter, "strike_slip", p) * (end * end - begin * begin) / deltaT[step];
      const double dip = value(Row::Parameter, "dip_slip", p) * (end - begin) / deltaT[step];
      const double cos = value(Row::Cos, "", p);
      const double sin = value(Row::Sin, "", p);
      CHECK(layer.rate(step, 0, p) == doctest::Approx(cos * strike + sin * dip).epsilon(1e-12));
      CHECK(layer.rate(step, 1, p) == doctest::Approx(cos * dip - sin * strike).epsilon(1e-12));
    }
    begin = end;
  }
}

TEST_CASE("The slip rate of a script is imposed as it is given at the end of a sub-step") {
  using namespace sliprate;
  const auto script = std::make_shared<const SlipRateScript>(
      expr::compileSderivModule("out def slip_rate_strike = 2.0 * t + dt\n"
                                "out def slip_rate_dip = x * sim\n"),
      "test");
  CHECK(!script->cumulative());
  CHECK(std::count(script->rowKinds().begin(), script->rowKinds().end(), Row::X) == 1);
  CHECK(std::count(script->rowKinds().begin(), script->rowKinds().end(), Row::Simulation) == 1);

  constexpr std::size_t Points = 10;
  const std::vector<double> deltaT = {0.5, 0.5};
  const auto value = [](Row kind, const std::string& /*name*/, std::size_t p) {
    switch (kind) {
    case Row::Cos:
      return 1.0;
    case Row::Sin:
      return 0.0;
    case Row::X:
      return static_cast<double>(p);
    default:
      return 3.0;
    }
  };
  Layer layer(*script, Points, deltaT.size(), value);
  evaluate(script, layer, 0.0, deltaT);
  for (std::size_t p = 0; p < Points; ++p) {
    CHECK(layer.rate(0, 0, p) == 2.0 * 0.5 + 0.5);
    CHECK(layer.rate(1, 0, p) == 2.0 * 1.0 + 0.5);
    CHECK(layer.rate(1, 1, p) == 3.0 * static_cast<double>(p));
  }
}

TEST_CASE("The slip of a script adds up over the sub-steps, and agrees with FL 34") {
  using namespace sliprate;
  // the cumulative slip of FL 34: the smooth step of the Gaussian source time function
  const auto script = std::make_shared<const SlipRateScript>(
      expr::compileSderivModule(
          "def tau = t - rupture_onset\n"
          "def ramp = exp((tau - rupture_rise_time) * (tau - rupture_rise_time) / "
          "(tau * (tau - 2.0 * rupture_rise_time)))\n"
          "def smooth = select(le(tau, 0.0), 0.0, select(lt(tau, rupture_rise_time), ramp, 1.0))\n"
          "out def slip_strike = strike_slip * smooth\n"
          "out def slip_dip = dip_slip * smooth\n"),
      "test");

  constexpr std::size_t Points = 64;
  const double riseTime = 0.8;
  const auto value = [&](Row kind, const std::string& name, std::size_t p) {
    switch (kind) {
    case Row::Cos:
      return 1.0;
    case Row::Sin:
      return 0.0;
    default:
      if (name == "rupture_onset") {
        return 0.01 * static_cast<double>(p);
      }
      if (name == "rupture_rise_time") {
        return riseTime;
      }
      return name == "strike_slip" ? 2.0 : 0.0;
    }
  };

  // steps of four sub-steps over [0, 2], as the friction law takes them
  const std::vector<double> deltaT = {0.01, 0.015, 0.0125, 0.0125};
  std::vector<double> slip(Points, 0.0);
  double begin = 0.0;
  double maxDifference = 0.0;
  for (std::size_t stepCount = 0; stepCount < 40; ++stepCount) {
    Layer layer(*script, Points, deltaT.size(), value);
    evaluate(script, layer, begin, deltaT);
    for (std::size_t step = 0; step < deltaT.size(); ++step) {
      const double end = begin + deltaT[step];
      for (std::size_t p = 0; p < Points; ++p) {
        slip[p] += layer.rate(step, 0, p) * deltaT[step];
        // the slip rate FL 34 imposes, from the same time function
        const double builtIn =
            2.0 *
            gaussianNucleationFunction::smoothStepIncrement(
                end - value(Row::Parameter, "rupture_onset", p), deltaT[step], riseTime) /
            deltaT[step];
        maxDifference = std::max(maxDifference, std::abs(layer.rate(step, 0, p) - builtIn));
        CHECK(layer.rate(step, 1, p) == 0.0);
      }
      begin = end;
    }
  }
  // the increments add up to the slip at the end
  for (std::size_t p = 0; p < Points; ++p) {
    const double onset = value(Row::Parameter, "rupture_onset", p);
    CHECK(slip[p] ==
          doctest::Approx(2.0 * gaussianNucleationFunction::smoothStep(begin - onset, riseTime))
              .epsilon(1e-13));
  }
  // rates of up to some 3 m/s, from differences of the slip over sub-steps of a hundredth of the
  // rise time
  CHECK(maxDifference < 1e-11);
}

TEST_CASE("A slip rate script that cannot be used says why") {
  using namespace sliprate;
  const auto rejects = [](const std::string& source) {
    CHECK_THROWS_AS(SlipRateScript(expr::compileSderivModule(source), "test"),
                    std::invalid_argument);
  };
  // one direction only
  rejects("out def slip_strike = t\n");
  // the slip along one direction, the rate along the other
  rejects("out def slip_strike = t\nout def slip_rate_dip = 1.0\n");
  // something else
  rejects("out def slip_strike = t\nout def slip_dip = t\nout def other = 1.0\n");
  // state
  rejects("state s = 0.0\nout def s = s + 1.0\nout def slip_strike = s\nout def slip_dip = s\n");
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_SLIPRATESCRIPT_T_H_
