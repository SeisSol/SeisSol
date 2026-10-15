// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Physics/Scenario/Scenarios.h"

#include "Common/ConfigLayout.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Typedefs.h"
#include "Model/CommonDatastructures.h"
#include "Model/MaterialType.h"
#include "Model/ScriptParameters.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/ScriptedScenario.h"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics {

namespace {

using scenario::concat;
using scenario::literal;

/// The expression of a field that is the same in the stresses s_xx, s_yy and s_zz.
std::string normalStresses(const std::string& value) {
  return concat(
      {"out def s_xx = ", value, "\nout def s_yy = ", value, "\nout def s_zz = ", value, "\n"});
}

} // namespace

std::vector<std::unique_ptr<InitialField>> zeroField(ConfigId config) {
  const auto names = scenario::stateNamesOf(config);
  const auto source =
      concat({"# no wave, written by SeisSol\nout def ", names.front(), " = 0.0\n"});
  return scenario::scriptedFields(source, config, "the zero field", false);
}

std::vector<std::unique_ptr<InitialField>> scholteWave(ConfigId config) {
  // the constants of the solution, as the closed form gives them; a cell is acoustic if its mu is
  // zero (below the epsilon of a double)
  const std::string source = concat({
      "# a Scholte wave along z = 0, written by SeisSol: acoustic above, elastic below\n",
      "def omega = ",
      literal(2.0 * std::acos(-1.0)),
      "\ndef omega2 = omega * omega\n",
      "def acoustic = lt(abs(mu), ",
      literal(std::numeric_limits<double>::epsilon()),
      ")\n",
      "def phase = omega * t - 1.406466352506808 * omega * x\n",
      "def s = sin(phase)\n",
      "def c = cos(phase)\n",
      "def ea = exp(-0.98901344820674908 * omega * z)\n",
      "def e1 = exp(0.98901344820674908 * omega * z)\n",
      "def e2 = exp(1.2825031256883821 * omega * z)\n",
      "def sa = 0.35944997730200889 * omega2 * ea * s\n",
      "out def s_xx = select(acoustic, sa, -2.7820282741590652 * omega2 * e1 * s + ",
      "3.5151973269883681 * omega2 * e2 * s)\n",
      "out def s_yy = select(acoustic, sa, -6.6613381477509402e-16 * omega2 * e1 * s + ",
      "0.27315475753283058 * omega2 * e2 * s)\n",
      "out def s_zz = select(acoustic, sa, 2.7820282741590621 * omega2 * e1 * s - ",
      "2.4225782968570462 * omega2 * e2 * s)\n",
      "out def s_xz = select(acoustic, 0.0, -2.956295201467618 * omega2 * e1 * c + ",
      "2.9562952014676029 * omega2 * e2 * c)\n",
      "out def v1 = select(acoustic, -0.50555429848461109 * omega2 * ea * s, ",
      "0.98901344820675241 * omega2 * e1 * s - 1.1525489264912381 * omega2 * e2 * s)\n",
      "out def v3 = select(acoustic, 0.35550086150929727 * omega2 * ea * c, ",
      "1.406466352506812 * omega2 * e1 * c - 1.050965490997515 * omega2 * e2 * c)\n",
  });
  return scenario::scriptedFields(source, config, "the Scholte wave", false);
}

std::vector<std::unique_ptr<InitialField>> snellsLaw(ConfigId config) {
  const std::string source = concat({
      "# a plane wave refracted at z = 0, written by SeisSol: acoustic above, elastic below\n",
      "def omega = ",
      literal(2.0 * std::acos(-1.0)),
      "\ndef acoustic = lt(abs(mu), ",
      literal(std::numeric_limits<double>::epsilon()),
      ")\n",
      // the incident and the reflected wave in the acoustic material
      "def ai = sin(omega * t - omega * (0.19866933079506119 * x + 0.98006657784124163 * z))\n",
      "def ar = sin(omega * t - omega * (0.19866933079506149 * x - 0.98006657784124152 * z))\n",
      // the refracted P and S wave in the elastic material
      "def ep = sin(omega * t - 1.0 / 2.0 * omega * (0.39733866159012299 * x + ",
      "0.91767204817721759 * z))\n",
      "def es = sin(omega * t - 1.0 / 3.0 * omega * (0.59600799238518454 * x + ",
      "0.8029785009656123 * z))\n",
      "def sa = 1.0 * omega * ai + 0.48055591432167399 * omega * ar\n",
      "out def s_xx = select(acoustic, sa, -0.59005639909185559 * omega * ep + ",
      "0.55554011463785213 * omega * es)\n",
      "out def s_yy = select(acoustic, sa, 0.14460396298676709 * omega * es)\n",
      "out def s_zz = select(acoustic, sa, 0.59005639909185559 * omega * ep + ",
      "0.89049951522981918 * omega * es)\n",
      "out def s_xz = select(acoustic, 0.0, -0.55363837274201066 * omega * ep + ",
      "0.55363837274201 * omega * es)\n",
      "out def v1 = select(acoustic, -0.19866933079506119 * omega * ai - ",
      "0.095471721907895893 * omega * ar, 0.37125533967075403 * omega * ep - ",
      "0.2585553530120539 * omega * es)\n",
      "out def v3 = select(acoustic, -0.98006657784124163 * omega * ai + ",
      "0.47097679041061191 * omega * ar, -0.16074816713222639 * omega * ep - ",
      "0.34834162029840349 * omega * es)\n",
  });
  return scenario::scriptedFields(source, config, "Snell's law", false);
}

std::vector<std::unique_ptr<InitialField>>
    ocean(int mode, double gravitationalAcceleration, ConfigId config) {
  if (mode < 0 || mode > 2) {
    throw std::runtime_error("Wave mode " + std::to_string(mode) + " is not supported.");
  }
  const auto g = gravitationalAcceleration;
  if (std::abs(g - 9.81e-3) > 10e-15) {
    logError() << "Ocean scenario only supports g=9.81e-3 currently!";
  }

  const double pi = std::acos(-1.0);
  const double lx = 10.0; // km
  const double ly = 10.0; // km
  const double kX = pi / lx;
  const double kY = pi / ly;
  // computed ahead, with more digits than they would have on the fly
  constexpr std::array<double, 3> KStars{0.4433813748841239, 1.5733628061766445, 4.713305873881573};
  constexpr std::array<double, 3> Omegas{
      0.0425599572628432, 2.4523337594491745, 7.1012991617572165};
  const auto kStar = KStars[mode];
  const auto omega = Omegas[mode];
  const auto b = g * kStar / (omega * omega);

  // the profile in depth: hyperbolic for the gravity mode, harmonic for the others
  const bool gravity = mode == 0;
  const std::string profile =
      gravity ? "(sinh(kz * z) + b * cosh(kz * z))" : "(sin(kz * z) + b * cos(kz * z))";
  const std::string profileZ =
      gravity ? "(cosh(kz * z) + b * sinh(kz * z))" : "(cos(kz * z) - b * sin(kz * z))";

  // an elastic material carries the pressure in its normal stresses, another one in its first
  // quantity
  const auto names = scenario::stateNamesOf(config);
  const auto velocityOffset = configLayout(config).velocityOffset;
  const bool elastic = configValue(config).materialType == model::MaterialType::Elastic;
  const std::string pressure = "pressure";
  const std::string stresses =
      elastic ? normalStresses(pressure) : concat({"out def ", names.front(), " = pressure\n"});

  const std::string source = concat({
      "# a standing wave in an ocean, written by SeisSol (mode ",
      std::to_string(mode),
      ")\n",
      "def kx = ",
      literal(kX),
      "\ndef ky = ",
      literal(kY),
      "\ndef kz = ",
      literal(kStar),
      "\ndef omega = ",
      literal(omega),
      "\ndef b = ",
      literal(b),
      "\ndef pressure = -sin(kx * x) * sin(ky * y) * sin(omega * t) * ",
      profile,
      "\n",
      stresses,
      "out def ",
      names.at(velocityOffset),
      " = kx / (omega * rho) * cos(kx * x) * sin(ky * y) * cos(omega * t) * ",
      profile,
      "\nout def ",
      names.at(velocityOffset + 1),
      " = ky / (omega * rho) * sin(kx * x) * cos(ky * y) * cos(omega * t) * ",
      profile,
      "\nout def ",
      names.at(velocityOffset + 2),
      " = kz / (omega * rho) * sin(kx * x) * sin(ky * y) * cos(omega * t) * ",
      profileZ,
      "\n",
  });

  const auto mu = model::parameterIndexOf(scenario::parameterNamesOf(config), "mu").value_or(0);
  return scenario::scriptedFields(
      source,
      config,
      "the ocean",
      false,
      [mu](const double* parameters) -> std::optional<std::string> {
        if (parameters[mu] > 10e-15) {
          return "Ocean scenario only works for acoustic material (mu = 0.0)!";
        }
        return std::nullopt;
      });
}

std::vector<std::unique_ptr<InitialField>>
    pressureInjection(const initializer::parameters::InitializationParameters& parameters,
                      ConfigId config) {
  const auto o1 = parameters.origin[0];
  const auto o2 = parameters.origin[1];
  const auto o3 = parameters.origin[2];
  logInfo() << "Prepare gaussian pressure perturbation with center at (" << o1 << ", " << o2 << ", "
            << o3 << "), magnitude = " << parameters.magnitude << ", width = " << parameters.width
            << ".";

  const std::string source = concat({
      "# a gaussian perturbation of the pressure, written by SeisSol\n",
      "def dx = x - ",
      literal(o1),
      "\ndef dy = y - ",
      literal(o2),
      "\ndef dz = z - ",
      literal(o3),
      "\nout def p = ",
      literal(parameters.magnitude),
      " * exp(-",
      literal(parameters.width),
      " * (dx * dx + dy * dy + dz * dz))\n",
  });
  return scenario::scriptedFields(source, config, "the pressure injection", false);
}

std::vector<std::unique_ptr<InitialField>>
    acousticTravellingWaveITM(const CellMaterialData& materialData,
                              const AcousticTravellingWaveParametersITM& parameters,
                              ConfigId config) {
  const double rho0 = materialData.local->getDensity();
  const double c0 = std::sqrt(materialData.local->getLambdaBar() / rho0);
  logInfo() << "Starting Test for Acoustic Travelling Wave with ITM";
  logInfo() << "rho0 = " << rho0;
  logInfo() << "c0 = " << c0;
  logInfo() << "k = " << parameters.k;

  // before, during and after the mirror, which starts at tm and lasts tau
  const std::string source = concat({
      "# an acoustic wave travelling through a time mirror, written by SeisSol\n",
      "def rho0 = ",
      literal(rho0),
      "\ndef c0 = ",
      literal(c0),
      "\ndef k = ",
      literal(parameters.k),
      "\ndef n = ",
      literal(parameters.itmVelocityScalingFactor),
      "\ndef tm = ",
      literal(parameters.itmStartingTime),
      "\ndef tau = ",
      literal(parameters.itmDuration),
      "\ndef tp = tm + tau\n",
      // before
      "def v0 = cos(k * x - c0 * k * t)\n",
      "def p0 = c0 * rho0 * v0\n",
      // during
      "def a1 = cos(k * x + c0 * k * n * t - (c0 * k * n + c0 * k) * tm)\n",
      "def b1 = cos(k * x - c0 * k * n * t + (c0 * k * n - c0 * k) * tm)\n",
      "def p1 = -0.5 * (n - 1.0) * c0 * rho0 * a1 + 0.5 * (n + 1.0) * c0 * rho0 * b1\n",
      "def v1m = 0.5 * (n - 1.0) * a1 + 0.5 * (n + 1.0) * b1\n",
      // after
      "def a2 = cos(k * x + c0 * k * t - 2.0 * c0 * k * tm - (c0 * k * n + c0 * k) * tau)\n",
      "def b2 = cos(k * x + c0 * k * t - 2.0 * c0 * k * tm + (c0 * k * n - c0 * k) * tau)\n",
      "def c2 = cos(k * x - c0 * k * t + (c0 * k * n + c0 * k) * tau)\n",
      "def d2 = cos(k * x - c0 * k * t - (c0 * k * n - c0 * k) * tau)\n",
      "def p2 = -0.25 * (1.0 / n) * c0 * rho0 * ((-n * n + 1.0) * a2 + (n * n - 1.0) * b2 + ",
      "(n * n - 2.0 * n + 1.0) * c2 + (-n * n - 2.0 * n - 1.0) * d2)\n",
      "def v2 = (-0.25 / n) * ((n * n - 1.0) * a2 + (-n * n + 1.0) * b2 + ",
      "(n * n - 2.0 * n + 1.0) * c2 + (-n * n - 2.0 * n - 1.0) * d2)\n",
      "def pressure = select(le(t, tm), p0, select(le(t, tp), p1, p2))\n",
      normalStresses("-pressure"),
      "out def v1 = select(le(t, tm), v0, select(le(t, tp), v1m, v2))\n",
  });
  return scenario::scriptedFields(source, config, "the acoustic travelling wave with ITM", false);
}

} // namespace seissol::physics
