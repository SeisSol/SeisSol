// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigRegistry.h"
#include "Common/Real.h"
#include "Equations/Datastructures.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/ParameterReader.h"
#include "Initializer/Typedefs.h"
#include "Model/MaterialType.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/PlaneWave.h"
#include "Physics/Scenario/Registry.h"
#include "Physics/Scenario/Scenarios.h"
#include "TestConfigs.h"
#include "TestHelper.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <math.h>
#include <memory>
#include <random>
#include <type_traits>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <yateto.h>

namespace seissol::unit_test::scenarios {

// The scenarios as SeisSol computed them before they became scripts, kept verbatim apart from the
// plumbing: one point into `state`, in the layout of an elastic material.

// NOLINTBEGIN(readability-identifier-naming,readability-function-cognitive-complexity)
inline std::array<double, 9> scholteReference(double time,
                                              const std::array<double, 3>& point,
                                              const CellMaterialData& materialData) {
  using RealT = double;
  std::array<double, 9> state{};
  const std::array<double, 3>* points = &point;
  const std::size_t i = 0;
  const RealT omega = 2.0 * std::acos(-1);
  const auto& x = points[i];
  const bool isAcousticPart =
      std::abs(materialData.local->getMuBar()) < std::numeric_limits<RealT>::epsilon();
  const auto x1 = x[0];
  const auto x3 = x[2];
  const auto t = time;
  if (isAcousticPart) {
    state[0] = 0.35944997730200889 * std::pow(omega, 2) *
               std::exp(-0.98901344820674908 * omega * x3) *
               std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_xx
    state[1] = 0.35944997730200889 * std::pow(omega, 2) *
               std::exp(-0.98901344820674908 * omega * x3) *
               std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_yy
    state[2] = 0.35944997730200889 * std::pow(omega, 2) *
               std::exp(-0.98901344820674908 * omega * x3) *
               std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_zz
    state[3] = 0;                                                    // sigma_xy
    state[4] = 0;                                                    // sigma_yz
    state[5] = 0;                                                    // sigma_xz
    state[6] = -0.50555429848461109 * std::pow(omega, 2) *
               std::exp(-0.98901344820674908 * omega * x3) *
               std::sin(omega * t - 1.406466352506808 * omega * x1); // u
    state[7] = 0;                                                    // v
    state[8] = 0.35550086150929727 * std::pow(omega, 2) *
               std::exp(-0.98901344820674908 * omega * x3) *
               std::cos(omega * t - 1.406466352506808 * omega * x1); // w
  } else {
    state[0] = -2.7820282741590652 * std::pow(omega, 2) *
                   std::exp(0.98901344820674908 * omega * x3) *
                   std::sin(omega * t - 1.406466352506808 * omega * x1) +
               3.5151973269883681 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
                   std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_xx
    state[1] =
        -6.6613381477509402e-16 * std::pow(omega, 2) * std::exp(0.98901344820674908 * omega * x3) *
            std::sin(omega * t - 1.406466352506808 * omega * x1) +
        0.27315475753283058 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
            std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_yy
    state[2] = 2.7820282741590621 * std::pow(omega, 2) *
                   std::exp(0.98901344820674908 * omega * x3) *
                   std::sin(omega * t - 1.406466352506808 * omega * x1) -
               2.4225782968570462 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
                   std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_zz
    state[3] = 0;                                                        // sigma_xy
    state[4] = 0;                                                        // sigma_yz
    state[5] = -2.956295201467618 * std::pow(omega, 2) *
                   std::exp(0.98901344820674908 * omega * x3) *
                   std::cos(omega * t - 1.406466352506808 * omega * x1) +
               2.9562952014676029 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
                   std::cos(omega * t - 1.406466352506808 * omega * x1); // sigma_xz
    state[6] = 0.98901344820675241 * std::pow(omega, 2) *
                   std::exp(0.98901344820674908 * omega * x3) *
                   std::sin(omega * t - 1.406466352506808 * omega * x1) -
               1.1525489264912381 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
                   std::sin(omega * t - 1.406466352506808 * omega * x1); // u
    state[7] = 0;                                                        // v
    state[8] = 1.406466352506812 * std::pow(omega, 2) * std::exp(0.98901344820674908 * omega * x3) *
                   std::cos(omega * t - 1.406466352506808 * omega * x1) -
               1.050965490997515 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
                   std::cos(omega * t - 1.406466352506808 * omega * x1); // w
  }
  return state;
}

inline std::array<double, 9> snellReference(double time,
                                            const std::array<double, 3>& point,
                                            const CellMaterialData& materialData) {
  using RealT = double;
  std::array<double, 9> state{};
  const std::array<double, 3>* points = &point;
  const std::size_t i = 0;
  const double pi = std::acos(-1);
  const double omega = 2.0 * pi;
  const auto& x = points[i];
  const bool isAcousticPart =
      std::abs(materialData.local->getMuBar()) < std::numeric_limits<RealT>::epsilon();

  const auto x1 = x[0];
  const auto x3 = x[2];
  const auto t = time;
  if (isAcousticPart) {
    state[0] =
        1.0 * omega *
            std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
        0.48055591432167399 * omega *
            std::sin(omega * t -
                     omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // sigma_xx
    state[1] =
        1.0 * omega *
            std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
        0.48055591432167399 * omega *
            std::sin(omega * t -
                     omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // sigma_yy
    state[2] =
        1.0 * omega *
            std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
        0.48055591432167399 * omega *
            std::sin(omega * t -
                     omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // sigma_zz
    state[3] = 0;                                                                    // sigma_xy
    state[4] = 0;                                                                    // sigma_yz
    state[5] = 0;                                                                    // sigma_xz
    state[6] =
        -0.19866933079506119 * omega *
            std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) -
        0.095471721907895893 * omega *
            std::sin(omega * t -
                     omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // u
    state[7] = 0;                                                                    // v
    state[8] =
        -0.98006657784124163 * omega *
            std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
        0.47097679041061191 * omega *
            std::sin(omega * t -
                     omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // w
  } else {
    state[0] = -0.59005639909185559 * omega *
                   std::sin(omega * t - 1.0 / 2.0 * omega *
                                            (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) +
               0.55554011463785213 * omega *
                   std::sin(omega * t -
                            1.0 / 3.0 * omega *
                                (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_xx
    state[1] =
        0.14460396298676709 * omega *
        std::sin(omega * t - 1.0 / 3.0 * omega *
                                 (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_yy
    state[2] = 0.59005639909185559 * omega *
                   std::sin(omega * t - 1.0 / 2.0 * omega *
                                            (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) +
               0.89049951522981918 * omega *
                   std::sin(omega * t -
                            1.0 / 3.0 * omega *
                                (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_zz
    state[3] = 0;                                                                      // sigma_xy
    state[4] = 0;                                                                      // sigma_yz
    state[5] = -0.55363837274201066 * omega *
                   std::sin(omega * t - 1.0 / 2.0 * omega *
                                            (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) +
               0.55363837274201 * omega *
                   std::sin(omega * t -
                            1.0 / 3.0 * omega *
                                (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_xz
    state[6] =
        0.37125533967075403 * omega *
            std::sin(omega * t -
                     1.0 / 2.0 * omega * (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) -
        0.2585553530120539 * omega *
            std::sin(omega * t -
                     1.0 / 3.0 * omega * (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // u
    state[7] = 0;                                                                               // v
    state[8] =
        -0.16074816713222639 * omega *
            std::sin(omega * t -
                     1.0 / 2.0 * omega * (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) -
        0.34834162029840349 * omega *
            std::sin(omega * t -
                     1.0 / 3.0 * omega * (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // w
  }
  return state;
}

inline std::array<double, 9> oceanReference(int mode_,
                                            double gravitationalAcceleration_,
                                            double time,
                                            const std::array<double, 3>& point,
                                            const CellMaterialData& materialData) {
  std::array<double, 9> state{};
  const std::array<double, 3>* points = &point;
  const std::size_t i = 0;
  const bool elastic_ = true;
  const std::size_t velocityOffset_ = 6;
  const auto x = points[i][0];
  const auto y = points[i][1];
  const auto z = points[i][2];
  const auto t = time;

  const auto g = gravitationalAcceleration_;
  const double pi = std::acos(-1);
  const double rho = materialData.local->getDensity();

  const double lx = 10.0;    // km
  const double ly = 10.0;    // km
  const double kX = pi / lx; // 1/km
  const double kY = pi / ly; // 1/km

  constexpr auto KStars =
      std::array<double, 3>{0.4433813748841239, 1.5733628061766445, 4.713305873881573};

  // Note: Could be computed on the fly but it's better to pre-compute them with higher precision!
  constexpr auto Omegas =
      std::array<double, 3>{0.0425599572628432, 2.4523337594491745, 7.1012991617572165};

  const auto kStar = KStars[mode_];
  const auto omega = Omegas[mode_];

  const auto b = g * kStar / (omega * omega);
  constexpr auto ScalingFactor = 1;

  const auto setStresses = [&](double value) {
    if (elastic_) {
      state[0] = value;
      state[1] = value;
      state[2] = value;

      // Shear stresses are zero for elastic
      state[3] = 0.0;
      state[4] = 0.0;
      state[5] = 0.0;
    } else {
      state[0] = value;
    }
  };

  const auto uIdx = velocityOffset_;
  const auto vIdx = velocityOffset_ + 1;
  const auto wIdx = velocityOffset_ + 2;

  if (mode_ == 0) {
    // Gravity mode
    const auto pressure = -std::sin(kX * x) * std::sin(kY * y) * std::sin(omega * t) *
                          (std::sinh(kStar * z) + b * std::cosh(kStar * z));

    setStresses(ScalingFactor * pressure);

    state[uIdx] = ScalingFactor * (kX / (omega * rho)) * std::cos(kX * x) * std::sin(kY * y) *
                  std::cos(omega * t) * (std::sinh(kStar * z) + b * std::cosh(kStar * z));
    state[vIdx] = ScalingFactor * (kY / (omega * rho)) * std::sin(kX * x) * std::cos(kY * y) *
                  std::cos(omega * t) * (std::sinh(kStar * z) + b * std::cosh(kStar * z));
    state[wIdx] = ScalingFactor * (kStar / (omega * rho)) * std::sin(kX * x) * std::sin(kY * y) *
                  std::cos(omega * t) * (std::cosh(kStar * z) + b * std::sinh(kStar * z));
  } else {
    // Elastic-acoustic mode
    const auto pressure = -std::sin(kX * x) * std::sin(kY * y) * std::sin(omega * t) *
                          (std::sin(kStar * z) + b * std::cos(kStar * z));

    setStresses(ScalingFactor * pressure);

    state[uIdx] = ScalingFactor * (kX / (omega * rho)) * std::cos(kX * x) * std::sin(kY * y) *
                  std::cos(omega * t) * (std::sin(kStar * z) + b * std::cos(kStar * z));
    state[vIdx] = ScalingFactor * (kY / (omega * rho)) * std::sin(kX * x) * std::cos(kY * y) *
                  std::cos(omega * t) * (std::sin(kStar * z) + b * std::cos(kStar * z));
    state[wIdx] = ScalingFactor * (kStar / (omega * rho)) * std::sin(kX * x) * std::sin(kY * y) *
                  std::cos(omega * t) * (std::cos(kStar * z) - b * std::sin(kStar * z));
  }
  return state;
}

inline std::array<double, 9> itmReference(double rho0_,
                                          double c0_,
                                          double k_,
                                          double tITMMinus_,
                                          double tau_,
                                          double n_,
                                          double time,
                                          const std::array<double, 3>& point) {
  std::array<double, 9> state{};
  const std::array<double, 3>* points = &point;
  const std::size_t i = 0;
  const double tITMPlus_ = tITMMinus_ + tau_;
  double pressure = 0.0;
  const auto& coordinates = points[i];
  const auto x = coordinates[0];
  const auto t = time;
  if (t <= tITMMinus_) {
    pressure = c0_ * rho0_ * std::cos(k_ * x - c0_ * k_ * t);
    state[0] = -pressure;                       // sigma_xx
    state[1] = -pressure;                       // sigma_yy
    state[2] = -pressure;                       // sigma_zz
    state[3] = 0.0;                             // sigma_xy
    state[4] = 0.0;                             // sigma_yz
    state[5] = 0.0;                             // sigma_xz
    state[6] = std::cos(k_ * x - c0_ * k_ * t); // u
    state[7] = 0.0;                             // v
    state[8] = 0.0;                             // w
  } else if (t <= tITMPlus_) {
    pressure = -0.5 * (n_ - 1) * c0_ * rho0_ *
                   std::cos(k_ * x + c0_ * k_ * n_ * t - (c0_ * k_ * n_ + c0_ * k_) * tITMMinus_) +
               0.5 * (n_ + 1) * c0_ * rho0_ *
                   std::cos(k_ * x - c0_ * k_ * n_ * t + (c0_ * k_ * n_ - c0_ * k_) * tITMMinus_);
    state[0] = -pressure; // sigma_xx
    state[1] = -pressure; // sigma_yy
    state[2] = -pressure; // sigma_zz
    state[3] = 0.0;       // sigma_xy
    state[4] = 0.0;       // sigma_yz
    state[5] = 0.0;       // sigma_xz
    state[6] =
        0.5 * (n_ - 1) *
            std::cos(k_ * x + c0_ * k_ * n_ * t - (c0_ * k_ * n_ + c0_ * k_) * tITMMinus_) +
        0.5 * (n_ + 1) *
            std::cos(k_ * x - c0_ * k_ * n_ * t + (c0_ * k_ * n_ - c0_ * k_) * tITMMinus_); // u
    state[7] = 0.0;                                                                         // v
    state[8] = 0.0;                                                                         // w
  } else {
    pressure = -0.25 * (1 / n_) * c0_ * rho0_ *
               ((-n_ * n_ + 1) * std::cos(k_ * x + c0_ * k_ * t - 2.0 * c0_ * k_ * tITMMinus_ -
                                          (c0_ * k_ * n_ + c0_ * k_) * tau_) +
                (n_ * n_ - 1) * std::cos(k_ * x + c0_ * k_ * t - 2.0 * c0_ * k_ * tITMMinus_ +
                                         (c0_ * k_ * n_ - c0_ * k_) * tau_) +
                (n_ * n_ - 2 * n_ + 1) *
                    std::cos(k_ * x - c0_ * k_ * t + (c0_ * k_ * n_ + c0_ * k_) * tau_) +
                (-n_ * n_ - 2 * n_ - 1) *
                    std::cos(k_ * x - c0_ * k_ * t - (c0_ * k_ * n_ - c0_ * k_) * tau_));
    state[0] = -pressure; // sigma_xx
    state[1] = -pressure; // sigma_yy
    state[2] = -pressure; // sigma_zz
    state[3] = 0.0;       // sigma_xy
    state[4] = 0.0;       // sigma_yz
    state[5] = 0.0;       // sigma_xz
    state[6] = (-0.25 / n_) *
               ((n_ * n_ - 1) * std::cos(k_ * x + c0_ * k_ * t - 2 * c0_ * k_ * tITMMinus_ -
                                         (c0_ * k_ * n_ + c0_ * k_) * tau_) +
                (-n_ * n_ + 1) * std::cos(k_ * x + c0_ * k_ * t - 2 * c0_ * k_ * tITMMinus_ +
                                          (c0_ * k_ * n_ - c0_ * k_) * tau_) +
                (n_ * n_ - 2 * n_ + 1) *
                    std::cos(k_ * x - c0_ * k_ * t + (c0_ * k_ * n_ + c0_ * k_) * tau_) +
                (-n_ * n_ - 2 * n_ - 1) *
                    std::cos(k_ * x - c0_ * k_ * t - (c0_ * k_ * n_ - c0_ * k_) * tau_)); // u
    state[7] = 0.0;                                                                       // v
    state[8] = 0.0;                                                                       // w
  }
  return state;
}
// NOLINTEND(readability-identifier-naming,readability-function-cognitive-complexity)

/// The state a field gives at `point` and `time` in the material `materialData`.
inline std::array<double, 9> stateOf(const physics::InitialField& field,
                                     const CellMaterialData& materialData,
                                     double time,
                                     const std::array<double, 3>& point) {
  std::array<double, 9> state{};
  state.fill(-1.0);
  yateto::DenseTensorView<2, double, unsigned> view(state.data(), {1, 9});
  field.evaluate(time, &point, 1, materialData, view);
  return state;
}

/// Points with x and y in [-1.5, 1.5] and z in [zMin, zMax].
inline std::vector<std::array<double, 3>>
    randomPoints(std::size_t count, unsigned seed, double zMin, double zMax) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> value(-1.5, 1.5);
  std::uniform_real_distribution<double> depth(zMin, zMax);
  std::vector<std::array<double, 3>> points(count);
  for (auto& point : points) {
    point = {value(rng), value(rng), depth(rng)};
  }
  return points;
}

/// The field against `reference(time, point)` at random points at the times `times`, with z in
/// [zMin, zMax].
template <typename ReferenceT>
void checkAgainst(const physics::InitialField& field,
                  const CellMaterialData& materialData,
                  const std::vector<double>& times,
                  const ReferenceT& reference,
                  double zMin = -1.5,
                  double zMax = 1.5) {
  const auto points = randomPoints(200, 11, zMin, zMax);
  for (const double time : times) {
    for (const auto& point : points) {
      const auto values = stateOf(field, materialData, time, point);
      const auto expected = reference(time, point);
      for (std::size_t j = 0; j < 9; ++j) {
        CAPTURE(time);
        CAPTURE(j);
        CHECK(values[j] == doctest::Approx(expected[j]).epsilon(1e-13).scale(1.0));
      }
    }
  }
}

TEST_CASE_TEMPLATE_DEFINE("Scenarios: the scripts SeisSol writes are the closed forms" *
                              doctest::test_suite("physics"),
                          Cfg,
                          ScriptsAreClosedForms) {
  const auto config = configIdOf<Cfg>();
  model::ElasticMaterial solid({2.7, 1.3, 2.1});
  model::ElasticMaterial fluid({1.0, 0.0, 2.25});
  CellMaterialData solidData;
  solidData.local = &solid;
  CellMaterialData fluidData;
  fluidData.local = &fluid;
  const std::vector<double> times{0.0, 0.37, 1.9};

  SUBCASE("Scholte wave, in both materials") {
    const auto fields = physics::scholteWave(config);
    REQUIRE(fields.size() == 1);
    // the elastic material below the interface z = 0, the acoustic one above, where the wave
    // decays away from it
    for (const auto* data : {&solidData, &fluidData}) {
      const bool below = data == &solidData;
      checkAgainst(
          *fields.front(),
          *data,
          times,
          [&](double time, const auto& point) { return scholteReference(time, point, *data); },
          below ? -1.5 : 0.0,
          below ? 0.0 : 1.5);
    }
  }

  SUBCASE("Snell's law, in both materials") {
    const auto fields = physics::snellsLaw(config);
    for (const auto* data : {&solidData, &fluidData}) {
      checkAgainst(*fields.front(), *data, times, [&](double time, const auto& point) {
        return snellReference(time, point, *data);
      });
    }
  }

  SUBCASE("ocean, in each mode") {
    model::ElasticMaterial denser({1.7, 0.0, 2.25});
    CellMaterialData denserData;
    denserData.local = &denser;
    for (const int mode : {0, 1, 2}) {
      CAPTURE(mode);
      const auto fields = physics::ocean(mode, 9.81e-3, config);
      for (const auto* data : {&fluidData, &denserData}) {
        checkAgainst(*fields.front(), *data, times, [&](double time, const auto& point) {
          return oceanReference(mode, 9.81e-3, time, point, *data);
        });
      }
    }
  }

  SUBCASE("acoustic travelling wave with ITM, before, during and after the mirror") {
    AcousticTravellingWaveParametersITM parameters;
    parameters.k = 2.0 * M_PI;
    parameters.itmStartingTime = 0.3;
    parameters.itmDuration = 0.25;
    parameters.itmVelocityScalingFactor = 1.7;
    const auto fields = physics::acousticTravellingWaveITM(fluidData, parameters, config);
    const double c0 = std::sqrt(fluid.getLambdaBar() / fluid.getDensity());
    checkAgainst(*fields.front(),
                 fluidData,
                 {0.0, 0.2, 0.3, 0.41, 0.55, 0.9},
                 [&](double time, const auto& point) {
                   return itmReference(fluid.getDensity(),
                                       c0,
                                       parameters.k,
                                       parameters.itmStartingTime,
                                       parameters.itmDuration,
                                       parameters.itmVelocityScalingFactor,
                                       time,
                                       point);
                 });
  }

  SUBCASE("zero") {
    const auto fields = physics::zeroField(config);
    checkAgainst(*fields.front(), solidData, times, [](double, const auto&) {
      return std::array<double, 9>{};
    });
  }
}

/// The configurations that hold a fluid and whose material is set up from rho, lambda and mu
/// alone, i.e. those the fluid scenarios are written for.
template <typename Cfg>
struct HoldsFluid : std::bool_constant<Cfg::MaterialType == model::MaterialType::Elastic ||
                                       Cfg::MaterialType == model::MaterialType::Acoustic> {};

/// The configurations whose material is set up from rho, lambda and mu alone. The attenuating ones
/// then carry no relaxation (omega = theta = 0), which leaves their memory variables decoupled.
template <typename Cfg>
struct SetUpByLame : std::bool_constant<HoldsFluid<Cfg>::value ||
                                        Cfg::MaterialType == model::MaterialType::Viscoelastic ||
                                        Cfg::MaterialType == model::MaterialType::Viscoacoustic> {};

/// Fills every entry before a scenario is evaluated.
constexpr double Sentinel = -1234.5;

/// The values of `pointCount` points of the configuration `Cfg`, followed by a guard as large as
/// the elastic quantities (nine) would need in addition. Nothing may be written into the guard.
template <typename Cfg>
std::vector<Real<Cfg>> guardedBuffer(std::size_t pointCount) {
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  return std::vector<Real<Cfg>>(pointCount * (NumQuantities + 9), Sentinel);
}

template <typename Cfg>
yateto::DenseTensorView<2, Real<Cfg>, unsigned> dofsView(std::vector<Real<Cfg>>& buffer,
                                                         std::size_t pointCount) {
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  return yateto::DenseTensorView<2, Real<Cfg>, unsigned>(
      buffer.data(), {static_cast<unsigned>(pointCount), static_cast<unsigned>(NumQuantities)});
}

template <typename Cfg>
void checkGuard(const std::vector<Real<Cfg>>& buffer, std::size_t pointCount) {
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  for (std::size_t i = pointCount * NumQuantities; i < buffer.size(); ++i) {
    CAPTURE(i);
    CHECK(buffer[i] == Sentinel);
  }
}

/// Checks the state of a fluid at `point`: `stress` (i.e. -p) on the normal stresses, no shear
/// stress, and `velocity` in the velocity of the configuration.
template <typename Cfg>
void checkFluidState(const yateto::DenseTensorView<2, Real<Cfg>, unsigned>& dofs,
                     std::size_t point,
                     double stress,
                     const std::array<double, 3>& velocity,
                     double tolerance) {
  using MaterialT = model::MaterialOf<Cfg>;
  CAPTURE(point);
  // a stress tensor has three normal components; the acoustic materials carry a single one
  constexpr std::size_t NormalComponents = std::min<std::size_t>(3, MaterialT::TractionComponents);
  for (std::size_t j = 0; j < MaterialT::TractionComponents; ++j) {
    CAPTURE(j);
    CHECK(dofs(point, j) == AbsApprox(j < NormalComponents ? stress : 0.0).epsilon(tolerance));
  }
  for (std::size_t d = 0; d < velocity.size(); ++d) {
    CAPTURE(d);
    CHECK(dofs(point, MaterialT::VelocityOffset + d) == AbsApprox(velocity[d]).epsilon(tolerance));
  }
}

template <typename Cfg>
model::MaterialOf<Cfg> fluidOf(double rho, double lambda) {
  model::MaterialOf<Cfg> material;
  material.rho = rho;
  material.setLameParameters(0.0, lambda);
  return material;
}

TEST_CASE_TEMPLATE_DEFINE(
    "AcousticTravellingWaveITM writes a plane wave into the pressure and velocity of the "
    "configuration" *
        doctest::test_suite("physics"),
    Cfg,
    TravellingWaveItmPlaneWave) {
  constexpr double Tolerance = 100 * std::numeric_limits<Real<Cfg>>::epsilon();

  constexpr double Rho = 2.0;
  constexpr double Lambda = 8.0;
  constexpr double K = 3.0;
  const double c = std::sqrt(Lambda / Rho);

  auto material = fluidOf<Cfg>(Rho, Lambda);
  const CellMaterialData materialData{&material, {}};

  // With a factor of one, the mirror leaves the wave as it is. So before, during and after it, the
  // plain travelling wave p = rho c cos(k x - c k t), u = cos(k x - c k t) has to come out.
  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = K;
  parameters.itmStartingTime = 1.0;
  parameters.itmDuration = 0.5;
  parameters.itmVelocityScalingFactor = 1.0;
  const auto fields =
      physics::acousticTravellingWaveITM(materialData, parameters, configIdOf<Cfg>());
  REQUIRE(fields.size() == 1);
  const auto& wave = *fields.front();

  const std::array<std::array<double, 3>, 3> points{
      {{0.1, 0.2, 0.3}, {0.45, -0.7, 1.1}, {0.9, 0.5, -0.4}}};
  for (const double time : {0.5, 1.25, 2.0}) {
    CAPTURE(time);
    auto buffer = guardedBuffer<Cfg>(points.size());
    auto dofs = dofsView<Cfg>(buffer, points.size());
    wave.evaluate(time, points.data(), points.size(), materialData, dofs);
    for (std::size_t i = 0; i < points.size(); ++i) {
      const double phase = K * points[i][0] - c * K * time;
      checkFluidState<Cfg>(
          dofs, i, -Rho * c * std::cos(phase), {std::cos(phase), 0.0, 0.0}, Tolerance);
    }
    checkGuard<Cfg>(buffer, points.size());
  }
}

TEST_CASE_TEMPLATE_DEFINE(
    "AcousticTravellingWaveITM keeps the pressure and scales the velocity at the mirror" *
        doctest::test_suite("physics"),
    Cfg,
    TravellingWaveItmMirror) {
  constexpr auto VelocityOffset = model::MaterialOf<Cfg>::VelocityOffset;

  constexpr double Start = 1.0;
  constexpr double Duration = 0.5;
  constexpr double N = 2.0;
  // the values right after a switch are taken one double later than the switch itself
  constexpr double SwitchTolerance =
      std::max(1.0e-10, 100 * static_cast<double>(std::numeric_limits<Real<Cfg>>::epsilon()));

  auto material = fluidOf<Cfg>(2.0, 8.0);
  const CellMaterialData materialData{&material, {}};

  AcousticTravellingWaveParametersITM parameters{};
  parameters.k = 3.0;
  parameters.itmStartingTime = Start;
  parameters.itmDuration = Duration;
  parameters.itmVelocityScalingFactor = N;
  const auto fields =
      physics::acousticTravellingWaveITM(materialData, parameters, configIdOf<Cfg>());
  REQUIRE(fields.size() == 1);
  const auto& wave = *fields.front();

  const std::array<std::array<double, 3>, 2> points{{{0.1, 0.2, 0.3}, {0.9, 0.5, -0.4}}};
  const auto evaluateAt = [&](double time) {
    auto buffer = guardedBuffer<Cfg>(points.size());
    auto dofs = dofsView<Cfg>(buffer, points.size());
    wave.evaluate(time, points.data(), points.size(), materialData, dofs);
    checkGuard<Cfg>(buffer, points.size());
    return buffer;
  };
  // The mirror scales the velocity by the factor n when it starts, and by 1 / n when it ends; the
  // pressure passes both switches continuously.
  const auto checkSwitch = [&](double switchTime, double velocityFactor) {
    auto before = evaluateAt(switchTime);
    auto after = evaluateAt(std::nextafter(switchTime, std::numeric_limits<double>::infinity()));
    const auto dofsBefore = dofsView<Cfg>(before, points.size());
    const auto dofsAfter = dofsView<Cfg>(after, points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
      CAPTURE(i);
      CHECK(dofsAfter(i, 0) == AbsApprox(dofsBefore(i, 0)).epsilon(SwitchTolerance));
      CHECK(dofsAfter(i, VelocityOffset) ==
            AbsApprox(velocityFactor * dofsBefore(i, VelocityOffset)).epsilon(SwitchTolerance));
    }
  };

  SUBCASE("start") { checkSwitch(Start, N); }
  SUBCASE("end") { checkSwitch(Start + Duration, 1.0 / N); }
}

TEST_CASE_TEMPLATE_DEFINE("Ocean writes a fluid stress and the velocity of the configuration" *
                              doctest::test_suite("physics"),
                          Cfg,
                          OceanLayout) {
  constexpr double Tolerance = 100 * std::numeric_limits<Real<Cfg>>::epsilon();

  // the medium the ocean modes are computed for: c = 1.5, g = 9.81e-3
  constexpr double Rho = 1.0;
  constexpr double G = 9.81e-3;
  constexpr double Time = 0.3;

  auto material = fluidOf<Cfg>(Rho, 2.25);
  const CellMaterialData materialData{&material, {}};

  const auto fields = physics::ocean(1, G, configIdOf<Cfg>());
  REQUIRE(fields.size() == 1);
  const auto& ocean = *fields.front();

  // no component vanishes here: x and y are no multiples of 5, kStar z is not near a multiple of
  // pi / 2
  const std::array<std::array<double, 3>, 3> points{
      {{1.3, 2.7, -0.4}, {3.1, 7.9, -1.2}, {8.2, 4.4, -2.5}}};
  auto buffer = guardedBuffer<Cfg>(points.size());
  auto dofs = dofsView<Cfg>(buffer, points.size());
  ocean.evaluate(Time, points.data(), points.size(), materialData, dofs);

  // mode 1 of the scenario
  const double kX = std::acos(-1.0) / 10.0;
  const double kY = kX;
  constexpr double KStar = 1.5733628061766445;
  constexpr double Omega = 2.4523337594491745;
  constexpr double B = G * KStar / (Omega * Omega);
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto [x, y, z] = points[i];
    const double depth = std::sin(KStar * z) + B * std::cos(KStar * z);
    const double stress = -std::sin(kX * x) * std::sin(kY * y) * std::sin(Omega * Time) * depth;
    const std::array<double, 3> velocity{
        kX / (Omega * Rho) * std::cos(kX * x) * std::sin(kY * y) * std::cos(Omega * Time) * depth,
        kY / (Omega * Rho) * std::sin(kX * x) * std::cos(kY * y) * std::cos(Omega * Time) * depth,
        KStar / (Omega * Rho) * std::sin(kX * x) * std::sin(kY * y) * std::cos(Omega * Time) *
            (std::cos(KStar * z) - B * std::sin(KStar * z))};
    checkFluidState<Cfg>(dofs, i, stress, velocity, Tolerance);
  }
  checkGuard<Cfg>(buffer, points.size());
}

TEST_CASE_TEMPLATE_DEFINE("SuperimposedPlanarwave is the sum of its three planar waves" *
                              doctest::test_suite("physics"),
                          Cfg,
                          SuperimposedPlanarwaveSum) {
  using real = Real<Cfg>;
  constexpr auto NumQuantities = model::MaterialOf<Cfg>::NumQuantities;
  constexpr double Tolerance = 100 * std::numeric_limits<real>::epsilon();

  model::MaterialOf<Cfg> material;
  material.rho = 1.0;
  material.setLameParameters(1.0, 2.0);
  const CellMaterialData materialData{&material, {}};

  // as many points as the initial field projection evaluates, i.e. more than basis functions
  constexpr std::size_t PointCount =
      (Cfg::ConvergenceOrder + 1) * (Cfg::ConvergenceOrder + 1) * (Cfg::ConvergenceOrder + 1);
  std::vector<std::array<double, 3>> points(PointCount);
  for (std::size_t i = 0; i < PointCount; ++i) {
    const auto s = static_cast<double>(i);
    points[i] = {std::fmod(0.37 * s, 1.0), std::fmod(0.61 * s, 1.0), std::fmod(0.83 * s, 1.0)};
  }

  constexpr double Time = 0.25;

  // the first simulation of each, which all have the phase zero
  std::vector<real> actual(PointCount * NumQuantities);
  auto actualDofs = dofsView<Cfg>(actual, PointCount);
  const auto superimposed = physics::superimposedPlanarwave(materialData, configIdOf<Cfg>());
  superimposed.front()->evaluate(Time, points.data(), PointCount, materialData, actualDofs);

  std::vector<real> expected(PointCount * NumQuantities, 0);
  std::vector<real> single(PointCount * NumQuantities);
  auto expectedDofs = dofsView<Cfg>(expected, PointCount);
  auto singleDofs = dofsView<Cfg>(single, PointCount);
  const std::array<Eigen::Vector3d, 3> kVecs{Eigen::Vector3d(M_PI, 0.0, 0.0),
                                             Eigen::Vector3d(0.0, M_PI, 0.0),
                                             Eigen::Vector3d(0.0, 0.0, M_PI)};
  for (const auto& kVec : kVecs) {
    const auto planarwave = physics::planarwave(materialData, configIdOf<Cfg>(), kVec);
    planarwave.front()->evaluate(Time, points.data(), PointCount, materialData, singleDofs);
    for (std::size_t j = 0; j < NumQuantities; ++j) {
      for (std::size_t i = 0; i < PointCount; ++i) {
        expectedDofs(i, j) += singleDofs(i, j);
      }
    }
  }

  std::size_t mismatches = 0;
  for (std::size_t i = 0; i < actual.size(); ++i) {
    if (!std::isfinite(actual[i]) || std::abs(actual[i] - expected[i]) > Tolerance) {
      ++mismatches;
    }
  }
  CHECK(mismatches == 0);
}

TEST_CASE_TEMPLATE("Scenarios are offered only for materials that can represent them" *
                       doctest::test_suite("physics"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using initializer::parameters::InitializationType;
  const auto available = [](InitializationType type) {
    return physics::scenario::availability(type, configIdOf<Cfg>()).available;
  };

  constexpr auto Type = Cfg::MaterialType;
  constexpr bool Elastic = Type == model::MaterialType::Elastic;
  constexpr bool Fluid = Elastic || Type == model::MaterialType::Acoustic;

  // built from the material of the configuration itself
  CHECK(available(InitializationType::Zero));
  CHECK(available(InitializationType::Planarwave));
  CHECK(available(InitializationType::SuperimposedPlanarwave));
  CHECK(available(InitializationType::Script));
  CHECK(available(InitializationType::Travelling) == (Cfg::RelaxationMechanisms == 0));

  // a solid next to a fluid
  CHECK(available(InitializationType::Scholte) == Elastic);
  CHECK(available(InitializationType::Snell) == Elastic);

  // a fluid
  CHECK(available(InitializationType::AcousticTravellingWithITM) == Fluid);
  CHECK(available(InitializationType::Ocean0) == Fluid);
  CHECK(available(InitializationType::Ocean1) == Fluid);
  CHECK(available(InitializationType::Ocean2) == Fluid);

  CHECK(available(InitializationType::PressureInjection) ==
        (Type == model::MaterialType::Poroelastic));
}

TEST_CASE("The acoustic travelling wave passes the time mirror only if it is enabled" *
          doctest::test_suite("physics")) {
  using namespace seissol::initializer::parameters;

  // mirror settings as in docs/parameters.par
  YAML::Node node = YAML::Load(R"(
    equations:
      itmenable: 1
      itmstartingtime: 2.0
      itmtime: 0.01
      itmvelocityscalingfactor: 2
    inicondition:
      cictype: AcousticTravellingwithITM
      k: 6.283
  )");
  const auto waveFor = [](const YAML::Node& parameterFile) {
    ParameterReader reader(parameterFile, "", false);
    const ITMParameters itm = readITMParameters(&reader);
    // the configuration only gives the default of the amplitudes, which the wave does not read
    const InitializationParameters initialization = readInitializationParameters(&reader, 0);
    return physics::scenario::acousticTravellingWaveParameters(initialization, itm);
  };

  SUBCASE("enabled") {
    const AcousticTravellingWaveParametersITM wave = waveFor(node);
    CHECK(wave.k == doctest::Approx(6.283));
    CHECK(wave.itmStartingTime == doctest::Approx(2.0));
    CHECK(wave.itmDuration == doctest::Approx(0.01));
    CHECK(wave.itmVelocityScalingFactor == doctest::Approx(2.0));
  }

  SUBCASE("disabled") {
    // The file keeps its mirror settings, but they are ignored now. The reference solution must
    // not model a mirror the run does not have.
    node["equations"]["itmenable"] = 0;
    const AcousticTravellingWaveParametersITM wave = waveFor(node);
    CHECK(wave.k == doctest::Approx(6.283));
    CHECK(wave.itmVelocityScalingFactor == doctest::Approx(1.0));
  }
}

TEST_CASE_TEMPLATE_APPLY(TravellingWaveItmPlaneWave, ConfigsWhere<HoldsFluid>);
TEST_CASE_TEMPLATE_APPLY(TravellingWaveItmMirror, ConfigsWhere<HoldsFluid>);
TEST_CASE_TEMPLATE_APPLY(OceanLayout, ConfigsWhere<HoldsFluid>);
TEST_CASE_TEMPLATE_APPLY(SuperimposedPlanarwaveSum, ConfigsWhere<SetUpByLame>);
TEST_CASE_TEMPLATE_APPLY(ScriptsAreClosedForms, ConfigsOfMaterial<model::MaterialType::Elastic>);

} // namespace seissol::unit_test::scenarios
