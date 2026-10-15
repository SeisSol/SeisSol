// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigRegistry.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "Initializer/Typedefs.h"
#include "Model/MaterialType.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/Scenarios.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <random>
#include <vector>

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

TEST_CASE("Scenarios: the scripts SeisSol writes are the closed forms" *
          doctest::test_suite("physics")) {
  if constexpr (Config::MaterialType == model::MaterialType::Elastic) {
    const auto config = configIdOf<Config>();
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
}

} // namespace seissol::unit_test::scenarios
