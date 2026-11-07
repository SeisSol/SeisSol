// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/ConfigRegistry.h"
#include "Equations/Datastructures.h"
#include "Initializer/Typedefs.h"
#include "Model/MaterialType.h"
#include "Model/StateNames.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/PlaneWave.h"
#include "TestConfigs.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <vector>

namespace seissol::unit_test::plane_wave {

/// The state of the modes at a point and time as the plane waves computed it before they became
/// a script: in complex arithmetic, a mode cut to one period about the origin if there is one.
inline std::vector<double>
    complexPlaneWave(const std::vector<std::vector<physics::PlaneWaveMode>>& waves,
                     double phase,
                     const std::optional<Eigen::Vector3d>& origin,
                     double time,
                     const std::array<double, 3>& point) {
  const auto count = waves.front().front().coefficients.size();
  std::vector<double> state(count, 0.0);
  const Eigen::Vector3d shift = origin.value_or(Eigen::Vector3d::Zero());
  for (const auto& wave : waves) {
    for (const auto& mode : wave) {
      const auto arg =
          std::complex<double>(0.0, 1.0) *
          (mode.omega * time - mode.k[0] * (point[0] - shift[0]) -
           mode.k[1] * (point[1] - shift[1]) - mode.k[2] * (point[2] - shift[2]) + phase);
      if (origin.has_value() && (arg.imag() <= -0.5 * M_PI || arg.imag() >= 1.5 * M_PI)) {
        continue;
      }
      for (std::size_t j = 0; j < count; ++j) {
        state[j] += (mode.coefficients[j] * std::exp(arg)).real();
      }
    }
  }
  return state;
}

/// The components a field gives at `points` at `time`, component j of point i at [j][i].
inline std::vector<std::vector<double>> evaluate(const physics::InitialField& field,
                                                 const CellMaterialData& materialData,
                                                 double time,
                                                 std::vector<std::array<double, 3>> points,
                                                 std::size_t count) {
  std::vector<double> data(points.size() * count, -1.0);
  yateto::DenseTensorView<2, double, unsigned> view(
      data.data(), {static_cast<unsigned>(points.size()), static_cast<unsigned>(count)});
  field.evaluate(time, points.data(), points.size(), materialData, view);
  std::vector<std::vector<double>> components(count, std::vector<double>(points.size()));
  for (std::size_t j = 0; j < count; ++j) {
    for (std::size_t i = 0; i < points.size(); ++i) {
      components[j][i] = view(i, j);
    }
  }
  return components;
}

inline std::vector<std::array<double, 3>> randomPoints(std::size_t count, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> value(-1.5, 1.5);
  std::vector<std::array<double, 3>> points(count);
  for (auto& point : points) {
    point = {value(rng), value(rng), value(rng)};
  }
  return points;
}

/// Checks the field against the complex plane wave at random points and times.
inline void checkAgainstComplex(const physics::InitialField& field,
                                const CellMaterialData& materialData,
                                const std::vector<std::vector<physics::PlaneWaveMode>>& waves,
                                double phase,
                                const std::optional<Eigen::Vector3d>& origin) {
  const auto count = waves.front().front().coefficients.size();
  const auto points = randomPoints(300, 7);
  for (const double time : {0.0, 0.37, 1.9}) {
    const auto values = evaluate(field, materialData, time, points, count);
    for (std::size_t i = 0; i < points.size(); ++i) {
      const auto expected = complexPlaneWave(waves, phase, origin, time, points[i]);
      for (std::size_t j = 0; j < count; ++j) {
        CAPTURE(i);
        CAPTURE(j);
        CHECK(values[j][i] == doctest::Approx(expected[j]).epsilon(1e-13).scale(1.0));
      }
    }
  }
}

TEST_CASE_TEMPLATE_DEFINE("PlaneWave: the generated script is the plane wave of the eigenpairs" *
                              doctest::test_suite("physics"),
                          Cfg,
                          GeneratedPlaneWave) {
  const auto config = configIdOf<Cfg>();
  model::ElasticMaterial material({2.7, 1.3, 2.1});
  CellMaterialData materialData;
  materialData.local = &material;
  const auto names = model::stateNames<model::MaterialOf<Cfg>>();

  SUBCASE("planar wave") {
    const auto [modes, amplitudes] = physics::defaultPlaneWaveModes(materialData, config);
    CHECK(modes == std::vector<int>{1, 8});
    const auto fields = physics::planarwave(materialData, config, Eigen::Vector3d::Zero());
    REQUIRE(fields.size() == Cfg::NumSimulations);
    const std::vector<std::vector<physics::PlaneWaveMode>> waves{physics::planeWaveModes(
        materialData, config, Eigen::Vector3d(M_PI, M_PI, M_PI), modes, amplitudes)};
    for (std::size_t sim = 0; sim < fields.size(); ++sim) {
      const double phase = (2.0 * M_PI * sim) / Cfg::NumSimulations;
      checkAgainstComplex(*fields[sim], materialData, waves, phase, std::nullopt);
    }
  }

  SUBCASE("planar wave along kVec") {
    const auto [modes, amplitudes] = physics::defaultPlaneWaveModes(materialData, config);
    const Eigen::Vector3d kVec(2.0 * M_PI, 0.0, -M_PI);
    const auto fields = physics::planarwave(materialData, config, kVec);
    const std::vector<std::vector<physics::PlaneWaveMode>> waves{
        physics::planeWaveModes(materialData, config, kVec, modes, amplitudes)};
    checkAgainstComplex(*fields.front(), materialData, waves, 0.0, std::nullopt);
  }

  SUBCASE("superimposed planar wave") {
    const auto [modes, amplitudes] = physics::defaultPlaneWaveModes(materialData, config);
    const auto fields = physics::superimposedPlanarwave(materialData, config);
    std::vector<std::vector<physics::PlaneWaveMode>> waves;
    for (const auto& k : {Eigen::Vector3d(M_PI, 0.0, 0.0),
                          Eigen::Vector3d(0.0, M_PI, 0.0),
                          Eigen::Vector3d(0.0, 0.0, M_PI)}) {
      waves.push_back(physics::planeWaveModes(materialData, config, k, modes, amplitudes));
    }
    checkAgainstComplex(*fields.front(), materialData, waves, 0.0, std::nullopt);
  }

  SUBCASE("travelling wave") {
    TravellingWaveParameters parameters;
    parameters.origin = Eigen::Vector3d(0.25, -0.5, 0.1);
    parameters.kVec = Eigen::Vector3d(2.0, -1.0, 0.5);
    parameters.varField = {2, 8};
    parameters.ampField = {0.5, std::complex<double>(1.0, -0.25)};
    const auto fields = physics::travellingWave(materialData, config, parameters);
    REQUIRE(fields.size() == 1);
    const std::vector<std::vector<physics::PlaneWaveMode>> waves{physics::planeWaveModes(
        materialData, config, parameters.kVec, parameters.varField, parameters.ampField)};
    checkAgainstComplex(*fields.front(), materialData, waves, 0.5 * M_PI, parameters.origin);
  }

  SUBCASE("an elastic material without shear carries the P mode alone") {
    model::ElasticMaterial fluid({1.0, 0.0, 2.25});
    CellMaterialData fluidData;
    fluidData.local = &fluid;
    const auto [modes, amplitudes] = physics::defaultPlaneWaveModes(fluidData, config);
    CHECK(modes == std::vector<int>{8});
    const auto fields = physics::planarwave(fluidData, config, Eigen::Vector3d::Zero());
    const std::vector<std::vector<physics::PlaneWaveMode>> waves{physics::planeWaveModes(
        fluidData, config, Eigen::Vector3d(M_PI, M_PI, M_PI), modes, amplitudes)};
    checkAgainstComplex(*fields.front(), fluidData, waves, 0.0, std::nullopt);
  }

  SUBCASE("the module is a readable script") {
    const auto [modes, amplitudes] = physics::defaultPlaneWaveModes(materialData, config);
    physics::PlaneWaveSum sum;
    sum.waves.push_back(physics::planeWaveModes(
        materialData, config, Eigen::Vector3d(M_PI, M_PI, M_PI), modes, amplitudes));
    const auto source = physics::planeWaveModule(sum, names);
    CHECK(source.find("def phase = 0.0") != std::string::npos);
    CHECK(source.find("def th1 = ") != std::string::npos);
    CHECK(source.find("def th2 = ") != std::string::npos);
    CHECK(source.find("out def v1 = ") != std::string::npos);
  }
}

TEST_CASE("PlaneWave: the eigenvectors do not depend on the eigensolver" *
          doctest::test_suite("physics")) {
  // an eigenvalue of its own and a repeated one, with an orthonormal basis of eigenvectors
  const std::vector<std::complex<double>> values{1.0, 2.0, 2.0};
  const Eigen::Matrix3cd basis =
      Eigen::HouseholderQR<Eigen::Matrix3cd>(Eigen::Matrix3cd::Random()).householderQ();
  // another solver may scale each eigenvector by a phase, and choose another basis of the
  // eigenspace of the repeated eigenvalue
  const std::complex<double> i(0.0, 1.0);
  Eigen::Matrix3cd other = basis;
  other.col(0) *= std::exp(0.7 * i);
  const double angle = 0.4;
  other.col(1) = std::exp(-1.1 * i) * (std::cos(angle) * basis.col(1) +
                                       std::sin(angle) * std::exp(0.3 * i) * basis.col(2));
  other.col(2) = std::exp(2.5 * i) * (-std::sin(angle) * std::exp(-0.3 * i) * basis.col(1) +
                                      std::cos(angle) * basis.col(2));

  const auto canonical = [&](const Eigen::Matrix3cd& vectors) {
    std::vector<std::complex<double>> columns(vectors.data(), vectors.data() + vectors.size());
    physics::canonicalizeEigenvectors(values, columns);
    return Eigen::Map<const Eigen::Matrix3cd>(columns.data()).eval();
  };
  const Eigen::Matrix3cd first = canonical(basis);
  const Eigen::Matrix3cd second = canonical(other);
  CHECK((first - second).norm() < 1e-12);
  // still eigenvectors of the same eigenvalues, of unit length
  for (Eigen::Index j = 0; j < 3; ++j) {
    CAPTURE(j);
    CHECK(first.col(j).norm() == doctest::Approx(1.0).epsilon(1e-14));
  }
  CHECK((first.col(0) - basis.col(0) * (basis.col(0).adjoint() * first.col(0))).norm() < 1e-12);
}

TEST_CASE_TEMPLATE_APPLY(GeneratedPlaneWave, ConfigsOfMaterial<model::MaterialType::Elastic>);

} // namespace seissol::unit_test::plane_wave
