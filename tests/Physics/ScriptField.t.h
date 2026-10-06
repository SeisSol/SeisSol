// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/elastic/Model/Datastructures.h"
#include "Initializer/Typedefs.h"
#include "Physics/ScriptField.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

namespace seissol::unit_test::script_field {

namespace {

/// A script written to a file of its own, removed again when done.
class ScriptFile {
  public:
  ScriptFile(const std::string& name, const std::string& content)
      : path_(std::filesystem::temp_directory_path() /
              ("seissol-scriptfield-" + std::to_string(std::random_device{}()) + "-" + name)) {
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

const std::vector<std::string> Quantities{
    "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

std::vector<std::array<double, 3>> randomPoints(std::size_t count, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> value(-2.0, 2.0);
  std::vector<std::array<double, 3>> points(count);
  for (auto& point : points) {
    point = {value(rng), value(rng), value(rng)};
  }
  return points;
}

} // namespace

TEST_CASE("ScriptField: a compiled script gives the quantities at any points and times" *
          doctest::test_suite("physics")) {
  const ScriptFile script("field.sderiv",
                          "out def v1 = sin(x) * cos(t)\n"
                          "out def s_xx = rho * x + mu * y + lambda * z + sim\n");
  const physics::ScriptField field("sderiv:" + script.path(), Quantities, 2, true);
  CHECK(field.compiled());

  model::ElasticMaterial material({2.5, 1.5, 3.0});
  CellMaterialData materialData;
  materialData.local = &material;

  // more points than a kernel is bound for, so the call runs in pieces
  const auto points = randomPoints(2500, 1);
  const double time = 0.75;
  std::vector<double> values(points.size() * Quantities.size(), -1.0);
  field.evaluateValues(time, points.data(), points.size(), materialData, values.data());

  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto& [x, y, z] = points[i];
    CHECK(values[6 * points.size() + i] ==
          doctest::Approx(std::sin(x) * std::cos(time)).epsilon(1e-14));
    CHECK(values[0 * points.size() + i] ==
          doctest::Approx(2.5 * x + 1.5 * y + 3.0 * z + 2.0).epsilon(1e-14));
    // what the script does not give is zero
    CHECK(values[7 * points.size() + i] == 0.0);
  }
}

TEST_CASE("ScriptField: many threads evaluate at once as one does alone" *
          doctest::test_suite("physics")) {
  const ScriptFile script("threads.sderiv",
                          "out def v2 = exp(-(x*x + y*y + z*z)) * sin(10.0 * t)\n"
                          "out def s_yz = x * y - z * t\n");
  const physics::ScriptField field("sderiv:" + script.path(), Quantities, 0, true);

  model::ElasticMaterial material({1.0, 1.0, 2.0});
  CellMaterialData materialData;
  materialData.local = &material;

  constexpr std::size_t Faces = 256;
  constexpr std::size_t PointsPerFace = 21;
  const auto points = randomPoints(Faces * PointsPerFace, 2);
  const std::size_t stride = PointsPerFace * Quantities.size();

  std::vector<double> serial(Faces * stride);
  for (std::size_t face = 0; face < Faces; ++face) {
    field.evaluateValues(0.1 * static_cast<double>(face),
                         points.data() + face * PointsPerFace,
                         PointsPerFace,
                         materialData,
                         serial.data() + face * stride);
  }
  std::vector<double> parallel(Faces * stride);
#pragma omp parallel for schedule(dynamic)
  for (std::size_t face = 0; face < Faces; ++face) {
    field.evaluateValues(0.1 * static_cast<double>(face),
                         points.data() + face * PointsPerFace,
                         PointsPerFace,
                         materialData,
                         parallel.data() + face * stride);
  }
  CHECK(std::memcmp(serial.data(), parallel.data(), serial.size() * sizeof(double)) == 0);
}

TEST_CASE("ScriptField: a script for the initial condition only is accepted until evaluated" *
          doctest::test_suite("physics")) {
  // the mesh group is known to the initial condition, not at a boundary node
  const ScriptFile script("group.sderiv", "out def v1 = group * x\n");
  CHECK_NOTHROW(physics::ScriptField("sderiv:" + script.path(), Quantities, 0, true));
}

TEST_CASE("ScriptField: an easi file is evaluated through its reader" *
          doctest::test_suite("physics")) {
  const ScriptFile script("field.yaml",
                          "!AffineMap\n"
                          "matrix:\n"
                          "  v3: [1.0, 2.0, -0.5]\n"
                          "translation:\n"
                          "  v3: 0.25\n");
  const physics::ScriptField field(script.path(), Quantities, 0, false);
  CHECK_FALSE(field.compiled());

  model::ElasticMaterial material({1.0, 1.0, 2.0});
  CellMaterialData materialData;
  materialData.local = &material;
  const auto points = randomPoints(10, 3);
  std::vector<double> values(points.size() * Quantities.size(), -1.0);
  field.evaluateValues(0.0, points.data(), points.size(), materialData, values.data());
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto& [x, y, z] = points[i];
    CHECK(values[8 * points.size() + i] ==
          doctest::Approx(x + 2.0 * y - 0.5 * z + 0.25).epsilon(1e-14));
    CHECK(values[0 * points.size() + i] == 0.0);
  }
}

} // namespace seissol::unit_test::script_field
