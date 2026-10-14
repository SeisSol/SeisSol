// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_MATERIALSAMPLING_T_H_
#define SEISSOL_TESTS_MODEL_MATERIALSAMPLING_T_H_

// Where a material varies inside a cell it is carried as samples, and some
// readers take it from there at points of their own: the points of the plastic
// strain and of the volume quadrature. Each of those is one matrix that takes
// the samples straight to the points, folded from the projection to the modal
// basis and the evaluation at the points.

#include <doctest.h>

#include "Common/Constants.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Numerical/Quadrature.h"

#include <array>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("Material at the points of the plastic strain and of the volume quadrature") {
  constexpr std::size_t Samples = tensor::materialNodes::Shape[0];
  constexpr double Tolerance = std::is_same_v<real, double> ? 1e-11 : 1e-5;

  // A field of degree two, or one where the basis does not reach that, which
  // every sample set carries exactly: sampled where the material is, it has to
  // arrive at the points it is read at as itself, in the order of those points.
  constexpr double Quadratic = ConvergenceOrder > 2 ? 1.0 : 0.0;
  const auto field = [&](const double* point) {
    return 1.0 + 0.3 * point[0] - 0.2 * point[1] + 0.5 * point[2] +
           Quadratic * (0.4 * point[0] * point[2] - 0.1 * point[1] * point[1]);
  };

  // A row that is zero, such as a vertex at the origin, may be left out of the
  // storage of a point set, so every coordinate is read through its range.
  const auto pointOf = [](const auto& view, std::size_t row) {
    std::array<double, 3> point{};
    for (std::size_t j = 0; j < 3; ++j) {
      if (view.isInRange(row, j)) {
        point[j] = view(row, j);
      }
    }
    return point;
  };

  const auto nodes = init::materialNodes::view::create(init::materialNodes::Values);
  std::vector<double> sampled(Samples);
  for (std::size_t sample = 0; sample < Samples; ++sample) {
    sampled[sample] = field(pointOf(nodes, sample).data());
  }

  const auto check = [&](const auto& interpolation, std::size_t points, const auto& pointAt) {
    for (std::size_t point = 0; point < points; ++point) {
      double value = 0.0;
      for (std::size_t sample = 0; sample < Samples; ++sample) {
        if (interpolation.isInRange(point, sample)) {
          value += interpolation(point, sample) * sampled[sample];
        }
      }
      const auto coordinates = pointAt(point);
      REQUIRE(value == doctest::Approx(field(coordinates.data())).epsilon(Tolerance));
    }
  };

  SUBCASE("where the plastic strain lives") {
    constexpr std::size_t Points = tensor::vNodes::Shape[0];
    static_assert(tensor::materialToPlasticity::Shape[0] == Points);
    static_assert(tensor::materialToPlasticity::Shape[1] == Samples);
    const auto plasticity = init::vNodes::view::create(init::vNodes::Values);
    check(init::materialToPlasticity::view::create(init::materialToPlasticity::Values),
          Points,
          [&](std::size_t point) { return pointOf(plasticity, point); });
  }

  SUBCASE("where the energies are integrated") {
    constexpr std::size_t PerDirection = ConvergenceOrder + 1;
    constexpr std::size_t Points = PerDirection * PerDirection * PerDirection;
    static_assert(tensor::materialToQuadrature::Shape[0] == Points);
    static_assert(tensor::materialToQuadrature::Shape[1] == Samples);
    std::array<std::array<double, 3>, Points> points{};
    double weights[Points]{};
    seissol::quadrature::TetrahedronQuadrature(points.data(), weights, PerDirection);
    check(init::materialToQuadrature::view::create(init::materialToQuadrature::Values),
          Points,
          [&](std::size_t point) {
            return std::array<double, 3>{points[point][0], points[point][1], points[point][2]};
          });
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_MATERIALSAMPLING_T_H_
