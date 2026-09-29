// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_MATERIALSAMPLING_T_H_
#define SEISSOL_TESTS_MODEL_MATERIALSAMPLING_T_H_

// Where a material varies inside a cell it is carried as samples, and whoever
// wants it somewhere else reads it there: the nodes of a face for the flux, the
// quadrature points of the dynamic rupture rule for a fault, the points of the
// plastic strain and of the volume quadrature. Each of those is one matrix that
// takes the samples straight to the points, folded from the projection to the
// modal basis and the evaluation at the points.

#include <doctest.h>

#include "Alignment.h"
#include "Common/Constants.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Numerical/Quadrature.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("Material at the points of a fault") {
  constexpr std::size_t Samples = tensor::materialNodes::Shape[0];
  constexpr std::size_t FaultPoints =
      tensor::materialAtFault::Shape[multisim::BasisFunctionDimension];
  constexpr std::size_t Modes = tensor::materialProject::Shape[0];
  // one matrix per side and face relation of a dynamic rupture face
  constexpr std::size_t Groups = Cell::NumFaces * dr::misc::NumFaceRelations;
  // what the folds may be off by, which follows the precision they are stored in
  constexpr double Tolerance = std::is_same_v<real, double> ? 1e-12 : 1e-5;

  // The folds below are stated in the orientation the mathematics has. A build
  // that bundles simulations stores whatever comes from a matrix file the other
  // way round, so those reads say where they come from.
  constexpr bool FileTransposed = multisim::NumSimulations > 1;
  const auto denseView =
      [](auto view, std::size_t rows, std::size_t columns, bool fromFile = false) {
        const bool flip = fromFile && FileTransposed;
        std::vector<double> dense(rows * columns, 0.0);
        for (std::size_t row = 0; row < rows; ++row) {
          for (std::size_t column = 0; column < columns; ++column) {
            const std::size_t first = flip ? column : row;
            const std::size_t second = flip ? row : column;
            if (view.isInRange(first, second)) {
              dense[row * columns + column] = view(first, second);
            }
          }
        }
        return dense;
      };

  // Every group of a family has the same shape and size here, so one group's
  // view reads another group's values; the loop checks that before it does.
  const auto sameLayout = [](const auto& sizes, const auto& shapes) {
    for (std::size_t group = 1; group < Groups; ++group) {
      if (sizes[group] != sizes[0] || shapes[group][0] != shapes[0][0] ||
          shapes[group][1] != shapes[0][1]) {
        return false;
      }
    }
    return true;
  };

  SUBCASE("the fold is the two steps it was folded from") {
    REQUIRE(sameLayout(tensor::materialToFault::Size, tensor::materialToFault::Shape));
    REQUIRE(sameLayout(tensor::V3mTo2n::Size, tensor::V3mTo2n::Shape));

    // materialToFault[side, relation] = V3mTo2n[side, relation] * materialProject:
    // the samples go to the fault points without the modal coefficients ever
    // being written, and that shortcut has to be the long way round.
    const auto project = denseView(
        init::materialProject::view::create(const_cast<real*>(init::materialProject::Values)),
        Modes,
        Samples);

    for (std::size_t group = 0; group < Groups; ++group) {
      const auto evaluate = denseView(
          init::V3mTo2n::view<0, 0>::create(const_cast<real*>(init::V3mTo2n::Values[group])),
          FaultPoints,
          Modes,
          true);
      const auto folded = denseView(init::materialToFault::view<0, 0>::create(
                                        const_cast<real*>(init::materialToFault::Values[group])),
                                    FaultPoints,
                                    Samples);

      for (std::size_t point = 0; point < FaultPoints; ++point) {
        for (std::size_t sample = 0; sample < Samples; ++sample) {
          double expected = 0.0;
          for (std::size_t mode = 0; mode < Modes; ++mode) {
            expected += evaluate[point * Modes + mode] * project[mode * Samples + sample];
          }
          const double scale = std::max(1.0, std::abs(expected));
          REQUIRE(folded[point * Samples + sample] ==
                  doctest::Approx(expected).epsilon(Tolerance).scale(scale));
        }
      }
    }
  }

  SUBCASE("a constant material reaches every point") {
    // Independent of any matrix product: reading a constant field at the points
    // of the fault has to give that constant, whichever side and whichever
    // face relation.
    alignas(Alignment) std::array<real, tensor::materialSamples::size()> samples{};
    alignas(Alignment) std::array<real, tensor::materialAtFault::size()> atFault{};

    constexpr double Value = 2.75;
    for (std::size_t sample = 0; sample < Samples; ++sample) {
      // the material is one field, so every fused simulation sees the same sample
      for (std::size_t sim = 0; sim < multisim::NumSimulations; ++sim) {
        samples[sim + multisim::NumSimulations * sample] = static_cast<real>(Value);
      }
    }

    dynamicRupture::kernel::projectMaterialToFault krnl{};
    krnl.materialSamples = samples.data();
    krnl.materialAtFault = atFault.data();
    krnl.bindGlobals(seissol::Pool::host());

    for (std::uint8_t side = 0; side < Cell::NumFaces; ++side) {
      for (std::uint8_t relation = 0; relation < dr::misc::NumFaceRelations; ++relation) {
        atFault.fill(0.0);
        krnl.execute(side, relation);
        for (std::size_t point = 0; point < FaultPoints * multisim::NumSimulations; ++point) {
          REQUIRE(atFault[point] == doctest::Approx(Value).epsilon(Tolerance));
        }
      }
    }
  }
}

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
    double points[Points][3]{};
    double weights[Points]{};
    seissol::quadrature::TetrahedronQuadrature(points, weights, PerDirection);
    check(init::materialToQuadrature::view::create(init::materialToQuadrature::Values),
          Points,
          [&](std::size_t point) {
            return std::array<double, 3>{points[point][0], points[point][1], points[point][2]};
          });
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_MATERIALSAMPLING_T_H_
