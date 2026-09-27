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
// quadrature points of the dynamic rupture rule for a fault. Each of those is
// one matrix that takes the samples straight to the points, folded from the
// projection to the modal basis and the evaluation at the points.

#include <doctest.h>

#include "Alignment.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Solver/MultipleSimulations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("Material at the points of a fault") {
  constexpr std::size_t Samples = tensor::materialNodes::Shape[0];
  constexpr std::size_t FaultPoints = tensor::materialAtFault::Shape[0];
  constexpr std::size_t Modes = tensor::materialProject::Shape[0];
  constexpr std::size_t Groups = 16;

  const auto denseView = [](auto view, std::size_t rows, std::size_t columns) {
    std::vector<double> dense(rows * columns, 0.0);
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t column = 0; column < columns; ++column) {
        if (view.isInRange(row, column)) {
          dense[row * columns + column] = view(row, column);
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
          Modes);
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
                  doctest::Approx(expected).epsilon(1e-12).scale(scale));
        }
      }
    }
  }

  SUBCASE("a constant material reaches every point") {
    // Independent of any matrix product: reading a constant field at the points
    // of the fault has to give that constant, whichever side and whichever
    // reparametrisation of the shared face.
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
    for (std::size_t side = 0; side < 4; ++side) {
      for (std::size_t relation = 0; relation < 4; ++relation) {
        krnl.materialToFault(side, relation) =
            init::materialToFault::Values[tensor::materialToFault::index(side, relation)];
      }
    }

    for (std::uint8_t side = 0; side < 4; ++side) {
      for (std::uint8_t relation = 0; relation < 4; ++relation) {
        atFault.fill(0.0);
        krnl.execute(side, relation);
        for (std::size_t point = 0; point < FaultPoints * multisim::NumSimulations; ++point) {
          REQUIRE(atFault[point] == doctest::Approx(Value).epsilon(1e-12));
        }
      }
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_MATERIALSAMPLING_T_H_
