// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Kernels/Receiver.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {

namespace receivertest {
/// d_j v_i of the velocity field the derived quantities are evaluated on. The entries are pairwise
/// distinct, so that mixing up a component and a direction changes the result, and they scale with
/// the simulation, so that mixing up simulations does as well.
inline real velocityGradient(std::size_t sim, std::size_t i, std::size_t j) {
  constexpr std::array<std::array<real, Cell::Dim>, Cell::Dim> Gradient{
      {{2, 3, 5}, {7, 11, 13}, {17, 19, 23}}};
  return static_cast<real>(sim + 1) * Gradient[i][j];
}

/// Evaluates a derived receiver quantity for every simulation on point values that carry the
/// velocity gradient above and a poison value everywhere else.
inline std::vector<real> evaluate(kernels::DerivedReceiverQuantity& derived) {
  // Larger than any entry of the velocity gradient, and finite, so that reading it does not depend
  // on how the floating-point mode treats NaN.
  constexpr real Poison = 1000;
  // Twice the size of the tensors: whatever is read besides the velocity gradient -- another
  // quantity, or something past the end of the tensor -- is the poison value and spoils the result,
  // rather than being whatever happens to follow in memory.
  std::vector<real> qAtPointData(2 * static_cast<std::size_t>(tensor::QAtPoint::size()), Poison);
  std::vector<real> qDerivativeAtPointData(
      2 * static_cast<std::size_t>(tensor::QDerivativeAtPoint::size()), Poison);
  auto qAtPoint = init::QAtPoint::view::create(qAtPointData.data());
  auto qDerivativeAtPoint = init::QDerivativeAtPoint::view::create(qDerivativeAtPointData.data());
  for (auto sim = multisim::MultisimStart; sim < multisim::MultisimEnd; ++sim) {
    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      for (std::size_t j = 0; j < Cell::Dim; ++j) {
        multisim::multisimWrap(qDerivativeAtPoint, sim, model::MaterialT::VelocityOffset + i, j) =
            velocityGradient(sim, i, j);
      }
    }
  }
  std::vector<real> output;
  for (auto sim = multisim::MultisimStart; sim < multisim::MultisimEnd; ++sim) {
    const auto before = output.size();
    derived.compute(sim, output, qAtPoint, qDerivativeAtPoint);
    // ReceiverCluster::ncols reserves exactly this many columns per simulation
    REQUIRE(output.size() - before == derived.quantities().size());
  }
  return output;
}
} // namespace receivertest

TEST_CASE("Receiver rotation is the curl of the velocity" * doctest::test_suite("kernel")) {
  kernels::ReceiverRotation rotation;
  const auto output = receivertest::evaluate(rotation);
  for (auto sim = multisim::MultisimStart; sim < multisim::MultisimEnd; ++sim) {
    const auto dv = [&](std::size_t i, std::size_t j) {
      return receivertest::velocityGradient(sim, i, j);
    };
    const std::array<real, 3> expected{
        dv(2, 1) - dv(1, 2), dv(0, 2) - dv(2, 0), dv(1, 0) - dv(0, 1)};
    for (std::size_t k = 0; k < expected.size(); ++k) {
      CHECK(output[(sim - multisim::MultisimStart) * expected.size() + k] ==
            doctest::Approx(expected[k]));
    }
  }
}

TEST_CASE("Receiver strain rate is the symmetric velocity gradient" *
          doctest::test_suite("kernel")) {
  kernels::ReceiverStrain strain;
  const auto output = receivertest::evaluate(strain);
  for (auto sim = multisim::MultisimStart; sim < multisim::MultisimEnd; ++sim) {
    const auto dv = [&](std::size_t i, std::size_t j) {
      return receivertest::velocityGradient(sim, i, j);
    };
    // xx, xy, xz, yy, yz, zz
    const std::array<real, 6> expected{dv(0, 0),
                                       (dv(0, 1) + dv(1, 0)) / 2,
                                       (dv(0, 2) + dv(2, 0)) / 2,
                                       dv(1, 1),
                                       (dv(1, 2) + dv(2, 1)) / 2,
                                       dv(2, 2)};
    for (std::size_t k = 0; k < expected.size(); ++k) {
      CHECK(output[(sim - multisim::MultisimStart) * expected.size() + k] ==
            doctest::Approx(expected[k]));
    }
  }
}

} // namespace seissol::unit_test
