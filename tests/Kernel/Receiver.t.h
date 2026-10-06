// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Receiver.h"
#include "Solver/MultipleSimulations.h"
#include "TestConfigs.h"

#include <array>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {

namespace receivertest {
/// d_j v_i of the velocity field the derived quantities are evaluated on. The entries are pairwise
/// distinct, so that mixing up a component and a direction changes the result, and they scale with
/// the simulation, so that mixing up simulations does as well.
template <typename Cfg>
Real<Cfg> velocityGradient(std::size_t sim, std::size_t i, std::size_t j) {
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  constexpr std::array<std::array<real, Cell::Dim>, Cell::Dim> Gradient{
      {{2, 3, 5}, {7, 11, 13}, {17, 19, 23}}};
  return static_cast<real>(sim + 1) * Gradient[i][j];
}

/// Evaluates a derived receiver quantity for every simulation on point values that carry the
/// velocity gradient above and a poison value everywhere else.
template <typename Cfg>
std::vector<double> evaluate(const kernels::DerivedReceiverQuantity& derived) {
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)
  using Multisim = multisim::MultisimHelperWrapper<Cfg>;
  // Larger than any entry of the velocity gradient, and finite, so that reading it does not depend
  // on how the floating-point mode treats NaN.
  constexpr real Poison = 1000;
  // Twice the size of the tensor: whatever is read besides the velocity gradient -- another
  // quantity, or something past the end of the tensor -- is the poison value and spoils the result,
  // rather than being whatever happens to follow in memory.
  std::vector<real> qDerivativeAtPointData(
      2 * static_cast<std::size_t>(tensor::QDerivativeAtPoint<Cfg>::size()), Poison);
  auto qDerivativeAtPoint =
      init::QDerivativeAtPoint<Cfg>::view::create(qDerivativeAtPointData.data());
  for (auto sim = Multisim::MultisimStart; sim < Multisim::MultisimEnd; ++sim) {
    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      for (std::size_t j = 0; j < Cell::Dim; ++j) {
        multisim::multisimWrap<Cfg>(
            qDerivativeAtPoint, sim, model::MaterialOf<Cfg>::VelocityOffset + i, j) =
            velocityGradient<Cfg>(sim, i, j);
      }
    }
  }
  std::vector<double> output;
  for (auto sim = Multisim::MultisimStart; sim < Multisim::MultisimEnd; ++sim) {
    const auto before = output.size();
    derived.compute(output, kernels::velocityGradient<Cfg>(qDerivativeAtPoint, sim));
    // ReceiverCluster::ncols reserves exactly this many columns per simulation
    REQUIRE(output.size() - before == derived.quantities().size());
  }
  return output;
}
} // namespace receivertest

TEST_CASE_TEMPLATE("Receiver rotation is the curl of the velocity" * doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using Multisim = multisim::MultisimHelperWrapper<Cfg>;
  const kernels::ReceiverRotation rotation;
  const auto output = receivertest::evaluate<Cfg>(rotation);
  for (auto sim = Multisim::MultisimStart; sim < Multisim::MultisimEnd; ++sim) {
    const auto dv = [&](std::size_t i, std::size_t j) {
      return receivertest::velocityGradient<Cfg>(sim, i, j);
    };
    const std::array<Real<Cfg>, 3> expected{
        dv(2, 1) - dv(1, 2), dv(0, 2) - dv(2, 0), dv(1, 0) - dv(0, 1)};
    for (std::size_t k = 0; k < expected.size(); ++k) {
      CHECK(output[(sim - Multisim::MultisimStart) * expected.size() + k] ==
            doctest::Approx(expected[k]));
    }
  }
}

TEST_CASE_TEMPLATE("Receiver strain rate is the symmetric velocity gradient" *
                       doctest::test_suite("kernel"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using Multisim = multisim::MultisimHelperWrapper<Cfg>;
  const kernels::ReceiverStrain strain;
  const auto output = receivertest::evaluate<Cfg>(strain);
  for (auto sim = Multisim::MultisimStart; sim < Multisim::MultisimEnd; ++sim) {
    const auto dv = [&](std::size_t i, std::size_t j) {
      return receivertest::velocityGradient<Cfg>(sim, i, j);
    };
    // xx, xy, xz, yy, yz, zz
    const std::array<Real<Cfg>, 6> expected{dv(0, 0),
                                            (dv(0, 1) + dv(1, 0)) / 2,
                                            (dv(0, 2) + dv(2, 0)) / 2,
                                            dv(1, 1),
                                            (dv(1, 2) + dv(2, 1)) / 2,
                                            dv(2, 2)};
    for (std::size_t k = 0; k < expected.size(); ++k) {
      CHECK(output[(sim - Multisim::MultisimStart) * expected.size() + k] ==
            doctest::Approx(expected[k]));
    }
  }
}

} // namespace seissol::unit_test
