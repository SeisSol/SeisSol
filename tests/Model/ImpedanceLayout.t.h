// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_IMPEDANCELAYOUT_T_H_
#define SEISSOL_TESTS_MODEL_IMPEDANCELAYOUT_T_H_

#include <doctest.h>

#include "Equations/Datastructures.h"
#include "Equations/Impedance.h" // IWYU pragma: keep
#include "Equations/ImpedanceBase.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "TestConfigs.h"

#include <cstddef>

namespace seissol::unit_test {

/**
 * The quantity indices of ImpedanceCompute decide which rows initializeDynamicRuptureMatrices fills
 * in the traction averaging matrices. The code generator defines the same layout again, in
 * extractTractions / extractVelocities (and in tractionMatrixSpp, which the pattern tests of the
 * traction averaging matrices tie to TractionIndices). If they disagreed, the averaging matrices
 * and the Riemann solver would not agree on which quantities are tractions; in a release build
 * nothing else would notice.
 */
template <typename Cfg>
void checkInterfaceQuantities() {
  using MaterialT = model::MaterialOf<Cfg>;
  using ImpedanceCompute = seissol::model::ImpedanceCompute<MaterialT>;

  const auto check = [](const auto* values,
                        const auto& shape,
                        const auto& start,
                        const auto& stop,
                        const auto& indices) {
    REQUIRE(shape[0] == ImpedanceCompute::Dim);
    // (the anelastic solver extracts from the elastic quantities only, so this may be less than
    // MaterialT::NumQuantities)
    const std::size_t quantities = shape[1];
    // the generated matrix only stores its bounding box, column major
    const std::size_t rows = stop[0] - start[0];
    for (std::size_t row = 0; row < ImpedanceCompute::Dim; ++row) {
      REQUIRE(indices[row] < quantities);
      for (std::size_t quantity = 0; quantity < quantities; ++quantity) {
        const bool inBox =
            row >= start[0] && row < stop[0] && quantity >= start[1] && quantity < stop[1];
        const double value = inBox ? values[(row - start[0]) + rows * (quantity - start[1])] : 0.0;
        CHECK(value == (quantity == indices[row] ? 1.0 : 0.0));
      }
    }
  };

  check(init::extractTractions<Cfg>::Values,
        tensor::extractTractions<Cfg>::Shape,
        init::extractTractions<Cfg>::Start,
        init::extractTractions<Cfg>::Stop,
        ImpedanceCompute::TractionIndices);
  check(init::extractVelocities<Cfg>::Values,
        tensor::extractVelocities<Cfg>::Shape,
        init::extractVelocities<Cfg>::Start,
        init::extractVelocities<Cfg>::Stop,
        ImpedanceCompute::VelocityIndices);
}

TEST_CASE_TEMPLATE_DEFINE(
    "Interface quantities of the impedance match the generated extraction matrices" *
        doctest::test_suite("dynamicrupture"),
    Cfg,
    InterfaceQuantities) {
  checkInterfaceQuantities<Cfg>();
}

TEST_CASE_TEMPLATE_APPLY(InterfaceQuantities, DynamicRuptureConfigs);

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_IMPEDANCELAYOUT_T_H_
