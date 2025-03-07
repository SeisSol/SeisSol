// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/Precision.h"
#include "Numerical/Eigenvalues.h"
#include "tests/TestHelper.h"

namespace seissol::unit_test {

template <unsigned Dim>
void testResidual(std::array<std::complex<double>, Dim * Dim>& m,
                  seissol::eigenvalues::Eigenpair<std::complex<double>, Dim>& eigenpair) {
  constexpr auto Epsilon = std::numeric_limits<double>::epsilon() * 1e3;

  // compute M*R
  std::array<std::complex<double>, Dim * Dim> mR{};
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      for (unsigned int k = 0; k < Dim; k++) {
        mR[i + Dim * j] += m[i + Dim * k] * eigenpair.vectors[k + Dim * j];
      }
    }
  }
  // compute R*L
  std::array<std::complex<double>, Dim * Dim> rL{};
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      rL[i + Dim * j] = eigenpair.vectors[i + Dim * j] * eigenpair.values[j];
    }
  }
  // compare residual
  for (unsigned i = 0; i < Dim * Dim; i++) {
    REQUIRE(std::abs(mR[i] - rL[i]) == AbsApprox(0.0).epsilon(Epsilon));
  }
}

TEST_CASE("Eigenvalues are correctly computed") {
  SUBCASE("Eigenvalues 3") {
    const unsigned dim = 3;

    std::array<std::array<std::complex<double>, dim * dim>, 2> matrices = {{
        {2.0, -3.0, -3.0, -2.0, 1.0, -2.0, -1.0, 1.0, 4.0},
        {-2.0, 3.0, 3.0, 2.0, -1.0, 2.0, 1.0, -1.0, -4.0},
    }};
    for (auto& m : matrices) {
      seissol::eigenvalues::Eigenpair<std::complex<double>, dim> eigenpair{};
      seissol::eigenvalues::computeEigenvalues(m, eigenpair);
      testResidual<dim>(m, eigenpair);
    }
  }
}

} // namespace seissol::unit_test
