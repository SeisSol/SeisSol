// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "Kernels/Precision.h"
#include "Numerical/Functions.h"
#include "Numerical/TimeBasis.h"
#include "TestHelper.h"

namespace seissol::unit_test {

TEST_CASE_TEMPLATE("Monomial Basis", RealT, float, double) {
  // needed due to absolute error
  constexpr double Delta = 10 * std::numeric_limits<double>::epsilon();

  // https://stackoverflow.com/a/50132998
  const auto taylorTerm = [](auto x, auto j) { return std::pow(x, j) / std::tgamma(j + 1); };
  for (std::size_t i = 1; i < 10; ++i) {
    const auto basis = numerical::MonomialBasis<RealT>(i);
    constexpr double Span = 15.0;
    for (const auto& point : {0.0, 1.0, 2.0, 0.123, 15.0}) {
      const auto exp = basis.point(point, Span);
      const auto expD = basis.derivative(point, Span);
      const auto expInt = basis.integrate(0, point, Span);
      for (std::size_t j = 0; j < i; ++j) {
        REQUIRE(exp[j] == AbsApprox(static_cast<RealT>(taylorTerm(point, j))).delta(Delta));
        REQUIRE(expInt[j] == AbsApprox(static_cast<RealT>(taylorTerm(point, j + 1))).delta(Delta));

        // needed due to size_t being unsigned
        if (j > 0) {
          REQUIRE(expD[j] == AbsApprox(static_cast<RealT>(taylorTerm(point, j - 1))).delta(Delta));
        } else {
          REQUIRE(expD[j] == AbsApprox(0.0));
        }
      }
    }
  }
}

TEST_CASE_TEMPLATE("Legendre Basis", RealT, float, double) {
  // needed due to absolute error
  constexpr double Delta = 10 * std::numeric_limits<RealT>::epsilon();
  constexpr double Epsilon = 1000 * std::numeric_limits<RealT>::epsilon();

  for (std::size_t i = 1; i < 10; ++i) {
    const auto basis = numerical::LegendreBasis<RealT>(i);
    constexpr double Span = 15.0;
    for (const auto& point : {0.0, 1.0, 2.0, 0.123, 15.0}) {
      const auto exp = basis.point(point, Span);
      const auto expD = basis.derivative(point, Span);
      const auto expInt = basis.integrate(0, point, Span);

      const auto tau = point / Span;
      for (std::size_t j = 0; j < i; ++j) {
        REQUIRE(exp[j] == AbsApprox(static_cast<RealT>(functions::shiftedLegendre(j, tau, 0)))
                              .delta(Delta)
                              .epsilon(Epsilon));
        REQUIRE(expD[j] * Span ==
                AbsApprox(static_cast<RealT>(functions::shiftedLegendre(j, tau, 1)))
                    .delta(Delta)
                    .epsilon(Epsilon));
        REQUIRE(expInt[j] / Span ==
                AbsApprox(static_cast<RealT>(functions::shiftedLegendre(j, tau, -1) -
                                             functions::shiftedLegendre(j, 0, -1)))
                    .delta(Delta)
                    .epsilon(Epsilon));
      }
    }
  }
}

} // namespace seissol::unit_test
