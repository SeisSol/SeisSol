// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_SYMMETRICTETRAHEDRONRULE_H_
#define SEISSOL_SRC_NUMERICAL_SYMMETRICTETRAHEDRONRULE_H_

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace seissol::quadrature {

//! @brief The highest degree of exactness symmetricTetrahedronRule provides.
constexpr std::size_t MaxSymmetricTetrahedronDegree = 14;

//! @brief The number of points of symmetricTetrahedronRule(degree); 0 above
//! MaxSymmetricTetrahedronDegree.
constexpr std::size_t symmetricTetrahedronRuleSize(std::size_t degree) {
  // (degree of exactness, number of points) of the tabulated rules
  constexpr std::array<std::array<std::size_t, 2>, 12> Sizes{{{1, 1},
                                                              {2, 4},
                                                              {3, 8},
                                                              {5, 14},
                                                              {6, 24},
                                                              {7, 35},
                                                              {8, 46},
                                                              {10, 79},
                                                              {11, 110},
                                                              {12, 128},
                                                              {13, 145},
                                                              {14, 188}}};
  for (const auto& [ruleDegree, size] : Sizes) {
    if (ruleDegree >= degree) {
      return size;
    }
  }
  return 0;
}

/**
 * @brief A fully symmetric quadrature rule on the reference tetrahedron (0,0,0), (1,0,0),
 * (0,1,0), (0,0,1) that is exact up to the given degree (at most
 * MaxSymmetricTetrahedronDegree).
 *
 * The rule is invariant under all permutations of the vertices, has positive weights (summing to
 * 1/6) and keeps all barycentric coordinates above 1e-3. It is the rule with the fewest points
 * among the tabulated ones that satisfies this, i.e. quad_rules.SymmetricTet of SeisSol/matrices,
 * with the points in the same order: orbit by orbit, each orbit a pattern of barycentric
 * coordinates (cf. SymmetricTetrahedronRule.cpp) whose distinct permutations follow in
 * lexicographic order.
 *
 * @return [points, weights]
 */
std::pair<std::vector<std::array<double, 3>>, std::vector<double>>
    symmetricTetrahedronRule(std::size_t degree);

} // namespace seissol::quadrature

#endif // SEISSOL_SRC_NUMERICAL_SYMMETRICTETRAHEDRONRULE_H_
