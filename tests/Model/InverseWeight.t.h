// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_INVERSEWEIGHT_T_H_
#define SEISSOL_TESTS_MODEL_INVERSEWEIGHT_T_H_

// The coefficient matrix of every direction is the inverse weight of the medium times one without
// material, A_d = W^-1 B_d (see WeightEntry): checked against the matrices the materials write
// themselves, for solids, a fluid in an elastic build and acoustic media.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "Model/Common.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cstddef>
#include <random>

namespace seissol::unit_test {

namespace inverseweight {

/// The transposed coefficient matrix of direction `dim` from the entries the material declares,
/// with the coefficients given.
template <typename MaterialT, std::size_t N, typename ValuesT>
Eigen::Matrix<double, N, N> fromEntries(const ValuesT& values, unsigned dim) {
  Eigen::Matrix<double, N, N> matrix = Eigen::Matrix<double, N, N>::Zero();
  for (const auto& entry : seissol::model::MaterialSetup<MaterialT>::CoefficientEntries) {
    if (entry.dim == dim) {
      matrix(entry.row, entry.column) += entry.factor * values.at(entry.coefficient);
    }
  }
  return matrix;
}

/// The inverse weight of a material, from the entries it declares.
template <typename MaterialT, std::size_t N>
Eigen::Matrix<double, N, N> inverseWeight(const MaterialT& material) {
  using Setup = seissol::model::MaterialSetup<MaterialT>;
  const auto values = Setup::getCoefficients(material);
  Eigen::Matrix<double, N, N> weight = Eigen::Matrix<double, N, N>::Zero();
  for (const auto& entry : Setup::WeightEntries) {
    weight(entry.row, entry.column) += entry.factor * values.at(entry.coefficient);
  }
  return weight;
}

/// A_d^T = B_d^T W^-1 for every direction, against the matrices the material writes itself.
template <typename MaterialT, std::size_t N>
void checkWeight(const MaterialT& material) {
  using Setup = seissol::model::MaterialSetup<MaterialT>;
  const auto weight = inverseWeight<MaterialT, N>(material);
  REQUIRE(weight == weight.transpose());
  for (unsigned dim = 0; dim < 3; ++dim) {
    Eigen::Matrix<double, N, N> reference = Eigen::Matrix<double, N, N>::Zero();
    Setup::getTransposedCoefficientMatrix(material, dim, reference);
    const Eigen::Matrix<double, N, N> candidate =
        fromEntries<MaterialT, N>(Setup::getUnweightedCoefficients(material), dim) * weight;
    const double scale = std::max(1.0, reference.cwiseAbs().maxCoeff());
    for (std::size_t row = 0; row < N; ++row) {
      for (std::size_t column = 0; column < N; ++column) {
        REQUIRE(candidate(row, column) ==
                doctest::Approx(reference(row, column)).epsilon(1e-13).scale(scale));
      }
    }
  }
}

} // namespace inverseweight

TEST_CASE("The operator of a medium is its inverse weight times the one without material") {
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937 rng(20261015);
  std::uniform_real_distribution<double> positive(0.5, 5.0);
  SUBCASE("elastic") {
    for (std::size_t sample = 0; sample < 16; ++sample) {
      inverseweight::checkWeight<seissol::model::ElasticMaterial, 9>(coefficients::elastic(
          positive(rng) * 1000.0, positive(rng) * 1e10, positive(rng) * 1e10));
    }
  }
  SUBCASE("a fluid in an elastic build") {
    inverseweight::checkWeight<seissol::model::ElasticMaterial, 9>(
        coefficients::elastic(1000.0, 0.0, 2.25e9));
  }
  SUBCASE("acoustic") {
    for (std::size_t sample = 0; sample < 16; ++sample) {
      inverseweight::checkWeight<seissol::model::AcousticMaterial, 4>(
          coefficients::acoustic(positive(rng) * 1000.0, positive(rng) * 1e10));
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_INVERSEWEIGHT_T_H_
