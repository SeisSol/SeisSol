// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/elastic/Model/Setup.h"
#include "GeneratedCode/init.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cstddef>
#include <doctest.h>
#include <random>

namespace seissol::unit_test {

namespace {

using Setup = seissol::model::MaterialSetup<seissol::model::ElasticMaterial>;

seissol::model::ElasticMaterial makeMaterial(double rho, double mu, double lambda) {
  seissol::model::ElasticMaterial material{};
  material.rho = rho;
  material.mu = mu;
  material.lambda = lambda;
  return material;
}

/// The coefficient matrix as the declared decomposition builds it.
Eigen::Matrix<double, 9, 9> assembled(const seissol::model::ElasticMaterial& material,
                                      unsigned dim) {
  const auto coefficients = Setup::getCoefficients(material);
  Eigen::Matrix<double, 9, 9> matrix = Eigen::Matrix<double, 9, 9>::Zero();
  for (const auto& entry : Setup::CoefficientEntries) {
    if (entry.dim == dim) {
      matrix(entry.row, entry.column) += entry.factor * coefficients.at(entry.coefficient);
    }
  }
  return matrix;
}

/// The packed star matrix as CellLocalMatrices builds it today: the three
/// directional matrices, each scaled by its row of the Jacobian.
std::array<double, seissol::tensor::star::size(0)>
    referenceStar(const seissol::model::ElasticMaterial& material, const double gradient[3]) {
  std::array<double, seissol::tensor::star::size(0)> matA{};
  std::array<double, seissol::tensor::star::size(0)> matB{};
  std::array<double, seissol::tensor::star::size(0)> matC{};
  auto viewA = seissol::init::star::view<0>::create(matA.data());
  auto viewB = seissol::init::star::view<0>::create(matB.data());
  auto viewC = seissol::init::star::view<0>::create(matC.data());
  Setup::getTransposedCoefficientMatrix(material, 0, viewA);
  Setup::getTransposedCoefficientMatrix(material, 1, viewB);
  Setup::getTransposedCoefficientMatrix(material, 2, viewC);

  std::array<double, seissol::tensor::star::size(0)> result{};
  for (std::size_t idx = 0; idx < result.size(); ++idx) {
    result[idx] = gradient[0] * matA[idx] + gradient[1] * matB[idx] + gradient[2] * matC[idx];
  }
  return result;
}

void compareAssembly(const seissol::model::ElasticMaterial& material, const double gradient[3]) {
  const auto reference = referenceStar(material, gradient);

  std::array<double, seissol::tensor::star::size(0)> assembledData{};
  auto view = seissol::init::star::view<0>::create(assembledData.data());
  seissol::model::assembleStarMatrix<seissol::model::ElasticMaterial>(
      Setup::getCoefficients(material), gradient, view);

  const double scale = std::max(1.0, *std::max_element(reference.begin(), reference.end()));
  for (std::size_t idx = 0; idx < reference.size(); ++idx) {
    REQUIRE(assembledData[idx] ==
            doctest::Approx(reference[idx]).epsilon(1e-13).scale(scale));
  }
}

void compareFor(const seissol::model::ElasticMaterial& material) {
  for (unsigned dim = 0; dim < 3; ++dim) {
    Eigen::Matrix<double, 9, 9> reference = Eigen::Matrix<double, 9, 9>::Zero();
    Setup::getTransposedCoefficientMatrix(material, dim, reference);

    const auto candidate = assembled(material, dim);
    for (std::size_t row = 0; row < 9; ++row) {
      for (std::size_t column = 0; column < 9; ++column) {
        REQUIRE(candidate(row, column) ==
                doctest::Approx(reference(row, column)).epsilon(1e-14));
      }
    }
  }
}

} // namespace

TEST_CASE("Elastic coefficient decomposition") {
  std::mt19937 rng(20260926);
  std::uniform_real_distribution<double> positive(0.5, 5.0);

  SUBCASE("elastic") {
    for (std::size_t sample = 0; sample < 64; ++sample) {
      compareFor(makeMaterial(positive(rng) * 1000.0, positive(rng) * 1e10, positive(rng) * 1e10));
    }
  }

  SUBCASE("acoustic") {
    for (std::size_t sample = 0; sample < 16; ++sample) {
      compareFor(makeMaterial(positive(rng) * 1000.0, 0.0, positive(rng) * 1e10));
    }
  }

  SUBCASE("assembly against the folded matrices") {
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (std::size_t sample = 0; sample < 64; ++sample) {
      const double gradient[3] = {gauss(rng), gauss(rng), gauss(rng)};
      compareAssembly(
          makeMaterial(positive(rng) * 1000.0, positive(rng) * 1e10, positive(rng) * 1e10),
          gradient);
      compareAssembly(makeMaterial(positive(rng) * 1000.0, 0.0, positive(rng) * 1e10), gradient);
    }
  }

  SUBCASE("one entry per star slot") {
    // the three directional matrices are disjointly occupied, so folding the
    // Jacobian costs exactly one multiplication per entry
    std::array<std::array<bool, 81>, 3> seen{};
    for (const auto& entry : Setup::CoefficientEntries) {
      auto& slot = seen.at(entry.dim).at(entry.row * 9 + entry.column);
      REQUIRE_FALSE(slot);
      slot = true;
    }
  }
}

} // namespace seissol::unit_test
