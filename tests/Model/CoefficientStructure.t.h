// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_COEFFICIENTSTRUCTURE_T_H_
#define SEISSOL_TESTS_MODEL_COEFFICIENTSTRUCTURE_T_H_

// A material's coefficient decomposition does not depend on the MaterialT of
// the build, so these tests run in every build. Only the assembly into the
// generated star layout needs a build whose quantities match.

#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "GeneratedCode/init.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <doctest.h>
#include <random>
#include <type_traits>

namespace seissol::unit_test {

namespace coefficients {

/// The coefficient matrix as the declared decomposition builds it.
template <typename MaterialT, std::size_t N>
Eigen::Matrix<double, N, N> assembled(const MaterialT& material, unsigned dim) {
  using Setup = seissol::model::MaterialSetup<MaterialT>;
  const auto coefficients = Setup::getCoefficients(material);
  Eigen::Matrix<double, N, N> matrix = Eigen::Matrix<double, N, N>::Zero();
  for (const auto& entry : Setup::CoefficientEntries) {
    if (entry.dim == dim) {
      matrix(entry.row, entry.column) += entry.factor * coefficients.at(entry.coefficient);
    }
  }
  return matrix;
}

/// The declaration has to reproduce what the material writes itself.
template <typename MaterialT, std::size_t N>
void checkDeclaration(const MaterialT& material) {
  for (unsigned dim = 0; dim < 3; ++dim) {
    Eigen::Matrix<double, N, N> reference = Eigen::Matrix<double, N, N>::Zero();
    seissol::model::MaterialSetup<MaterialT>::getTransposedCoefficientMatrix(
        material, dim, reference);

    const auto candidate = assembled<MaterialT, N>(material, dim);
    const double scale = std::max(1.0, reference.cwiseAbs().maxCoeff());
    for (std::size_t row = 0; row < N; ++row) {
      for (std::size_t column = 0; column < N; ++column) {
        REQUIRE(candidate(row, column) ==
                doctest::Approx(reference(row, column)).epsilon(1e-13).scale(scale));
      }
    }
  }
}

/// No entry may be claimed twice within one direction, or the decomposition
/// would silently add two contributions into the same slot.
template <typename MaterialT, std::size_t N>
void checkUnique() {
  std::array<std::array<bool, N * N>, 3> seen{};
  for (const auto& entry : seissol::model::MaterialSetup<MaterialT>::CoefficientEntries) {
    auto& slot = seen.at(entry.dim).at(entry.row * N + entry.column);
    REQUIRE_FALSE(slot);
    slot = true;
  }
}

inline seissol::model::ElasticMaterial elastic(double rho, double mu, double lambda) {
  seissol::model::ElasticMaterial material{};
  material.rho = rho;
  material.mu = mu;
  material.lambda = lambda;
  return material;
}

inline seissol::model::AcousticMaterial acoustic(double rho, double lambda) {
  seissol::model::AcousticMaterial material{};
  material.rho = rho;
  material.lambda = lambda;
  return material;
}

inline seissol::model::AnisotropicMaterial anisotropic(std::mt19937& rng) {
  std::uniform_real_distribution<double> modulus(1e9, 1e11);
  seissol::model::AnisotropicMaterial material{};
  material.rho = 2500.0;
  for (auto* component :
       {&material.c11, &material.c12, &material.c13, &material.c14, &material.c15,
        &material.c16, &material.c22, &material.c23, &material.c24, &material.c25,
        &material.c26, &material.c33, &material.c34, &material.c35, &material.c36,
        &material.c44, &material.c45, &material.c46, &material.c55, &material.c56,
        &material.c66}) {
    *component = modulus(rng);
  }
  return material;
}

/// A source term's entries, as the declared decomposition builds them.
template <typename MaterialT, std::size_t N>
Eigen::Matrix<double, N, N> assembledSource(const MaterialT& material, std::size_t mech) {
  using Setup = seissol::model::MaterialSetup<MaterialT>;
  const auto coefficients = Setup::getSourceCoefficients(material, mech);
  Eigen::Matrix<double, N, N> matrix = Eigen::Matrix<double, N, N>::Zero();
  for (const auto& entry : Setup::SourceEntries) {
    matrix(entry.row, entry.column) += entry.factor * coefficients.at(entry.coefficient);
  }
  return matrix;
}

template <typename MaterialT, std::size_t N>
void checkSourceDeclaration(const MaterialT& material, std::size_t mech) {
  Eigen::Matrix<double, N, N> reference = Eigen::Matrix<double, N, N>::Zero();
  seissol::model::MaterialSetup<MaterialT>::forEachSourceEntry(
      material, mech, [&reference](std::size_t row, std::size_t column, double value) {
        reference(row, column) += value;
      });

  const auto candidate = assembledSource<MaterialT, N>(material, mech);
  const double scale = std::max(1.0, reference.cwiseAbs().maxCoeff());
  for (std::size_t row = 0; row < N; ++row) {
    for (std::size_t column = 0; column < N; ++column) {
      REQUIRE(candidate(row, column) ==
              doctest::Approx(reference(row, column)).epsilon(1e-13).scale(scale));
    }
  }
}

template <std::size_t Mechanisms>
seissol::model::ViscoElasticMaterial<Mechanisms> viscoelastic(std::mt19937& rng) {
  std::uniform_real_distribution<double> value(-1e11, -1e9);
  seissol::model::ViscoElasticMaterial<Mechanisms> material{};
  material.rho = 2500.0;
  material.mu = 3e10;
  material.lambda = 2e10;
  for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
    for (std::size_t component = 0; component < 3; ++component) {
      material.theta[mech][component] = value(rng);
    }
  }
  return material;
}

inline seissol::model::PoroElasticMaterial poroelastic(std::mt19937& rng) {
  std::uniform_real_distribution<double> unit(0.1, 0.9);
  seissol::model::PoroElasticMaterial material{};
  material.rho = 2500.0;
  material.mu = 1e10;
  material.lambda = 1.2e10;
  material.bulkSolid = 4e10;
  material.porosity = unit(rng) * 0.3;
  material.permeability = 1e-12;
  material.tortuosity = 1.0 + unit(rng);
  material.bulkFluid = 2.2e9;
  material.rhoFluid = 1000.0;
  material.viscosity = 1e-3;
  return material;
}

} // namespace coefficients

TEST_CASE("Coefficient decomposition") {
  std::mt19937 rng(20260926);
  std::uniform_real_distribution<double> positive(0.5, 5.0);

  SUBCASE("elastic") {
    for (std::size_t sample = 0; sample < 64; ++sample) {
      coefficients::checkDeclaration<seissol::model::ElasticMaterial, 9>(coefficients::elastic(
          positive(rng) * 1000.0, positive(rng) * 1e10, positive(rng) * 1e10));
    }
    // acoustic material suppresses the shear-coupled 1/rho, which the
    // decomposition carries as a coefficient of its own
    for (std::size_t sample = 0; sample < 16; ++sample) {
      coefficients::checkDeclaration<seissol::model::ElasticMaterial, 9>(
          coefficients::elastic(positive(rng) * 1000.0, 0.0, positive(rng) * 1e10));
    }
    coefficients::checkUnique<seissol::model::ElasticMaterial, 9>();
  }

  SUBCASE("acoustic") {
    for (std::size_t sample = 0; sample < 16; ++sample) {
      coefficients::checkDeclaration<seissol::model::AcousticMaterial, 4>(
          coefficients::acoustic(positive(rng) * 1000.0, positive(rng) * 1e10));
    }
    coefficients::checkUnique<seissol::model::AcousticMaterial, 4>();
  }

  SUBCASE("anisotropic") {
    for (std::size_t sample = 0; sample < 32; ++sample) {
      coefficients::checkDeclaration<seissol::model::AnisotropicMaterial, 9>(
          coefficients::anisotropic(rng));
    }
    // no uniqueness check here: the three directional matrices share their
    // stress block, so a slot legitimately takes one entry per direction
  }

  SUBCASE("poroelastic") {
    for (std::size_t sample = 0; sample < 32; ++sample) {
      coefficients::checkDeclaration<seissol::model::PoroElasticMaterial, 13>(
          coefficients::poroelastic(rng));
    }
    // the declaration assumes an isotropic frame, which is what makes cBar
    // three values rather than twenty-one; the comparison above is what
    // catches an anisotropic one
    coefficients::checkUnique<seissol::model::PoroElasticMaterial, 13>();
  }

  SUBCASE("solver operator including the relaxation blocks") {
    constexpr std::size_t Mechanisms = 3;
    using Material = seissol::model::ViscoElasticMaterial<Mechanisms>;
    constexpr std::size_t N = Material::NumQuantities;

    for (std::size_t sample = 0; sample < 16; ++sample) {
      auto material = coefficients::viscoelastic<Mechanisms>(rng);
      for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
        material.omega[mech] = positive(rng);
      }

      for (unsigned dim = 0; dim < 3; ++dim) {
        Eigen::Matrix<double, N, N> reference = Eigen::Matrix<double, N, N>::Zero();
        seissol::model::SolverSetup<typename Material::Solver, Material>::
            getTransposedCoefficientMatrix(material, dim, reference);

        using Setup = seissol::model::SolverSetup<typename Material::Solver, Material>;
        const auto solverCoefficients = Setup::getCoefficients(material);
        Eigen::Matrix<double, N, N> candidate = Eigen::Matrix<double, N, N>::Zero();
        Setup::forEachCoefficientEntry([&](std::size_t coefficient,
                                           std::size_t entryDim,
                                           std::size_t row,
                                           std::size_t column,
                                           double factor) {
          if (entryDim == dim) {
            candidate(row, column) += factor * solverCoefficients.at(coefficient);
          }
        });

        const double scale = std::max(1.0, reference.cwiseAbs().maxCoeff());
        for (std::size_t row = 0; row < N; ++row) {
          for (std::size_t column = 0; column < N; ++column) {
            REQUIRE(candidate(row, column) ==
                    doctest::Approx(reference(row, column)).epsilon(1e-13).scale(scale));
          }
        }
      }
    }
  }

  SUBCASE("viscoelastic source") {
    constexpr std::size_t Mechanisms = 3;
    for (std::size_t sample = 0; sample < 16; ++sample) {
      const auto material = coefficients::viscoelastic<Mechanisms>(rng);
      // the flux is the base material's, so its decomposition has to be too
      coefficients::checkDeclaration<seissol::model::ViscoElasticMaterial<Mechanisms>, 9>(material);
      for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
        coefficients::checkSourceDeclaration<seissol::model::ViscoElasticMaterial<Mechanisms>, 6>(
            material, mech);
      }
    }
  }
}

TEST_CASE("Star assembly from coefficients") {
  if constexpr (std::is_same_v<seissol::model::MaterialT, seissol::model::ElasticMaterial>) {
    std::mt19937 rng(20260927);
    std::uniform_real_distribution<double> positive(0.5, 5.0);
    std::normal_distribution<double> gauss(0.0, 1.0);

    for (std::size_t sample = 0; sample < 64; ++sample) {
      const auto material = coefficients::elastic(
          positive(rng) * 1000.0, positive(rng) * 1e10, positive(rng) * 1e10);
      const double gradient[3] = {gauss(rng), gauss(rng), gauss(rng)};

      // what CellLocalMatrices builds today: the three directional matrices,
      // each scaled by its row of the Jacobian
      std::array<std::array<double, seissol::tensor::star::size(0)>, 3> directional{};
      for (unsigned dim = 0; dim < 3; ++dim) {
        auto view = seissol::init::star::view<0>::create(directional.at(dim).data());
        seissol::model::MaterialSetup<seissol::model::ElasticMaterial>::
            getTransposedCoefficientMatrix(material, dim, view);
      }

      std::array<double, seissol::tensor::star::size(0)> assembledData{};
      auto view = seissol::init::star::view<0>::create(assembledData.data());
      seissol::model::assembleStarMatrix<seissol::model::ElasticMaterial>(
          seissol::model::MaterialSetup<seissol::model::ElasticMaterial>::getCoefficients(material),
          gradient,
          view);

      for (std::size_t idx = 0; idx < assembledData.size(); ++idx) {
        const double reference = gradient[0] * directional[0].at(idx) +
                                 gradient[1] * directional[1].at(idx) +
                                 gradient[2] * directional[2].at(idx);
        const double scale = std::max(1.0, std::abs(reference));
        REQUIRE(assembledData.at(idx) == doctest::Approx(reference).epsilon(1e-13).scale(scale));
      }
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_COEFFICIENTSTRUCTURE_T_H_
