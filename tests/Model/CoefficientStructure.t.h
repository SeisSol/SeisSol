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

#include <doctest.h>

#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/init.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Model/Common.h"
#include "Model/CommonDatastructures.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
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
       {&material.c11, &material.c12, &material.c13, &material.c14, &material.c15, &material.c16,
        &material.c22, &material.c23, &material.c24, &material.c25, &material.c26, &material.c33,
        &material.c34, &material.c35, &material.c36, &material.c44, &material.c45, &material.c46,
        &material.c55, &material.c56, &material.c66}) {
    *component = modulus(rng);
  }
  return material;
}

/// Not every viscous material attenuates shear: the acoustic one has no shear
/// to attenuate, so it carries one quality factor instead of two.
template <typename T, typename = void>
struct HasShearQuality : std::false_type {};
template <typename T>
struct HasShearQuality<T, std::void_t<decltype(std::declval<T&>().qs)>> : std::true_type {};

template <typename T>
void setShearQuality(T& material, double value) {
  if constexpr (HasShearQuality<T>::value) {
    material.qs = value;
  }
}

/// The tensor the configured solver states its source term in.
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
using SourceTensor = seissol::tensor::E;
using SourceInit = seissol::init::E;
#else
using SourceTensor = seissol::tensor::ET;
using SourceInit = seissol::init::ET;
#endif

/// A material of whatever type the build is configured for, with enough in it
/// that its source term is non-trivial.
template <typename MaterialT>
MaterialT configuredMaterial(std::mt19937& rng) {
  const auto parameters = [] {
    seissol::initializer::parameters::ModelParameters p{};
    p.freqCentral = 1.0;
    p.freqRatio = 100.0;
    return p;
  }();
  std::uniform_real_distribution<double> unit(0.3, 0.9);

  MaterialT material{};
  if constexpr (std::is_base_of_v<seissol::model::PoroElasticMaterial, MaterialT>) {
    material.rho = 2500.0;
    material.lambda = 1.2e10;
    material.mu = 1.0e10;
    material.bulkSolid = 4.0e10;
    material.porosity = unit(rng) * 0.3;
    material.permeability = 6.0e-13 * unit(rng);
    material.tortuosity = 1.0 + 2.0 * unit(rng);
    material.bulkFluid = 2.5e9;
    material.rhoFluid = 1040.0;
    material.viscosity = 1.0e-3 * unit(rng);
  } else {
    material.rho = 2500.0 * unit(rng);
    material.lambda = 2e10 * unit(rng);
    material.mu = 3e10 * unit(rng);
    material.qp = 100.0 * unit(rng);
    if constexpr (HasShearQuality<MaterialT>::value) {
      setShearQuality(material, 50.0 * unit(rng));
    }
    material.initialize(parameters);
  }
  return material;
}

template <typename ArrayT>
double largest(const ArrayT& values) {
  double result = 0.0;
  for (const auto value : values) {
    result = std::max(result, std::abs(static_cast<double>(value)));
  }
  return result;
}

/// The source term the solver declares, written out entry by entry.
template <typename ViewT, typename ScalarsT>
void writeSource(ViewT& view, const ScalarsT& scalars) {
  using Material = seissol::model::MaterialT;
  using Setup = seissol::model::SolverSetup<typename Material::Solver, Material>;
  Setup::forEachSourceCoefficientEntry([&](std::size_t coefficient, auto... rest) {
    const std::array<double, sizeof...(rest)> arguments{static_cast<double>(rest)...};
    constexpr std::size_t Indices = sizeof...(rest) - 1;
    const double weight = arguments[Indices] * scalars.at(coefficient);
    if constexpr (Indices == 2) {
      view(static_cast<std::size_t>(arguments[0]), static_cast<std::size_t>(arguments[1])) +=
          weight;
    } else {
      view(static_cast<std::size_t>(arguments[0]),
           static_cast<std::size_t>(arguments[1]),
           static_cast<std::size_t>(arguments[2])) += weight;
    }
  });
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
  // the quality factors are what the attenuation fit turns into theta, so a
  // material that leaves them at zero has no relaxation to speak of
  material.qp = 100.0;
  material.qs = 50.0;
  for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
    for (std::size_t component = 0; component < 3; ++component) {
      material.theta[mech][component] = value(rng);
    }
  }
  return material;
}

#ifdef SEISSOL_KERNELS_STP
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
#endif

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

#ifdef SEISSOL_KERNELS_STP
  // the poroelastic declaration only exists where its solver is built
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

  SUBCASE("poroelastic source") {
    // the Biot drag, which a cell carries as two scalars where the material
    // varies inside it
    for (std::size_t sample = 0; sample < 32; ++sample) {
      coefficients::checkSourceDeclaration<seissol::model::PoroElasticMaterial, 13>(
          coefficients::poroelastic(rng), 0);
    }
  }
#endif

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
        seissol::model::SolverSetup<typename Material::Solver,
                                    Material>::getTransposedCoefficientMatrix(material,
                                                                              dim,
                                                                              reference);

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

/// The declaration says a coefficient is either a field or one number for the
/// whole run. That is a claim about the model, and this holds it: a coefficient
/// marked Global must not move when the material does, and one marked Material
/// has to be reachable from the material at all.
TEST_CASE("Coefficient origins") {
  auto rng = std::mt19937(20260926);

  const auto parameters = [] {
    seissol::initializer::parameters::ModelParameters p{};
    p.freqCentral = 1.0;
    p.freqRatio = 100.0;
    return p;
  }();

  const auto check = [&](auto prototype) {
    using Material = decltype(prototype);
    using Setup = seissol::model::SolverSetup<typename Material::Solver, Material>;
    constexpr auto Origins = Setup::CoefficientOrigins;

    auto base = prototype;
    base.initialize(parameters);
    const auto reference = Setup::getCoefficients(base);
    REQUIRE(Origins.size() == reference.size());

    // how far a coefficient moves over the whole sweep, so that a Material one
    // can be shown to depend on the material at all
    std::array<double, Setup::NumCoefficients> movement{};

    for (const auto& [name, member] : Material::ParameterMap) {
      auto perturbed = prototype;
      // a factor rather than an offset, so every parameter is probed on its
      // own scale
      perturbed.*member *= 1.5;
      perturbed.initialize(parameters);
      const auto moved = Setup::getCoefficients(perturbed);

      for (std::size_t i = 0; i < reference.size(); ++i) {
        const double scale = std::max(std::abs(reference[i]), std::abs(moved[i]));
        const double relative = scale > 0.0 ? std::abs(moved[i] - reference[i]) / scale : 0.0;
        movement[i] = std::max(movement[i], relative);

        if (Origins[i] == seissol::model::CoefficientOrigin::Global) {
          INFO("coefficient " << i << " is declared Global but " << name << " moves it");
          REQUIRE(relative <= 1e-14);
        }
      }
    }

    for (std::size_t i = 0; i < reference.size(); ++i) {
      if (Origins[i] == seissol::model::CoefficientOrigin::Material) {
        INFO("coefficient " << i << " is declared Material but no parameter reaches it");
        REQUIRE(movement[i] > 0.0);
      }
    }
  };

  /// The same question for the source term, where a material has one. Every
  /// scalar it is stated in has to be reachable from the material, or a cell
  /// that samples the material inside itself would be sampling a constant.
  /// The relaxation frequencies are the counterpart: they follow the frequency
  /// band alone, which is why they stay one number for the whole domain.
  const auto checkSource = [&](auto prototype) {
    using Material = decltype(prototype);
    using Setup = seissol::model::MaterialSetup<Material>;
    constexpr std::size_t Mechanisms = std::max<std::size_t>(Material::Mechanisms, 1);

    auto base = prototype;
    base.initialize(parameters);

    std::array<std::array<double, Setup::NumSourceCoefficients>, Mechanisms> movement{};
    double frequencyMovement = 0.0;

    for (const auto& [name, member] : Material::ParameterMap) {
      auto perturbed = prototype;
      perturbed.*member *= 1.5;
      perturbed.initialize(parameters);

      for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
        const auto reference = Setup::getSourceCoefficients(base, mech);
        const auto moved = Setup::getSourceCoefficients(perturbed, mech);
        for (std::size_t i = 0; i < reference.size(); ++i) {
          const double scale = std::max(std::abs(reference[i]), std::abs(moved[i]));
          if (scale > 0.0) {
            movement[mech][i] =
                std::max(movement[mech][i], std::abs(moved[i] - reference[i]) / scale);
          }
        }
        if constexpr (Material::Mechanisms > 0) {
          const double scale =
              std::max(std::abs(base.omega[mech]), std::abs(perturbed.omega[mech]));
          if (scale > 0.0) {
            frequencyMovement = std::max(
                frequencyMovement, std::abs(perturbed.omega[mech] - base.omega[mech]) / scale);
          }
        }
      }
    }

    for (std::size_t mech = 0; mech < Mechanisms; ++mech) {
      for (std::size_t i = 0; i < Setup::NumSourceCoefficients; ++i) {
        INFO("source coefficient " << i << " of mechanism " << mech
                                   << " is not reached by any material parameter");
        REQUIRE(movement[mech][i] > 0.0);
      }
    }
    INFO("a relaxation frequency moves with the material");
    REQUIRE(frequencyMovement <= 1e-14);
  };

  SUBCASE("elastic") { check(coefficients::elastic(2700.0, 3.24e10, 3.24e10)); }
  SUBCASE("acoustic") { check(coefficients::acoustic(1000.0, 2.25e9)); }
  SUBCASE("viscoelastic") {
    check(coefficients::viscoelastic<3>(rng));
    checkSource(coefficients::viscoelastic<3>(rng));
  }
#ifdef SEISSOL_KERNELS_STP
  SUBCASE("poroelastic") { checkSource(coefficients::poroelastic(rng)); }
#endif
}

/// The flux operator is not linear in either material -- the Riemann solver is
/// not -- but it occupies few entries and spans few dimensions, so it decomposes
/// the same way the star does: scalars read off a computed operator, times fixed
/// entries. This holds that, including where one side is acoustic and where
/// there is no other side at all.
TEST_CASE("Flux decomposition") {
  using Material = seissol::model::ElasticMaterial;
  // the operator is stated over the quantities of the Riemann problem and the
  // columns the star writes, which are not the same count where a solver
  // carries the relaxation in the same matrix
  constexpr std::size_t N = seissol::tensor::QgodLocal::Shape[0];
  constexpr std::size_t Columns = seissol::tensor::star::Shape[0][1];
  using Matrix = Eigen::Matrix<double, N, Columns>;
  using Square = Eigen::Matrix<double, N, N>;

  auto rng = std::mt19937(20260926);
  auto positive = std::uniform_real_distribution<double>(0.4, 2.5);
  const auto draw = [&](bool acoustic) {
    Material material{};
    material.rho = positive(rng);
    material.lambda = positive(rng);
    material.mu = acoustic ? 0.0 : positive(rng);
    return material;
  };

  const auto fluxOperator = [](const Material& local,
                               const Material& neighbor,
                               bool plus,
                               seissol::FaceType faceType) {
    alignas(Alignment) std::array<double, seissol::tensor::QgodLocal::size()> localData{};
    alignas(Alignment) std::array<double, seissol::tensor::QgodNeighbor::size()> neighborData{};
    auto godLocal = seissol::init::QgodLocal::view::create(localData.data());
    auto godNeighbor = seissol::init::QgodNeighbor::view::create(neighborData.data());
    seissol::model::getTransposedGodunovState(local, neighbor, faceType, godLocal, godNeighbor);

    Matrix coefficientMatrix = Matrix::Zero();
    seissol::model::getTransposedCoefficientMatrix(plus ? local : neighbor, 0, coefficientMatrix);
    Square godunov = Square::Zero();
    auto& view = plus ? godLocal : godNeighbor;
    for (std::size_t row = 0; row < N; ++row) {
      for (std::size_t column = 0; column < N; ++column) {
        if (view.isInRange(row, column)) {
          godunov(row, column) = view(row, column);
        }
      }
    }
    return Matrix(godunov * coefficientMatrix);
  };

  const auto check = [&](bool acousticLocal, bool acousticNeighbor, seissol::FaceType faceType) {
    for (std::size_t sample = 0; sample < 64; ++sample) {
      const auto local = draw(acousticLocal);
      // a face without another cell behind it is handed the cell's own
      // material, the way the initialization does it
      const auto neighbor = faceType == seissol::FaceType::Regular ? draw(acousticNeighbor) : local;
      // a face without another cell behind it has its neighbor operator
      // poisoned with a signalling NaN on purpose, so there is nothing to
      // decompose there
      const bool hasNeighbor = faceType == seissol::FaceType::Regular;
      for (const bool plus :
           hasNeighbor ? std::vector<bool>{true, false} : std::vector<bool>{true}) {
        const Matrix reference = fluxOperator(local, neighbor, plus, faceType);

        std::array<double, seissol::generated::FluxNumCoefficients> coefficients{};
        for (std::size_t a = 0; a < coefficients.size(); ++a) {
          const auto& source = seissol::generated::FluxCoefficientSources[a];
          coefficients[a] = reference(source.row, source.column);
        }

        Matrix candidate = Matrix::Zero();
        for (const auto& entry : seissol::generated::FluxCoefficientEntries) {
          candidate(entry.row, entry.column) += entry.factor * coefficients[entry.coefficient];
        }

        const double scale = std::max(1.0, reference.cwiseAbs().maxCoeff());
        for (std::size_t row = 0; row < N; ++row) {
          for (std::size_t column = 0; column < Columns; ++column) {
            REQUIRE(candidate(row, column) ==
                    doctest::Approx(reference(row, column)).epsilon(1e-14).scale(scale));
          }
        }
      }
    }
  };

  SUBCASE("elastic against elastic") { check(false, false, seissol::FaceType::Regular); }
  SUBCASE("acoustic against elastic") { check(true, false, seissol::FaceType::Regular); }
  SUBCASE("elastic against acoustic") { check(false, true, seissol::FaceType::Regular); }
  SUBCASE("acoustic against acoustic") { check(true, true, seissol::FaceType::Regular); }
  SUBCASE("free surface") { check(false, false, seissol::FaceType::FreeSurface); }
  SUBCASE("free surface, acoustic") { check(true, true, seissol::FaceType::FreeSurface); }
}

/// The two cells sharing a face parametrise it differently, and the generated
/// fP carries the map between the two in the modal basis with a mass factor.
/// Strip the factor and go to the nodes and it is a renumbering, nothing more --
/// which is what a value given at the nodes needs, and what a flux built from
/// the material of both sides is made of.
TEST_CASE("Face orientation renumbering") {
  constexpr auto Nodes = seissol::generated::FaceNodes;
  using Matrix = Eigen::Matrix<double, Nodes, Nodes>;

  // Every matrix here comes from a matrix file, and a build that bundles
  // simulations stores those the other way round; the identity below is the
  // mathematical one.
  const auto dense = [](auto view, std::size_t rows, std::size_t columns) {
    constexpr bool Flip = multisim::NumSimulations > 1;
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(rows, columns);
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t column = 0; column < columns; ++column) {
        const std::size_t first = Flip ? column : row;
        const std::size_t second = Flip ? row : column;
        if (view.isInRange(first, second)) {
          matrix(row, column) = view(first, second);
        }
      }
    }
    return matrix;
  };

  const auto nodalToModal = dense(seissol::nodal::init::MV2nTo2m::view::create(
                                      const_cast<double*>(seissol::nodal::init::MV2nTo2m::Values)),
                                  Nodes,
                                  Nodes);
  // the way back is its inverse; only the one direction is generated
  const Eigen::MatrixXd modalToNodal = nodalToModal.inverse();
  const auto massInverse =
      dense(seissol::init::M2inv::view::create(const_cast<double*>(seissol::init::M2inv::Values)),
            Nodes,
            Nodes);

  REQUIRE(seissol::generated::FaceOrientations == 3);
  const auto check = [&](auto tag) {
    constexpr unsigned Orientation = decltype(tag)::value;
    const auto facePermutation =
        dense(seissol::init::fP::view<Orientation>::create(
                  const_cast<double*>(seissol::init::fP::Values[Orientation])),
              Nodes,
              Nodes);
    const std::size_t orientation = Orientation;

    // the same map without the mass factor, taken to the nodes
    const Eigen::MatrixXd atNodes = modalToNodal * (massInverse * facePermutation) * nodalToModal;

    Matrix expected = Matrix::Zero();
    const auto& renumbering = seissol::generated::FaceOrientationPermutations[orientation];
    for (std::size_t node = 0; node < Nodes; ++node) {
      expected(node, renumbering[node]) = 1.0;
    }

    for (std::size_t row = 0; row < Nodes; ++row) {
      for (std::size_t column = 0; column < Nodes; ++column) {
        INFO("orientation " << orientation << " at (" << row << "," << column << ")");
        REQUIRE(atNodes(row, column) == doctest::Approx(expected(row, column)).epsilon(1e-10));
      }
    }
  };
  check(std::integral_constant<unsigned, 0>{});
  check(std::integral_constant<unsigned, 1>{});
  check(std::integral_constant<unsigned, 2>{});
}

/// The source term of the configured solver, put together from the scalars a
/// cell carries and the entries the generator wrote into the kernel, against
/// the tensor the solver builds for itself. The two are composed the same way
/// -- one block per relaxation mechanism, wherever that solver puts it -- so a
/// block placed at the wrong offset, or a relaxation frequency counted twice,
/// shows up here.
TEST_CASE("Source assembly from coefficients") {
  using Material = seissol::model::MaterialT;
  using Setup = seissol::model::SolverSetup<typename Material::Solver, Material>;

  if constexpr (Setup::NumSourceCoefficients > 0) {
    std::mt19937 rng(20260927);

    for (std::size_t sample = 0; sample < 16; ++sample) {
      const auto material = coefficients::configuredMaterial<Material>(rng);

      std::array<real, coefficients::SourceTensor::size()> referenceData{};
      auto reference = coefficients::SourceInit::view::create(referenceData.data());
      Setup::getTransposedSourceCoefficientTensor(material, reference);

      std::array<real, coefficients::SourceTensor::size()> candidateData{};
      auto candidate = coefficients::SourceInit::view::create(candidateData.data());
      candidate.setZero();
      const auto scalars = Setup::getSourceCoefficients(material);
      coefficients::writeSource(candidate, scalars);

      const double scale = std::max(1.0, coefficients::largest(referenceData));
      for (std::size_t i = 0; i < referenceData.size(); ++i) {
        REQUIRE(candidateData.at(i) ==
                doctest::Approx(referenceData.at(i)).epsilon(1e-13).scale(scale));
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
      const auto material =
          coefficients::elastic(positive(rng) * 1000.0, positive(rng) * 1e10, positive(rng) * 1e10);
      const double gradient[3] = {gauss(rng), gauss(rng), gauss(rng)};

      // what CellLocalMatrices builds today: the three directional matrices,
      // each scaled by its row of the Jacobian
      std::array<std::array<double, seissol::tensor::star::size(0)>, 3> directional{};
      for (unsigned dim = 0; dim < 3; ++dim) {
        auto view = seissol::init::star::view<0>::create(directional.at(dim).data());
        seissol::model::MaterialSetup<
            seissol::model::ElasticMaterial>::getTransposedCoefficientMatrix(material, dim, view);
      }

      std::array<double, seissol::tensor::star::size(0)> assembledData{};
      auto view = seissol::init::star::view<0>::create(assembledData.data());
      const auto coefficients = seissol::model::getStarCoefficients(material);
      seissol::model::assembleStarMatrix<seissol::model::ElasticMaterial>(
          coefficients.data(), gradient, view);

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

TEST_CASE("Face rotation structure") {
  // The face rotation is stored by its pattern, and its inverse is a function of
  // the forward matrix rather than an independent one. Both are properties of
  // how the rotation is built, so they are checked here: a flux that wants the
  // rotation per face node can rely on them and keep one matrix instead of two.
  using Material = seissol::model::MaterialT;
  constexpr std::size_t N = Material::NumQuantities;

  // A face takes an arbitrary orthonormal frame -- the normal is a face normal
  // and the first tangent one of its edges -- so nothing may lean on a
  // particular choice of tangents.
  std::mt19937 rng(20260927);
  std::normal_distribution<double> gauss(0.0, 1.0);

  for (std::size_t sample = 0; sample < 64; ++sample) {
    Eigen::Matrix3d gaussian;
    for (unsigned row = 0; row < 3; ++row) {
      for (unsigned column = 0; column < 3; ++column) {
        gaussian(row, column) = gauss(rng);
      }
    }
    Eigen::Matrix3d frame = Eigen::HouseholderQR<Eigen::Matrix3d>(gaussian).householderQ();
    if (frame.determinant() < 0.0) {
      frame.col(2) *= -1.0;
    }
    const VrtxCoords normal{frame(0, 0), frame(1, 0), frame(2, 0)};
    const VrtxCoords tangent1{frame(0, 1), frame(1, 1), frame(2, 1)};
    const VrtxCoords tangent2{frame(0, 2), frame(1, 2), frame(2, 2)};

    alignas(Alignment) std::array<double, seissol::tensor::T::size()> forwardData{};
    alignas(Alignment) std::array<double, seissol::tensor::Tinv::size()> inverseData{};
    auto forwardView = seissol::init::T::view::create(forwardData.data());
    auto inverseView = seissol::init::Tinv::view::create(inverseData.data());
    seissol::model::getFaceRotationMatrix<Material>(
        normal, tangent1, tangent2, forwardView, inverseView);

    const auto densify = [](auto& view, std::size_t rows, std::size_t columns) {
      Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(rows, columns);
      for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < columns; ++column) {
          if (view.isInRange(row, column)) {
            matrix(row, column) = view(row, column);
          }
        }
      }
      return matrix;
    };
    const Eigen::MatrixXd forward =
        densify(forwardView, forwardView.shape(0), forwardView.shape(1));
    const Eigen::MatrixXd inverse =
        densify(inverseView, inverseView.shape(0), inverseView.shape(1));

    // One dense block per quantity group, and nothing outside them -- not
    // merely numerically zero, but absent from the storage.
    std::size_t offset = 0;
    for (const auto& group : Material::RotationGroups) {
      const std::size_t extent = group.extent();
      for (std::size_t row = 0; row < forwardView.shape(0); ++row) {
        for (std::size_t column = offset; column < offset + extent; ++column) {
          const bool inBlock = row >= offset && row < offset + extent;
          REQUIRE(forwardView.isInRange(row, column) == inBlock);
        }
      }
      offset += extent;
    }

    // The inverse undoes the forward rotation on the quantities it spans. Where
    // the two span different sets -- a solver that rotates one anelastic block
    // forwards and none back -- that is the leading square of the forward one.
    const auto inverted = static_cast<std::size_t>(inverseView.shape(0));
    const Eigen::MatrixXd product = inverse * forward.topLeftCorner(inverted, inverted);
    for (std::size_t row = 0; row < inverted; ++row) {
      for (std::size_t column = 0; column < inverted; ++column) {
        const double expected = (row == column) ? 1.0 : 0.0;
        REQUIRE(std::abs(product(row, column) - expected) < 1e-12);
      }
    }

    // Per group, the inverse follows from the forward matrix alone: a rotated
    // vector is orthogonal, so its inverse is its transpose, and a symmetric
    // second-order tensor in Voigt form differs from the transpose by the
    // weights of its shear components.
    offset = 0;
    for (const auto& group : Material::InverseRotationGroups) {
      const std::size_t extent = group.extent();
      const Eigen::MatrixXd block = forward.block(offset, offset, extent, extent);
      const Eigen::MatrixXd inverseBlock = inverse.block(offset, offset, extent, extent);
      Eigen::VectorXd weights = Eigen::VectorXd::Ones(extent);
      if (group.kind == seissol::model::QuantityKind::SymTensor2) {
        for (std::size_t component = 3; component < extent; ++component) {
          weights(component) = 0.5;
        }
      }
      const Eigen::MatrixXd derived =
          weights.asDiagonal() * block.transpose() * weights.cwiseInverse().asDiagonal();
      for (std::size_t row = 0; row < extent; ++row) {
        for (std::size_t column = 0; column < extent; ++column) {
          // the two are built from the same products of the same frame, so the
          // agreement is exact rather than merely close
          REQUIRE(derived(row, column) == inverseBlock(row, column));
        }
      }
      offset += extent;
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_COEFFICIENTSTRUCTURE_T_H_
