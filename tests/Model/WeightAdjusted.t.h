// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_WEIGHTADJUSTED_T_H_
#define SEISSOL_TESTS_MODEL_WEIGHTADJUSTED_T_H_

// A weight-adjusted cell lifts its face terms without material into a buffer of its own and
// applies the inverse weight of its medium at the sample points once they are all in
// (WeightAdjustable). Checked here, besides the inverse weight itself (InverseWeight.t.h), which
// makes the two lifts the same where the material does not vary:
//
// - the flux without material decomposes into the scalars of a face;
// - the kernel that applies the weight, against the matrices it is made of, where the material
//   differs from sample to sample;
// - with one material for the cell, what a weight-adjusted cell adds -- the volume term, the
//   correctors of its faces, the lift of a fault face -- is what a cell that lifts with the
//   material adds, to the field or to the extended field of the solver that keeps the memory
//   variables apart;
// - which cells are weight-adjusted, and how much shorter their time step is.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Common/Constants.h"
#include "DynamicRupture/Typedefs.h"
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "FaceTraces.h"
#include "FaultFlux.t.h" // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Model/FaultFlux.h"
#include "Initializer/Model/WeightAdjusted.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Typedefs.h"
#include "InverseWeight.t.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"
#include "NodalCorrector.t.h" // IWYU pragma: keep
#include "Solver/MultipleSimulations.h"
#include "VolumeTerm.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

namespace weightadjusted {

/// The operator of a face without material -- the Godunov state times the coefficient matrix
/// without material -- is made of the scalars per point the flux of a face carries, as the one
/// with the material is (see "Flux decomposition"); a weight-adjusted face fills fewer of them.
template <bool Enabled>
void unweightedFluxDecomposes() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t N = tensor::QgodLocal::Shape[0];
    constexpr std::size_t Columns = tensor::star::Shape[0][1];
    using Matrix = Eigen::Matrix<double, N, Columns>;

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261019);
    const auto check = [&](bool acousticLocal, bool acousticNeighbor, FaceType faceType) {
      for (std::size_t sample = 0; sample < 32; ++sample) {
        const auto local = coefficients::configuredMaterial<Material>(rng, acousticLocal);
        const auto neighbor =
            faceType == FaceType::Regular
                ? coefficients::configuredMaterial<Material>(rng, acousticNeighbor)
                : local;
        alignas(Alignment) std::array<real, tensor::QgodLocal::size()> localData{};
        alignas(Alignment) std::array<real, tensor::QgodNeighbor::size()> neighborData{};
        auto godLocal = init::QgodLocal::view::create(localData.data());
        auto godNeighbor = init::QgodNeighbor::view::create(neighborData.data());
        seissol::model::getTransposedGodunovState(local, neighbor, faceType, godLocal, godNeighbor);
        // both sides contract with the matrix of the local material, see fluxScalarsOfNode
        Matrix unweighted = Matrix::Zero();
        seissol::model::getTransposedUnweightedCoefficientMatrix(local, 0, unweighted);

        // a face without another cell behind it has no neighbour operator to decompose: the
        // setup poisons it on purpose
        const bool hasNeighbor = faceType == FaceType::Regular;
        for (const bool plus :
             hasNeighbor ? std::vector<bool>{true, false} : std::vector<bool>{true}) {
          Matrix reference = Matrix::Zero();
          for (std::size_t row = 0; row < N; ++row) {
            for (std::size_t k = 0; k < N; ++k) {
              const auto& view = plus ? godLocal : godNeighbor;
              if (view.isInRange(row, k)) {
                reference.row(row) += view(row, k) * unweighted.row(k);
              }
            }
          }
          std::array<double, generated::FluxNumCoefficients> scalars{};
          for (std::size_t a = 0; a < scalars.size(); ++a) {
            const auto& source = generated::FluxCoefficientSources[a];
            scalars[a] = reference(source.row, source.column);
          }
          Matrix candidate = Matrix::Zero();
          for (const auto& entry : generated::FluxCoefficientEntries) {
            candidate(entry.row, entry.column) += entry.factor * scalars[entry.coefficient];
          }
          const double scale = std::max(1.0, reference.cwiseAbs().maxCoeff());
          for (std::size_t row = 0; row < N; ++row) {
            for (std::size_t column = 0; column < Columns; ++column) {
              REQUIRE(candidate(row, column) == doctest::Approx(reference(row, column))
                                                    .epsilon(coefficients::tolerance(1e-14))
                                                    .scale(scale));
            }
          }
        }
      }
    };
    check(false, false, FaceType::Regular);
    check(true, false, FaceType::Regular);
    check(false, true, FaceType::Regular);
    check(true, true, FaceType::Regular);
    check(false, false, FaceType::FreeSurface);
    check(true, true, FaceType::FreeSurface);
  }
}

} // namespace weightadjusted

TEST_CASE("The flux without material decomposes into the scalars of a face") {
  weightadjusted::unweightedFluxDecomposes<WeightAdjustable && NodalFlux>();
}

TEST_CASE("Which cells are weight-adjusted, and their time step") {
  using seissol::initializer::ImpedanceSpread;
  using seissol::initializer::isWeightAdjusted;
  using seissol::initializer::weightAdjustedTimestepScale;
  const auto spreadOf = [](const auto& materials) {
    ImpedanceSpread spread;
    for (const auto& material : materials) {
      spread.add(material);
    }
    return spread;
  };

  SUBCASE("one material") {
    const auto material = coefficients::elastic(2500.0, 3e10, 2e10);
    const auto spread = spreadOf(std::vector{material, material, material});
    REQUIRE(spread.ratio() == doctest::Approx(1.0));
    REQUIRE(spread.timestepFactor() == doctest::Approx(1.0));
    REQUIRE_FALSE(isWeightAdjusted(spread.ratio(), 1.0, 2.0, 1.5));
  }
  SUBCASE("a jump in the bulk modulus by 25") {
    // the acoustic case of the prototype: an impedance ratio of 5, the factor c / 5
    const auto spread =
        spreadOf(std::vector{coefficients::acoustic(1.0, 1.0), coefficients::acoustic(1.0, 25.0)});
    REQUIRE(spread.ratio() == doctest::Approx(5.0));
    for (std::size_t order = 2; order <= 8; ++order) {
      REQUIRE(spread.timestepFactor(order) ==
              doctest::Approx(weightAdjustedTimestepScale(order) / 5.0));
    }
    REQUIRE(spread.timestepFactor() == spread.timestepFactor(ConvergenceOrder));
    // a jump that accounts for a ratio of 3 of it, and a smooth spread that leaves 1.2 to the jump
    REQUIRE(isWeightAdjusted(spread.ratio(), 3.0, 2.0, 1.5));
    REQUIRE_FALSE(isWeightAdjusted(spread.ratio(), 3.0, 5.0, 1.5));
    REQUIRE_FALSE(isWeightAdjusted(spread.ratio(), 3.0, 0.0, 1.5));
    REQUIRE_FALSE(isWeightAdjusted(spread.ratio(), 1.2, 2.0, 1.5));
    // without the jump test, every cell that spreads by more than the ratio
    REQUIRE(isWeightAdjusted(spread.ratio(), 1.2, 2.0, 1.0));
    REQUIRE(isWeightAdjusted(spread.ratio(), 1.0, 2.0, 0.0));
  }
  SUBCASE("what the orders permit") {
    // orders 4 and 6 permit the least of those measured, and an order not measured no more
    for (std::size_t order = 2; order <= 8; ++order) {
      REQUIRE(weightAdjustedTimestepScale(order) >= weightAdjustedTimestepScale(4));
    }
    REQUIRE(weightAdjustedTimestepScale(4) == 1.0);
    REQUIRE(weightAdjustedTimestepScale(6) == weightAdjustedTimestepScale(4));
    for (const std::size_t order : {2, 7, 8}) {
      REQUIRE(weightAdjustedTimestepScale(order) == weightAdjustedTimestepScale(4));
    }
    // never longer than the one the wave speeds permit: a spread of 1.5 at order 3
    const auto spread =
        spreadOf(std::vector{coefficients::acoustic(1.0, 1.0), coefficients::acoustic(1.0, 2.25)});
    REQUIRE(spread.ratio() == doctest::Approx(1.5));
    REQUIRE(spread.timestepFactor(3) == 1.0);
  }
  SUBCASE("the shear waves count with their speed") {
    // the P impedances agree and the S impedances differ by 5, the faster S wave half as fast as
    // the P waves: c / 5 for the S waves, weighed with that half
    const auto soft = coefficients::elastic(1.0, 1.0, 98.0);
    const auto hard = coefficients::elastic(1.0, 25.0, 50.0);
    const auto spread = spreadOf(std::vector{soft, hard});
    REQUIRE(spread.ratio() == doctest::Approx(5.0));
    REQUIRE(spread.timestepFactor(4) == doctest::Approx(1.0 / (0.5 * 5.0)));
  }
  SUBCASE("a point without shear has no shear impedance") {
    const auto solid = coefficients::elastic(2500.0, 3e10, 2e10);
    const auto fluid = coefficients::elastic(1000.0, 0.0, 2.25e9);
    const auto spread = spreadOf(std::vector{solid, fluid});
    const double pRatio =
        std::sqrt((solid.lambda + 2 * solid.mu) * solid.rho / (fluid.lambda * fluid.rho));
    REQUIRE(spread.ratio() == doctest::Approx(pRatio));
  }
}

TEST_CASE("A jump inside a cell, against a material that varies smoothly") {
  using seissol::initializer::FacePointGenerator;
  using seissol::initializer::ImpedanceJump;
  const auto nodes = init::materialNodes::view::create(init::materialNodes::Values);
  // the points of the faces the material is read at where the flux varies along a face, in
  // reference coordinates
  const auto facePoints = FacePointGenerator::fluxPoints();
  const auto facePointAt = [&](std::size_t point) {
    return FacePointGenerator::referencePoint(
        facePoints, point / facePoints.count, point % facePoints.count);
  };
  const auto faceBasis = ImpedanceJump::faceBasis(facePointAt);
  // the material `materialAt` gives at every sample of the cell and at those points, in reference
  // coordinates
  const auto jumpOf = [&](const auto& materialAt) {
    ImpedanceJump jump(faceBasis);
    for (std::size_t sample = 0; sample < ImpedanceJump::Samples; ++sample) {
      jump.add(sample, materialAt(nodes(sample, 0), nodes(sample, 1), nodes(sample, 2)));
    }
    for (std::size_t point = 0; point < ImpedanceJump::FacePoints; ++point) {
      const auto xi = facePointAt(point);
      jump.addFace(point, materialAt(xi[0], xi[1], xi[2]));
    }
    return jump.ratio();
  };
  // the impedances scaled by `factor`: the moduli by its square, the density kept
  const auto scaled = [](double factor) {
    return coefficients::elastic(2500.0, 3e10 * factor * factor, 2e10 * factor * factor);
  };

  SUBCASE("the modes at the points of the faces") {
    // the basis the samples are fitted in, evaluated at the samples, is the one the projection
    // reads them with (init::materialEval)
    const auto evaluate = init::materialEval::view::create(init::materialEval::Values);
    for (std::size_t sample = 0; sample < ImpedanceJump::Samples; ++sample) {
      const auto basis = ImpedanceJump::faceBasis([&](std::size_t /*point*/) {
        return std::array<double, Cell::Dim>{nodes(sample, 0), nodes(sample, 1), nodes(sample, 2)};
      });
      for (std::size_t point = 0; point < ImpedanceJump::FacePoints; ++point) {
        for (std::size_t mode = 0; mode < ImpedanceJump::FitModes; ++mode) {
          const double expected = evaluate.isInRange(sample, mode) ? evaluate(sample, mode) : 0.0;
          REQUIRE(basis[point][mode] == doctest::Approx(expected));
        }
      }
    }
  }
  SUBCASE("a smooth material") {
    // impedances that grow exponentially along a direction, by a factor of 49 over the cell, and
    // with the square of a coordinate, by 55: their logarithms are a polynomial of degree one,
    // which the projection follows exactly from order 3 on, and of degree two, from order 4 on --
    // at the samples and at the points of the faces alike
    if constexpr (ConvergenceOrder >= 3) {
      REQUIRE(jumpOf([&](double x, double y, double z) {
                return scaled(std::exp(2.3 * x - 1.6 * y + 0.9 * z));
              }) == doctest::Approx(1.0));
    }
    if constexpr (ConvergenceOrder >= 4) {
      REQUIRE(jumpOf([&](double x, double y, double z) {
                return scaled(std::exp(4.0 * (x + y + z) * (x + y + z)));
              }) == doctest::Approx(1.0));
    }
  }
  SUBCASE("a jump") {
    // by 10 across a plane through the cell, of which every order and both sets of samples keep
    // more than a factor of five, with the points of the faces or without
    REQUIRE(jumpOf([&](double x, double y, double z) {
              return scaled(x + y + z < 0.5 ? 1.0 : 10.0);
            }) > 5.0);
    // the same in the shear modulus alone, the P impedance kept
    REQUIRE(jumpOf([&](double x, double y, double z) {
              const double mu = x + y + z < 0.5 ? 4e8 : 4e10;
              return coefficients::elastic(2500.0, mu, 1.2e11 - 2.0 * mu);
            }) > 5.0);
  }
  SUBCASE("a jump next to a face") {
    // between the face zeta = 0 and the samples nearest to it, which only the points of that face
    // see, pulled into the cell by far less: what the samples give misses it, the points of the
    // face do not, where the flux reads the material there
    double nearest = 1.0;
    for (std::size_t sample = 0; sample < ImpedanceJump::Samples; ++sample) {
      nearest = std::min(nearest, static_cast<double>(nodes(sample, 2)));
    }
    const double ratio = jumpOf([&](double /*x*/, double /*y*/, double z) {
      return scaled(z < 0.5 * nearest ? 10.0 : 1.0);
    });
    if constexpr (ImpedanceJump::FacePoints > 0) {
      REQUIRE(ratio > 5.0);
    } else {
      REQUIRE(ratio == doctest::Approx(1.0));
    }
  }
  SUBCASE("a fluid among the samples") {
    // where the shear modulus vanishes at some samples only, it jumps, by any ratio
    REQUIRE(std::isinf(jumpOf([&](double x, double y, double z) {
      return x + y + z < 0.5 ? coefficients::elastic(1000.0, 0.0, 2.25e9)
                             : coefficients::elastic(2500.0, 3e10, 2e10);
    })));
    // and a fluid throughout has no shear impedance to fit
    if constexpr (ConvergenceOrder >= 3) {
      REQUIRE(jumpOf([&](double x, double y, double /*z*/) {
                return coefficients::elastic(1000.0, 0.0, 2.25e9 * std::exp(x - y));
              }) == doctest::Approx(1.0));
    }
  }
}

namespace weightadjusted {

/// What the volume term and the faces of a cell go to, and what the weight of a weight-adjusted
/// cell comes back to: the field, or the extended field of the solver that keeps the memory
/// variables apart, which hands both on to the field and the memory variables alike.
struct Lifted {
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
  using Tensor = tensor::Qext;
  using Init = init::Qext;
  template <typename KernelT>
  static void bind(KernelT& krnl, real* target) {
    krnl.Qext = target;
  }
  template <typename LocalT>
  static void addVolume(real* timeIntegrated, real* target, const LocalT& local) {
    kernel::volumeExt volume{};
    volume.bindGlobals(seissol::Pool::host());
    volume.I = timeIntegrated;
    volume.Qext = target;
    kernels::bindStarOperands(volume, local);
    volume.execute();
  }
#else
  using Tensor = tensor::Q;
  using Init = init::Q;
  template <typename KernelT>
  static void bind(KernelT& krnl, real* target) {
    krnl.Q = target;
  }
  template <typename LocalT>
  static void addVolume(real* timeIntegrated, real* target, const LocalT& local) {
    addVolumeTerm<kernel::volume>(timeIntegrated, target, local);
  }
#endif
  static constexpr std::size_t Basis = Tensor::Shape[multisim::BasisFunctionDimension];
  static constexpr std::size_t Quantities = Tensor::Shape[multisim::BasisFunctionDimension + 1];
};

/// what a lifted field holds for one fused simulation, as a matrix of modes by quantities
template <typename ViewT>
Eigen::MatrixXd fieldOf(ViewT view, std::size_t simulation) {
  auto sliced = multisim::simtensor(view, simulation);
  Eigen::MatrixXd field = Eigen::MatrixXd::Zero(Lifted::Basis, Lifted::Quantities);
  for (std::size_t row = 0; row < Lifted::Basis; ++row) {
    for (std::size_t column = 0; column < Lifted::Quantities; ++column) {
      field(row, column) = sliced(row, column);
    }
  }
  return field;
}

/// Two updates of a cell agree, every quantity measured against the size of its group: the
/// quantities differ by orders of magnitude, a stress by the moduli and a velocity by the inverse
/// density, and within a group they mix as a face turns them. A fluid shows why the group and not
/// the column: the isotropic stress a face lifts in its own coordinates comes back with a shear
/// component of roundoff, which the weight, without shear modulus, makes zero.
inline void requireSameUpdate(real* candidate, real* reference) {
  using Material = seissol::model::MaterialT;
  constexpr std::size_t Quantities = Lifted::Quantities;
  std::array<std::size_t, Quantities> groupOf{};
  std::size_t groups = 0;
  std::size_t covered = 0;
  for (const auto& group : Material::RotationGroups) {
    for (std::size_t i = 0; i < group.extent() && covered < Quantities; ++i) {
      groupOf[covered++] = groups;
    }
    ++groups;
  }
  REQUIRE(covered == Quantities);

  auto viewCandidate = Lifted::Init::view::create(candidate);
  auto viewReference = Lifted::Init::view::create(reference);
  constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
  for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
    const auto got = fieldOf(viewCandidate, simulation);
    const auto expected = fieldOf(viewReference, simulation);
    std::vector<double> scaleOfGroup(groups, 1.0);
    for (std::size_t column = 0; column < Quantities; ++column) {
      scaleOfGroup[groupOf[column]] =
          std::max(scaleOfGroup[groupOf[column]], expected.col(column).cwiseAbs().maxCoeff());
    }
    for (std::size_t column = 0; column < Quantities; ++column) {
      for (Eigen::Index row = 0; row < expected.rows(); ++row) {
        REQUIRE(got(row, column) == doctest::Approx(expected(row, column))
                                        .epsilon(Tolerance)
                                        .scale(scaleOfGroup[groupOf[column]]));
      }
    }
  }
}

/// Adds the inverse weight applied to `buffer` to `target`, with the coefficients `local` carries
/// at the sample points, as a weight-adjusted cell does at the end of its neighbour integral.
template <typename KernelT>
void applyWeight(real* buffer, real* target, const LocalIntegrationData& local) {
  KernelT krnl{};
  krnl.bindGlobals(seissol::Pool::host());
  Lifted::bind(krnl, target);
  krnl.fluxBuffer = buffer;
  kernels::bindWeightOperands(krnl, local);
  krnl.execute();
}

template <typename T>
void fillRandom(T& values, std::mt19937& rng) {
  std::normal_distribution<double> gauss(0.0, 1.0);
  for (auto& value : values) {
    value = static_cast<real>(gauss(rng));
  }
}

/// The kernel against the matrices it is made of: the buffer read at the sample points, times the
/// inverse weight of the material of every point, projected back and added. The rows beyond the
/// quantities of the medium -- what a face lifts for the memory variables -- carry no material and
/// are added as they are.
template <bool Enabled>
void kernelAgainstMatrices() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    using KernelT = std::conditional_t<Enabled, kernel::weightAdjust, void>;
    constexpr std::size_t Points = tensor::materialEval::Shape[0];
    constexpr std::size_t Basis = tensor::materialEval::Shape[1];
    constexpr std::size_t Quantities = Lifted::Quantities;
    constexpr std::size_t Weighted = seissol::model::totalExtent(Material::PrimaryGroups);
    constexpr std::size_t Size = Lifted::Tensor::size();
    static_assert(Points == MaterialSampleCount);
    static_assert(!Enabled || Weighted <= Quantities);
    static_assert(!Enabled || kernels::size<tensor::fluxBuffer>() == Size);

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261016);
    const Eigen::MatrixXd evaluate = nodalcorrector::denseOf(
        init::materialEval::view::create(const_cast<real*>(init::materialEval::Values)),
        Points,
        Basis);
    const Eigen::MatrixXd project = nodalcorrector::denseOf(
        init::materialProject::view::create(const_cast<real*>(init::materialProject::Values)),
        Basis,
        Points);

    for (std::size_t sample = 0; sample < 3; ++sample) {
      // a material of its own at every point, a fluid among them where the medium has shear; and
      // where a cell may be curved, a 1 / J of its own at every point, which the weight carries
      LocalIntegrationData local{};
      std::vector<Eigen::MatrixXd> weights;
      std::uniform_real_distribution<double> inverseDeterminant(0.3, 3.0);
      for (std::size_t point = 0; point < Points; ++point) {
        const auto material = coefficients::configuredMaterial<Material>(rng, point % 7 == 3);
        const auto values = seissol::model::getStarCoefficients(material);
        for (std::size_t a = 0; a < values.size(); ++a) {
          local.materialCoefficients[a][point] = values[a];
        }
        weights.emplace_back(inverseweight::inverseWeight<Material, Weighted>(material));
        if constexpr (Curvilinear) {
          local.weightScale[point] = static_cast<real>(inverseDeterminant(rng));
          weights.back() *= static_cast<double>(local.weightScale[point]);
        }
      }

      alignas(Alignment) std::array<real, Size> buffer{};
      alignas(Alignment) std::array<real, Size> start{};
      fillRandom(buffer, rng);
      fillRandom(start, rng);
      alignas(Alignment) std::array<real, Size> update = start;
      applyWeight<KernelT>(buffer.data(), update.data(), local);

      alignas(Alignment) std::array<real, Size> expectedData{};
      auto expectedView = Lifted::Init::view::create(expectedData.data());
      for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
        const Eigen::MatrixXd lifted =
            fieldOf(Lifted::Init::view::create(buffer.data()), simulation);
        const Eigen::MatrixXd values = evaluate * lifted.leftCols(Weighted);
        Eigen::MatrixXd product = Eigen::MatrixXd::Zero(Points, Weighted);
        for (std::size_t point = 0; point < Points; ++point) {
          product.row(point) = values.row(point) * weights[point];
        }
        Eigen::MatrixXd expected = fieldOf(Lifted::Init::view::create(start.data()), simulation);
        expected.leftCols(Weighted) += project * product;
        expected.rightCols(Quantities - Weighted) += lifted.rightCols(Quantities - Weighted);
        auto sliced = multisim::simtensor(expectedView, simulation);
        for (std::size_t row = 0; row < Basis; ++row) {
          for (std::size_t column = 0; column < Quantities; ++column) {
            sliced(row, column) = expected(row, column);
          }
        }
      }
      requireSameUpdate(update.data(), expectedData.data());
    }
  }
}

/// With one material for the cell, what a weight-adjusted cell adds -- its volume term, and the
/// correctors its faces lift without material into its buffer, with the weight applied to them --
/// is what a cell adds that lifts with the material. Every face takes a path of its own: a regular
/// face both traces, a free surface its mirrored state, a fault face the subtraction of its own
/// trace.
template <bool Enabled>
void sameWhereConstant() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    using KernelT = std::conditional_t<Enabled, kernel::weightAdjust, void>;
    constexpr std::size_t Size = Lifted::Tensor::size();

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261017);
    const auto cell = nodalcorrector::testCell();
    const std::array<FaceType, Cell::NumFaces> faceTypes{
        FaceType::Regular, FaceType::FreeSurface, FaceType::DynamicRupture, FaceType::Regular};

    for (std::size_t sample = 0; sample < 3; ++sample) {
      // a fluid as well, where the medium has shear: its shear stresses reach no velocity
      const auto material = coefficients::configuredMaterial<Material>(rng, sample == 2);
      const auto field = [&](const nodalcorrector::VectorT& /*x*/) { return material; };

      LocalIntegrationData plainLocal{};
      LocalIntegrationData adjustedLocal{};
      std::array<std::array<real, tensor::AplusT::size()>, Cell::NumFaces> weakAplus{};
      nodalcorrector::setupCell<Material>(cell, field, faceTypes, plainLocal, weakAplus, false);
      nodalcorrector::setupCell<Material>(cell, field, faceTypes, adjustedLocal, weakAplus, true);

      alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
      fillRandom(dofs, rng);

      // The traces of the cell on its faces, each face in another parametrisation than the one
      // before, and what the neighbours write there: the trace of a field of their own on another
      // face of theirs, in the same parametrisation. A face without neighbour reads zero there.
      alignas(Alignment) std::array<FaceTraceT, Cell::NumFaces> own{};
      alignas(Alignment) std::array<FaceTraceT, Cell::NumFaces> across{};
      for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
        writeFaceTrace(dofs.data(), side, side % 2, own[side].data());
        if (faceTypes[side] == FaceType::Regular) {
          alignas(Alignment) std::array<real, tensor::I::size()> neighbor{};
          fillRandom(neighbor, rng);
          writeFaceTrace(
              neighbor.data(), (side + 1) % Cell::NumFaces, 1 - side % 2, across[side].data());
        }
      }

      // the faces, lifted with the operators the setup gave them into `target`
      const auto addFaces = [&](const LocalIntegrationData& local, real* target) {
        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          if (appliesFaceCorrector(faceTypes[side])) {
            applyFaceCorrector(
                own[side].data(), across[side].data(), local, side, side % 2, target);
          }
        }
      };

      alignas(Alignment) std::array<real, Size> plain{};
      Lifted::addVolume(dofs.data(), plain.data(), plainLocal);
      addFaces(plainLocal, plain.data());

      alignas(Alignment) std::array<real, Size> adjusted{};
      alignas(Alignment) std::array<real, Size> buffer{};
      Lifted::addVolume(dofs.data(), adjusted.data(), adjustedLocal);
      addFaces(adjustedLocal, buffer.data());
      applyWeight<KernelT>(buffer.data(), adjusted.data(), adjustedLocal);

      requireSameUpdate(adjusted.data(), plain.data());
    }
  }
}

/// The same for the lift of a fault face: the side of a weight-adjusted cell lifts the imposed
/// state without material, and with the weight applied that is the lift with the material.
template <bool Enabled>
void faultSameWhereConstant() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    using KernelT = std::conditional_t<Enabled, kernel::weightAdjust, void>;
    constexpr std::size_t Points = dr::misc::NumBoundaryGaussPoints;
    constexpr std::size_t Size = Lifted::Tensor::size();

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261018);
    std::uniform_real_distribution<double> positive(0.4, 2.5);

    for (std::size_t sample = 0; sample < 3; ++sample) {
      const auto material = coefficients::configuredMaterial<Material>(rng);
      const auto frame = faultflux::randomFrame(rng);
      alignas(Alignment) std::array<real, tensor::T::size()> matT{};
      alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
      auto viewT = init::T::view::create(matT.data());
      auto viewTinv = init::Tinv::view::create(matTinv.data());
      seissol::model::getFaceRotationMatrix(
          frame.normal, frame.tangent1, frame.tangent2, viewT, viewTinv);
      std::array<double, 36> bond{};
      seissol::model::getBondMatrix(frame.normal, frame.tangent1, frame.tangent2, bond);

      // the rotation the face stores: the one of the face, at every point where it may be curved
      alignas(Alignment) std::array<real, dr::FaultFluxLayout::RotationSize> rotation{};
      if constexpr (Curvilinear) {
        init::TPoints::view::create(rotation.data()).forall([&](const auto* entry, auto& value) {
          value = viewT(entry[1], entry[2]);
        });
      } else {
        std::copy(matT.begin(), matT.end(), rotation.begin());
      }
      std::array<double, Points> fluxScale{};
      fluxScale.fill(-0.37 * positive(rng));
      // where a cell may be curved, the side of a weight-adjusted cell lifts with the surface
      // Jacobian alone, and its weight carries the 1 / J of the scale (see
      // initializer::CurvedCell::setWeightScale)
      const double determinant = Curvilinear ? positive(rng) : 1.0;
      std::array<double, Points> adjustedScale{};
      for (std::size_t point = 0; point < Points; ++point) {
        adjustedScale[point] = fluxScale[point] * determinant;
      }
      std::array<Material, dr::ImpedancePoints> atPoints{};
      atPoints.fill(material);

      alignas(Alignment) std::array<real, dr::FaultFluxLayout::Size> plainLift{};
      alignas(Alignment) std::array<real, dr::FaultFluxLayout::Size> adjustedLift{};
      seissol::initializer::setPointwiseFaultFlux(
          plainLift.data(), rotation.data(), fluxScale, atPoints, material, bond, false);
      seissol::initializer::setPointwiseFaultFlux(
          adjustedLift.data(), rotation.data(), adjustedScale, atPoints, material, bond, true);

      alignas(Alignment) std::array<real, tensor::QInterpolated::size()> imposed{};
      fillRandom(imposed, rng);

      LocalIntegrationData local{};
      const auto values = seissol::model::getStarCoefficients(material);
      for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
        for (std::size_t a = 0; a < values.size(); ++a) {
          local.materialCoefficients[a][point] = values[a];
        }
        if constexpr (Curvilinear) {
          local.weightScale[point] = static_cast<real>(1.0 / determinant);
        }
      }

      const auto lift = [&](const real* scalars, real* target, unsigned side, unsigned relation) {
        dynamicRupture::kernel::nodalFlux krnl{};
        krnl.bindGlobals(seissol::Pool::host());
        kernels::bindFaultFluxOperands(krnl, scalars);
        krnl.QInterpolated = imposed.data();
        Lifted::bind(krnl, target);
        krnl.execute(side, relation);
      };
      for (unsigned side = 0; side < Cell::NumFaces; ++side) {
        for (unsigned relation = 0; relation < 2; ++relation) {
          alignas(Alignment) std::array<real, Size> plain{};
          lift(plainLift.data(), plain.data(), side, relation);
          alignas(Alignment) std::array<real, Size> buffer{};
          alignas(Alignment) std::array<real, Size> adjusted{};
          lift(adjustedLift.data(), buffer.data(), side, relation);
          applyWeight<KernelT>(buffer.data(), adjusted.data(), local);
          requireSameUpdate(adjusted.data(), plain.data());
        }
      }
    }
  }
}

} // namespace weightadjusted

TEST_CASE("Applying the inverse weight against the matrices it is made of") {
  weightadjusted::kernelAgainstMatrices<WeightAdjustable>();
}

TEST_CASE("A weight-adjusted cell adds what a plain one does where the material does not vary") {
  weightadjusted::sameWhereConstant<WeightAdjustable>();
}

TEST_CASE("A weight-adjusted fault side lifts what a plain one does where the material does not "
          "vary") {
  weightadjusted::faultSameWhereConstant<WeightAdjustable && NodalFaultFlux>();
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_WEIGHTADJUSTED_T_H_
