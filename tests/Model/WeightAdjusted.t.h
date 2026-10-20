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
// - with one material for the cell, what a weight-adjusted cell adds -- the volume term, the local
//   and the neighbour fluxes, the lift of a fault face -- is what a cell that lifts with the
//   material adds;
// - which cells are weight-adjusted, and how much shorter their time step is.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Common/Constants.h"
#include "DynamicRupture/Typedefs.h"
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "FaultFlux.t.h" // IWYU pragma: keep
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Model/FaultFlux.h"
#include "Initializer/Model/WeightAdjusted.h"
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
  weightadjusted::unweightedFluxDecomposes<WeightAdjustable>();
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
    REQUIRE_FALSE(isWeightAdjusted(spread.ratio(), 2.0));
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
    REQUIRE(isWeightAdjusted(spread.ratio(), 2.0));
    REQUIRE_FALSE(isWeightAdjusted(spread.ratio(), 5.0));
    REQUIRE_FALSE(isWeightAdjusted(spread.ratio(), 0.0));
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

#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

// The anelastic solver writes its corrector into the extended quantities and adds the relaxation
// in a kernel of its own; it has no single kernel sequence to compare, as for NodalCorrector.t.h.
namespace weightadjusted {

/// what a field holds for one fused simulation, as a matrix of modes by quantities
template <typename ViewT>
Eigen::MatrixXd fieldOf(ViewT view, std::size_t simulation) {
  constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
  constexpr std::size_t Quantities = tensor::Q::Shape[multisim::BasisFunctionDimension + 1];
  auto sliced = multisim::simtensor(view, simulation);
  Eigen::MatrixXd field = Eigen::MatrixXd::Zero(Basis, Quantities);
  for (std::size_t row = 0; row < Basis; ++row) {
    for (std::size_t column = 0; column < Quantities; ++column) {
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
  constexpr std::size_t Quantities = tensor::Q::Shape[multisim::BasisFunctionDimension + 1];
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

  auto viewCandidate = init::Q::view::create(candidate);
  auto viewReference = init::Q::view::create(reference);
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
  krnl.Q = target;
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
/// inverse weight of the material of every point, projected back and added.
template <bool Enabled>
void kernelAgainstMatrices() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    using KernelT = std::conditional_t<Enabled, kernel::weightAdjust, void>;
    constexpr std::size_t Points = tensor::materialEval::Shape[0];
    constexpr std::size_t Basis = tensor::materialEval::Shape[1];
    constexpr std::size_t Quantities = Material::NumQuantities;
    static_assert(Points == MaterialSampleCount);

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
      // a material of its own at every point, a fluid among them where the medium has shear
      LocalIntegrationData local{};
      std::vector<Eigen::MatrixXd> weights;
      for (std::size_t point = 0; point < Points; ++point) {
        const auto material = coefficients::configuredMaterial<Material>(rng, point % 7 == 3);
        const auto values = seissol::model::getStarCoefficients(material);
        for (std::size_t a = 0; a < values.size(); ++a) {
          local.materialCoefficients[a][point] = values[a];
        }
        weights.emplace_back(inverseweight::inverseWeight<Material, Quantities>(material));
      }

      alignas(Alignment) std::array<real, tensor::Q::size()> buffer{};
      alignas(Alignment) std::array<real, tensor::Q::size()> start{};
      fillRandom(buffer, rng);
      fillRandom(start, rng);
      alignas(Alignment) std::array<real, tensor::Q::size()> update = start;
      applyWeight<KernelT>(buffer.data(), update.data(), local);

      alignas(Alignment) std::array<real, tensor::Q::size()> expectedData{};
      auto expectedView = init::Q::view::create(expectedData.data());
      for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
        const Eigen::MatrixXd values =
            evaluate * fieldOf(init::Q::view::create(buffer.data()), simulation);
        Eigen::MatrixXd product = Eigen::MatrixXd::Zero(Points, Quantities);
        for (std::size_t point = 0; point < Points; ++point) {
          product.row(point) = values.row(point) * weights[point];
        }
        const Eigen::MatrixXd expected =
            fieldOf(init::Q::view::create(start.data()), simulation) + project * product;
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

/// With one material for the cell, the corrector of a weight-adjusted cell -- the volume term into
/// the field, the local and the neighbour fluxes without material into the buffer, the weight
/// applied to it -- is the one of a cell that lifts with the material. Every face takes a path of
/// its own: a regular face both fluxes, a free surface its mirrored state, a fault face the
/// subtraction of its own trace.
template <bool Enabled>
void sameWhereConstant() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    using KernelT = std::conditional_t<Enabled, kernel::weightAdjust, void>;
    using Volume = std::conditional_t<Enabled, kernel::volume, void>;

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
      NeighboringIntegrationData plainNeighboring{};
      LocalIntegrationData adjustedLocal{};
      NeighboringIntegrationData adjustedNeighboring{};
      std::array<std::array<real, tensor::AplusT::size()>, Cell::NumFaces> weakAplus{};
      nodalcorrector::setupCell<Material>(
          cell, field, faceTypes, plainLocal, plainNeighboring, weakAplus, false);
      nodalcorrector::setupCell<Material>(
          cell, field, faceTypes, adjustedLocal, adjustedNeighboring, weakAplus, true);

      alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
      fillRandom(dofs, rng);
      std::array<std::array<real, tensor::I::size()>, Cell::NumFaces> neighbors{};
      for (auto& neighbor : neighbors) {
        fillRandom(neighbor, rng);
      }

      // the faces, lifted with the scalars the setup gave them into `target`
      const auto addFaces = [&](const LocalIntegrationData& local,
                                const NeighboringIntegrationData& neighboring,
                                real* target) {
        kernel::localFlux flux{};
        flux.bindGlobals(seissol::Pool::host());
        flux.I = dofs.data();
        flux.Q = target;
        kernel::neighboringFlux neighbor{};
        neighbor.bindGlobals(seissol::Pool::host());
        neighbor.Q = target;
        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          if (appliesLocalFlux(faceTypes[side])) {
            kernels::bindLocalFluxOperands(flux, local, side);
            flux.execute(side);
          }
          if (faceTypes[side] == FaceType::Regular) {
            neighbor.I = neighbors[side].data();
            kernels::bindNeighborFluxOperands(neighbor, local, neighboring, side);
            neighbor.execute(0, side);
          }
        }
      };

      alignas(Alignment) std::array<real, tensor::Q::size()> plain{};
      addVolumeTerm<Volume>(dofs.data(), plain.data(), plainLocal);
      addFaces(plainLocal, plainNeighboring, plain.data());

      alignas(Alignment) std::array<real, tensor::Q::size()> adjusted{};
      alignas(Alignment) std::array<real, tensor::Q::size()> buffer{};
      addVolumeTerm<Volume>(dofs.data(), adjusted.data(), adjustedLocal);
      addFaces(adjustedLocal, adjustedNeighboring, buffer.data());
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
      std::array<Material, dr::ImpedancePoints> atPoints{};
      atPoints.fill(material);

      alignas(Alignment) std::array<real, dr::FaultFluxLayout::Size> plainLift{};
      alignas(Alignment) std::array<real, dr::FaultFluxLayout::Size> adjustedLift{};
      seissol::initializer::setPointwiseFaultFlux(
          plainLift.data(), rotation.data(), fluxScale, atPoints, material, bond, false);
      seissol::initializer::setPointwiseFaultFlux(
          adjustedLift.data(), rotation.data(), fluxScale, atPoints, material, bond, true);

      alignas(Alignment) std::array<real, tensor::QInterpolated::size()> imposed{};
      fillRandom(imposed, rng);

      LocalIntegrationData local{};
      const auto values = seissol::model::getStarCoefficients(material);
      for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
        for (std::size_t a = 0; a < values.size(); ++a) {
          local.materialCoefficients[a][point] = values[a];
        }
      }

      const auto lift = [&](const real* scalars, real* target, unsigned side, unsigned relation) {
        dynamicRupture::kernel::nodalFlux krnl{};
        krnl.bindGlobals(seissol::Pool::host());
        kernels::bindFaultFluxOperands(krnl, scalars);
        krnl.QInterpolated = imposed.data();
        krnl.Q = target;
        krnl.execute(side, relation);
      };
      for (unsigned side = 0; side < Cell::NumFaces; ++side) {
        for (unsigned relation = 0; relation < 2; ++relation) {
          alignas(Alignment) std::array<real, tensor::Q::size()> plain{};
          lift(plainLift.data(), plain.data(), side, relation);
          alignas(Alignment) std::array<real, tensor::Q::size()> buffer{};
          alignas(Alignment) std::array<real, tensor::Q::size()> adjusted{};
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

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_MODEL_WEIGHTADJUSTED_T_H_
