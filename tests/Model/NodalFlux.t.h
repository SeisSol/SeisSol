// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_NODALFLUX_T_H_
#define SEISSOL_TESTS_MODEL_NODALFLUX_T_H_

// Where the material varies along a face the flux is applied at the nodes of
// that face, from a few scalars per node instead of one matrix. A material that
// does not vary has to give back exactly what the matrix form gives, and the
// reference here is the real one: the operator that computeFluxSolverLocal
// builds, contracted with the matrices the modal kernel uses. A face taken for
// another, a transposition, or a rotation left out shows up in that comparison
// and in nothing weaker.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshTools.h"
#include "Initializer/Model/CurvedCell.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <type_traits>
#include <vector>

// The anelastic solver writes its face contribution into the extended
// quantities and has no localFlux of its own, so the comparison below has no
// kernel to run there. Its decomposition is a separate table and is checked
// where that table is.
#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

namespace seissol::unit_test {

namespace nodalflux {

/// Makes the kernel names dependent on a template parameter, so that a build
/// without a varying material never instantiates a body that names operands it
/// does not have.
template <bool Enabled>
struct Kernels {
  using LocalFlux = std::conditional_t<Enabled, kernel::localFlux, kernel::localFlux>;
  using FluxSolver = std::conditional_t<Enabled, kernel::computeFluxSolverLocal, void>;
};

/// The boundary conditions whose operator the flux itself carries. The others
/// -- gravity, Dirichlet, analytical -- state their ghost cell at the nodes of
/// the face first and apply the neighbour matrix to it, so they are not a
/// question about this decomposition.
constexpr FaceType LinearFaceTypes[] = {
    FaceType::Regular, FaceType::FreeSurface, FaceType::Outflow};

template <bool Enabled>
void compareAgainstModal(FaceType faceType) {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t NQ = tensor::QgodLocal::Shape[0];
    constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
    constexpr std::size_t Nodes = generated::FaceNodes;
    constexpr std::size_t Coefficients = generated::FluxNumCoefficients;

    std::mt19937 rng(20260927);
    std::uniform_real_distribution<double> positive(0.4, 2.5);
    std::normal_distribution<double> gauss(0.0, 1.0);

    for (std::size_t sample = 0; sample < 8; ++sample) {
      const auto local = coefficients::configuredMaterial<Material>(rng);
      const auto neighbor = coefficients::configuredMaterial<Material>(rng);

      // a general orthonormal frame: the first tangent of a face is one of its
      // edges, so nothing may lean on a particular choice
      Eigen::Matrix3d noise;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          noise(i, j) = gauss(rng);
        }
      }
      Eigen::Matrix3d frame = Eigen::HouseholderQR<Eigen::Matrix3d>(noise).householderQ();
      if (frame.determinant() < 0) {
        frame.col(2) *= -1.0;
      }
      const CoordinateT normal{frame(0, 0), frame(1, 0), frame(2, 0)};
      const CoordinateT tangent1{frame(0, 1), frame(1, 1), frame(2, 1)};
      const CoordinateT tangent2{frame(0, 2), frame(1, 2), frame(2, 2)};

      alignas(Alignment) std::array<real, tensor::T::size()> matT{};
      alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
      auto viewT = init::T::view::create(matT.data());
      auto viewTinv = init::Tinv::view::create(matTinv.data());
      seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, viewT, viewTinv);

      const double fluxScale = -0.37 * positive(rng);

      // the operator in face coordinates, and the star the setup kernel folds in
      alignas(Alignment) std::array<real, tensor::QgodLocal::size()> godLocal{};
      alignas(Alignment) std::array<real, tensor::QgodNeighbor::size()> godNeighbor{};
      auto viewGodLocal = init::QgodLocal::view::create(godLocal.data());
      auto viewGodNeighbor = init::QgodNeighbor::view::create(godNeighbor.data());
      seissol::model::getTransposedGodunovState(
          local, neighbor, faceType, viewGodLocal, viewGodNeighbor);

      alignas(Alignment) std::array<real, tensor::star::size(0)> star{};
      auto viewStar = init::star::view<0>::create(star.data());
      seissol::model::getTransposedCoefficientMatrix(local, 0, viewStar);

      // 1. the reference operator, from the kernel the setup actually runs
      alignas(Alignment) std::array<real, tensor::AplusT::size()> aplus{};
      alignas(Alignment) std::array<real, tensor::QcorrLocal::size()> qcorr{};
      typename Kernels<Enabled>::FluxSolver solver{};
      solver.fluxScale = fluxScale;
      solver.AplusT = aplus.data();
      solver.QgodLocal = godLocal.data();
      solver.QcorrLocal = qcorr.data();
      solver.T = matT.data();
      solver.Tinv = matTinv.data();
      solver.star(0) = star.data();
      solver.execute();

      // 2. the scalars, read off the face-coordinate operator
      Eigen::MatrixXd faceOperator = Eigen::MatrixXd::Zero(NQ, NQ);
      for (std::size_t row = 0; row < NQ; ++row) {
        for (std::size_t column = 0; column < NQ; ++column) {
          double value = 0.0;
          for (std::size_t inner = 0; inner < NQ; ++inner) {
            const double god = viewGodLocal.isInRange(row, inner) ? viewGodLocal(row, inner) : 0.0;
            const double a = viewStar.isInRange(inner, column) ? viewStar(inner, column) : 0.0;
            value += god * a;
          }
          faceOperator(row, column) = value;
        }
      }

      std::array<std::array<real, Nodes>, Coefficients> scalars{};
      for (std::size_t c = 0; c < Coefficients; ++c) {
        const auto& source = generated::FluxCoefficientSources[c];
        // the scale the matrix form carries in AplusT rides on the scalars here
        scalars[c].fill(static_cast<real>(fluxScale * faceOperator(source.row, source.column)));
      }

      // a field to apply both forms to
      alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
      for (auto& value : dofs) {
        value = static_cast<real>(gauss(rng));
      }

      // the side is a template parameter so that every member of a family is
      // read through its own view: they need not share a layout
      const auto checkSide = [&](auto tag) {
        constexpr std::uint8_t side = decltype(tag)::value;
        // 3. the nodal form
        alignas(Alignment) std::array<real, tensor::Q::size()> nodal{};
        typename Kernels<Enabled>::LocalFlux krnl{};
        krnl.I = dofs.data();
        krnl.Q = nodal.data();
        // the rotation of the face, at each of its nodes where a face may be curved
        alignas(Alignment)
            std::array<real, sizeof(LocalIntegrationData::faceRotation[0]) / sizeof(real)>
                rotation{};
        initializer::CurvedCell::setConstantRotation(matT.data(), rotation.data());
        kernels::bindFaceRotation(krnl, rotation.data());
        // the constants in the arrangement the kernel was generated against,
        // which only the pool holds
        krnl.bindGlobals(seissol::Pool::host());
        for (std::size_t c = 0; c < Coefficients; ++c) {
          krnl.fluxCoefficientsLocal(c) = scalars[c].data();
        }
        krnl.execute(side);

        // 4. the modal form, by hand out of the matrices it is made of
        const auto denseOf = [](auto view, std::size_t rows, std::size_t columns) {
          Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(rows, columns);
          for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t column = 0; column < columns; ++column) {
              if (view.isInRange(row, column)) {
                dense(row, column) = view(row, column);
              }
            }
          }
          return dense;
        };
        // a build that bundles simulations stores the global matrices the other
        // way round, so read them by their own extents and put them back
        const auto mathMatrix = [&denseOf](auto view, std::size_t rows, std::size_t columns) {
          const Eigen::MatrixXd stored = denseOf(view, rows, columns);
          return multisim::NumSimulations > 1 ? Eigen::MatrixXd(stored.transpose()) : stored;
        };
        const auto rDivM = mathMatrix(
            init::rDivM::view<side>::create(const_cast<real*>(init::rDivM::Values[side])),
            tensor::rDivM::Shape[side][0],
            tensor::rDivM::Shape[side][1]);
        const auto fMrT =
            mathMatrix(init::fMrT::view<side>::create(const_cast<real*>(init::fMrT::Values[side])),
                       tensor::fMrT::Shape[side][0],
                       tensor::fMrT::Shape[side][1]);
        constexpr std::size_t Written = tensor::AplusT::Shape[1];
        const auto operatorT = denseOf(init::AplusT::view::create(aplus.data()), NQ, Written);
        auto viewI = init::I::view::create(dofs.data());
        auto viewQ = init::Q::view::create(nodal.data());

        // the operator is shared by the simulations a build bundles, the field
        // is not, so every one of them has to come back right
        for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
          auto slicedI = multisim::simtensor(viewI, simulation);
          auto slicedQ = multisim::simtensor(viewQ, simulation);

          Eigen::MatrixXd field = Eigen::MatrixXd::Zero(Basis, NQ);
          for (std::size_t row = 0; row < Basis; ++row) {
            for (std::size_t column = 0; column < NQ; ++column) {
              field(row, column) = slicedI(row, column);
            }
          }
          const Eigen::MatrixXd expected = rDivM * fMrT * field * operatorT;

          // Entries of an output quantity cancel against each other, so what one
          // of them may be off by follows the size of that quantity in the
          // update, not its own. The quantities themselves differ by orders of
          // magnitude -- a stress by the moduli, a velocity by the inverse
          // density -- so each is measured against its own column.
          constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
          for (std::size_t column = 0; column < Written; ++column) {
            const double scale = std::max(1.0, expected.col(column).cwiseAbs().maxCoeff());
            for (std::size_t row = 0; row < Basis; ++row) {
              REQUIRE(slicedQ(row, column) ==
                      doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
            }
          }
        }
      };
      checkSide(std::integral_constant<std::uint8_t, 0>{});
      checkSide(std::integral_constant<std::uint8_t, 1>{});
      checkSide(std::integral_constant<std::uint8_t, 2>{});
      checkSide(std::integral_constant<std::uint8_t, 3>{});
    }
  }
}

} // namespace nodalflux

TEST_CASE("Nodal flux against the matrix form") {
  for (const auto faceType : nodalflux::LinearFaceTypes) {
    nodalflux::compareAgainstModal<NodalFlux>(faceType);
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_NODALFLUX_T_H_
