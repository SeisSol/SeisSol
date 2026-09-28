// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_FAULTFLUX_T_H_
#define SEISSOL_TESTS_MODEL_FAULTFLUX_T_H_

// The lift of a fault face takes the imposed state of a side, given in the coordinates of the
// face, into the side's cell. Two things are checked here: that the operator it applies is the
// flux of the fault normal, whatever the orientation of the face and whatever the material; and
// that where a face carries it per point, a material that does not vary gives back exactly what
// the one matrix per side gives.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "DynamicRupture/Typedefs.h"
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshTools.h"
#include "Initializer/Model/FaultFlux.h"
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

// The anelastic solver lifts into the extended quantities, which the checks below do not model;
// its lift is the same code with a wider target.
#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

namespace seissol::unit_test {

namespace faultflux {

/// A random right-handed orthonormal frame: the fault normal and the two tangents.
struct Frame {
  VrtxCoords normal;
  VrtxCoords tangent1;
  VrtxCoords tangent2;
};

inline Frame randomFrame(std::mt19937& rng) {
  std::normal_distribution<double> gauss(0.0, 1.0);
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
  return Frame{{frame(0, 0), frame(1, 0), frame(2, 0)},
               {frame(0, 1), frame(1, 1), frame(2, 1)},
               {frame(0, 2), frame(1, 2), frame(2, 2)}};
}

/// The matrix of one side, as a build that keeps one per side forms it.
template <typename MaterialT>
Eigen::MatrixXd matrixForm(const MaterialT& material,
                           const std::array<double, 36>& bond,
                           real* matT,
                           double fluxScale) {
  alignas(Alignment) std::array<real, tensor::star::size(0)> star{};
  auto viewStar = init::star::view<0>::create(star.data());
  seissol::model::getTransposedCoefficientMatrix(
      seissol::model::getRotatedMaterialCoefficients(bond, material), 0, viewStar);

  alignas(Alignment) std::array<real, tensor::fluxSolver::size()> fluxSolver{};
  dynamicRupture::kernel::rotateFluxMatrix krnl;
  krnl.T = matT;
  krnl.fluxSolver = fluxSolver.data();
  krnl.fluxScaleDR = fluxScale;
  krnl.star(0) = star.data();
  krnl.execute();

  const auto view = init::fluxSolver::view::create(fluxSolver.data());
  Eigen::MatrixXd dense =
      Eigen::MatrixXd::Zero(tensor::fluxSolver::Shape[0], tensor::fluxSolver::Shape[1]);
  for (std::size_t row = 0; row < tensor::fluxSolver::Shape[0]; ++row) {
    for (std::size_t column = 0; column < tensor::fluxSolver::Shape[1]; ++column) {
      if (view.isInRange(row, column)) {
        dense(row, column) = view(row, column);
      }
    }
  }
  return dense;
}

/// The lift is the coefficient matrix of the normal direction, in global coordinates, applied to
/// the imposed state turned back into them: T A' = (n_x A + n_y B + n_z C) T, with A' the matrix
/// of the first direction of the material seen in face coordinates. Checked on the matrix the
/// setup folds; a material taken in the wrong coordinates -- which an isotropic one cannot show
/// -- fails it.
inline void liftIsNormalFlux() {
  using Material = seissol::model::MaterialT;
  constexpr std::size_t Rows = tensor::fluxSolver::Shape[0];
  constexpr std::size_t Columns = tensor::fluxSolver::Shape[1];

  std::mt19937 rng(20260928);
  std::uniform_real_distribution<double> positive(0.4, 2.5);

  for (std::size_t sample = 0; sample < 8; ++sample) {
    const auto material = coefficients::configuredMaterial<Material>(rng);
    const auto frame = randomFrame(rng);

    alignas(Alignment) std::array<real, tensor::T::size()> matT{};
    alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
    auto viewT = init::T::view::create(matT.data());
    auto viewTinv = init::Tinv::view::create(matTinv.data());
    seissol::model::getFaceRotationMatrix(
        frame.normal, frame.tangent1, frame.tangent2, viewT, viewTinv);
    std::array<double, 36> bond{};
    seissol::model::getBondMatrix(frame.normal, frame.tangent1, frame.tangent2, bond);

    const double fluxScale = -0.37 * positive(rng);
    const Eigen::MatrixXd lift = matrixForm(material, bond, matT.data(), fluxScale);

    // the coefficient matrices of the three global directions, transposed as SeisSol keeps them
    std::array<Eigen::MatrixXd, 3> transposed{};
    for (std::size_t dim = 0; dim < 3; ++dim) {
      std::vector<double> data(Rows * Columns, 0.0);
      yateto::DenseTensorView<2, double> view(data.data(), {Rows, Columns});
      seissol::model::getTransposedCoefficientMatrix(material, dim, view);
      transposed[dim] = Eigen::Map<Eigen::MatrixXd>(data.data(), Rows, Columns);
    }
    Eigen::MatrixXd rotation = Eigen::MatrixXd::Zero(Columns, Columns);
    for (std::size_t row = 0; row < Columns; ++row) {
      for (std::size_t column = 0; column < Columns; ++column) {
        if (viewT.isInRange(row, column)) {
          rotation(row, column) = viewT(row, column);
        }
      }
    }

    // lift(q, p) is what face quantity q gives global quantity p
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(Rows, Columns);
    for (std::size_t dim = 0; dim < 3; ++dim) {
      expected += fluxScale * frame.normal[dim] *
                  (transposed[dim].transpose() * rotation.topLeftCorner(Rows, Rows)).transpose();
    }

    constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
    const double scale = expected.cwiseAbs().maxCoeff();
    for (std::size_t row = 0; row < Rows; ++row) {
      for (std::size_t column = 0; column < Columns; ++column) {
        REQUIRE(lift(row, column) ==
                doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
      }
    }
  }
}

/// Makes the kernel operands dependent on a template parameter, so that a build whose face keeps
/// one matrix per side never instantiates a body that names the operands of the pointwise form.
template <bool Enabled>
void pointwiseAgainstMatrix() {
  if constexpr (Enabled && multisim::NumSimulations == 1) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t Basis = tensor::Q::Shape[0];
    constexpr std::size_t Points = dr::misc::NumBoundaryGaussPoints;
    constexpr std::size_t Quantities = tensor::QInterpolated::Shape[1];
    constexpr std::size_t Written = tensor::fluxSolver::Shape[1];

    std::mt19937 rng(20260929);
    std::uniform_real_distribution<double> positive(0.4, 2.5);
    std::normal_distribution<double> gauss(0.0, 1.0);

    for (std::size_t sample = 0; sample < 4; ++sample) {
      const auto material = coefficients::configuredMaterial<Material>(rng);
      const auto frame = randomFrame(rng);

      alignas(Alignment) std::array<real, tensor::T::size()> matT{};
      alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
      auto viewT = init::T::view::create(matT.data());
      auto viewTinv = init::Tinv::view::create(matTinv.data());
      seissol::model::getFaceRotationMatrix(
          frame.normal, frame.tangent1, frame.tangent2, viewT, viewTinv);
      std::array<double, 36> bond{};
      seissol::model::getBondMatrix(frame.normal, frame.tangent1, frame.tangent2, bond);

      const double fluxScale = -0.37 * positive(rng);
      const Eigen::MatrixXd lift = matrixForm(material, bond, matT.data(), fluxScale);

      // the same material at every point
      std::array<Material, dr::ImpedancePoints> atPoints{};
      atPoints.fill(material);
      alignas(Alignment) std::array<real, dr::FaultFluxLayout::Size> pointwise{};
      seissol::initializer::setPointwiseFaultFlux(
          pointwise.data(), matT.data(), fluxScale, atPoints, material, bond);

      alignas(Alignment) std::array<real, tensor::QInterpolated::size()> imposed{};
      for (auto& value : imposed) {
        value = static_cast<real>(gauss(rng));
      }
      const auto viewImposed = init::QInterpolated::view::create(imposed.data());
      Eigen::MatrixXd state = Eigen::MatrixXd::Zero(Points, Quantities);
      for (std::size_t point = 0; point < Points; ++point) {
        for (std::size_t quantity = 0; quantity < Quantities; ++quantity) {
          state(point, quantity) = viewImposed(point, quantity);
        }
      }

      const auto checkSide = [&](auto sideTag, auto relationTag) {
        constexpr unsigned Side = decltype(sideTag)::value;
        constexpr unsigned Relation = decltype(relationTag)::value;

        alignas(Alignment) std::array<real, tensor::Q::size()> dofs{};
        dynamicRupture::kernel::nodalFlux krnl{};
        krnl.bindGlobals(seissol::Pool::host());
        kernels::bindFaultFluxOperands(krnl, pointwise.data());
        krnl.QInterpolated = imposed.data();
        krnl.Q = dofs.data();
        krnl.execute(Side, Relation);

        constexpr auto Index = tensor::V3mTo2nTWDivM::index(Side, Relation);
        const auto viewLift =
            init::V3mTo2nTWDivM::view<Side, Relation>::create(init::V3mTo2nTWDivM::Values[Index]);
        Eigen::MatrixXd project = Eigen::MatrixXd::Zero(Basis, Points);
        for (std::size_t row = 0; row < Basis; ++row) {
          for (std::size_t column = 0; column < Points; ++column) {
            project(row, column) = viewLift(row, column);
          }
        }
        const Eigen::MatrixXd expected = project * state * lift;

        const auto viewDofs = init::Q::view::create(dofs.data());
        constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
        for (std::size_t column = 0; column < Written; ++column) {
          const double scale = std::max(1.0, expected.col(column).cwiseAbs().maxCoeff());
          for (std::size_t row = 0; row < Basis; ++row) {
            REQUIRE(viewDofs(row, column) ==
                    doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
          }
        }
      };
      const auto forRelations = [&](auto sideTag) {
        checkSide(sideTag, std::integral_constant<unsigned, 0>{});
        checkSide(sideTag, std::integral_constant<unsigned, 1>{});
      };
      forRelations(std::integral_constant<unsigned, 0>{});
      forRelations(std::integral_constant<unsigned, 1>{});
      forRelations(std::integral_constant<unsigned, 2>{});
      forRelations(std::integral_constant<unsigned, 3>{});
    }
  }
}

} // namespace faultflux

TEST_CASE("The lift of a fault face is the flux of its normal") { faultflux::liftIsNormalFlux(); }

TEST_CASE("Pointwise fault lift against the matrix form") {
  faultflux::pointwiseAgainstMatrix<NodalFaultFlux>();
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_FAULTFLUX_T_H_
