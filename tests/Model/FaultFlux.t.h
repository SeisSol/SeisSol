// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_FAULTFLUX_T_H_
#define SEISSOL_TESTS_MODEL_FAULTFLUX_T_H_

// The lift of a fault face takes the imposed state of a side, given in the coordinates of the
// face, into the side's cell. Two things are checked here: that the matrix a face keeps per side
// is the flux of the fault normal, whatever the orientation of the face and whatever the material;
// and that where a face carries the lift per point, every point lifts with the material there,
// just as that matrix does for the material of the point.

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
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <type_traits>
#include <vector>

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

/// The entries a view stores, as a dense matrix.
template <typename ViewT>
Eigen::MatrixXd denseOf(const ViewT& view, std::size_t rows, std::size_t columns) {
  Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(rows, columns);
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t column = 0; column < columns; ++column) {
      if (view.isInRange(row, column)) {
        dense(row, column) = view(row, column);
      }
    }
  }
  return dense;
}

/// The matrix of one side, formed by what a build that keeps one per side runs.
template <typename MaterialT>
Eigen::MatrixXd matrixForm(const MaterialT& material,
                           const std::array<double, 36>& bond,
                           const real* matT,
                           double fluxScale) {
  alignas(Alignment) std::array<real, tensor::fluxSolver::size()> fluxSolver{};
  seissol::initializer::setMatrixFaultFlux(fluxSolver.data(), matT, fluxScale, material, bond);
  return denseOf(init::fluxSolver::view::create(fluxSolver.data()),
                 tensor::fluxSolver::Shape[0],
                 tensor::fluxSolver::Shape[1]);
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
    const Eigen::MatrixXd rotation = denseOf(viewT, Columns, Columns).topLeftCorner(Rows, Rows);

    // lift(q, p) is what face quantity q gives global quantity p; `magnitude` adds up the size of
    // the terms each entry of it is a sum of
    Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(Rows, Columns);
    Eigen::MatrixXd magnitude = Eigen::MatrixXd::Zero(Rows, Columns);
    for (std::size_t dim = 0; dim < 3; ++dim) {
      const double weight = fluxScale * frame.normal[dim];
      expected += weight * (transposed[dim].transpose() * rotation).transpose();
      magnitude += std::abs(weight) *
                   (transposed[dim].cwiseAbs().transpose() * rotation.cwiseAbs()).transpose();
    }

    // The entries are orders of magnitude apart -- a face stress enters the velocities with the
    // inverse density, a face velocity the stresses with the moduli and the memory variables with
    // a weight of one -- so every entry is measured against the size of the terms of its row
    // within the quantity group of its column, and a small entry is not hidden behind a large
    // one. The terms rather than the entries, since a row may cancel to zero exactly in one form
    // and to roundoff in the other; an entry without any terms has to vanish exactly.
    std::array<std::size_t, Columns> groupOf{};
    std::size_t groups = 0;
    std::size_t covered = 0;
    for (const auto& group : Material::RotationGroups) {
      for (std::size_t i = 0; i < group.extent() && covered < Columns; ++i) {
        groupOf[covered++] = groups;
      }
      ++groups;
    }
    REQUIRE(covered == Columns);

    constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
    for (std::size_t row = 0; row < Rows; ++row) {
      std::vector<double> scaleOfGroup(groups, 0.0);
      for (std::size_t column = 0; column < Columns; ++column) {
        scaleOfGroup[groupOf[column]] =
            std::max(scaleOfGroup[groupOf[column]], magnitude(row, column));
      }
      for (std::size_t column = 0; column < Columns; ++column) {
        const double scale = scaleOfGroup[groupOf[column]];
        if (scale == 0.0) {
          REQUIRE(lift(row, column) == expected(row, column));
        } else {
          REQUIRE(lift(row, column) ==
                  doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
        }
      }
    }
  }
}

/// Where the lift writes: the quantities of the cell, or, for the solver that keeps its memory
/// variables apart from them, the extended quantities that include those.
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
using LiftTarget = tensor::Qext;
using LiftTargetInit = init::Qext;
#else
using LiftTarget = tensor::Q;
using LiftTargetInit = init::Q;
#endif

template <typename KernelT>
void bindLiftTarget(KernelT& krnl, real* target) {
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
  krnl.Qext = target;
#else
  krnl.Q = target;
#endif
}

/// Makes the kernel operands dependent on a template parameter, so that a build whose face keeps
/// one matrix per side never instantiates a body that names the operands of the pointwise form.
template <bool Enabled>
void pointwiseAgainstMatrix() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t Basis = LiftTarget::Shape[multisim::BasisFunctionDimension];
    constexpr std::size_t Points = dr::misc::NumBoundaryGaussPoints;
    constexpr std::size_t Quantities =
        tensor::QInterpolated::Shape[multisim::BasisFunctionDimension + 1];
    constexpr std::size_t Written = tensor::fluxSolver::Shape[1];
    constexpr std::size_t Simulations = multisim::NumSimulations;

    std::mt19937 rng(20260929);
    std::uniform_real_distribution<double> positive(0.4, 2.5);
    std::normal_distribution<double> gauss(0.0, 1.0);

    for (std::size_t sample = 0; sample < 4; ++sample) {
      const auto cell = coefficients::configuredMaterial<Material>(rng);
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

      // A material of its own at every point, set up the way the fault sets up the material at
      // its points: default constructed, with only the fields the material declares taken from
      // the samples, the same for every fused simulation, and the padding given the material of
      // the cell. What is derived rather than sampled -- the relaxation frequencies -- is left at
      // its default there, so a lift that took it from a point instead of the cell would show.
      // The reference of a point is the matrix form of the material that point stands for: the
      // cell's, with the sampled fields of the point.
      std::array<Material, dr::ImpedancePoints> atPoints{};
      std::vector<Eigen::MatrixXd> lifts(Points);
      for (std::size_t point = 0; point < Points; ++point) {
        const auto drawn = coefficients::configuredMaterial<Material>(rng);
        Material sampled{};
        Material reference = cell;
        for (const auto& [name, member] : Material::ParameterMap) {
          sampled.*member = drawn.*member;
          reference.*member = drawn.*member;
        }
        for (std::size_t simulation = 0; simulation < Simulations; ++simulation) {
          atPoints[point * Simulations + simulation] = sampled;
        }
        lifts[point] = matrixForm(reference, bond, matT.data(), fluxScale);
      }
      for (std::size_t point = Points * Simulations; point < atPoints.size(); ++point) {
        atPoints[point] = cell;
      }

      alignas(Alignment) std::array<real, dr::FaultFluxLayout::Size> pointwise{};
      seissol::initializer::setPointwiseFaultFlux(
          pointwise.data(), matT.data(), fluxScale, atPoints, cell, bond);

      alignas(Alignment) std::array<real, tensor::QInterpolated::size()> imposed{};
      for (auto& value : imposed) {
        value = static_cast<real>(gauss(rng));
      }

      const auto checkSide = [&](auto sideTag, auto relationTag) {
        constexpr unsigned Side = decltype(sideTag)::value;
        constexpr unsigned Relation = decltype(relationTag)::value;

        alignas(Alignment) std::array<real, LiftTarget::size()> dofs{};
        dynamicRupture::kernel::nodalFlux krnl{};
        krnl.bindGlobals(seissol::Pool::host());
        kernels::bindFaultFluxOperands(krnl, pointwise.data());
        krnl.QInterpolated = imposed.data();
        bindLiftTarget(krnl, dofs.data());
        krnl.execute(Side, Relation);

        // a build that bundles simulations stores the global matrices the other way round
        constexpr auto Index = tensor::V3mTo2nTWDivM::index(Side, Relation);
        const auto viewLift =
            init::V3mTo2nTWDivM::view<Side, Relation>::create(init::V3mTo2nTWDivM::Values[Index]);
        const Eigen::MatrixXd project =
            Simulations > 1 ? Eigen::MatrixXd(denseOf(viewLift, Points, Basis).transpose())
                            : denseOf(viewLift, Basis, Points);

        auto viewImposed = init::QInterpolated::view::create(imposed.data());
        auto viewDofs = LiftTargetInit::view::create(dofs.data());
        // the material is shared by the simulations a build bundles, the state is not
        for (std::size_t simulation = 0; simulation < Simulations; ++simulation) {
          auto stateOfSim = multisim::simtensor(viewImposed, simulation);
          auto dofsOfSim = multisim::simtensor(viewDofs, simulation);

          // the lift of every point applied to the state there, projected into the cell; and
          // the size of the terms that adds up
          Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(Basis, Written);
          Eigen::MatrixXd magnitude = Eigen::MatrixXd::Zero(Basis, Written);
          for (std::size_t point = 0; point < Points; ++point) {
            Eigen::RowVectorXd state = Eigen::RowVectorXd::Zero(Quantities);
            for (std::size_t quantity = 0; quantity < Quantities; ++quantity) {
              state(quantity) = stateOfSim(point, quantity);
            }
            expected += project.col(point) * (state * lifts[point]);
            magnitude +=
                project.col(point).cwiseAbs() * (state.cwiseAbs() * lifts[point].cwiseAbs());
          }

          // every quantity against its own size, as in the matrix check above
          constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
          for (std::size_t column = 0; column < Written; ++column) {
            const double scale = magnitude.col(column).maxCoeff();
            for (std::size_t row = 0; row < Basis; ++row) {
              if (scale == 0.0) {
                REQUIRE(dofsOfSim(row, column) == expected(row, column));
              } else {
                REQUIRE(dofsOfSim(row, column) ==
                        doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
              }
            }
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

#endif // SEISSOL_TESTS_MODEL_FAULTFLUX_T_H_
