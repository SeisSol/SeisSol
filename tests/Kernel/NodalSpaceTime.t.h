// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_KERNEL_NODALSPACETIME_T_H_
#define SEISSOL_TESTS_KERNEL_NODALSPACETIME_T_H_

// The space-time predictor for a material that varies inside the cell. Its
// block sweep is replaced by a fixed point there, and the solve that fixed
// point is built on is factorised once for the cell -- source term and all --
// so what the samples ask for beyond that rides along in the iteration.
//
// What holds all of it is the residual of the system the predictor is supposed
// to solve, with the operator assembled at the sample points from the material
// itself. A material that does not vary leaves that residual at zero to the
// last bit; one that varies leaves what the iteration has not caught up with,
// which carries a factor of the timestep per step.

#include <doctest.h>

#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Helper.h"
#include "Equations/poroelastic/Model/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Kernels/STP/Setup.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"
#include "Numerical/Transformation.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

namespace nodalspacetime {

constexpr double Dt = 1.05109e-06;

/// A poroelastic material whose permeability and viscosity are scaled, so that
/// the Biot drag -- the source term -- moves with the sample point while the
/// wave speeds barely do.
inline seissol::model::PoroElasticMaterial drawMaterial(double scale) {
  const auto values = std::vector<double>{
      {40.0e9, 2500, 12.0e9, 10.0e9, 0.2, 600.0e-15 * scale, 3, 2.5e9, 1040, 0.001 / scale}};
  return seissol::model::PoroElasticMaterial(values);
}

template <typename ViewT>
Eigen::MatrixXd denseOf(ViewT view, std::size_t rows, std::size_t columns) {
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

/// The star matrix of one reference direction for one material.
inline Eigen::MatrixXd starOf(const seissol::model::PoroElasticMaterial& material,
                              const double gradient[3]) {
  constexpr std::size_t NQ = seissol::model::PoroElasticMaterial::NumQuantities;
  Eigen::MatrixXd star = Eigen::MatrixXd::Zero(NQ, NQ);
  for (std::size_t dim = 0; dim < 3; ++dim) {
    alignas(Alignment) std::array<real, seissol::tensor::star::size(0)> data{};
    auto view = seissol::init::star::view<0>::create(data.data());
    seissol::model::getTransposedCoefficientMatrix(material, dim, view);
    star += gradient[dim] * denseOf(view, NQ, NQ);
  }
  return star;
}

inline Eigen::MatrixXd sourceOf(const seissol::model::PoroElasticMaterial& material) {
  constexpr std::size_t NQ = seissol::model::PoroElasticMaterial::NumQuantities;
  alignas(Alignment) std::array<real, seissol::tensor::ET::size()> data{};
  auto view = seissol::init::ET::view::create(data.data());
  seissol::model::getTransposedSourceCoefficientTensor(material, view);
  return denseOf(view, NQ, NQ);
}

} // namespace nodalspacetime

TEST_CASE("Space time predictor at the material samples" * doctest::test_suite("kernel")) {
  using namespace nodalspacetime;
  using Material = seissol::model::PoroElasticMaterial;

  if constexpr (NodalMaterial) {
    constexpr std::size_t NQ = Material::NumQuantities;
    constexpr std::size_t Basis = tensor::Q::Shape[0];
    constexpr std::size_t Order = tensor::spaceTimePredictor::Shape[2];
    constexpr std::size_t Points = MaterialSampleCount;

    // no variation at all, then a tenth of the permeability across the cell
    for (const double variation : {0.0, 0.1}) {
      std::mt19937 rng(20260927);
      std::uniform_real_distribution<real> unit(0.0, 1.0);

      double x[Cell::NumVertices];
      double y[Cell::NumVertices];
      double z[Cell::NumVertices];
      for (std::size_t vertex = 0; vertex < Cell::NumVertices; ++vertex) {
        x[vertex] = unit(rng);
        y[vertex] = unit(rng);
        z[vertex] = unit(rng);
      }
      double gradXi[3];
      double gradEta[3];
      double gradZeta[3];
      transformations::tetrahedronGlobalToReferenceJacobian(x, y, z, gradXi, gradEta, gradZeta);
      const double* const gradients[3] = {gradXi, gradEta, gradZeta};

      // the cell's own material, and one per sample point around it
      const auto cellMaterial = drawMaterial(1.0);
      std::vector<Material> sampled;
      for (std::size_t point = 0; point < Points; ++point) {
        const double offset =
            variation * (2.0 * static_cast<double>(point) / static_cast<double>(Points) - 1.0);
        sampled.push_back(drawMaterial(1.0 + offset));
      }

      LocalIntegrationData localIntegration{};
      for (std::size_t dim = 0; dim < 3; ++dim) {
        for (std::size_t component = 0; component < 3; ++component) {
          localIntegration.referenceGradients[dim][component] = gradients[dim][component];
        }
      }
      const auto meanSource = seissol::model::getSourceCoefficients(cellMaterial);
      for (std::size_t point = 0; point < Points; ++point) {
        const auto coefficients = seissol::model::getStarCoefficients(sampled[point]);
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
          localIntegration.materialCoefficients[i][point] = coefficients[i];
        }
        const auto source = seissol::model::getSourceCoefficients(sampled[point]);
        for (std::size_t i = 0; i < source.size(); ++i) {
          localIntegration.sourceCoefficients[i][point] = source[i];
          if constexpr (NodalSourceDeviation) {
            localIntegration.sourceDeviation[i][point] = source[i] - meanSource[i];
          }
        }
      }

      // the time system, factorised for the cell
      alignas(Alignment) std::array<real, tensor::ET::size()> sourceData{};
      auto sourceView = init::ET::view::create(sourceData.data());
      seissol::model::getTransposedSourceCoefficientTensor(cellMaterial, sourceView);
      std::array<std::array<real, tensor::Zinv::size(0)>, NQ> zMatrix{};

      alignas(PagesizeStack) std::array<real, tensor::spaceTimePredictor::size()> stp{};
      alignas(PagesizeStack) std::array<real, tensor::I::size()> integrated{};
      alignas(PagesizeStack) std::array<real, tensor::Q::size()> dofs{};

      // scaled the way the quantities are, so that the residual is not read off
      // a field the material barely moves
      const std::array<real, NQ> factor{{1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
      auto viewQ = init::Q::view::create(dofs.data());
      for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
        for (std::size_t mode = 0; mode < Basis; ++mode) {
          viewQ(mode, quantity) = unit(rng) * factor[quantity];
        }
      }

      kernel::spaceTimePredictor krnl;
      krnl.bindGlobals(seissol::Pool::host());
      kernels::bindStarOperands(krnl, localIntegration);
      kernels::bindSourceDeviationOperands(krnl, localIntegration);
      for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
        auto zinv = init::Zinv::view<0>::create(zMatrix[quantity].data());
        seissol::model::calcZinv(
            zinv, sourceView, quantity, seissol::model::isStiffRow<Material>(quantity), Dt);
        krnl.Zinv(quantity) = zMatrix[quantity].data();
      }
      for (std::size_t i = 0; i < Material::StiffSourceRows.size(); ++i) {
        const auto& row = Material::StiffSourceRows[i];
        krnl.G(i) = sourceView(row.quantity, row.target) * Dt;
      }
      krnl.Q = dofs.data();
      krnl.I = integrated.data();
      krnl.timestep = Dt;
      krnl.spaceTimePredictor = stp.data();
      krnl.execute();

      // The residual of the system, with the operator assembled where the
      // material is sampled: project a mode onto the points, apply the star and
      // the deviation of the source there, and project back.
      const auto evaluate =
          denseOf(init::materialEval::view::create(const_cast<real*>(init::materialEval::Values)),
                  Points,
                  Basis);
      const auto project = denseOf(
          init::materialProject::view::create(const_cast<real*>(init::materialProject::Values)),
          Basis,
          Points);
      const auto timeOperator =
          denseOf(init::Z::view::create(const_cast<real*>(init::Z::Values)), Order, Order);
      const auto wHat = init::wHat::Values;
      const auto cellSource = sourceOf(cellMaterial);

      auto viewStp = init::spaceTimePredictor::view::create(stp.data());
      Eigen::MatrixXd field = Eigen::MatrixXd::Zero(Basis, NQ * Order);
      for (std::size_t mode = 0; mode < Basis; ++mode) {
        for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
          for (std::size_t step = 0; step < Order; ++step) {
            field(mode, quantity * Order + step) = viewStp(mode, quantity, step);
          }
        }
      }

      Eigen::MatrixXd deviationTerm = Eigen::MatrixXd::Zero(Basis, NQ * Order);
      // the left side: the time operator, and the source the cell carries
      Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(Basis, NQ * Order);
      for (std::size_t mode = 0; mode < Basis; ++mode) {
        for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
          for (std::size_t step = 0; step < Order; ++step) {
            double value = 0.0;
            for (std::size_t inner = 0; inner < Order; ++inner) {
              value += timeOperator(step, inner) * field(mode, quantity * Order + inner);
            }
            for (std::size_t other = 0; other < NQ; ++other) {
              value -= Dt * cellSource(other, quantity) * field(mode, other * Order + step);
            }
            residual(mode, quantity * Order + step) = value;
          }
        }
      }

      // the right side: the load, the flux at the samples, and what the samples
      // ask of the source beyond the cell
      for (std::size_t mode = 0; mode < Basis; ++mode) {
        for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
          for (std::size_t step = 0; step < Order; ++step) {
            residual(mode, quantity * Order + step) -= wHat[step] * viewQ(mode, quantity);
          }
        }
      }
      // every member of a family has a layout of its own
      const auto addDirection = [&](auto tag) {
        constexpr std::size_t dim = decltype(tag)::value;
        const auto stiffness =
            denseOf(init::kDivMT::view<dim>::create(const_cast<real*>(init::kDivMT::Values[dim])),
                    Basis,
                    Basis);
        const Eigen::MatrixXd atPoints = evaluate * stiffness * field;
        Eigen::MatrixXd applied = Eigen::MatrixXd::Zero(Points, NQ * Order);
        for (std::size_t point = 0; point < Points; ++point) {
          const auto star = starOf(sampled[point], gradients[dim]);
          for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
            for (std::size_t step = 0; step < Order; ++step) {
              double value = 0.0;
              for (std::size_t other = 0; other < NQ; ++other) {
                value += atPoints(point, other * Order + step) * star(other, quantity);
              }
              applied(point, quantity * Order + step) = value;
            }
          }
        }
        residual -= Dt * project * applied;
      };
      addDirection(std::integral_constant<std::size_t, 0>{});
      addDirection(std::integral_constant<std::size_t, 1>{});
      addDirection(std::integral_constant<std::size_t, 2>{});
      {
        const Eigen::MatrixXd atPoints = evaluate * field;
        Eigen::MatrixXd applied = Eigen::MatrixXd::Zero(Points, NQ * Order);
        for (std::size_t point = 0; point < Points; ++point) {
          const Eigen::MatrixXd deviation = sourceOf(sampled[point]) - cellSource;
          for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
            for (std::size_t step = 0; step < Order; ++step) {
              double value = 0.0;
              for (std::size_t other = 0; other < NQ; ++other) {
                value += atPoints(point, other * Order + step) * deviation(other, quantity);
              }
              applied(point, quantity * Order + step) = value;
            }
          }
        }
        deviationTerm = Dt * project * applied;
        residual -= deviationTerm;
      }

      // The load is what the residual is measured against, quantity by
      // quantity: they differ by nine orders of magnitude here, so one number
      // over all of them would only ever report the stresses.
      double relative = 0.0;
      double deviationRelative = 0.0;
      for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
        double load = 0.0;
        double worst = 0.0;
        double term = 0.0;
        for (std::size_t mode = 0; mode < Basis; ++mode) {
          for (std::size_t step = 0; step < Order; ++step) {
            load = std::max(load, std::abs(wHat[step] * viewQ(mode, quantity)));
            worst = std::max(worst, std::abs(residual(mode, quantity * Order + step)));
            term = std::max(term, std::abs(deviationTerm(mode, quantity * Order + step)));
          }
        }
        if (load > 0.0) {
          relative = std::max(relative, worst / load);
          deviationRelative = std::max(deviationRelative, term / load);
        }
      }
      // A material that does not vary leaves nothing: the iteration is then the
      // block sweep it replaces and it is exact. One that varies leaves what it
      // has not caught up with, and every step of the iteration carries one
      // more factor of how far the samples are from the cell.
      const double bar = variation > 0.0 ? 1e-4 : (std::is_same_v<real, double> ? 1e-12 : 1e-5);
      INFO("variation " << variation << ", relative residual " << relative << ", deviation term "
                        << deviationRelative);
      REQUIRE(relative < bar);
      // and the term is worth carrying: without it the residual would be its
      // own size, which is far above what the iteration leaves behind
      if (variation > 0.0) {
        REQUIRE(deviationRelative > 10.0 * bar);
      }
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_KERNEL_NODALSPACETIME_T_H_
