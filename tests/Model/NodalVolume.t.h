// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_NODALVOLUME_T_H_
#define SEISSOL_TESTS_MODEL_NODALVOLUME_T_H_

// The volume kernel where the material varies inside a cell: the operator and
// the source term are both formed where the material is sampled. A material
// that does not vary has to give back exactly what the matrices give, and the
// reference here is those matrices, contracted by hand. A sample point read in
// the wrong order, a source placed at the wrong rows, or a projection left out
// shows up in that comparison.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <random>
#include <vector>

// The anelastic solver writes its volume contribution into the extended
// quantities and carries its source term in a tensor of its own, so this
// comparison has no kernel to run there.
#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

namespace seissol::unit_test {

namespace nodalvolume {

template <bool Enabled>
struct Kernels {
  using Volume = std::conditional_t<Enabled, kernel::volume, kernel::volume>;
};

template <bool Enabled>
void compareAgainstModal() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t NQ = seissol::tensor::star::Shape[0][0];
    constexpr std::size_t Columns = seissol::tensor::star::Shape[0][1];
    constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];

    std::mt19937 rng(20260927);
    std::normal_distribution<double> gauss(0.0, 1.0);

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
    // a build that bundles simulations stores the matrices from the matrix
    // files the other way round
    const auto mathMatrix = [&denseOf](auto view, std::size_t rows, std::size_t columns) {
      const Eigen::MatrixXd stored = denseOf(view, rows, columns);
      return multisim::NumSimulations > 1 ? Eigen::MatrixXd(stored.transpose()) : stored;
    };

    for (std::size_t sample = 0; sample < 4; ++sample) {
      const auto material = coefficients::configuredMaterial<Material>(rng);

      double gradients[3][3];
      for (auto& row : gradients) {
        for (auto& value : row) {
          value = gauss(rng);
        }
      }

      // 1. what the cell carries: the same coefficients at every sample point,
      //    since the material does not vary here
      LocalIntegrationData localIntegration{};
      for (std::size_t dim = 0; dim < 3; ++dim) {
        for (std::size_t component = 0; component < 3; ++component) {
          localIntegration.referenceGradients[dim][component] = gradients[dim][component];
        }
      }
      const auto starCoefficients = seissol::model::getStarCoefficients(material);
      for (std::size_t i = 0; i < starCoefficients.size(); ++i) {
        for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
          localIntegration.materialCoefficients[i][point] = starCoefficients[i];
        }
      }
      if constexpr (NodalSource) {
        const auto sourceCoefficients = seissol::model::getSourceCoefficients(material);
        for (std::size_t i = 0; i < sourceCoefficients.size(); ++i) {
          for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
            localIntegration.sourceCoefficients[i][point] = sourceCoefficients[i];
          }
        }
      }

      alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
      for (auto& value : dofs) {
        value = static_cast<real>(gauss(rng));
      }
      alignas(Alignment) std::array<real, tensor::Q::size()> update{};

      typename Kernels<Enabled>::Volume krnl{};
      krnl.bindGlobals(seissol::Pool::host());
      krnl.I = dofs.data();
      krnl.Q = update.data();
      kernels::bindStarOperands(krnl, localIntegration);
      kernels::bindSourceOperands(krnl, localIntegration);
      krnl.execute();

      // 2. the same thing out of the matrices it is made of
      std::array<Eigen::MatrixXd, 3> directional;
      for (std::size_t dim = 0; dim < 3; ++dim) {
        alignas(Alignment) std::array<real, tensor::star::size(0)> data{};
        auto view = init::star::view<0>::create(data.data());
        seissol::model::getTransposedCoefficientMatrix(material, dim, view);
        directional[dim] = denseOf(view, NQ, Columns);
      }
      alignas(Alignment) std::array<real, tensor::ET::size()> sourceData{};
      auto sourceView = init::ET::view::create(sourceData.data());
      seissol::model::getTransposedSourceCoefficientTensor(material, sourceView);
      const auto source = denseOf(sourceView, NQ, Columns);

      auto viewI = init::I::view::create(dofs.data());
      auto viewQ = init::Q::view::create(update.data());
      for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
        auto slicedI = multisim::simtensor(viewI, simulation);
        auto slicedQ = multisim::simtensor(viewQ, simulation);

        Eigen::MatrixXd field = Eigen::MatrixXd::Zero(Basis, NQ);
        for (std::size_t row = 0; row < Basis; ++row) {
          for (std::size_t column = 0; column < NQ; ++column) {
            field(row, column) = slicedI(row, column);
          }
        }

        Eigen::MatrixXd expected = field * source;
        // every member of a family has a layout of its own, so it has to be
        // read through its own view
        const auto addDirection = [&](auto tag) {
          constexpr std::size_t Dim = decltype(tag)::value;
          Eigen::MatrixXd star = Eigen::MatrixXd::Zero(NQ, Columns);
          for (std::size_t component = 0; component < 3; ++component) {
            star += gradients[Dim][component] * directional[component];
          }
          const auto stiffness = mathMatrix(
              init::kDivM::view<Dim>::create(const_cast<real*>(init::kDivM::Values[Dim])),
              tensor::kDivM::Shape[Dim][0],
              tensor::kDivM::Shape[Dim][1]);
          expected += stiffness * field * star;
        };
        addDirection(std::integral_constant<std::size_t, 0>{});
        addDirection(std::integral_constant<std::size_t, 1>{});
        addDirection(std::integral_constant<std::size_t, 2>{});

        for (std::size_t row = 0; row < Basis; ++row) {
          for (std::size_t column = 0; column < Columns; ++column) {
            // entries of the update cancel against each other, so what an
            // entry may be off by follows the size of the update, not its own
            const double scale = std::max(1.0, expected.cwiseAbs().maxCoeff());
            REQUIRE(slicedQ(row, column) ==
                    doctest::Approx(expected(row, column)).epsilon(1e-11).scale(scale));
          }
        }
      }
    }
  }
}

} // namespace nodalvolume

TEST_CASE("Nodal volume against the matrix form") {
  nodalvolume::compareAgainstModal<NodalMaterial>();
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_NODALVOLUME_T_H_
