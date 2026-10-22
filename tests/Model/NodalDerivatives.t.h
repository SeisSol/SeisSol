// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_NODALDERIVATIVES_T_H_
#define SEISSOL_TESTS_MODEL_NODALDERIVATIVES_T_H_

// The time derivatives of the predictor where the material varies inside a cell, as the kernel a
// cell runs computes them:
//
// - with one material for the cell they are what the matrices give, contracted by hand, also where
//   the material is one per cell and the chain goes through the basis of a lower degree. Where the
//   build lets the k-th derivative keep only the polynomials up to the degree
//   ConvergenceOrder - 1 - k (generated::NarrowedDerivatives), this is where one dropped that a
//   constant operator gives would show up;
// - with a material that differs from point to point they are the chain formed by hand at the
//   points: the derivative read at the points, multiplied by the operator there, and projected
//   back onto the polynomials the derivative keeps.

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
#include "Initializer/Model/CurvedCell.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Functions.h"
#include "Numerical/Projection.h"
#include "Numerical/Quadrature.h"
#include "Solver/MultipleSimulations.h"
#include "VolumeTerm.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cstddef>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>
#include <yateto.h>

// The anelastic solver keeps its own chain, with the relaxation in tensors of
// its own; there is no single derivative to compare here.
#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

namespace seissol::unit_test {

namespace nodalderivatives {

/// Makes the kernel name dependent on a template parameter, so that a build
/// without a varying material never instantiates a body that names operands it
/// does not have.
template <bool Enabled>
struct Kernels {
  using Derivative = std::conditional_t<Enabled, kernel::derivative, kernel::derivative>;
};

/// Whether a kernel reads the source term as one matrix of the cell (ET), where the build carries
/// it that way.
template <typename KernelT, typename = void>
struct TakesSourceMatrix : std::false_type {};

template <typename KernelT>
struct TakesSourceMatrix<KernelT, std::void_t<decltype(std::declval<KernelT>().ET)>>
    : std::true_type {};

template <typename ViewT>
auto denseOf(const ViewT& view, std::size_t rows, std::size_t columns) -> Eigen::MatrixXd {
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

/// a global matrix as the mathematics has it: a build that bundles simulations
/// stores the matrices from the matrix files the other way round
template <typename ViewT>
auto mathMatrix(const ViewT& view, std::size_t rows, std::size_t columns) -> Eigen::MatrixXd {
  const Eigen::MatrixXd stored = denseOf(view, rows, columns);
  return multisim::NumSimulations > 1 ? Eigen::MatrixXd(stored.transpose()) : stored;
}

/// calls `function` with every member index of a family, as a constant
template <typename F, std::size_t... Ks>
void forEachMember(F&& function, std::index_sequence<Ks...> /*members*/) {
  (function(std::integral_constant<std::size_t, Ks>{}), ...);
}

/// The projection of what is given at the samples onto the polynomials up to `degree`, in the
/// basis of the build: the L2 projection with the weights of the rule of the samples, of which
/// `project` is the one onto every function of the cell. It is formed from the Dubiner functions up
/// to that degree, which span those polynomials whatever the basis, and the mass matrix of the
/// cell.
inline auto projectionOntoDegree(std::size_t degree, const Eigen::MatrixXd& project)
    -> Eigen::MatrixXd {
  constexpr std::size_t Points = tensor::materialEval::Shape[0];
  constexpr std::size_t Basis = tensor::materialEval::Shape[1];
  const auto nodes =
      init::materialNodes::view::create(const_cast<real*>(init::materialNodes::Values));
  const auto indices = numerical::projection::modalIndices<Cell::Dim>(degree + 1);
  Eigen::MatrixXd lower(Points, indices.size());
  for (std::size_t point = 0; point < Points; ++point) {
    const std::array<double, Cell::Dim> xi{nodes(point, 0), nodes(point, 1), nodes(point, 2)};
    for (std::size_t function = 0; function < indices.size(); ++function) {
      lower(point, function) = functions::DubinerP<Cell::Dim>(indices[function], xi);
    }
  }
  // those functions in the basis of the cell, which the projection onto all of it keeps
  const Eigen::MatrixXd inCell = project * lower;
  // a rule of ConvergenceOrder points per direction integrates the product of two functions
  const auto rule = seissol::quadrature::simplexRule<Cell::Dim>(ConvergenceOrder);
  Eigen::MatrixXd mass = Eigen::MatrixXd::Zero(Basis, Basis);
  for (std::size_t q = 0; q < rule.first.size(); ++q) {
    const basisFunction::SampledBasisFunctions<double> sampled(
        ConvergenceOrder, rule.first[q][0], rule.first[q][1], rule.first[q][2]);
    const Eigen::Map<const Eigen::VectorXd> values(sampled.data().data(), Basis);
    mass += rule.second[q] * values * values.transpose();
  }
  return inCell * (inCell.transpose() * mass * inCell).inverse() * inCell.transpose() * mass *
         project;
}

/// The basis of a degree in the basis of the cell, column by column: its functions sampled at the
/// samples of the material and projected with `project`, which keeps a polynomial of the degree of
/// the cell. What raises a derivative the build keeps in the basis of its degree
/// (generated::DerivativesByDegree).
inline auto embeddingOfDegree(std::size_t degree, const Eigen::MatrixXd& project)
    -> Eigen::MatrixXd {
  constexpr std::size_t Points = tensor::materialEval::Shape[0];
  const auto nodes =
      init::materialNodes::view::create(const_cast<real*>(init::materialNodes::Values));
  const auto indices = numerical::projection::modalIndices<Cell::Dim>(degree + 1);
  Eigen::MatrixXd sampled(Points, indices.size());
  for (std::size_t point = 0; point < Points; ++point) {
    const std::array<double, Cell::Dim> xi{nodes(point, 0), nodes(point, 1), nodes(point, 2)};
    for (std::size_t function = 0; function < indices.size(); ++function) {
      sampled(point, function) =
          functions::BasisP<Config::Basis, Cell::Dim>(degree, indices[function], xi);
    }
  }
  return project * sampled;
}

/// What the kernel gives for one material per sample point, against `reference`, which forms
/// the derivatives of the field `dofs` from the operators per point (the star matrices of the
/// three directions, and the source) and is handed the field of one simulation.
template <bool Enabled, typename MaterialT, typename ReferenceT>
void compareChain(const std::vector<MaterialT>& materials,
                  const Eigen::Matrix3d& gradients,
                  std::mt19937& rng,
                  const ReferenceT& reference) {
  constexpr std::size_t NQ = tensor::star::Shape[0][0];
  constexpr std::size_t Columns = tensor::star::Shape[0][1];
  constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
  constexpr std::size_t Derivatives = yateto::numFamilyMembers<tensor::dQ>();
  static_assert(Derivatives == ConvergenceOrder);
  std::normal_distribution<double> gauss(0.0, 1.0);

  LocalIntegrationData localIntegration{};
  initializer::CurvedCell::setConstantGradients(gradients, localIntegration.referenceGradients);
  for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
    const auto starCoefficients = seissol::model::getStarCoefficients(materials[point]);
    for (std::size_t i = 0; i < starCoefficients.size(); ++i) {
      localIntegration.materialCoefficients[i][point] = starCoefficients[i];
    }
    if constexpr (NodalSource) {
      const auto sourceCoefficients = seissol::model::getSourceCoefficients(materials[point]);
      for (std::size_t i = 0; i < sourceCoefficients.size(); ++i) {
        localIntegration.sourceCoefficients[i][point] = sourceCoefficients[i];
      }
    }
  }

  alignas(Alignment) std::array<real, tensor::Q::size()> dofs{};
  for (auto& value : dofs) {
    value = static_cast<real>(gauss(rng));
  }
  alignas(Alignment) std::array<real, yateto::computeFamilySize<tensor::dQ>()> derivatives{};
  alignas(Alignment) std::array<real, tensor::I::size()> integrated{};
  std::array<double, Derivatives> powers{};
  for (auto& power : powers) {
    power = gauss(rng);
  }

  typename Kernels<Enabled>::Derivative krnl{};
  krnl.bindGlobals(seissol::Pool::host());
  kernels::bindStarOperands(krnl, localIntegration);
  kernels::bindSourceOperands(krnl, localIntegration);
  // a kernel that reads the source term as one matrix reads the one of the material of the cell,
  // which is one material for the cell there
  alignas(Alignment) std::array<real, tensor::ET::size()> sourceMatrix{};
  if constexpr (TakesSourceMatrix<typename Kernels<Enabled>::Derivative>::value) {
    auto sourceView = init::ET::view::create(sourceMatrix.data());
    seissol::model::getTransposedSourceCoefficientTensor(materials[0], sourceView);
    krnl.ET = sourceMatrix.data();
  }
  krnl.dQ(0) = dofs.data();
  for (std::size_t k = 1; k < Derivatives; ++k) {
    krnl.dQ(k) = derivatives.data() + yateto::computeFamilySize<tensor::dQ>(1, k);
  }
  krnl.I = integrated.data();
  for (std::size_t k = 0; k < Derivatives; ++k) {
    krnl.power(k) = powers[k];
  }
  // where the build forms the volume term of the corrector in the chain, the chain writes it too
  kernels::LocalTmp tmp(0.0);
  bindVolumeFused(krnl, tmp);
  krnl.execute();

  auto viewQ = init::Q::view::create(dofs.data());
  auto viewI = init::I::view::create(integrated.data());
  // the projection from the samples of the material, which raises a derivative kept in the basis of
  // its degree into the basis of the cell
  const Eigen::MatrixXd project =
      denseOf(init::materialProject::view::create(const_cast<real*>(init::materialProject::Values)),
              Basis,
              tensor::materialEval::Shape[0]);
  for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
    auto slicedQ = multisim::simtensor(viewQ, simulation);
    Eigen::MatrixXd field = Eigen::MatrixXd::Zero(Basis, NQ);
    for (std::size_t row = 0; row < Basis; ++row) {
      for (std::size_t column = 0; column < NQ; ++column) {
        field(row, column) = slicedQ(row, column);
      }
    }
    const std::vector<Eigen::MatrixXd> expected = reference(field);
    REQUIRE(expected.size() == Derivatives);

    // entries of a derivative cancel against each other, so what an entry may be off by follows
    // the size of its column, which differ by orders of magnitude between stresses and velocities
    // -- and a column that vanishes is held to the size of the whole derivative
    const auto check = [](const Eigen::MatrixXd& computed, const Eigen::MatrixXd& wanted) {
      const double whole = wanted.cwiseAbs().maxCoeff();
      for (std::size_t column = 0; column < Columns; ++column) {
        const double scale = std::max(1e-10 * whole, wanted.col(column).cwiseAbs().maxCoeff());
        for (std::size_t row = 0; row < Basis; ++row) {
          REQUIRE(computed(row, column) == doctest::Approx(wanted(row, column))
                                               .epsilon(std::is_same_v<real, double> ? 1e-10 : 1e-3)
                                               .scale(scale));
        }
      }
    };

    Eigen::MatrixXd expectedIntegral = powers[0] * expected[0];
    // every member of a family has a layout of its own, so it has to be read
    // through its own view
    const auto checkDerivative = [&](auto tag) {
      constexpr std::size_t K = decltype(tag)::value;
      if constexpr (K > 0 && K < Derivatives) {
        auto view = init::dQ::view<K>::create(derivatives.data() +
                                              yateto::computeFamilySize<tensor::dQ>(1, K));
        auto sliced = multisim::simtensor(view, simulation);
        Eigen::MatrixXd computed = denseOf(sliced, Basis, Columns);
        if constexpr (generated::DerivativesByDegree) {
          // the derivative in the basis of its degree, in the leading rows
          const auto embedding = embeddingOfDegree(ConvergenceOrder - 1 - K, project);
          const Eigen::MatrixXd raised = embedding * computed.topRows(embedding.cols());
          computed = raised;
        }
        const std::size_t derivative = K;
        CAPTURE(derivative);
        check(computed, expected[K]);
        expectedIntegral += powers[K] * expected[K];
      }
    };
    forEachMember(checkDerivative, std::make_index_sequence<Derivatives>{});

    auto slicedI = multisim::simtensor(viewI, simulation);
    check(denseOf(slicedI, Basis, Columns), expectedIntegral);
  }
}

/// the operators of one material: the star matrices of the three directions for the rows of the
/// inverse Jacobian `gradients`, and the source
template <typename MaterialT>
auto operatorsOf(const MaterialT& material, const Eigen::Matrix3d& gradients)
    -> std::array<Eigen::MatrixXd, Cell::Dim + 1> {
  constexpr std::size_t NQ = tensor::star::Shape[0][0];
  constexpr std::size_t Columns = tensor::star::Shape[0][1];
  std::array<Eigen::MatrixXd, Cell::Dim> directional;
  for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
    alignas(Alignment) std::array<real, tensor::star::size(0)> data{};
    auto view = init::star::view<0>::create(data.data());
    seissol::model::getTransposedCoefficientMatrix(material, dim, view);
    directional[dim] = denseOf(view, NQ, Columns);
  }
  std::array<Eigen::MatrixXd, Cell::Dim + 1> operators;
  for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
    operators[dim] = Eigen::MatrixXd::Zero(NQ, Columns);
    for (std::size_t component = 0; component < Cell::Dim; ++component) {
      operators[dim] += gradients(dim, component) * directional[component];
    }
  }
  alignas(Alignment) std::array<real, tensor::ET::size()> sourceData{};
  auto sourceView = init::ET::view::create(sourceData.data());
  seissol::model::getTransposedSourceCoefficientTensor(material, sourceView);
  operators[Cell::Dim] = denseOf(sourceView, NQ, Columns);
  return operators;
}

/// the derivative matrices the chain applies, as the mathematics has them
inline auto stiffnessMatrices() -> std::array<Eigen::MatrixXd, Cell::Dim> {
  std::array<Eigen::MatrixXd, Cell::Dim> stiffness;
  const auto read = [&](auto tag) {
    constexpr std::size_t Dim = decltype(tag)::value;
    stiffness[Dim] =
        mathMatrix(init::kDivMT::view<Dim>::create(const_cast<real*>(init::kDivMT::Values[Dim])),
                   tensor::kDivMT::Shape[Dim][0],
                   tensor::kDivMT::Shape[Dim][1]);
  };
  read(std::integral_constant<std::size_t, 0>{});
  read(std::integral_constant<std::size_t, 1>{});
  read(std::integral_constant<std::size_t, 2>{});
  return stiffness;
}

inline auto randomGradients(std::mt19937& rng) -> Eigen::Matrix3d {
  std::normal_distribution<double> gauss(0.0, 1.0);
  Eigen::Matrix3d gradients;
  for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
    for (std::size_t component = 0; component < Cell::Dim; ++component) {
      gradients(dim, component) = gauss(rng);
    }
  }
  return gradients;
}

template <bool Enabled>
void constantMaterial() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261013);
    const auto stiffness = stiffnessMatrices();
    for (std::size_t sample = 0; sample < 4; ++sample) {
      const auto material = coefficients::configuredMaterial<Material>(rng);
      const auto gradients = randomGradients(rng);
      const auto operators = operatorsOf(material, gradients);
      // the chain of the matrix form, every mode of every derivative: a derivative the kernel
      // keeps fewer modes of has to be zero beyond them
      const auto reference = [&](const Eigen::MatrixXd& field) {
        std::vector<Eigen::MatrixXd> chain{field};
        for (std::size_t k = 1; k < ConvergenceOrder; ++k) {
          Eigen::MatrixXd next = chain.back() * operators[Cell::Dim];
          for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
            next += stiffness[dim] * chain.back() * operators[dim];
          }
          chain.push_back(next);
        }
        return chain;
      };
      compareChain<Enabled>(
          std::vector<Material>(MaterialSampleCount, material), gradients, rng, reference);
    }
  }
}

template <bool Enabled>
void varyingMaterial() {
  if constexpr (Enabled) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t Points = tensor::materialEval::Shape[0];
    constexpr std::size_t Basis = tensor::materialEval::Shape[1];
    static_assert(Points == MaterialSampleCount);
    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261014);
    const auto stiffness = stiffnessMatrices();
    // the material points and their weights: read a field there, project back from there
    const Eigen::MatrixXd evaluate =
        denseOf(init::materialEval::view::create(const_cast<real*>(init::materialEval::Values)),
                Points,
                Basis);
    const Eigen::MatrixXd project = denseOf(
        init::materialProject::view::create(const_cast<real*>(init::materialProject::Values)),
        Basis,
        Points);
    // what the k-th derivative is projected onto
    std::vector<Eigen::MatrixXd> projections{project};
    for (std::size_t k = 1; k < ConvergenceOrder; ++k) {
      projections.push_back(generated::NarrowedDerivatives
                                ? projectionOntoDegree(ConvergenceOrder - 1 - k, project)
                                : project);
    }

    for (std::size_t sample = 0; sample < 3; ++sample) {
      std::vector<Material> materials;
      for (std::size_t point = 0; point < Points; ++point) {
        materials.push_back(coefficients::configuredMaterial<Material>(rng));
      }
      const auto gradients = randomGradients(rng);
      std::vector<std::array<Eigen::MatrixXd, Cell::Dim + 1>> operators;
      for (const auto& material : materials) {
        operators.push_back(operatorsOf(material, gradients));
      }
      const auto reference = [&](const Eigen::MatrixXd& field) {
        std::vector<Eigen::MatrixXd> chain{field};
        for (std::size_t k = 1; k < ConvergenceOrder; ++k) {
          const Eigen::MatrixXd& previous = chain.back();
          Eigen::MatrixXd product = Eigen::MatrixXd::Zero(Points, previous.cols());
          const Eigen::MatrixXd atPoints = evaluate * previous;
          std::array<Eigen::MatrixXd, Cell::Dim> derivativeAtPoints;
          for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
            derivativeAtPoints[dim] = evaluate * stiffness[dim] * previous;
          }
          for (std::size_t point = 0; point < Points; ++point) {
            product.row(point) = atPoints.row(point) * operators[point][Cell::Dim];
            for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
              product.row(point) += derivativeAtPoints[dim].row(point) * operators[point][dim];
            }
          }
          chain.push_back(projections[k] * product);
        }
        return chain;
      };
      compareChain<Enabled>(materials, gradients, rng, reference);
    }
  }
}

} // namespace nodalderivatives

TEST_CASE("Time derivatives against the matrix form") {
  // every build that carries the coefficients of its operator binds the chain the same way
  nodalderivatives::constantMaterial<FactoredStar>();
}

TEST_CASE("Nodal time derivatives where the material varies inside the cell") {
  nodalderivatives::varyingMaterial<NodalMaterial>();
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_NODALDERIVATIVES_T_H_
