// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_CURVEDCELL_T_H_
#define SEISSOL_TESTS_MODEL_CURVEDCELL_T_H_

// What a cell that may be curved carries, and what the kernels make of it.
//
// - A straight-sided cell evaluated at every point carries what it carries as one cell: the
//   metric at every operator point, the rotation and the scale at every node of a face.
// - On a curved face, the scale is the surface Jacobian over the Jacobian determinant and the
//   rotation the frame of the face at the node -- here measured by differencing the map.
// - A field that is linear in space is a polynomial on a curved cell of order two, and the
//   volume term has to give its derivative exactly. With the metric of one point for the whole
//   cell it does not, so this is what checks that the metric is the one at the point.
// - A constant state stays constant on a curved cell whose material varies.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/IsoparametricTransform.h"
#include "Initializer/Model/CellFlux.h"
#include "Initializer/Model/CurvedCell.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"
#include "Numerical/Quadrature.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <type_traits>
#include <vector>

// the anelastic solver writes into the extended quantities, see NodalCorrector.t.h
#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

namespace seissol::unit_test {

namespace curvedcell {

/// Makes the kernel names dependent on a template parameter, so that a build whose cells are
/// straight-sided never instantiates a body that names operands it does not have.
template <bool Enabled>
struct Kernels {
  using Volume = std::conditional_t<Enabled, kernel::volume, kernel::volume>;
  using LocalFlux = std::conditional_t<Enabled, kernel::localFlux, kernel::localFlux>;
};

using initializer::CurvedCell;
using VectorT = seissol::geometry::CellTransform::VectorEigenT;
using MatrixT = seissol::geometry::CellTransform::MatrixEigenT;
using FaceVectorT = seissol::geometry::FaceTransform::FaceVectorT;

inline auto straightVertices() -> std::array<VectorT, Cell::NumVertices> {
  return {VectorT(0.1, 0.2, -0.3),
          VectorT(2.1, 0.0, 0.4),
          VectorT(-0.2, 1.7, 0.1),
          VectorT(0.3, 0.1, 2.2)};
}

/// the cell with every edge bent off the straight line by a few percent of its length
inline auto curvedCell() -> seissol::geometry::IsoparametricTransform {
  const auto vertices = straightVertices();
  constexpr std::array<std::array<std::size_t, 2>, 6> Edges{
      {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
  std::array<VectorT, 6> midpoints{};
  for (std::size_t edge = 0; edge < Edges.size(); ++edge) {
    midpoints[edge] = 0.5 * (vertices[Edges[edge][0]] + vertices[Edges[edge][1]]) +
                      0.06 * VectorT(1.0 - 0.3 * edge, 0.2 * edge - 0.5, 0.1 * edge * edge - 0.4);
  }
  return seissol::geometry::IsoparametricTransform::fromEdgeMidpoints(vertices, midpoints);
}

template <bool Enabled>
void straightCell() {
  if constexpr (Enabled) {
    const auto cell = seissol::geometry::AffineTransform(straightVertices());
    LocalIntegrationData pointwise{};
    LocalIntegrationData constant{};
    CurvedCell::setGradients(cell, pointwise.referenceGradients);
    CurvedCell::setConstantGradients(
        cell.refToSpaceJacobianInverse(VectorT(Cell::ReferenceBarycenter.data())),
        constant.referenceGradients);
    for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
      for (std::size_t i = 0; i < sizeof(pointwise.referenceGradients[dim]) / sizeof(real); ++i) {
        REQUIRE(pointwise.referenceGradients[dim][i] == constant.referenceGradients[dim][i]);
      }
    }

    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const seissol::geometry::AffineFaceTransform face(cell,
                                                        seissol::geometry::ReferenceFaceMap(side));
      std::array<double, FluxFaceNodes> scale{};
      CurvedCell::setFace(cell, face, side, pointwise.faceRotation[side], scale);
      for (const auto value : scale) {
        REQUIRE(value == doctest::Approx(-2.0 * face.area() / cell.determinant()).epsilon(1e-13));
      }

      const auto basis = face.faceAlignedBasis();
      alignas(Alignment) std::array<real, tensor::T::size()> matT{};
      alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
      auto viewT = init::T::view::create(matT.data());
      auto viewTinv = init::Tinv::view::create(matTinv.data());
      seissol::model::getFaceRotationMatrix(
          basis[0].normalized(), basis[1].normalized(), basis[2].normalized(), viewT, viewTinv);
      CurvedCell::setConstantRotation(matT.data(), constant.faceRotation[side]);
      for (std::size_t i = 0; i < sizeof(pointwise.faceRotation[side]) / sizeof(real); ++i) {
        REQUIRE(pointwise.faceRotation[side][i] ==
                AbsApprox(constant.faceRotation[side][i]).epsilon(1e-13));
      }
    }
  }
}

template <bool Enabled>
void curvedFace() {
  if constexpr (Enabled) {
    const auto cell = curvedCell();
    LocalIntegrationData local{};
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const seissol::geometry::IsoparametricFaceTransform face(
          cell, seissol::geometry::ReferenceFaceMap(side));
      std::array<double, FluxFaceNodes> scale{};
      CurvedCell::setFace(cell, face, side, local.faceRotation[side], scale);
      auto rotations = init::TNodes::view::create(local.faceRotation[side]);

      const VectorT center = cell.refToSpace(VectorT(Cell::ReferenceBarycenter.data()));
      const auto corners = face.vertices();
      for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
        const VectorT reference = CurvedCell::faceNode(side, node);
        const FaceVectorT point = face.map().cellToFace(reference);

        // the normal and the surface Jacobian by differencing the map of the face
        constexpr double Step = 1e-6;
        const VectorT alongChi = (face.refToSpace(FaceVectorT(point + FaceVectorT(Step, 0))) -
                                  face.refToSpace(FaceVectorT(point - FaceVectorT(Step, 0)))) /
                                 (2 * Step);
        const VectorT alongTau = (face.refToSpace(FaceVectorT(point + FaceVectorT(0, Step))) -
                                  face.refToSpace(FaceVectorT(point - FaceVectorT(0, Step)))) /
                                 (2 * Step);
        VectorT normal = alongChi.cross(alongTau);
        const double determinant = cell.refToSpaceJacobian(reference).determinant();
        REQUIRE(determinant > 0);
        REQUIRE(scale[node] == doctest::Approx(-normal.norm() / determinant).epsilon(1e-8));

        // the face points away from the cell
        REQUIRE(normal.dot(face.refToSpace(point) - center) > 0);

        // and the rotation is the one into the frame the face has at the node
        normal.normalize();
        const VectorT edge = corners[1] - corners[0];
        const VectorT tangent1 = (edge - edge.dot(normal) * normal).normalized();
        const VectorT tangent2 = normal.cross(tangent1);
        alignas(Alignment) std::array<real, tensor::T::size()> matT{};
        alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
        auto viewT = init::T::view::create(matT.data());
        auto viewTinv = init::Tinv::view::create(matTinv.data());
        seissol::model::getFaceRotationMatrix(normal, tangent1, tangent2, viewT, viewTinv);
        for (std::size_t row = 0; row < tensor::T::Shape[0]; ++row) {
          for (std::size_t column = 0; column < tensor::T::Shape[1]; ++column) {
            if (viewT.isInRange(row, column)) {
              REQUIRE(rotations(node, row, column) == AbsApprox(viewT(row, column)).epsilon(1e-7));
            }
          }
        }
      }
    }
  }
}

/// the values of a field at the quadrature points the initial state is read at, back to the modes
template <typename FieldT>
void project(const seissol::geometry::CellTransform& cell, const FieldT& field, real* dofs) {
  constexpr std::size_t NQ = tensor::star::Shape[0][0];
  constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
  const auto rule = seissol::quadrature::simplexRule<Cell::Dim>(ConvergenceOrder + 1);
  std::vector<Eigen::VectorXd> values;
  values.reserve(rule.first.size());
  for (const auto& point : rule.first) {
    values.emplace_back(field(cell.refToSpace(VectorT(point.data()))));
  }
  const auto projection = init::projectQP::view::create(const_cast<real*>(init::projectQP::Values));
  auto view = init::I::view::create(dofs);
  for (std::size_t mode = 0; mode < Basis; ++mode) {
    for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
      double value = 0;
      for (std::size_t point = 0; point < values.size(); ++point) {
        if (projection.isInRange(mode, point)) {
          value += projection(mode, point) * values[point](quantity);
        }
      }
      view(mode, quantity) = static_cast<real>(value);
    }
  }
}

template <bool Enabled>
void linearField() {
  // the global matrices of a build that bundles simulations are stored the other way round, and
  // the projection here reads them as they are
  if constexpr (Enabled && multisim::NumSimulations == 1 && !NodalSource) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t NQ = tensor::star::Shape[0][0];
    constexpr std::size_t Columns = tensor::star::Shape[0][1];
    constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
    // the field has the degree of the map, which a basis of this order has to hold
    static_assert(ConvergenceOrder >= 3, "The basis does not hold a field of degree two.");

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261008);
    std::normal_distribution<double> gauss(0.0, 1.0);
    const auto cell = curvedCell();
    const auto material = coefficients::configuredMaterial<Material>(rng);

    LocalIntegrationData local{};
    CurvedCell::setGradients(cell, local.referenceGradients);
    const auto coefficients = seissol::model::getStarCoefficients(material);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
      for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
        local.materialCoefficients[i][point] = coefficients[i];
      }
    }

    // q(x) = a + B x, quantity by quantity
    Eigen::VectorXd offset(NQ);
    Eigen::MatrixXd slope(NQ, Cell::Dim);
    for (std::size_t quantity = 0; quantity < NQ; ++quantity) {
      offset(quantity) = gauss(rng);
      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        slope(quantity, d) = gauss(rng);
      }
    }
    const auto field = [&](const VectorT& x) -> Eigen::VectorXd { return offset + slope * x; };

    alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
    project(cell, field, dofs.data());
    alignas(Alignment) std::array<real, tensor::Q::size()> update{};

    typename Kernels<Enabled>::Volume volume{};
    volume.bindGlobals(seissol::Pool::host());
    volume.I = dofs.data();
    volume.Q = update.data();
    kernels::bindStarOperands(volume, local);
    volume.execute();

    // what the strong form gives for it: the derivative applied with the star of every
    // direction, which is a constant, with the sign of the term
    Eigen::RowVectorXd constant = Eigen::RowVectorXd::Zero(Columns);
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      alignas(Alignment) std::array<real, tensor::star::size(0)> data{};
      auto star = init::star::view<0>::create(data.data());
      seissol::model::getTransposedCoefficientMatrix(material, d, star);
      for (std::size_t row = 0; row < NQ; ++row) {
        for (std::size_t column = 0; column < Columns; ++column) {
          if (star.isInRange(row, column)) {
            constant(column) -= slope(row, d) * star(row, column);
          }
        }
      }
    }
    alignas(Alignment) std::array<real, tensor::I::size()> expected{};
    project(
        cell,
        [&](const VectorT& /*x*/) -> Eigen::VectorXd {
          Eigen::VectorXd value = Eigen::VectorXd::Zero(NQ);
          for (std::size_t quantity = 0; quantity < std::min(NQ, Columns); ++quantity) {
            value(quantity) = constant(quantity);
          }
          return value;
        },
        expected.data());

    auto viewQ = init::Q::view::create(update.data());
    auto viewExpected = init::I::view::create(expected.data());
    const double scale = std::max(1.0, constant.cwiseAbs().maxCoeff());
    for (std::size_t mode = 0; mode < Basis; ++mode) {
      for (std::size_t quantity = 0; quantity < std::min(NQ, Columns); ++quantity) {
        REQUIRE(viewQ(mode, quantity) ==
                doctest::Approx(viewExpected(mode, quantity)).epsilon(1e-9).scale(scale));
      }
    }
  }
}

template <bool Enabled>
void constantState() {
  if constexpr (Enabled && !NodalSource) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t Columns = tensor::star::Shape[0][1];
    constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261009);
    std::normal_distribution<double> gauss(0.0, 1.0);
    const auto cell = curvedCell();
    const auto base = coefficients::configuredMaterial<Material>(rng);
    const auto materialAt = [&](const VectorT& x) {
      Material material = base;
      double phase = 0.0;
      for (const auto& [name, member] : Material::ParameterMap) {
        material.*member =
            base.*member * (1.0 + 0.2 * std::sin(0.9 * x(0) - 0.7 * x(1) + 1.1 * x(2) + phase));
        phase += 1.3;
      }
      return material;
    };

    LocalIntegrationData local{};
    NeighboringIntegrationData neighboring{};
    CurvedCell::setGradients(cell, local.referenceGradients);
    const auto samples =
        init::materialNodes::view::create(const_cast<real*>(init::materialNodes::Values));
    for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
      VectorT reference = VectorT::Zero();
      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        reference(d) = samples.isInRange(point, d) ? samples(point, d) : 0;
      }
      const auto coefficients =
          seissol::model::getStarCoefficients(materialAt(cell.refToSpace(reference)));
      for (std::size_t i = 0; i < coefficients.size(); ++i) {
        local.materialCoefficients[i][point] = coefficients[i];
      }
    }
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const seissol::geometry::IsoparametricFaceTransform face(
          cell, seissol::geometry::ReferenceFaceMap(side));
      std::array<double, FluxFaceNodes> scale{};
      CurvedCell::setFace(cell, face, side, local.faceRotation[side], scale);
      for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
        const auto material = materialAt(cell.refToSpace(CurvedCell::faceNode(side, node)));
        std::array<double, FluxCoefficientCount> plus{};
        std::array<double, FluxCoefficientCount> minus{};
        initializer::fluxScalarsOfNode(material,
                                       material,
                                       FaceType::Regular,
                                       initializer::parameters::NumericalFlux::Godunov,
                                       scale[node],
                                       plus,
                                       minus);
        for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
          local.fluxCoefficients[side][c][node] = static_cast<real>(plus[c]);
          neighboring.fluxCoefficients[side][c][node] = static_cast<real>(minus[c]);
        }
      }
    }

    alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
    auto viewI = init::I::view::create(dofs.data());
    for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
      auto slicedI = multisim::simtensor(viewI, simulation);
      for (std::size_t quantity = 0; quantity < tensor::star::Shape[0][0]; ++quantity) {
        slicedI(0, quantity) = static_cast<real>(gauss(rng));
      }
    }

    alignas(Alignment) std::array<real, tensor::Q::size()> localOnly{};
    alignas(Alignment) std::array<real, tensor::Q::size()> update{};
    typename Kernels<Enabled>::LocalFlux flux{};
    flux.bindGlobals(seissol::Pool::host());
    flux.I = dofs.data();
    flux.Q = localOnly.data();
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      kernels::bindLocalFluxOperands(flux, local, side);
      flux.execute(side);
    }

    typename Kernels<Enabled>::Volume volume{};
    volume.bindGlobals(seissol::Pool::host());
    volume.I = dofs.data();
    volume.Q = update.data();
    kernels::bindStarOperands(volume, local);
    volume.execute();
    // the neighbour's trace is the cell's own where the state is continuous, so its contribution
    // is the local kernel applied with the neighbour's scalars
    flux.Q = update.data();
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      kernels::bindLocalFluxOperands(flux, local, side);
      flux.execute(side);
      for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
        flux.fluxCoefficientsLocal(c) = neighboring.fluxCoefficients[side][c];
      }
      flux.execute(side);
    }

    auto viewLocal = init::Q::view::create(localOnly.data());
    auto viewQ = init::Q::view::create(update.data());
    for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
      auto slicedLocal = multisim::simtensor(viewLocal, simulation);
      auto slicedQ = multisim::simtensor(viewQ, simulation);
      constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
      for (std::size_t column = 0; column < Columns; ++column) {
        double scale = 0.0;
        for (std::size_t row = 0; row < Basis; ++row) {
          scale = std::max(scale, std::abs(static_cast<double>(slicedLocal(row, column))));
        }
        for (std::size_t row = 0; row < Basis; ++row) {
          REQUIRE(std::abs(static_cast<double>(slicedQ(row, column))) <=
                  Tolerance * std::max(scale, 1.0));
        }
      }
    }
  }
}

} // namespace curvedcell

TEST_CASE("A straight-sided cell carries the same geometry at every point") {
  curvedcell::straightCell<Curvilinear>();
}

TEST_CASE("A curved face carries its scale and its frame at every node") {
  curvedcell::curvedFace<Curvilinear>();
}

TEST_CASE("The volume term differentiates a linear field exactly on a curved cell") {
  curvedcell::linearField<Curvilinear>();
}

TEST_CASE("A curved cell keeps a constant state constant") {
  curvedcell::constantState<Curvilinear>();
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_CURVEDCELL_T_H_
