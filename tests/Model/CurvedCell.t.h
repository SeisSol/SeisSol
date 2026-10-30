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
// - On a curved face, the scale is the surface Jacobian over the Jacobian determinant, or the
//   surface Jacobian alone where the cell is weight-adjusted, and the rotation the frame of the
//   face at the node -- here measured by differencing the map.
// - A field that is linear in space is a polynomial on a curved cell of order two, and the
//   volume term of the predictor has to give its derivative exactly. With the metric of one point
//   for the whole cell it does not, so this is what checks that the metric is the one at the
//   point.
// - The corrector of a curved cell (kernels::CurvedSplit): on a straight-sided cell, the split
//   form with the weight applied is the strong form; on a curved cell, the split volume term does
//   no work in the energy, and a constant state stays constant where the material varies.

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "FaceTraces.h"
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
#include "Kernels/CurvedSplit.h"
#include "Kernels/Precision.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"
#include "Numerical/BasisFunction.h"
#include "Numerical/Quadrature.h"
#include "Solver/MultipleSimulations.h"
#include "VolumeTerm.h"

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
  using SplitVolume = std::conditional_t<Enabled, kernel::splitVolume, void>;
  using WeightAdjust = std::conditional_t<Enabled, kernel::weightAdjust, void>;
  using Split = kernels::CurvedSplit<Enabled, kernel::splitVolume, kernel::ownTraction>;
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

    // a weight-adjusted cell carries 1 / J in its weight, and its faces the surface Jacobian alone
    CurvedCell::setWeightScale(cell, pointwise.weightScale);
    for (std::size_t point = 0; point < OperatorPointCount; ++point) {
      REQUIRE(pointwise.weightScale[point] ==
              doctest::Approx(1.0 / std::abs(cell.determinant())).epsilon(1e-13));
    }

    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const seissol::geometry::AffineFaceTransform face(cell,
                                                        seissol::geometry::ReferenceFaceMap(side));
      std::array<double, FluxFaceNodes> scale{};
      CurvedCell::setFace(cell, face, side, false, pointwise.faceRotation[side], scale);
      for (const auto value : scale) {
        REQUIRE(value == doctest::Approx(-2.0 * face.area() / cell.determinant()).epsilon(1e-13));
      }
      CurvedCell::setFace(cell, face, side, true, pointwise.faceRotation[side], scale);
      for (const auto value : scale) {
        REQUIRE(value == doctest::Approx(-2.0 * face.area()).epsilon(1e-13));
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
      std::array<double, FluxFaceNodes> adjustedScale{};
      CurvedCell::setFace(cell, face, side, true, local.faceRotation[side], adjustedScale);
      CurvedCell::setFace(cell, face, side, false, local.faceRotation[side], scale);
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
        REQUIRE(adjustedScale[node] == doctest::Approx(-normal.norm()).epsilon(1e-8));

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

    addVolumeTerm<typename Kernels<Enabled>::Volume>(dofs.data(), update.data(), local);

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

/// A material that varies smoothly in space around `base`, every parameter with a phase of its own.
template <typename MaterialT>
auto smoothMaterial(const MaterialT& base) {
  return [base](const VectorT& x) {
    MaterialT material = base;
    double phase = 0.0;
    for (const auto& [name, member] : MaterialT::ParameterMap) {
      material.*member =
          base.*member * (1.0 + 0.2 * std::sin(0.9 * x(0) - 0.7 * x(1) + 1.1 * x(2) + phase));
      phase += 1.3;
    }
    return material;
  };
}

/// The data the corrector of a curved cell reads, as the setup gives it to every curved cell
/// (initializer::initializeCellLocalMatrices): weight-adjusted, with the material `materialAt`
/// gives at a position, and every face regular. `FaceT` is the map of a face of `cell`.
template <typename FaceT, typename MaterialT, typename CellT, typename FieldT>
void setupSplitCell(const CellT& cell, const FieldT& materialAt, LocalIntegrationData& local) {
  std::array<MaterialT, MaterialSampleCount> samples{};
  for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
    samples[point] = materialAt(cell.refToSpace(CurvedCell::operatorPoint(point)));
    const auto coefficients = seissol::model::getStarCoefficients(samples[point]);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
      local.materialCoefficients[i][point] = coefficients[i];
    }
  }
  CurvedCell::setGradients(cell, local.referenceGradients);
  CurvedCell::setWeightScale(cell, local.weightScale);
  CurvedCell::setSplitCoefficients(cell, samples.data(), local.splitCoefficients);

  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    const FaceT face(cell, seissol::geometry::ReferenceFaceMap(side));
    std::array<double, FluxFaceNodes> scale{};
    CurvedCell::setFace(cell, face, side, true, local.faceRotation[side], scale);
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
                                     minus,
                                     true);
      std::array<double, FluxCoefficientCount> scalars{};
      initializer::correctorScalars(FaceType::Regular, plus, minus, scalars);
      for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
        local.fluxCoefficients[side][c][node] = static_cast<real>(scalars[c]);
      }
      CurvedCell::setOwnTraction(material, scale[node], node, local.ownTractionCoefficients[side]);
    }
  }
}

/// A cell whose faces are all regular, each written in another parametrisation than the one before.
inline auto regularCell() -> CellLocalInformation {
  CellLocalInformation information{};
  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    information.faceTypes[side] = FaceType::Regular;
    information.traceSides[side] = side % 2;
  }
  return information;
}

/// The volume term and the own traction of a cell in the split form, from the time integral
/// `dofs`, into `buffer`, as kernels::CurvedSplit adds them; the traces it reads are written into
/// `traces` on the way.
template <bool Enabled>
void addSplit(const real* dofs,
              const LocalIntegrationData& local,
              const CellLocalInformation& information,
              std::array<FaceTraceT, Cell::NumFaces>& traces,
              real* buffer) {
  std::array<real*, Cell::NumFaces> pointers{};
  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    writeFaceTrace(dofs, side, information.traceSides[side], traces[side].data());
    pointers[side] = traces[side].data();
  }
  typename Kernels<Enabled>::Split split{};
  split.bindGlobals(seissol::Pool::host());
  split.apply(dofs, local, information, pointers, buffer);
}

/// Fills a time integral with random values, simulation by simulation.
inline void randomField(real* dofs, std::mt19937& rng) {
  std::normal_distribution<double> gauss(0.0, 1.0);
  for (std::size_t i = 0; i < tensor::I::size(); ++i) {
    dofs[i] = static_cast<real>(gauss(rng));
  }
}

/// Requires `got` to be `expected`, quantity by quantity against the largest value of the quantity
/// in `scaleOf`; a quantity that is zero there against the smallest of the others, since the
/// quantities of a medium differ by the orders of magnitude of its material.
inline void requireColumnsClose(const real* got,
                                const real* expected,
                                const real* scaleOf,
                                double tolerance) {
  constexpr std::size_t Columns = tensor::star::Shape[0][1];
  constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
  auto viewGot = init::Q::view::create(const_cast<real*>(got));
  auto viewExpected = init::Q::view::create(const_cast<real*>(expected));
  auto viewScale = init::Q::view::create(const_cast<real*>(scaleOf));
  for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
    auto slicedGot = multisim::simtensor(viewGot, simulation);
    auto slicedExpected = multisim::simtensor(viewExpected, simulation);
    auto slicedScale = multisim::simtensor(viewScale, simulation);
    std::array<double, Columns> scales{};
    double smallest = 0.0;
    for (std::size_t column = 0; column < Columns; ++column) {
      for (std::size_t row = 0; row < Basis; ++row) {
        scales[column] =
            std::max(scales[column], std::abs(static_cast<double>(slicedScale(row, column))));
      }
      if (scales[column] > 0.0 && (smallest == 0.0 || scales[column] < smallest)) {
        smallest = scales[column];
      }
    }
    REQUIRE(smallest > 0.0);
    for (std::size_t column = 0; column < Columns; ++column) {
      const double scale = scales[column] > 0.0 ? scales[column] : smallest;
      for (std::size_t row = 0; row < Basis; ++row) {
        REQUIRE(std::abs(static_cast<double>(slicedGot(row, column)) -
                         static_cast<double>(slicedExpected(row, column))) <= tolerance * scale);
      }
    }
  }
}

/// On a straight-sided cell, the split form is the strong form: its volume term and the own
/// traction of its faces, through its weight, are the volume term of the predictor, here with a
/// material that varies inside the cell.
template <bool Enabled>
void splitIsStrongOnStraightCell() {
  if constexpr (Enabled && !NodalSource) {
    using Material = seissol::model::MaterialT;
    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261029);
    const auto cell = seissol::geometry::AffineTransform(straightVertices());
    const auto materialAt = smoothMaterial(coefficients::configuredMaterial<Material>(rng));
    LocalIntegrationData local{};
    setupSplitCell<seissol::geometry::AffineFaceTransform, Material>(cell, materialAt, local);

    alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
    randomField(dofs.data(), rng);

    alignas(Alignment) std::array<real, tensor::Q::size()> expected{};
    addVolumeTerm<typename Kernels<Enabled>::Volume>(dofs.data(), expected.data(), local);

    alignas(Alignment) std::array<FaceTraceT, Cell::NumFaces> traces{};
    alignas(Alignment) std::array<real, tensor::Q::size()> buffer{};
    addSplit<Enabled>(dofs.data(), local, regularCell(), traces, buffer.data());
    alignas(Alignment) std::array<real, tensor::Q::size()> got{};
    typename Kernels<Enabled>::WeightAdjust weight{};
    weight.bindGlobals(seissol::Pool::host());
    weight.Q = got.data();
    weight.fluxBuffer = buffer.data();
    kernels::bindWeightOperands(weight, local);
    weight.execute();

    constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
    requireColumnsClose(got.data(), expected.data(), expected.data(), Tolerance);
  }
}

/// The weights of the energy of the medium quantity by quantity, which make every operator without
/// material symmetric: rate p of quantity q gives w_p B_d(q, p) = w_q B_d(p, q). Twice for a shear
/// stress, whose rate is a component of the tensor of the strain rate.
template <typename MaterialT>
auto energyPairing(const MaterialT& material) -> std::vector<double> {
  constexpr std::size_t NQ = tensor::star::Shape[0][0];
  std::array<std::array<real, tensor::star::size(0)>, Cell::Dim> data{};
  for (std::size_t d = 0; d < Cell::Dim; ++d) {
    auto view = init::star::view<0>::create(data[d].data());
    seissol::model::getTransposedUnweightedCoefficientMatrix(material, d, view);
  }
  const auto entry = [&](std::size_t d, std::size_t row, std::size_t column) -> double {
    auto view = init::star::view<0>::create(data[d].data());
    return view.isInRange(row, column) ? static_cast<double>(view(row, column)) : 0.0;
  };
  std::vector<double> pairing(NQ, 0.0);
  pairing[0] = 1.0;
  for (bool changed = true; changed;) {
    changed = false;
    for (std::size_t d = 0; d < Cell::Dim; ++d) {
      for (std::size_t q = 0; q < NQ; ++q) {
        for (std::size_t p = 0; p < NQ; ++p) {
          const double forward = entry(d, q, p);
          const double backward = entry(d, p, q);
          if (forward != 0.0 && backward != 0.0 && pairing[q] != 0.0 && pairing[p] == 0.0) {
            pairing[p] = pairing[q] * backward / forward;
            changed = true;
          }
        }
      }
    }
  }
  for (std::size_t q = 0; q < NQ; ++q) {
    REQUIRE(pairing[q] > 0.0);
  }
  return pairing;
}

/// On a curved cell, the two parts of the split volume term cancel in the energy at every point:
/// the field paired with the moments of the term with the weights of the energy gives nothing.
template <bool Enabled>
void splitVolumeDoesNoWork() {
  if constexpr (Enabled && !NodalSource && multisim::NumSimulations == 1) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t NQ = tensor::star::Shape[0][0];
    constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261030);
    const auto cell = curvedCell();
    const auto base = coefficients::configuredMaterial<Material>(rng);
    const auto materialAt = smoothMaterial(base);
    LocalIntegrationData local{};
    setupSplitCell<seissol::geometry::IsoparametricFaceTransform, Material>(
        cell, materialAt, local);
    const auto pairing = energyPairing(base);

    const auto mass = init::M3::view::create(const_cast<real*>(init::M3::Values));
    for (std::size_t sample = 0; sample < 3; ++sample) {
      alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
      randomField(dofs.data(), rng);
      alignas(Alignment) std::array<real, tensor::Q::size()> buffer{};
      typename Kernels<Enabled>::SplitVolume volume{};
      volume.bindGlobals(seissol::Pool::host());
      volume.I = dofs.data();
      volume.fluxBuffer = buffer.data();
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
        volume.referenceGradients(dim) = local.referenceGradients[dim];
      }
      for (std::size_t a = 0; a < StarCoefficientCount; ++a) {
        volume.splitCoefficients(a) = local.splitCoefficients[a];
      }
      volume.execute();

      // the buffer is the inverse mass matrix applied to the moments
      auto viewI = init::I::view::create(dofs.data());
      auto viewBuffer = init::Q::view::create(buffer.data());
      double work = 0.0;
      double size = 0.0;
      for (std::size_t k = 0; k < Basis; ++k) {
        for (std::size_t q = 0; q < NQ; ++q) {
          double moment = 0.0;
          for (std::size_t l = 0; l < Basis; ++l) {
            if (mass.isInRange(k, l)) {
              moment += static_cast<double>(mass(k, l)) * static_cast<double>(viewBuffer(l, q));
            }
          }
          const double term = static_cast<double>(viewI(k, q)) * pairing[q] * moment;
          work += term;
          size += std::abs(term);
        }
      }
      REQUIRE(size > 0.0);
      constexpr double Tolerance = std::is_same_v<real, double> ? 1e-12 : 1e-5;
      REQUIRE(std::abs(work) <= Tolerance * size);
    }
  }
}

/// A constant state stays constant on a curved cell whose material varies: the weak form of the
/// momentum rows against the own traction the faces take back, and the correctors of a state that
/// is continuous across every face.
template <bool Enabled>
void constantState() {
  if constexpr (Enabled && !NodalSource) {
    using Material = seissol::model::MaterialT;
    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261009);
    std::normal_distribution<double> gauss(0.0, 1.0);
    const auto cell = curvedCell();
    const auto materialAt = smoothMaterial(coefficients::configuredMaterial<Material>(rng));
    LocalIntegrationData local{};
    setupSplitCell<seissol::geometry::IsoparametricFaceTransform, Material>(
        cell, materialAt, local);

    alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
    auto viewI = init::I::view::create(dofs.data());
    for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
      auto slicedI = multisim::simtensor(viewI, simulation);
      for (std::size_t quantity = 0; quantity < tensor::star::Shape[0][0]; ++quantity) {
        const double value = gauss(rng);
        for (std::size_t mode = 0; mode < tensor::Q::Shape[multisim::BasisFunctionDimension];
             ++mode) {
          const double coefficient = basisFunction::constantCoefficient(mode);
          if (coefficient != 0.0) {
            slicedI(mode, quantity) = static_cast<real>(value * coefficient);
          }
        }
      }
    }

    // the size of what cancels: the weak form of the momentum rows alone
    const auto information = regularCell();
    alignas(Alignment) std::array<FaceTraceT, Cell::NumFaces> traces{};
    alignas(Alignment) std::array<real, tensor::Q::size()> volumeOnly{};
    {
      auto onlyVolume = information;
      for (auto& faceType : onlyVolume.faceTypes) {
        // a fault face takes no own traction from the cell
        faceType = FaceType::DynamicRupture;
      }
      addSplit<Enabled>(dofs.data(), local, onlyVolume, traces, volumeOnly.data());
    }

    alignas(Alignment) std::array<real, tensor::Q::size()> buffer{};
    addSplit<Enabled>(dofs.data(), local, information, traces, buffer.data());
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      // the neighbour's trace is the cell's own where the state is continuous
      applyFaceCorrector(traces[side].data(),
                         traces[side].data(),
                         local,
                         side,
                         information.traceSides[side],
                         buffer.data());
    }

    alignas(Alignment) const std::array<real, tensor::Q::size()> zero{};
    constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
    requireColumnsClose(buffer.data(), zero.data(), volumeOnly.data(), Tolerance);
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

TEST_CASE("On a straight-sided cell the split form is the strong form") {
  curvedcell::splitIsStrongOnStraightCell<Curvilinear>();
}

TEST_CASE("The split volume term of a curved cell does no work") {
  curvedcell::splitVolumeDoesNoWork<Curvilinear>();
}

TEST_CASE("A curved cell keeps a constant state constant") {
  curvedcell::constantState<Curvilinear>();
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_CURVEDCELL_T_H_
