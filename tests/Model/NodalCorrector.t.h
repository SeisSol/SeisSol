// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_NODALCORRECTOR_T_H_
#define SEISSOL_TESTS_MODEL_NODALCORRECTOR_T_H_

// Where the material varies inside a cell, the corrector is the strong form: the volume term
// differentiates the field, and the local flux of every face takes off the normal flux of the
// cell's own trace, a fault face included (StrongCorrector). Both checks below run the kernels a
// cell runs, with the flux scalars the setup computes:
//
// - with one material for the cell, volume term and local fluxes together are the weak form, as
//   the matrices give it contracted by hand. A sign, a scale or a face left out of the
//   subtraction shows up there;
// - a constant state stays constant where the material varies inside the cell and along its
//   faces. The weak form does not do that: it lifts the cell's own trace with the material inside
//   the cell, and the flux takes it off with the material at the face.

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
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Initializer/Model/CellFlux.h"
#include "Initializer/Model/CurvedCell.h"
#include "Initializer/Parameters/ModelParameters.h"
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

// The anelastic solver writes its corrector into the extended quantities and
// adds the relaxation in a kernel of its own; the comparison has no single
// kernel sequence to run there.
#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

namespace seissol::unit_test {

namespace nodalcorrector {

/// Makes the kernel names dependent on a template parameter, so that a build
/// without a flux at the nodes never instantiates a body that names operands
/// it does not have.
template <bool Enabled>
struct Kernels {
  using Volume = std::conditional_t<Enabled, kernel::volume, kernel::volume>;
  using LocalFlux = std::conditional_t<Enabled, kernel::localFlux, kernel::localFlux>;
  using NeighborFlux = std::conditional_t<Enabled, kernel::neighboringFlux, void>;
  using FluxSolver = std::conditional_t<Enabled, kernel::computeFluxSolverLocal, void>;
};

using VectorT = seissol::geometry::CellTransform::VectorEigenT;
using FaceVectorT = seissol::geometry::FaceTransform::FaceVectorT;

/// a well-shaped cell that is not aligned with anything
inline auto testCell() -> seissol::geometry::AffineTransform {
  return seissol::geometry::AffineTransform(
      std::array<VectorT, Cell::NumVertices>{VectorT(0.1, 0.2, -0.3),
                                             VectorT(2.1, 0.0, 0.4),
                                             VectorT(-0.2, 1.7, 0.1),
                                             VectorT(0.3, 0.1, 2.2)});
}

/// the face node `node` of `side`, in the reference coordinates of the cell
inline auto faceNode(std::size_t side, std::size_t node) -> VectorT {
  const auto nodes =
      nodal::init::nodes2D::view::create(const_cast<real*>(nodal::init::nodes2D::Values));
  return seissol::geometry::ReferenceFaceMap(side).faceToCell(
      FaceVectorT(nodes(node, 0), nodes(node, 1)));
}

/// the sample point `point` of the material, in the reference coordinates of the cell
inline auto samplePoint(std::size_t point) -> VectorT {
  const auto nodes =
      init::materialNodes::view::create(const_cast<real*>(init::materialNodes::Values));
  VectorT result = VectorT::Zero();
  for (std::size_t d = 0; d < Cell::Dim; ++d) {
    if (nodes.isInRange(point, d)) {
      result(d) = nodes(point, d);
    }
  }
  return result;
}

/// `base`, with every parameter it binds varying smoothly in space by up to a fifth
template <typename MaterialT>
auto varying(const MaterialT& base, const VectorT& x) -> MaterialT {
  MaterialT material = base;
  double phase = 0.0;
  for (const auto& [name, member] : MaterialT::ParameterMap) {
    material.*member =
        base.*member * (1.0 + 0.2 * std::sin(0.9 * x(0) - 0.7 * x(1) + 1.1 * x(2) + phase));
    phase += 1.3;
  }
  return material;
}

/// The cell data the volume and the local flux read, for a material given as a function of the
/// position; `faceTypes` decides what the local flux of a face is.
template <typename MaterialT, typename FieldT>
void setupCell(const seissol::geometry::AffineTransform& cell,
               const FieldT& field,
               const std::array<FaceType, Cell::NumFaces>& faceTypes,
               LocalIntegrationData& local,
               NeighboringIntegrationData& neighboring,
               std::array<std::array<real, tensor::AplusT::size()>, Cell::NumFaces>& weakAplus) {
  const auto gradients = cell.refToSpaceJacobianInverse(VectorT(Cell::ReferenceBarycenter.data()));
  initializer::CurvedCell::setConstantGradients(gradients, local.referenceGradients);
  for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
    const auto material = field(cell.refToSpace(samplePoint(point)));
    const auto coefficients = seissol::model::getStarCoefficients(material);
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
      local.materialCoefficients[i][point] = coefficients[i];
    }
    if constexpr (NodalSource) {
      const auto source = seissol::model::getSourceCoefficients(material);
      for (std::size_t i = 0; i < source.size(); ++i) {
        local.sourceCoefficients[i][point] = source[i];
      }
    }
  }

  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    const seissol::geometry::AffineFaceTransform face(cell,
                                                      seissol::geometry::ReferenceFaceMap(side));
    const auto basis = face.faceAlignedBasis();
    alignas(Alignment) std::array<real, tensor::T::size()> matT{};
    alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
    auto viewT = init::T::view::create(matT.data());
    auto viewTinv = init::Tinv::view::create(matTinv.data());
    seissol::model::getFaceRotationMatrix(
        basis[0].normalized(), basis[1].normalized(), basis[2].normalized(), viewT, viewTinv);
    initializer::CurvedCell::setConstantRotation(matT.data(), local.faceRotation[side]);

    // |S| / |J| with the sign of a subtracted flux, as the setup takes it
    const double fluxScale = -2.0 * face.area() / std::abs(cell.determinant());

    for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
      const auto material = field(cell.refToSpace(faceNode(side, node)));
      std::array<double, FluxCoefficientCount> plus{};
      std::array<double, FluxCoefficientCount> minus{};
      seissol::initializer::fluxScalarsOfNode(material,
                                              material,
                                              faceTypes[side],
                                              initializer::parameters::NumericalFlux::Godunov,
                                              fluxScale,
                                              plus,
                                              minus);
      for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
        local.fluxCoefficients[side][c][node] = static_cast<real>(plus[c]);
        neighboring.fluxCoefficients[side][c][node] = static_cast<real>(minus[c]);
      }
    }

    // the weak form's local operator of the face, from the material at its barycentre: what the
    // matrix form of a build without a varying material applies. A fault face has none.
    weakAplus[side].fill(0);
    if (faceTypes[side] != FaceType::DynamicRupture) {
      const auto material = field(face.center());
      alignas(Alignment) std::array<real, tensor::QgodLocal::size()> godLocal{};
      alignas(Alignment) std::array<real, tensor::QgodNeighbor::size()> godNeighbor{};
      auto viewGodLocal = init::QgodLocal::view::create(godLocal.data());
      auto viewGodNeighbor = init::QgodNeighbor::view::create(godNeighbor.data());
      seissol::model::getTransposedGodunovState(
          material, material, faceTypes[side], viewGodLocal, viewGodNeighbor);
      alignas(Alignment) std::array<real, tensor::star::size(0)> star{};
      auto viewStar = init::star::view<0>::create(star.data());
      seissol::model::getTransposedCoefficientMatrix(material, 0, viewStar);
      alignas(Alignment) std::array<real, tensor::QcorrLocal::size()> qcorr{};

      kernel::computeFluxSolverLocal solver{};
      solver.fluxScale = fluxScale;
      solver.AplusT = weakAplus[side].data();
      solver.QgodLocal = godLocal.data();
      solver.QcorrLocal = qcorr.data();
      solver.T = matT.data();
      solver.Tinv = matTinv.data();
      solver.star(0) = star.data();
      solver.execute();
    }
  }
}

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

/// a global matrix as the mathematics has it: a build that bundles simulations stores the
/// matrices from the matrix files the other way round
template <typename ViewT>
auto mathMatrix(const ViewT& view, std::size_t rows, std::size_t columns) -> Eigen::MatrixXd {
  const Eigen::MatrixXd stored = denseOf(view, rows, columns);
  return multisim::NumSimulations > 1 ? Eigen::MatrixXd(stored.transpose()) : stored;
}

template <bool Enabled>
void weakEquivalence() {
  if constexpr (Enabled && StrongCorrector) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t NQ = tensor::star::Shape[0][0];
    constexpr std::size_t Columns = tensor::star::Shape[0][1];
    constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261006);
    std::normal_distribution<double> gauss(0.0, 1.0);
    const auto cell = testCell();

    // every face takes a different path through the subtraction: a fault face has nothing but it
    const std::array<FaceType, Cell::NumFaces> faceTypes{
        FaceType::Regular, FaceType::FreeSurface, FaceType::DynamicRupture, FaceType::Outflow};

    for (std::size_t sample = 0; sample < 3; ++sample) {
      const auto material = coefficients::configuredMaterial<Material>(rng);
      const auto field = [&](const VectorT& /*x*/) { return material; };

      LocalIntegrationData local{};
      NeighboringIntegrationData neighboring{};
      std::array<std::array<real, tensor::AplusT::size()>, Cell::NumFaces> weakAplus{};
      setupCell<Material>(cell, field, faceTypes, local, neighboring, weakAplus);

      alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
      for (auto& value : dofs) {
        value = static_cast<real>(gauss(rng));
      }
      alignas(Alignment) std::array<real, tensor::Q::size()> update{};

      typename Kernels<Enabled>::Volume volume{};
      volume.bindGlobals(seissol::Pool::host());
      volume.I = dofs.data();
      volume.Q = update.data();
      kernels::bindStarOperands(volume, local);
      kernels::bindSourceOperands(volume, local);
      volume.execute();

      typename Kernels<Enabled>::LocalFlux flux{};
      flux.bindGlobals(seissol::Pool::host());
      flux.I = dofs.data();
      flux.Q = update.data();
      for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
        // A fault that takes the subtraction of the strong form on itself (FaultSubtractsOwnTrace)
        // leaves its face out of the local flux, and subtracts the same at its own points; the
        // identity is the one of the strong form, so the local flux stands in for it here.
        REQUIRE((appliesLocalFlux(faceTypes[side]) ||
                 (FaultSubtractsOwnTrace && faceTypes[side] == FaceType::DynamicRupture)));
        kernels::bindLocalFluxOperands(flux, local, side);
        flux.execute(side);
      }

      // the weak form, by hand out of the matrices it is made of
      std::array<Eigen::MatrixXd, Cell::Dim> directional;
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
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
        const auto inverse =
            cell.refToSpaceJacobianInverse(VectorT(Cell::ReferenceBarycenter.data()));
        const auto addDirection = [&](auto tag) {
          constexpr std::size_t Dim = decltype(tag)::value;
          Eigen::MatrixXd star = Eigen::MatrixXd::Zero(NQ, Columns);
          for (std::size_t component = 0; component < Cell::Dim; ++component) {
            star += inverse(Dim, component) * directional[component];
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

        const auto addFace = [&](auto tag) {
          constexpr std::size_t Side = decltype(tag)::value;
          const auto rDivM = mathMatrix(
              init::rDivM::view<Side>::create(const_cast<real*>(init::rDivM::Values[Side])),
              tensor::rDivM::Shape[Side][0],
              tensor::rDivM::Shape[Side][1]);
          const auto fMrT = mathMatrix(
              init::fMrT::view<Side>::create(const_cast<real*>(init::fMrT::Values[Side])),
              tensor::fMrT::Shape[Side][0],
              tensor::fMrT::Shape[Side][1]);
          const auto aplus =
              denseOf(init::AplusT::view::create(weakAplus[Side].data()), NQ, Columns);
          expected += rDivM * fMrT * field * aplus;
        };
        addFace(std::integral_constant<std::size_t, 0>{});
        addFace(std::integral_constant<std::size_t, 1>{});
        addFace(std::integral_constant<std::size_t, 2>{});
        addFace(std::integral_constant<std::size_t, 3>{});

        // the quantities differ by orders of magnitude -- a stress by the moduli, a velocity by
        // the inverse density -- and entries of one quantity cancel against each other, so each
        // is measured against the size of its own column
        constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
        for (std::size_t column = 0; column < Columns; ++column) {
          const double scale = std::max(1.0, expected.col(column).cwiseAbs().maxCoeff());
          for (std::size_t row = 0; row < Basis; ++row) {
            REQUIRE(slicedQ(row, column) ==
                    doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
          }
        }
      }
    }
  }
}

template <bool Enabled>
void constantState() {
  if constexpr (Enabled && StrongCorrector) {
    using Material = seissol::model::MaterialT;
    constexpr std::size_t Columns = tensor::star::Shape[0][1];
    constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];

    // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
    std::mt19937 rng(20261007);
    std::normal_distribution<double> gauss(0.0, 1.0);
    const auto cell = testCell();
    const std::array<FaceType, Cell::NumFaces> faceTypes{
        FaceType::Regular, FaceType::Regular, FaceType::Regular, FaceType::Regular};

    for (std::size_t sample = 0; sample < 3; ++sample) {
      const auto base = coefficients::configuredMaterial<Material>(rng);
      // one field for the cell and its neighbours, so that the material is continuous across the
      // faces and varies along them as well as inside the cell
      const auto field = [&](const VectorT& x) { return varying(base, x); };

      LocalIntegrationData local{};
      NeighboringIntegrationData neighboring{};
      std::array<std::array<real, tensor::AplusT::size()>, Cell::NumFaces> weakAplus{};
      setupCell<Material>(cell, field, faceTypes, local, neighboring, weakAplus);

      // a constant state: the first mode, which the basis keeps constant, and nothing else
      alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
      auto viewI = init::I::view::create(dofs.data());
      for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
        auto slicedI = multisim::simtensor(viewI, simulation);
        for (std::size_t quantity = 0; quantity < tensor::star::Shape[0][0]; ++quantity) {
          slicedI(0, quantity) = static_cast<real>(gauss(rng));
        }
      }

      // the size of what cancels: the local flux of the faces on its own
      alignas(Alignment) std::array<real, tensor::Q::size()> localOnly{};
      // and everything the corrector adds for a regular face
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
      kernels::bindSourceOperands(volume, local);
      volume.execute();

      flux.Q = update.data();
      typename Kernels<Enabled>::NeighborFlux neighbor{};
      neighbor.bindGlobals(seissol::Pool::host());
      neighbor.I = dofs.data();
      neighbor.Q = update.data();
      for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
        kernels::bindLocalFluxOperands(flux, local, side);
        flux.execute(side);
        // the neighbour's trace of a constant state is that state whichever side it is read from
        kernels::bindNeighborFluxOperands(neighbor, local, neighboring, side);
        neighbor.execute(0, side);
      }

      auto viewLocal = init::Q::view::create(localOnly.data());
      auto viewQ = init::Q::view::create(update.data());
      for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
        auto slicedLocal = multisim::simtensor(viewLocal, simulation);
        auto slicedQ = multisim::simtensor(viewQ, simulation);
        // a source term acts on a constant state as on any other, so there is nothing that has to
        // vanish where the solver has one
        if constexpr (!NodalSource) {
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
  }
}

} // namespace nodalcorrector

TEST_CASE("Nodal corrector is the weak form where the material does not vary") {
  nodalcorrector::weakEquivalence<NodalFlux>();
}

TEST_CASE("Nodal corrector keeps a constant state constant") {
  nodalcorrector::constantState<NodalFlux>();
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_NODALCORRECTOR_T_H_
