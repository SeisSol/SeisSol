// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_MODEL_FACECORRECTOR_T_H_
#define SEISSOL_TESTS_MODEL_FACECORRECTOR_T_H_

// The corrector of a face works from the traces the two cells sharing it write there (see
// _configureTraces in the generator). Checked here against what it stands for: the flux of the
// cell's own field and of the neighbour's, contracted by hand out of the matrices of the volume --
// rDivM(i) fMrT(i) I A+ for the cell's own side of its face i, and rDivM(i) fP(0) rT(j) I' A- for
// the neighbour, whose face j it is -- with the operators of the face the setup forms. The two
// cells write their traces in the parametrisation of one side of the face, so a cell is checked on
// either side. A face taken for another, a reflection applied twice or not at all, or a rotation
// left out shows up in that comparison and in nothing weaker.
//
// Where the flux is the scalars at the points of a face, the corrector applies one set of them to
// the jump of the two traces, the neighbour's less the cell's own; the part of either side is then
// checked with the other trace zero, with the scalars of that part.
//
// A face without a neighbour reads a trace of zeros as the neighbour's, against the operator the
// setup stores for the side outside the cell, which is zero where there is none
// (initializer::hasExteriorOperator); it is checked that way, and has to come back with the cell's
// own part alone. A condition that folds a map into the operator of the cell's own side has to
// come back with that operator, where the scalars of the face leave the folded part to a kernel of
// its own (appliesFoldApart).

#include <doctest.h>

// the material builders live with the decomposition they were written for
#include "CoefficientStructure.t.h" // IWYU pragma: keep
#include "Equations/Datastructures.h"
#include "Equations/Setup.h"
#include "FaceTraces.h"
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BoundarySetup.h"
#include "Initializer/Model/CellFlux.h"
#include "Initializer/Model/CurvedCell.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
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

// The anelastic solver lifts into the extended quantities, and its decomposition is a separate
// table; the comparison here is for the field.
#ifndef SEISSOL_KERNELS_LINEARCKANELASTIC

namespace seissol::unit_test {

namespace facecorrectortest {

/// The boundary conditions whose operator the flux itself carries, besides a face with a
/// neighbour. Gravity and Dirichlet fold the map of their ghost state into the operator of the
/// cell's own side, checked below (compareFoldedBoundary); what is left of them, and the analytical
/// condition, is a ghost state at the nodes of the face that the neighbour matrix lifts, which is
/// not a question about the traces.
constexpr FaceType LinearFaceTypes[] = {
    FaceType::Regular, FaceType::FreeSurface, FaceType::Outflow};

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

/// The member `Face` of a family of global matrices with one per face, as the mathematics has it.
#define SEISSOL_FACE_MATRIX(NAME, FACE)                                                            \
  mathMatrix(init::NAME::view<FACE>::create(const_cast<real*>(init::NAME::Values[FACE])),          \
             tensor::NAME::Shape[FACE][0],                                                         \
             tensor::NAME::Shape[FACE][1])

/// rDivM(i) and fMrT(i) of every face i, and fP(0) M2^-1, which takes the trace a cell writes on
/// one of its faces, M2 rT(j), to the one its neighbour reads there, fP(0) rT(j).
struct VolumeForm {
  std::array<Eigen::MatrixXd, Cell::NumFaces> rDivM;
  std::array<Eigen::MatrixXd, Cell::NumFaces> fMrT;
  Eigen::MatrixXd neighborRead;

  VolumeForm() {
    const auto setFace = [&](auto tag) {
      constexpr std::size_t Face = decltype(tag)::value;
      rDivM[Face] = SEISSOL_FACE_MATRIX(rDivM, Face);
      fMrT[Face] = SEISSOL_FACE_MATRIX(fMrT, Face);
    };
    setFace(std::integral_constant<std::size_t, 0>{});
    setFace(std::integral_constant<std::size_t, 1>{});
    setFace(std::integral_constant<std::size_t, 2>{});
    setFace(std::integral_constant<std::size_t, 3>{});
    const auto reflected = SEISSOL_FACE_MATRIX(fP, 0);
    const auto inverseMass =
        mathMatrix(init::M2inv::view::create(const_cast<real*>(init::M2inv::Values)),
                   tensor::M2inv::Shape[0],
                   tensor::M2inv::Shape[1]);
    neighborRead = reflected * inverseMass;
  }
};

#undef SEISSOL_FACE_MATRIX

/// what a field holds for one fused simulation, as a matrix of modes by quantities
template <typename ViewT>
auto fieldOf(ViewT view, std::size_t simulation, std::size_t quantities) -> Eigen::MatrixXd {
  constexpr std::size_t Basis = tensor::Q::Shape[multisim::BasisFunctionDimension];
  auto sliced = multisim::simtensor(view, simulation);
  Eigen::MatrixXd field = Eigen::MatrixXd::Zero(Basis, quantities);
  for (std::size_t row = 0; row < Basis; ++row) {
    for (std::size_t column = 0; column < quantities; ++column) {
      field(row, column) = sliced(row, column);
    }
  }
  return field;
}

/// The scalars of the flux of a face, read off its operator in face coordinates, `god` times the
/// coefficient matrix of the normal, scaled as the setup scales them. A material that does not
/// vary along the face has the same at every point.
template <typename GodViewT, typename StarViewT>
void setScalars(const GodViewT& god,
                const StarViewT& star,
                double scale,
                std::size_t face,
                LocalIntegrationData& local) {
  constexpr std::size_t NQ = tensor::QgodLocal::Shape[0];
  for (std::size_t c = 0; c < FluxCoefficientCount; ++c) {
    const auto& source = generated::FluxCoefficientSources[c];
    double value = 0.0;
    for (std::size_t inner = 0; inner < NQ; ++inner) {
      const double g = god.isInRange(source.row, inner) ? god(source.row, inner) : 0.0;
      const double a = star.isInRange(inner, source.column) ? star(inner, source.column) : 0.0;
      value += g * a;
    }
    for (std::size_t node = 0; node < FluxFaceNodes; ++node) {
      local.fluxCoefficients[face][c][node] = static_cast<real>(scale * value);
    }
  }
}

inline void compareAgainstVolumeForm(FaceType faceType) {
  using Material = seissol::model::MaterialT;
  constexpr std::size_t NQ = tensor::QgodLocal::Shape[0];
  constexpr std::size_t Written = tensor::AplusT::Shape[1];

  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937 rng(20261020);
  std::uniform_real_distribution<double> positive(0.4, 2.5);
  std::normal_distribution<double> gauss(0.0, 1.0);
  const VolumeForm volume;

  for (std::size_t sample = 0; sample < 4; ++sample) {
    const auto local = coefficients::configuredMaterial<Material>(rng);
    const auto neighbor = coefficients::configuredMaterial<Material>(rng);

    // a general orthonormal frame: the first tangent of a face is one of its edges, so nothing may
    // lean on a particular choice
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
    alignas(Alignment) std::array<real, tensor::T::size()> matT{};
    alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
    auto viewT = init::T::view::create(matT.data());
    auto viewTinv = init::Tinv::view::create(matTinv.data());
    seissol::model::getFaceRotationMatrix(CoordinateT{frame(0, 0), frame(1, 0), frame(2, 0)},
                                          CoordinateT{frame(0, 1), frame(1, 1), frame(2, 1)},
                                          CoordinateT{frame(0, 2), frame(1, 2), frame(2, 2)},
                                          viewT,
                                          viewTinv);

    const double fluxScale = -0.37 * positive(rng);

    // the operator of either side in face coordinates, and the star the setup folds in
    alignas(Alignment) std::array<real, tensor::QgodLocal::size()> godLocal{};
    alignas(Alignment) std::array<real, tensor::QgodNeighbor::size()> godNeighbor{};
    auto viewGodLocal = init::QgodLocal::view::create(godLocal.data());
    auto viewGodNeighbor = init::QgodNeighbor::view::create(godNeighbor.data());
    seissol::model::getTransposedGodunovState(
        local, neighbor, faceType, viewGodLocal, viewGodNeighbor);
    alignas(Alignment) std::array<real, tensor::star::size(0)> star{};
    auto viewStar = init::star::view<0>::create(star.data());
    seissol::model::getTransposedCoefficientMatrix(local, 0, viewStar);

    // the operator of either side with the rotation folded in, from the kernels the setup runs
    alignas(Alignment) std::array<real, tensor::AplusT::size()> aplus{};
    alignas(Alignment) std::array<real, tensor::AminusT::size()> aminus{};
    alignas(Alignment) std::array<real, tensor::QcorrLocal::size()> qcorrLocal{};
    alignas(Alignment) std::array<real, tensor::QcorrNeighbor::size()> qcorrNeighbor{};
    kernel::computeFluxSolverLocal localSolver{};
    localSolver.fluxScale = fluxScale;
    localSolver.AplusT = aplus.data();
    localSolver.QgodLocal = godLocal.data();
    localSolver.QcorrLocal = qcorrLocal.data();
    localSolver.T = matT.data();
    localSolver.Tinv = matTinv.data();
    localSolver.star(0) = star.data();
    localSolver.execute();
    kernel::computeFluxSolverNeighbor neighborSolver{};
    neighborSolver.fluxScale = fluxScale;
    neighborSolver.AminusT = aminus.data();
    neighborSolver.QgodNeighbor = godNeighbor.data();
    neighborSolver.QcorrNeighbor = qcorrNeighbor.data();
    neighborSolver.T = matT.data();
    neighborSolver.Tinv = matTinv.data();
    neighborSolver.star(0) = star.data();
    if (initializer::hasExteriorOperator(faceType)) {
      neighborSolver.execute();
    }
    const auto aplusDense = denseOf(init::AplusT::view::create(aplus.data()), NQ, Written);
    const auto aminusDense = denseOf(init::AminusT::view::create(aminus.data()), NQ, Written);

    // the fields of the cell and of its neighbour
    alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
    alignas(Alignment) std::array<real, tensor::I::size()> neighborDofs{};
    for (auto& value : dofs) {
      value = static_cast<real>(gauss(rng));
    }
    for (auto& value : neighborDofs) {
      value = static_cast<real>(gauss(rng));
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      LocalIntegrationData localData{};
      NeighboringIntegrationData neighboringData{};
      if constexpr (NodalFlux) {
        initializer::CurvedCell::setConstantRotation(matT.data(), localData.faceRotation[face]);
      } else {
        std::copy(aplus.begin(), aplus.end(), localData.nApNm1[face]);
        std::copy(aminus.begin(), aminus.end(), neighboringData.nAmNm1[face]);
      }

      // only a regular face has a neighbour; the others keep the trace of zeros as its trace
      const bool hasNeighbor = faceType == FaceType::Regular;
      for (std::size_t side = 0; side < 2; ++side) {
        for (std::size_t neighborFace = 0; neighborFace < (hasNeighbor ? Cell::NumFaces : 1);
             ++neighborFace) {
          // the traces the two cells write on the face, both in the parametrisation of `side`
          alignas(Alignment) FaceTraceT own{};
          alignas(Alignment) FaceTraceT across{};
          alignas(Alignment) const FaceTraceT zero{};
          writeFaceTrace(dofs.data(), face, side, own.data());
          if (hasNeighbor) {
            writeFaceTrace(neighborDofs.data(), neighborFace, 1 - side, across.data());
          }

          alignas(Alignment) std::array<real, tensor::Q::size()> update{};
          if constexpr (NodalFlux) {
            // the jump of the traces is the neighbour's less the cell's own, so the cell's own part
            // takes its scalars negated
            setScalars(viewGodLocal, viewStar, -fluxScale, face, localData);
            applyFaceCorrector(
                own.data(), zero.data(), localData, neighboringData, face, side, update.data());
            if (hasNeighbor) {
              setScalars(viewGodNeighbor, viewStar, fluxScale, face, localData);
              applyFaceCorrector(zero.data(),
                                 across.data(),
                                 localData,
                                 neighboringData,
                                 face,
                                 side,
                                 update.data());
            }
          } else {
            applyFaceCorrector(
                own.data(), across.data(), localData, neighboringData, face, side, update.data());
          }

          auto viewI = init::I::view::create(dofs.data());
          auto viewNeighbor = init::I::view::create(neighborDofs.data());
          auto viewQ = init::Q::view::create(update.data());
          // the operator is shared by the simulations a build bundles, the field is not, so every
          // one of them has to come back right
          for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
            const auto field = fieldOf(viewI, simulation, NQ);
            const auto fieldAcross = fieldOf(viewNeighbor, simulation, NQ);
            const auto got = fieldOf(viewQ, simulation, Written);
            Eigen::MatrixXd expected = volume.rDivM[face] * volume.fMrT[face] * field * aplusDense;
            if (hasNeighbor) {
              expected += volume.rDivM[face] * volume.neighborRead * volume.fMrT[neighborFace] *
                          fieldAcross * aminusDense;
            }

            // Entries of an output quantity cancel against each other, so what one of them may
            // be off by follows the size of that quantity in the update, not its own. The
            // quantities themselves differ by orders of magnitude -- a stress by the moduli, a
            // velocity by the inverse density -- so each is measured against its own column.
            constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
            for (std::size_t column = 0; column < Written; ++column) {
              const double scale = std::max(1.0, expected.col(column).cwiseAbs().maxCoeff());
              for (Eigen::Index row = 0; row < expected.rows(); ++row) {
                REQUIRE(got(row, column) ==
                        doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
              }
            }
          }
        }
      }
    }
  }
}

/// A Dirichlet face takes a ghost state that is a map of the cell's own state, constant over the
/// face, which the setup folds into the operator of the cell's own side (foldDirichlet). Where the
/// corrector reads the scalars of the face for that side instead, the operator keeps the folded
/// part alone, which the local integral lifts from the time integral of the cell
/// (appliesFoldApart, foldedBoundaryFlux). Either way, what the cell's own side adds is the flux of
/// the operator with the map folded in.
inline void compareFoldedBoundary() {
  using Material = seissol::model::MaterialT;
  constexpr std::size_t NQ = tensor::QgodLocal::Shape[0];
  constexpr std::size_t Written = tensor::AplusT::Shape[1];

  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937 rng(20261021);
  std::uniform_real_distribution<double> positive(0.4, 2.5);
  std::normal_distribution<double> gauss(0.0, 1.0);
  const VolumeForm volume;

  const auto material = coefficients::configuredMaterial<Material>(rng);

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
  alignas(Alignment) std::array<real, tensor::T::size()> matT{};
  alignas(Alignment) std::array<real, tensor::Tinv::size()> matTinv{};
  auto viewT = init::T::view::create(matT.data());
  auto viewTinv = init::Tinv::view::create(matTinv.data());
  seissol::model::getFaceRotationMatrix(CoordinateT{frame(0, 0), frame(1, 0), frame(2, 0)},
                                        CoordinateT{frame(0, 1), frame(1, 1), frame(2, 1)},
                                        CoordinateT{frame(0, 2), frame(1, 2), frame(2, 2)},
                                        viewT,
                                        viewTinv);
  // the ghost state is given in the coordinates of the face already
  alignas(Alignment) std::array<real, tensor::Tinv::size()> identityTinv{};
  {
    const auto identity = init::identityT::view::create(init::identityT::Values);
    init::Tinv::view::create(identityTinv.data()).forall([&](const auto* entry, real& value) {
      value = identity.isInRange(entry[0], entry[1]) ? identity(entry[0], entry[1]) : 0;
    });
  }

  const double fluxScale = -0.37 * positive(rng);

  // the face has the cell's own material on both sides, as a boundary face has
  alignas(Alignment) std::array<real, tensor::QgodLocal::size()> godLocal{};
  alignas(Alignment) std::array<real, tensor::QgodNeighbor::size()> godNeighbor{};
  auto viewGodLocal = init::QgodLocal::view::create(godLocal.data());
  auto viewGodNeighbor = init::QgodNeighbor::view::create(godNeighbor.data());
  seissol::model::getTransposedGodunovState(
      material, material, FaceType::Dirichlet, viewGodLocal, viewGodNeighbor);
  alignas(Alignment) std::array<real, tensor::star::size(0)> star{};
  auto viewStar = init::star::view<0>::create(star.data());
  seissol::model::getTransposedCoefficientMatrix(material, 0, viewStar);

  alignas(Alignment) std::array<real, tensor::AplusT::size()> aplus{};
  alignas(Alignment) std::array<real, tensor::AminusT::size()> aminus{};
  alignas(Alignment) std::array<real, tensor::QcorrLocal::size()> qcorrLocal{};
  alignas(Alignment) std::array<real, tensor::QcorrNeighbor::size()> qcorrNeighbor{};
  kernel::computeFluxSolverLocal localSolver{};
  localSolver.fluxScale = fluxScale;
  localSolver.AplusT = aplus.data();
  localSolver.QgodLocal = godLocal.data();
  localSolver.QcorrLocal = qcorrLocal.data();
  localSolver.T = matT.data();
  localSolver.Tinv = matTinv.data();
  localSolver.star(0) = star.data();
  localSolver.execute();
  kernel::computeFluxSolverNeighbor neighborSolver{};
  neighborSolver.fluxScale = fluxScale;
  neighborSolver.AminusT = aminus.data();
  neighborSolver.QgodNeighbor = godNeighbor.data();
  neighborSolver.QcorrNeighbor = qcorrNeighbor.data();
  neighborSolver.T = matT.data();
  neighborSolver.Tinv = identityTinv.data();
  neighborSolver.star(0) = star.data();
  neighborSolver.execute();

  // a map that mixes every quantity of the Riemann problem with every other one
  alignas(Alignment) std::array<real, tensor::dirichletMap::size()> dirichletMap{};
  for (auto& value : dirichletMap) {
    value = static_cast<real>(gauss(rng));
  }
  // the part the fold adds, from an operator of zeros
  alignas(Alignment) std::array<real, tensor::AplusT::size()> folded{};
  kernel::foldDirichlet foldSolver{};
  foldSolver.AplusT = folded.data();
  foldSolver.AminusT = aminus.data();
  foldSolver.Tinv = matTinv.data();
  foldSolver.dirichletMap = dirichletMap.data();
  foldSolver.execute();

  const auto aplusDense = denseOf(init::AplusT::view::create(aplus.data()), NQ, Written);
  const auto foldedDense = denseOf(init::AplusT::view::create(folded.data()), NQ, Written);

  alignas(Alignment) std::array<real, tensor::I::size()> dofs{};
  for (auto& value : dofs) {
    value = static_cast<real>(gauss(rng));
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    LocalIntegrationData localData{};
    NeighboringIntegrationData neighboringData{};
    std::copy(aminus.begin(), aminus.end(), neighboringData.nAmNm1[face]);
    if constexpr (appliesFoldApart(FaceType::Dirichlet)) {
      // the setup keeps the folded part alone
      std::copy(folded.begin(), folded.end(), localData.nApNm1[face]);
      initializer::CurvedCell::setConstantRotation(matT.data(), localData.faceRotation[face]);
      // a face without a neighbour takes its own part negated against the trace of zeros
      setScalars(viewGodLocal, viewStar, -fluxScale, face, localData);
    } else {
      for (std::size_t i = 0; i < aplus.size(); ++i) {
        localData.nApNm1[face][i] = aplus[i] + folded[i];
      }
    }

    for (std::size_t side = 0; side < 2; ++side) {
      alignas(Alignment) FaceTraceT own{};
      alignas(Alignment) const FaceTraceT zero{};
      writeFaceTrace(dofs.data(), face, side, own.data());

      alignas(Alignment) std::array<real, tensor::Q::size()> update{};
      applyFaceCorrector(
          own.data(), zero.data(), localData, neighboringData, face, side, update.data());
      if constexpr (appliesFoldApart(FaceType::Dirichlet)) {
        kernel::foldedBoundaryFlux foldedFlux{};
        foldedFlux.bindGlobals(seissol::Pool::host());
        bindLiftTarget(foldedFlux, update.data());
        foldedFlux.I = dofs.data();
        foldedFlux.AplusT = localData.nApNm1[face];
        foldedFlux.execute(face);
      }

      auto viewI = init::I::view::create(dofs.data());
      auto viewQ = init::Q::view::create(update.data());
      for (std::size_t simulation = 0; simulation < multisim::NumSimulations; ++simulation) {
        const auto field = fieldOf(viewI, simulation, NQ);
        const auto got = fieldOf(viewQ, simulation, Written);
        const Eigen::MatrixXd expected =
            volume.rDivM[face] * volume.fMrT[face] * field * (aplusDense + foldedDense);

        constexpr double Tolerance = std::is_same_v<real, double> ? 1e-10 : 1e-4;
        for (std::size_t column = 0; column < Written; ++column) {
          const double scale = std::max(1.0, expected.col(column).cwiseAbs().maxCoeff());
          for (Eigen::Index row = 0; row < expected.rows(); ++row) {
            REQUIRE(got(row, column) ==
                    doctest::Approx(expected(row, column)).epsilon(Tolerance).scale(scale));
          }
        }
      }
    }
  }
}

} // namespace facecorrectortest

TEST_CASE("The face corrector is the flux of both sides of a face") {
  for (const auto faceType : facecorrectortest::LinearFaceTypes) {
    facecorrectortest::compareAgainstVolumeForm(faceType);
  }
}

TEST_CASE("A boundary condition folded into the flux of the cell's own side") {
  facecorrectortest::compareFoldedBoundary();
}

} // namespace seissol::unit_test

#endif // SEISSOL_KERNELS_LINEARCKANELASTIC

#endif // SEISSOL_TESTS_MODEL_FACECORRECTOR_T_H_
