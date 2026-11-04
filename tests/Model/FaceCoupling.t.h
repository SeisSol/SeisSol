// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Equations/Datastructures.h"
#include "Geometry/MeshDefinition.h"
#include "Initializer/BasicTypedefs.h"
#include "Model/Common.h"
#include "Model/FaceCoupling.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>
#include <yateto.h>

namespace seissol::unit_test {

namespace facecouplingtest {

inline model::ElasticMaterial rock() { return model::ElasticMaterial({2700, 3.2e10, 3.4e10}); }

inline model::ElasticMaterial softRock() { return model::ElasticMaterial({2400, 1.1e10, 1.9e10}); }

inline model::AcousticMaterial water() { return model::AcousticMaterial({1000, 2.25e9}); }

/// The Voigt stiffness of an isotropic solid with the Lame parameters `lambda` and `mu`, with a
/// perturbation of all its entries by up to a tenth of `mu` times `amount`, which keeps it positive
/// definite.
inline Eigen::Matrix<double, 6, 6> perturbedStiffness(double lambda, double mu, double amount) {
  Eigen::Matrix<double, 6, 6> stiffness = Eigen::Matrix<double, 6, 6>::Zero();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      stiffness(i, j) = lambda + (i == j ? 2 * mu : 0);
    }
    stiffness(3 + i, 3 + i) = mu;
  }
  for (int i = 0; i < 6; ++i) {
    for (int j = i; j < 6; ++j) {
      const double perturbation = amount * 0.1 * mu * std::sin(1.0 + 3.0 * i + 0.7 * j);
      stiffness(i, j) += perturbation;
      stiffness(j, i) = stiffness(i, j);
    }
  }
  return stiffness;
}

/// An anisotropic solid of the density `rho` and the Voigt stiffness `stiffness` (in the order xx,
/// yy, zz, yz, xz, xy of its parameters).
inline model::AnisotropicMaterial anisotropic(double rho,
                                              const Eigen::Matrix<double, 6, 6>& stiffness) {
  std::vector<double> values{rho};
  for (int i = 0; i < 6; ++i) {
    for (int j = i; j < 6; ++j) {
      values.push_back(stiffness(i, j));
    }
  }
  return model::AnisotropicMaterial(values);
}

/// An anisotropic solid without any symmetry, of the magnitudes of rock.
inline model::AnisotropicMaterial triclinic() {
  return anisotropic(2700, perturbedStiffness(3.4e10, 3.2e10, 1));
}

/// The impedance of a half space of an anisotropic solid of the density `rho` and the Voigt
/// stiffness `stiffness` towards a face with the unit normal `normal`, in the global coordinates:
/// the square root of rho times the Christoffel matrix of the normal.
inline Eigen::Matrix3d christoffelImpedance(double rho,
                                            const Eigen::Matrix<double, 6, 6>& stiffness,
                                            const Eigen::Vector3d& normal) {
  // the Voigt index of a pair of directions, as the parameters number them
  constexpr std::array<std::array<int, 3>, 3> Voigt{{{0, 5, 4}, {5, 1, 3}, {4, 3, 2}}};
  Eigen::Matrix3d christoffel = Eigen::Matrix3d::Zero();
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 3; ++k) {
      for (int j = 0; j < 3; ++j) {
        for (int l = 0; l < 3; ++l) {
          christoffel(i, k) += stiffness(Voigt[i][j], Voigt[k][l]) * normal(j) * normal(l);
        }
      }
    }
  }
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(rho * christoffel);
  return solver.eigenvectors() * solver.eigenvalues().cwiseSqrt().asDiagonal() *
         solver.eigenvectors().transpose();
}

/// The scale of the impedances in the tests.
constexpr double ImpedanceScale = 1e7;

/// The quantity `quantity` of the Riemann problem of the material `MaterialT`, without units: a
/// traction, divided by the scale of the impedances; a velocity, as it is.
template <typename MaterialT>
double quantityScale(Eigen::Index quantity) {
  return quantity < static_cast<Eigen::Index>(model::RiemannQuantitiesOf<MaterialT>::VelocityOffset)
             ? 1 / ImpedanceScale
             : 1.0;
}

/// Checks the rows `rows` of the matrix `computed` (from the quantities of the Riemann problem of
/// `ColumnT` to those of `RowT`) against `reference`, without units.
template <typename RowT, typename ColumnT>
void checkRows(const Eigen::MatrixXd& computed,
               const Eigen::MatrixXd& reference,
               const std::vector<Eigen::Index>& rows) {
  REQUIRE(computed.rows() == reference.rows());
  REQUIRE(computed.cols() == reference.cols());
  for (const auto row : rows) {
    for (Eigen::Index column = 0; column < computed.cols(); ++column) {
      CAPTURE(row);
      CAPTURE(column);
      const double scale = quantityScale<RowT>(row) / quantityScale<ColumnT>(column);
      CHECK(computed(row, column) * scale ==
            doctest::Approx(reference(row, column) * scale).epsilon(1e-10).scale(1.0));
    }
  }
}

/// The Godunov state of the cell of the material `local` at a face towards `neighbor` of the same
/// material, with the normal along x: the parts from the cell and from the neighbor.
template <typename MaterialT>
std::array<Eigen::MatrixXd, 2> godunovState(const MaterialT& local, const MaterialT& neighbor) {
  constexpr auto Count = static_cast<unsigned>(MaterialT::NumQuantities);
  std::vector<double> localData(Count * Count);
  std::vector<double> neighborData(Count * Count);
  yateto::DenseTensorView<2, double, unsigned> localView(localData.data(), {Count, Count});
  yateto::DenseTensorView<2, double, unsigned> neighborView(neighborData.data(), {Count, Count});
  model::getTransposedGodunovState(local, neighbor, FaceType::Regular, localView, neighborView);
  // stored transposed
  return {Eigen::Map<const Eigen::MatrixXd>(localData.data(), Count, Count).transpose(),
          Eigen::Map<const Eigen::MatrixXd>(neighborData.data(), Count, Count).transpose()};
}

/// The rows of the traction and the velocity among the quantities of the Riemann problem of the
/// material `MaterialT`, which its flux reads in the coordinates of the face.
template <typename MaterialT>
std::vector<Eigen::Index> interfaceRowIndices() {
  const Eigen::MatrixXd rows = model::interfaceRows<MaterialT>();
  std::vector<Eigen::Index> indices;
  for (Eigen::Index quantity = 0; quantity < rows.cols(); ++quantity) {
    if (rows.col(quantity).cwiseAbs().maxCoeff() > 0) {
      indices.push_back(quantity);
    }
  }
  return indices;
}

/// The unit normal and tangents of a face in general position.
struct Face {
  CoordinateT normal{};
  CoordinateT tangent1{};
  CoordinateT tangent2{};

  [[nodiscard]] Eigen::Vector3d n() const { return {normal[0], normal[1], normal[2]}; }

  /// The rotation from the coordinates of the face into the global ones.
  [[nodiscard]] Eigen::Matrix3d rotation() const {
    Eigen::Matrix3d matrix;
    matrix << normal[0], tangent1[0], tangent2[0], normal[1], tangent1[1], tangent2[1], normal[2],
        tangent1[2], tangent2[2];
    return matrix;
  }
};

inline Face generalFace() {
  const Eigen::Vector3d n = Eigen::Vector3d(1, -2, 3).normalized();
  const Eigen::Vector3d t1 = n.cross(Eigen::Vector3d(0, 0, 1)).normalized();
  const Eigen::Vector3d t2 = n.cross(t1);
  return {{n(0), n(1), n(2)}, {t1(0), t1(1), t1(2)}, {t2(0), t2(1), t2(2)}};
}

/// The stress tensor of the Voigt vector `voigt` (xx, yy, zz, xy, yz, xz).
inline Eigen::Matrix3d stressTensor(const Eigen::VectorXd& voigt) {
  Eigen::Matrix3d stress;
  stress << voigt(0), voigt(3), voigt(5), voigt(3), voigt(1), voigt(4), voigt(5), voigt(4),
      voigt(2);
  return stress;
}

/// The traction on the face and the velocity of the state `state` of the quantities of the Riemann
/// problem of `MaterialT`, in the global coordinates.
template <typename MaterialT>
std::array<Eigen::Vector3d, 2> tractionAndVelocity(const Eigen::VectorXd& state,
                                                   const Eigen::Vector3d& normal) {
  constexpr auto Velocity = model::RiemannQuantitiesOf<MaterialT>::VelocityOffset;
  const Eigen::Vector3d velocity = state.segment<3>(Velocity);
  if constexpr (Velocity == 1) {
    // the mean normal stress of a fluid
    return {state(0) * normal, velocity};
  } else {
    return {stressTensor(state.head<6>()) * normal, velocity};
  }
}

/// A state of the quantities of the Riemann problem of `MaterialT`, with stresses of the order of
/// MPa and velocities of the order of m/s.
template <typename MaterialT>
Eigen::VectorXd sampleState(double seed) {
  constexpr auto Count = model::RiemannQuantitiesOf<MaterialT>::NumQuantities;
  constexpr auto Velocity = model::RiemannQuantitiesOf<MaterialT>::VelocityOffset;
  Eigen::VectorXd state(Count);
  for (Eigen::Index i = 0; i < state.size(); ++i) {
    const double magnitude = i < static_cast<Eigen::Index>(Velocity) ? 1e6 : 1.0;
    state(i) = magnitude * std::sin(seed + 1.9 * static_cast<double>(i));
  }
  return state;
}

/// The state at a face with the face `face` between a cell of the material `local` in the state
/// `localState` and a neighbor `neighbor` in the state `neighborState`, all in the global
/// coordinates and in the quantities of the cell: the part of the cell through its Godunov state
/// rotated as T and Tinv of getFaceRotationMatrix rotate it, the part of the neighbor through the
/// coupling matrix of the face.
template <typename MaterialT, typename NeighborT>
Eigen::VectorXd stateAtFace(const MaterialT& local,
                            const NeighborT& neighbor,
                            const Face& face,
                            const Eigen::VectorXd& localState,
                            const Eigen::VectorXd& neighborState) {
  std::array<double, 36> bond{};
  model::getBondMatrix(face.normal, face.tangent1, face.tangent2, bond);
  const auto godunov = model::coupledGodunov(model::getRotatedMaterialCoefficients(bond, local),
                                             model::getRotatedMaterialCoefficients(bond, neighbor));
  constexpr auto Groups = model::RiemannQuantitiesOf<MaterialT>::PrimaryGroups;
  const Eigen::MatrixXd toGlobal =
      model::quantityRotation<false>(Groups, face.normal, face.tangent1, face.tangent2);
  const Eigen::MatrixXd toFace =
      model::quantityRotation<true>(Groups, face.normal, face.tangent1, face.tangent2);
  const Eigen::MatrixXd coupling = model::coupledNeighborMap<MaterialT, NeighborT>(
      godunov, face.normal, face.tangent1, face.tangent2);
  REQUIRE(coupling.rows() == neighborState.size());
  REQUIRE(coupling.cols() == localState.size());
  return toGlobal * godunov.local * toFace * localState + coupling.transpose() * neighborState;
}

/// Checks that the traction `traction` and the velocity `velocity` at a face with the unit normal
/// `normal`, in the global coordinates, solve the Riemann problem between a cell with the traction
/// `localTraction`, the velocity `localVelocity` and the impedance `localImpedance` and a neighbor
/// with `neighborTraction`, `neighborVelocity` and `neighborImpedance`: the waves into the cell and
/// into the neighbor carry the jumps along their impedances. A fluid has no shear impedance, and
/// the shear traction on its side of the face is zero.
inline void checkRiemannSolution(const Eigen::Vector3d& traction,
                                 const Eigen::Vector3d& velocity,
                                 const Eigen::Vector3d& localTraction,
                                 const Eigen::Vector3d& localVelocity,
                                 const Eigen::Matrix3d& localImpedance,
                                 const Eigen::Vector3d& neighborTraction,
                                 const Eigen::Vector3d& neighborVelocity,
                                 const Eigen::Matrix3d& neighborImpedance) {
  const auto checkVector = [](const Eigen::Vector3d& actual, const Eigen::Vector3d& expected) {
    CAPTURE(actual.transpose());
    CAPTURE(expected.transpose());
    // the tractions are of the order of MPa
    CHECK((actual - expected).norm() <= 1e-10 * (1e6 + expected.norm()));
  };
  checkVector(traction - localTraction, localImpedance * (velocity - localVelocity));
  checkVector(traction - neighborTraction, -neighborImpedance * (velocity - neighborVelocity));
}

/// The impedance of a fluid of the material `fluid` towards a face with the unit normal `normal`,
/// in the global coordinates.
inline Eigen::Matrix3d fluidImpedance(const model::AcousticMaterial& fluid,
                                      const Eigen::Vector3d& normal) {
  return std::sqrt(fluid.rho * fluid.lambda) * normal * normal.transpose();
}

/// Checks the state at the face `face` between the solid `solid` of the impedance `solidImpedance`
/// (in the global coordinates) and the fluid `fluid`, from both sides.
template <typename SolidT>
void checkSolidAndFluid(const SolidT& solid,
                        const Eigen::Matrix3d& solidImpedance,
                        const model::AcousticMaterial& fluid,
                        const Face& face) {
  using FluidT = model::AcousticMaterial;
  const Eigen::Vector3d normal = face.n();
  const auto fluidZ = fluidImpedance(fluid, normal);

  const auto solidState = sampleState<SolidT>(0.4);
  const auto fluidState = sampleState<FluidT>(2.3);
  const auto [solidTraction, solidVelocity] = tractionAndVelocity<SolidT>(solidState, normal);
  const auto [fluidTraction, fluidVelocity] = tractionAndVelocity<FluidT>(fluidState, normal);

  SUBCASE("A solid with a fluid neighbor") {
    const auto state = stateAtFace(solid, fluid, face, solidState, fluidState);
    const auto [traction, velocity] = tractionAndVelocity<SolidT>(state, normal);
    checkRiemannSolution(traction,
                         velocity,
                         solidTraction,
                         solidVelocity,
                         solidImpedance,
                         fluidTraction,
                         fluidVelocity,
                         fluidZ);
    // no shear traction on the face of a fluid
    const Eigen::Vector3d shear = traction - normal.dot(traction) * normal;
    CHECK(shear.norm() <= 1e-10 * 1e6);
  }

  SUBCASE("A fluid with a solid neighbor") {
    const auto state = stateAtFace(fluid, solid, face, fluidState, solidState);
    const auto [traction, velocity] = tractionAndVelocity<FluidT>(state, normal);
    // the tangential velocity at the face is the one that the wave into the solid leaves behind,
    // as the fluid exerts no shear traction on it
    checkRiemannSolution(traction,
                         velocity,
                         fluidTraction,
                         fluidVelocity,
                         fluidZ,
                         solidTraction,
                         solidVelocity,
                         solidImpedance);
  }
}

} // namespace facecouplingtest

TEST_CASE("The state at a face carries the jumps of the traction and the velocity along the "
          "impedances of both sides" *
          doctest::test_suite("model")) {
  using namespace facecouplingtest;
  // symmetric and positive definite, with impedances of rock
  Eigen::Matrix3d zLocal;
  zLocal << 1.8e7, 0.2e7, -0.1e7, 0.2e7, 1.1e7, 0.05e7, -0.1e7, 0.05e7, 0.9e7;
  Eigen::Matrix3d zNeighbor;
  zNeighbor << 1.2e7, -0.1e7, 0.3e7, -0.1e7, 0.7e7, 0.1e7, 0.3e7, 0.1e7, 0.8e7;

  const auto check = [](const Eigen::Matrix3d& zCell, const Eigen::Matrix3d& zOther) {
    const auto godunov = model::interfaceGodunov(zCell, zOther);
    Eigen::Matrix<double, 6, 1> local;
    Eigen::Matrix<double, 6, 1> neighbor;
    for (int i = 0; i < 6; ++i) {
      local(i) = (i < 3 ? 1e6 : 1.0) * std::sin(0.5 + i);
      neighbor(i) = (i < 3 ? 1e6 : 1.0) * std::cos(1.5 + 2 * i);
    }
    const Eigen::Matrix<double, 6, 1> state = godunov.local * local + godunov.neighbor * neighbor;
    checkRiemannSolution(state.head<3>(),
                         state.tail<3>(),
                         local.head<3>(),
                         local.tail<3>(),
                         zCell,
                         neighbor.head<3>(),
                         neighbor.tail<3>(),
                         zOther);

    // a state that is the same on both sides is the one at the face
    const model::InterfaceMap sum = godunov.local + godunov.neighbor;
    CHECK((sum - model::InterfaceMap::Identity()).norm() <= 1e-12);
  };

  SUBCASE("Two solids") { check(zLocal, zNeighbor); }

  // a fluid only resists a normal motion
  const Eigen::Matrix3d zFluid = Eigen::Vector3d(1.5e6, 0, 0).asDiagonal();
  SUBCASE("A solid and a fluid") { check(zLocal, zFluid); }
  SUBCASE("A fluid and a solid") { check(zFluid, zNeighbor); }
}

TEST_CASE("Between solids, the state at a face is the one of their Godunov state where their "
          "flux reads it" *
          doctest::test_suite("model")) {
  using namespace facecouplingtest;

  SUBCASE("Two isotropic solids") {
    const auto reference = godunovState(rock(), softRock());
    const auto godunov = model::coupledGodunov(rock(), softRock());
    const auto rows = interfaceRowIndices<model::ElasticMaterial>();
    checkRows<model::ElasticMaterial, model::ElasticMaterial>(godunov.local, reference[0], rows);
    checkRows<model::ElasticMaterial, model::ElasticMaterial>(godunov.neighbor, reference[1], rows);
  }

  SUBCASE("Two anisotropic solids") {
    const auto local = triclinic();
    const auto neighbor = anisotropic(2400, perturbedStiffness(1.9e10, 1.1e10, -0.8));
    const auto reference = godunovState(local, neighbor);
    const auto godunov = model::coupledGodunov(local, neighbor);
    const auto rows = interfaceRowIndices<model::AnisotropicMaterial>();
    checkRows<model::AnisotropicMaterial, model::AnisotropicMaterial>(
        godunov.local, reference[0], rows);
    checkRows<model::AnisotropicMaterial, model::AnisotropicMaterial>(
        godunov.neighbor, reference[1], rows);
  }

  SUBCASE("A solid with a fluid neighbor, which is a solid without shear modulus") {
    // the fluid with its pressure as the isotropic stress of a solid
    Eigen::MatrixXd asSolid = Eigen::MatrixXd::Zero(9, 4);
    asSolid(0, 0) = asSolid(1, 0) = asSolid(2, 0) = 1;
    asSolid.bottomRightCorner<3, 3>().setIdentity();

    const auto reference = godunovState(rock(), model::ElasticMaterial(water()));
    const auto godunov = model::coupledGodunov(rock(), water());
    const auto rows = interfaceRowIndices<model::ElasticMaterial>();
    checkRows<model::ElasticMaterial, model::ElasticMaterial>(godunov.local, reference[0], rows);
    checkRows<model::ElasticMaterial, model::AcousticMaterial>(
        godunov.neighbor, reference[1] * asSolid, rows);
  }

  SUBCASE("A fluid with a solid neighbor, which is a fluid of its P-wave modulus") {
    // the normal stress of the solid as the pressure of a fluid
    Eigen::MatrixXd asFluid = Eigen::MatrixXd::Zero(4, 9);
    asFluid(0, 0) = 1;
    asFluid.bottomRightCorner<3, 3>().setIdentity();

    const auto solid = rock();
    const auto reference =
        godunovState(water(), model::AcousticMaterial({solid.rho, solid.lambda + 2 * solid.mu}));
    const auto godunov = model::coupledGodunov(water(), solid);
    // the pressure and the normal velocity, which the flux of a fluid reads
    const std::vector<Eigen::Index> rows{0, 1};
    checkRows<model::AcousticMaterial, model::AcousticMaterial>(godunov.local, reference[0], rows);
    checkRows<model::AcousticMaterial, model::ElasticMaterial>(
        godunov.neighbor, reference[1] * asFluid, rows);
  }
}

TEST_CASE("Between a solid and a fluid, the state at a face in general position solves their "
          "Riemann problem" *
          doctest::test_suite("model")) {
  using namespace facecouplingtest;
  const auto face = generalFace();

  SUBCASE("An isotropic solid") {
    const auto solid = rock();
    const Eigen::Vector3d normal = face.n();
    const double zp = std::sqrt(solid.rho * (solid.lambda + 2 * solid.mu));
    const double zs = std::sqrt(solid.rho * solid.mu);
    const Eigen::Matrix3d impedance =
        zp * normal * normal.transpose() +
        zs * (Eigen::Matrix3d::Identity() - normal * normal.transpose());
    checkSolidAndFluid(solid, impedance, water(), face);
  }

  SUBCASE("An anisotropic solid") {
    const auto stiffness = perturbedStiffness(3.4e10, 3.2e10, 1);
    checkSolidAndFluid(anisotropic(2700, stiffness),
                       christoffelImpedance(2700, stiffness, face.n()),
                       water(),
                       face);
  }
}

} // namespace seissol::unit_test
