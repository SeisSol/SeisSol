// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Equations/Datastructures.h"
#include "Equations/poroelastic/Model/Impedance.h"
#include "Geometry/MeshDefinition.h"
#include "Initializer/BasicTypedefs.h"
#include "Model/Common.h"
#include "Model/FaceCoupling.h"
#include "Model/Quantities.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>
#include <yateto.h>

namespace seissol::unit_test {

namespace facecouplingtest {

constexpr double Unbounded = std::numeric_limits<double>::infinity();

inline model::ElasticMaterial rock() { return model::ElasticMaterial({2700, 3.2e10, 3.4e10}); }

inline model::ElasticMaterial softRock() { return model::ElasticMaterial({2400, 1.1e10, 1.9e10}); }

inline model::AcousticMaterial water() { return model::AcousticMaterial({1000, 2.25e9}); }

/// A sandstone and a more porous one, saturated with water.
inline model::PoroElasticMaterial sandstone() {
  return model::PoroElasticMaterial({40e9, 2500, 12e9, 10e9, 0.2, 6e-13, 3, 2.5e9, 1040, 0.001});
}

inline model::PoroElasticMaterial porousSandstone() {
  return model::PoroElasticMaterial({30e9, 2600, 10e9, 10e9, 0.3, 1e-13, 3, 2.5e9, 1040, 0.001});
}

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

/// The Voigt stiffness of the anisotropic solid `solid`.
inline Eigen::Matrix<double, 6, 6> stiffnessOf(const model::AnisotropicMaterial& solid) {
  const std::array<double, 21> upper{
      solid.c11, solid.c12, solid.c13, solid.c14, solid.c15, solid.c16, solid.c22,
      solid.c23, solid.c24, solid.c25, solid.c26, solid.c33, solid.c34, solid.c35,
      solid.c36, solid.c44, solid.c45, solid.c46, solid.c55, solid.c56, solid.c66};
  Eigen::Matrix<double, 6, 6> stiffness;
  std::size_t next = 0;
  for (int i = 0; i < 6; ++i) {
    for (int j = i; j < 6; ++j) {
      stiffness(i, j) = upper[next++];
      stiffness(j, i) = stiffness(i, j);
    }
  }
  return stiffness;
}

/// The unit normal and tangents of a face in general position.
struct Face {
  CoordinateT normal{};
  CoordinateT tangent1{};
  CoordinateT tangent2{};

  [[nodiscard]] Eigen::Vector3d n() const { return {normal[0], normal[1], normal[2]}; }

  /// The rotation from the coordinates of the face into the global ones, with the pore fluid as
  /// the fourth interface variable.
  [[nodiscard]] Eigen::Matrix4d rotation() const {
    Eigen::Matrix4d matrix = Eigen::Matrix4d::Identity();
    for (int i = 0; i < 3; ++i) {
      matrix(i, 0) = normal[i];
      matrix(i, 1) = tangent1[i];
      matrix(i, 2) = tangent2[i];
    }
    return matrix;
  }

  /// The same face, seen from the other side.
  [[nodiscard]] Face flipped() const {
    Face face;
    for (int i = 0; i < 3; ++i) {
      face.normal[i] = -normal[i];
      face.tangent1[i] = tangent1[i];
      face.tangent2[i] = -tangent2[i];
    }
    return face;
  }
};

inline Face generalFace() {
  const Eigen::Vector3d n = Eigen::Vector3d(1, -2, 3).normalized();
  const Eigen::Vector3d t1 = n.cross(Eigen::Vector3d(0, 0, 1)).normalized();
  const Eigen::Vector3d t2 = n.cross(t1);
  return {{n(0), n(1), n(2)}, {t1(0), t1(1), t1(2)}, {t2(0), t2(1), t2(2)}};
}

/// The impedance of a solid without pores of the impedance `impedance` (3x3), in the interface
/// variables.
inline Eigen::Matrix4d solidImpedance(const Eigen::Matrix3d& impedance) {
  Eigen::Matrix4d result = Eigen::Matrix4d::Zero();
  result.topLeftCorner<3, 3>() = impedance;
  return result;
}

/// The impedance of a half space of the material `material` towards the face `face`, in the
/// interface variables of the global coordinates (whose fourth one is the pore fluid along the
/// normal of the face), computed without the face coordinates: for an anisotropic solid, the square
/// root of rho times the Christoffel matrix of the normal; for a fluid, its impedance along the
/// normal both for the normal stress and the pore pressure; for a poroelastic solid, whose frame is
/// isotropic, the one of a face with the normal along x, rotated.
template <typename MaterialT>
Eigen::Matrix4d globalImpedance(const MaterialT& material, const Face& face) {
  const Eigen::Vector3d normal = face.n();
  if constexpr (std::is_same_v<MaterialT, model::AnisotropicMaterial>) {
    // the Voigt index of a pair of directions, as the parameters number them
    constexpr std::array<std::array<int, 3>, 3> Voigt{{{0, 5, 4}, {5, 1, 3}, {4, 3, 2}}};
    const auto stiffness = stiffnessOf(material);
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
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(material.rho * christoffel);
    return solidImpedance(solver.eigenvectors() * solver.eigenvalues().cwiseSqrt().asDiagonal() *
                          solver.eigenvectors().transpose());
  } else if constexpr (std::is_same_v<MaterialT, model::ElasticMaterial>) {
    const double zp = std::sqrt(material.rho * (material.lambda + 2 * material.mu));
    const double zs = std::sqrt(material.rho * material.mu);
    return solidImpedance(zp * normal * normal.transpose() +
                          zs * (Eigen::Matrix3d::Identity() - normal * normal.transpose()));
  } else if constexpr (std::is_same_v<MaterialT, model::AcousticMaterial>) {
    Eigen::Vector4d direction;
    direction << normal, 1;
    return std::sqrt(material.rho * material.lambda) * direction * direction.transpose();
  } else {
    static_assert(std::is_same_v<MaterialT, model::PoroElasticMaterial>);
    return face.rotation() *
           model::ImpedanceCompute<model::PoroElasticMaterial>::impedance(material) *
           face.rotation().transpose();
  }
}

/// The scale of the impedances in the tests.
constexpr double ImpedanceScale = 1e7;

/// The quantity `quantity` of the Riemann problem of the material `MaterialT`, without units: a
/// traction or a pressure, divided by the scale of the impedances; a velocity, as it is.
template <typename MaterialT>
double quantityScale(Eigen::Index quantity) {
  std::size_t offset = 0;
  for (const auto& group : model::RiemannQuantitiesOf<MaterialT>::PrimaryGroups) {
    if (static_cast<std::size_t>(quantity) < offset + group.extent()) {
      const bool traction =
          group.role == model::FaceRole::Traction || group.role == model::FaceRole::ExtraTraction;
      return traction ? 1 / ImpedanceScale : 1.0;
    }
    offset += group.extent();
  }
  return 1.0;
}

/// Checks the rows `rows` of the matrix `computed` (from the quantities of the Riemann problem of
/// `ColumnT` to those of `RowT`) against `reference`, without units, to `epsilon`.
template <typename RowT, typename ColumnT>
void checkRows(const Eigen::MatrixXd& computed,
               const Eigen::MatrixXd& reference,
               const std::vector<Eigen::Index>& rows,
               double epsilon = 1e-10) {
  REQUIRE(computed.rows() == reference.rows());
  REQUIRE(computed.cols() == reference.cols());
  for (const auto row : rows) {
    for (Eigen::Index column = 0; column < computed.cols(); ++column) {
      CAPTURE(row);
      CAPTURE(column);
      const double scale = quantityScale<RowT>(row) / quantityScale<ColumnT>(column);
      CHECK(computed(row, column) * scale ==
            doctest::Approx(reference(row, column) * scale).epsilon(epsilon).scale(1.0));
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

/// The quantities of the Riemann problem of the material `MaterialT` that its flux reads in the
/// coordinates of the face: those of the interface variables.
template <typename MaterialT>
std::vector<Eigen::Index> interfaceQuantities() {
  const Eigen::MatrixXd in = model::interfaceIn<MaterialT>();
  std::vector<Eigen::Index> indices;
  for (Eigen::Index quantity = 0; quantity < in.cols(); ++quantity) {
    if (in.col(quantity).cwiseAbs().maxCoeff() > 0) {
      indices.push_back(quantity);
    }
  }
  return indices;
}

/// A traction and a velocity in the interface variables, with tractions of the order of MPa and
/// velocities of the order of m/s.
inline Eigen::Matrix<double, 8, 1> sampleInterfaceState(double seed) {
  Eigen::Matrix<double, 8, 1> state;
  for (int i = 0; i < 8; ++i) {
    state(i) = (i < 4 ? 1e6 : 1.0) * std::sin(seed + 1.9 * i);
  }
  return state;
}

/// Checks that the traction `traction` on the side of the cell and the velocity `velocity` at a
/// face solve the Riemann problem between the cell with the traction `localTraction`, the velocity
/// `localVelocity` and the impedance `localImpedance` and the neighbor with `neighborTraction`,
/// `neighborVelocity` and `neighborImpedance`, all in the interface variables, with a skin of the
/// conductance `conductance` between them: the waves into the cell and into the neighbor carry the
/// jumps along their impedances, and the pore pressure jumps across the skin with the filtration
/// through it, which a sealed face has none of.
inline void checkRiemannSolution(const Eigen::Vector4d& traction,
                                 const Eigen::Vector4d& velocity,
                                 const Eigen::Vector4d& localTraction,
                                 const Eigen::Vector4d& localVelocity,
                                 const Eigen::Matrix4d& localImpedance,
                                 const Eigen::Vector4d& neighborTraction,
                                 const Eigen::Vector4d& neighborVelocity,
                                 const Eigen::Matrix4d& neighborImpedance,
                                 double conductance) {
  const auto checkVector = [](const Eigen::VectorXd& actual, const Eigen::VectorXd& expected) {
    CAPTURE(actual.transpose());
    CAPTURE(expected.transpose());
    // the tractions are of the order of MPa
    CHECK((actual - expected).norm() <= 1e-10 * (1e6 + expected.norm()));
  };
  checkVector(traction - localTraction, localImpedance * (velocity - localVelocity));

  // the traction on the side of the neighbor
  Eigen::Vector4d neighborSide = traction;
  if (conductance > 0 && !std::isinf(conductance)) {
    neighborSide(3) += velocity(3) / conductance;
  }
  const Eigen::Vector4d expected =
      neighborTraction - neighborImpedance * (velocity - neighborVelocity);
  if (conductance == 0) {
    checkVector(neighborSide.head<3>(), expected.head<3>());
    CHECK(std::abs(velocity(3)) <= 1e-12);
  } else {
    checkVector(neighborSide, expected);
  }
}

/// The stress tensor of the Voigt vector `voigt` (xx, yy, zz, xy, yz, xz).
inline Eigen::Matrix3d stressTensor(const Eigen::VectorXd& voigt) {
  Eigen::Matrix3d stress;
  stress << voigt(0), voigt(3), voigt(5), voigt(3), voigt(1), voigt(4), voigt(5), voigt(4),
      voigt(2);
  return stress;
}

/// The interface variables of the state `state` of the quantities of the Riemann problem of
/// `MaterialT` at a face with the unit normal `normal`, in the global coordinates: the traction on
/// the face and the pore pressure (-p), and the velocity and the filtration through the face. A
/// fluid takes its pressure as both the normal stress and the pore pressure.
template <typename MaterialT>
std::array<Eigen::Vector4d, 2> interfaceVariables(const Eigen::VectorXd& state,
                                                  const Eigen::Vector3d& normal) {
  using QuantitiesT = model::RiemannQuantitiesOf<MaterialT>;
  constexpr auto Velocity = QuantitiesT::VelocityOffset;
  Eigen::Vector4d traction = Eigen::Vector4d::Zero();
  Eigen::Vector4d velocity = Eigen::Vector4d::Zero();
  velocity.head<3>() = state.segment<3>(Velocity);
  if constexpr (std::is_same_v<QuantitiesT, model::AcousticMaterial>) {
    traction.head<3>() = state(0) * normal;
    traction(3) = state(0);
  } else {
    traction.head<3>() = stressTensor(state.head<6>()) * normal;
    if constexpr (std::is_same_v<QuantitiesT, model::PoroElasticMaterial>) {
      traction(3) = -state(9);
      velocity(3) = state.segment<3>(10).dot(normal);
    }
  }
  return {traction, velocity};
}

/// A state of the quantities of the Riemann problem of `MaterialT`, with stresses and pressures of
/// the order of MPa and velocities of the order of m/s.
template <typename MaterialT>
Eigen::VectorXd sampleState(double seed) {
  Eigen::VectorXd state(model::RiemannQuantitiesOf<MaterialT>::NumQuantities);
  for (Eigen::Index i = 0; i < state.size(); ++i) {
    const double magnitude = quantityScale<MaterialT>(i) < 1 ? 1e6 : 1.0;
    state(i) = magnitude * std::sin(seed + 1.9 * static_cast<double>(i));
  }
  return state;
}

/// The conductance of a face at the boundary of a region between cells of the materials `local`
/// and `neighbor`.
template <typename MaterialT, typename NeighborT>
double regionConductance(const MaterialT& local, const NeighborT& neighbor) {
  return model::seriesConductance(model::sideConductance(local, true),
                                  model::sideConductance(neighbor, true));
}

/// The state at the face `face` between a cell of the material `local` in the state `localState`
/// and a neighbor `neighbor` in the state `neighborState`, all in the global coordinates and in
/// the quantities of the cell, at the boundary of a region: the part of the cell through its
/// Godunov state rotated as T and Tinv of getFaceRotationMatrix rotate it, the part of the
/// neighbor through the coupling matrix of the face.
template <typename MaterialT, typename NeighborT>
Eigen::VectorXd stateAtFace(const MaterialT& local,
                            const NeighborT& neighbor,
                            const Face& face,
                            const Eigen::VectorXd& localState,
                            const Eigen::VectorXd& neighborState) {
  std::array<double, 36> bond{};
  model::getBondMatrix(face.normal, face.tangent1, face.tangent2, bond);
  const auto godunov = model::coupledGodunov(model::getRotatedMaterialCoefficients(bond, local),
                                             model::getRotatedMaterialCoefficients(bond, neighbor),
                                             regionConductance(local, neighbor));
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

/// Checks the state at the face `face` from a cell of the material `local` in the state
/// `localState` towards a neighbor `neighbor` in the state `neighborState` against the Riemann
/// problem between them in the global coordinates, and returns it. The cell has to be a solid or a
/// poroelastic solid, or the face sealed, so that its quantities hold all interface variables.
template <typename MaterialT, typename NeighborT>
Eigen::VectorXd checkFace(const MaterialT& local,
                          const NeighborT& neighbor,
                          const Face& face,
                          const Eigen::VectorXd& localState,
                          const Eigen::VectorXd& neighborState) {
  const Eigen::Vector3d normal = face.n();
  const auto state = stateAtFace(local, neighbor, face, localState, neighborState);
  const auto [localTraction, localVelocity] = interfaceVariables<MaterialT>(localState, normal);
  const auto [neighborTraction, neighborVelocity] =
      interfaceVariables<NeighborT>(neighborState, normal);
  const auto [traction, velocity] = interfaceVariables<MaterialT>(state, normal);
  checkRiemannSolution(traction,
                       velocity,
                       localTraction,
                       localVelocity,
                       globalImpedance(local, face),
                       neighborTraction,
                       neighborVelocity,
                       globalImpedance(neighbor, face),
                       regionConductance(local, neighbor));
  return state;
}

} // namespace facecouplingtest

TEST_CASE("The state at a face carries the jumps of the traction and the velocity along the "
          "impedances of both sides" *
          doctest::test_suite("model")) {
  using namespace facecouplingtest;
  const Eigen::Matrix4d zSolid = model::faceImpedance(rock());
  Eigen::Matrix3d zAnisotropic;
  zAnisotropic << 1.8e7, 0.2e7, -0.1e7, 0.2e7, 1.1e7, 0.05e7, -0.1e7, 0.05e7, 0.9e7;
  const Eigen::Matrix4d zFluid = model::faceImpedance(water());
  const Eigen::Matrix4d zPoro = model::faceImpedance(sandstone());
  const Eigen::Matrix4d zPorous = model::faceImpedance(porousSandstone());

  // `filtration`: whether both sides can have a pore pressure and a filtration
  const auto check = [](const Eigen::Matrix4d& zCell,
                        const Eigen::Matrix4d& zOther,
                        double conductance,
                        bool filtration) {
    CAPTURE(conductance);
    const auto godunov = model::interfaceGodunov(zCell, zOther, conductance);
    auto local = sampleInterfaceState(0.5);
    auto neighbor = sampleInterfaceState(1.7);
    if (!filtration) {
      local(3) = local(7) = neighbor(3) = neighbor(7) = 0;
    }
    const Eigen::Matrix<double, 8, 1> state = godunov.local * local + godunov.neighbor * neighbor;
    checkRiemannSolution(state.head<4>(),
                         state.tail<4>(),
                         local.head<4>(),
                         local.tail<4>(),
                         zCell,
                         neighbor.head<4>(),
                         neighbor.tail<4>(),
                         zOther,
                         conductance);

    // a state that is the same on both sides, without a filtration through the face, is the one
    // at the face
    auto uniform = sampleInterfaceState(2.9);
    uniform(7) = 0;
    if (!filtration) {
      uniform(3) = 0;
    }
    const Eigen::Matrix<double, 8, 1> same = (godunov.local + godunov.neighbor) * uniform;
    CHECK((same - uniform).head<4>().norm() <= 1e-12 * uniform.head<4>().norm());
    CHECK((same - uniform).tail<4>().norm() <= 1e-12 * uniform.tail<4>().norm());
  };

  SUBCASE("Two isotropic solids") { check(zSolid, zSolid * 0.7, 0, false); }
  SUBCASE("An isotropic and an anisotropic solid") {
    check(zSolid, solidImpedance(zAnisotropic), 0, false);
  }
  SUBCASE("A solid and a fluid") { check(zSolid, zFluid, 0, false); }
  SUBCASE("A fluid and an anisotropic solid") {
    check(zFluid, solidImpedance(zAnisotropic), 0, false);
  }
  SUBCASE("A poroelastic solid and a solid, sealed") { check(zPoro, zSolid, 0, true); }
  SUBCASE("A poroelastic solid and a fluid, open") { check(zPoro, zFluid, Unbounded, true); }
  SUBCASE("A fluid and a poroelastic solid, with a skin") {
    check(zFluid, zPoro, 1 / zPoro(3, 3), true);
  }
  SUBCASE("Two poroelastic solids, open") { check(zPoro, zPorous, Unbounded, true); }
  SUBCASE("Two poroelastic solids, with a skin") {
    check(zPoro, zPorous, 0.3 / (zPoro(3, 3) + zPorous(3, 3)), true);
  }
  SUBCASE("Two poroelastic solids, sealed") { check(zPoro, zPorous, 0, true); }
}

TEST_CASE("The state at a face is the one of the Godunov state of the materials where their flux "
          "reads it" *
          doctest::test_suite("model")) {
  using namespace facecouplingtest;

  SUBCASE("Two isotropic solids") {
    const auto reference = godunovState(rock(), softRock());
    const auto godunov = model::coupledGodunov(rock(), softRock(), 0);
    const auto rows = interfaceQuantities<model::ElasticMaterial>();
    checkRows<model::ElasticMaterial, model::ElasticMaterial>(godunov.local, reference[0], rows);
    checkRows<model::ElasticMaterial, model::ElasticMaterial>(godunov.neighbor, reference[1], rows);
  }

  SUBCASE("Two anisotropic solids") {
    const auto local = triclinic();
    const auto neighbor = anisotropic(2400, perturbedStiffness(1.9e10, 1.1e10, -0.8));
    const auto reference = godunovState(local, neighbor);
    const auto godunov = model::coupledGodunov(local, neighbor, 0);
    const auto rows = interfaceQuantities<model::AnisotropicMaterial>();
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
    const auto godunov = model::coupledGodunov(rock(), water(), 0);
    const auto rows = interfaceQuantities<model::ElasticMaterial>();
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
    const auto godunov = model::coupledGodunov(water(), solid, 0);
    // the pressure and the normal velocity, which the flux of a fluid reads
    const std::vector<Eigen::Index> rows{0, 1};
    checkRows<model::AcousticMaterial, model::AcousticMaterial>(godunov.local, reference[0], rows);
    checkRows<model::AcousticMaterial, model::ElasticMaterial>(
        godunov.neighbor, reference[1] * asFluid, rows);
  }

// the Godunov state of a poroelastic solid needs its setup, which only a build with a configuration
// of the STP kernels has (Equations/poroelastic/Model/Setup.h)
#ifdef SEISSOL_KERNELS_STP
  SUBCASE("Two poroelastic solids with open pores") {
    const auto reference = godunovState(sandstone(), porousSandstone());
    const auto godunov = model::coupledGodunov(sandstone(), porousSandstone(), Unbounded);
    const auto rows = interfaceQuantities<model::PoroElasticMaterial>();
    // the reference goes through a numerical eigendecomposition (cf. the Godunov state tests)
    constexpr double Epsilon = 1e-8;
    checkRows<model::PoroElasticMaterial, model::PoroElasticMaterial>(
        godunov.local, reference[0], rows, Epsilon);
    checkRows<model::PoroElasticMaterial, model::PoroElasticMaterial>(
        godunov.neighbor, reference[1], rows, Epsilon);
  }
#endif
}

TEST_CASE("Between different materials, the state at a face in general position solves their "
          "Riemann problem" *
          doctest::test_suite("model")) {
  using namespace facecouplingtest;
  const auto face = generalFace();
  const auto fluid = water();

  SUBCASE("An isotropic solid and a fluid") {
    const auto solid = rock();
    checkFace(solid,
              fluid,
              face,
              sampleState<model::ElasticMaterial>(0.4),
              sampleState<model::AcousticMaterial>(2.3));
    checkFace(fluid,
              solid,
              face,
              sampleState<model::AcousticMaterial>(1.1),
              sampleState<model::ElasticMaterial>(0.8));
  }

  SUBCASE("An anisotropic solid and a fluid") {
    const auto solid = triclinic();
    checkFace(solid,
              fluid,
              face,
              sampleState<model::AnisotropicMaterial>(0.4),
              sampleState<model::AcousticMaterial>(2.3));
    checkFace(fluid,
              solid,
              face,
              sampleState<model::AcousticMaterial>(1.1),
              sampleState<model::AnisotropicMaterial>(0.8));
  }

#ifdef SEISSOL_KERNELS_STP
  const Eigen::Vector3d normal = face.n();

  SUBCASE("A poroelastic and an anisotropic solid") {
    const auto solid = triclinic();
    const auto poro = sandstone();
    const auto poroState = sampleState<model::PoroElasticMaterial>(0.4);
    const auto solidState = sampleState<model::AnisotropicMaterial>(2.3);
    const auto fromPoro = checkFace(poro, solid, face, poroState, solidState);
    // from the other side, the same traction and velocity at the face
    const auto fromSolid = checkFace(solid, poro, face.flipped(), solidState, poroState);
    const auto [poroTraction, poroVelocity] =
        interfaceVariables<model::PoroElasticMaterial>(fromPoro, normal);
    const auto [solidTraction, solidVelocity] =
        interfaceVariables<model::AnisotropicMaterial>(fromSolid, normal);
    CHECK((poroTraction.head<3>() - solidTraction.head<3>()).norm() <= 1e-10 * poroTraction.norm());
    CHECK((poroVelocity.head<3>() - solidVelocity.head<3>()).norm() <= 1e-10 * poroVelocity.norm());
  }

  const auto checkPoroelasticAndFluid = [&](double conductance) {
    auto poro = sandstone();
    poro.interfaceConductance = conductance;
    const auto poroState = sampleState<model::PoroElasticMaterial>(0.4);
    const auto fluidState = sampleState<model::AcousticMaterial>(2.3);
    const auto fromPoro = checkFace(poro, fluid, face, poroState, fluidState);
    // from the side of the fluid: its pressure is the normal stress at the face, and it moves with
    // the frame of the poroelastic solid plus the filtration
    const auto fromFluid = stateAtFace(fluid, poro, face.flipped(), fluidState, poroState);
    const auto [traction, velocity] =
        interfaceVariables<model::PoroElasticMaterial>(fromPoro, normal);
    CHECK(fromFluid(0) == doctest::Approx(traction.head<3>().dot(normal)).epsilon(1e-10));
    CHECK(fromFluid.tail<3>().dot(normal) ==
          doctest::Approx(velocity.head<3>().dot(normal) + velocity(3)).epsilon(1e-10));
  };
  SUBCASE("A poroelastic solid and a fluid, open") { checkPoroelasticAndFluid(Unbounded); }
  SUBCASE("A poroelastic solid and a fluid, with a skin") { checkPoroelasticAndFluid(3e-8); }

  SUBCASE("Two poroelastic solids with a skin") {
    auto local = sandstone();
    auto neighbor = porousSandstone();
    local.interfaceConductance = 2e-7;
    neighbor.interfaceConductance = 1e-7;
    checkFace(local,
              neighbor,
              face,
              sampleState<model::PoroElasticMaterial>(0.4),
              sampleState<model::PoroElasticMaterial>(2.3));
  }
#endif
}

TEST_CASE("The conductances of both sides of a face add up in series" *
          doctest::test_suite("model")) {
  using namespace facecouplingtest;
  CHECK(model::seriesConductance(2, 3) == doctest::Approx(1.2));
  CHECK(model::seriesConductance(2, Unbounded) == 2);
  CHECK(model::seriesConductance(Unbounded, Unbounded) == Unbounded);
  CHECK(model::seriesConductance(0, Unbounded) == 0);

  auto poro = sandstone();
  poro.interfaceConductance = 4;
  // a solid without pores seals a face, a fluid takes no part in a skin
  CHECK(model::sideConductance(rock(), true) == 0);
  CHECK(model::sideConductance(water(), true) == Unbounded);
  CHECK(model::sideConductance(poro, true) == 4);
  // within a region, the pore pressure is continuous
  CHECK(model::sideConductance(poro, false) == Unbounded);
  // the default, if a material file leaves it out
  CHECK(sandstone().interfaceConductance == Unbounded);
}

} // namespace seissol::unit_test
