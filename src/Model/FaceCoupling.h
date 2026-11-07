// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_FACECOUPLING_H_
#define SEISSOL_SRC_MODEL_FACECOUPLING_H_

#include "Equations/Datastructures.h"
#include "Equations/anisotropic/Model/Impedance.h"
#include "Equations/poroelastic/Model/Impedance.h"
#include "Geometry/MeshDefinition.h"
#include "Model/Common.h"
#include "Model/Quantities.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>
#include <yateto.h>

/**
 * The Riemann problem at a face between cells of different families (a solid, a fluid and a
 * poroelastic solid), which pose it in different quantities, and at a face between poroelastic
 * cells with a skin between them.
 *
 * All of them meet in the interface variables of the face, in the coordinates of the face (the
 * normal first, pointing from the cell to its neighbor): the traction T = (s_nn, s_ns, s_nd, -p)
 * and the velocity U = (v_n, v_s, v_d, q_n), with the pore pressure p and the normal filtration
 * velocity q_n of a poroelastic solid (-p q_n is the power that the pore fluid carries through the
 * face). A wave that leaves the face into a half space carries the traction T and the velocity U
 * of the half space to the state at the face along its impedance Z, a symmetric 4x4 matrix:
 * towards the cell, T* - T = Z (U* - U); towards the neighbor, T*' - T = -Z (U* - U). The velocity
 * U* at the face is the same on both sides, and so are the tractions T* and T*' but for the pore
 * pressure, which jumps across a skin of the hydraulic conductance k: T*' - T* = K U* with
 * K = e4 e4^T / k, so that the filtration through the face is k (p* - p*'). With
 * W = (Z_cell + Z_neighbor + K)^-1,
 *
 *   U* = W (Z_cell U_cell - T_cell) + W (Z_neighbor U_neighbor + T_neighbor)
 *   T* = T_cell + Z_cell (U* - U_cell).
 *
 * A solid without pores has neither a pore pressure nor a filtration: its impedance has an empty
 * fourth row and column, and its conductance is zero, a sealed face, at which q_n* vanishes and the
 * pore pressure of the other side is free. A fluid has no shear impedance, so that the shear
 * traction at the face vanishes; it moves with the frame plus the filtration of the other side,
 * v_n + q_n, and its pressure is both the normal stress and the pore pressure there: a fluid of the
 * pressure (mean stress) p' and the impedance Z_f has the traction T = p' c and the impedance
 * Z_f c c^T with c = (1, 0, 0, 1), and an unbounded conductance.
 *
 * The flux of a cell reads the state at the face through the interface variables alone, so the
 * other quantities in the Godunov state of a cell are those of the cell.
 */
namespace seissol::model {

/// The interface variables of a face, the traction T and the velocity U one after the other.
using InterfaceMap = Eigen::Matrix<double, 8, 8>;

/// The impedance of a half space in the interface variables.
using InterfaceImpedance = Eigen::Matrix4d;

/// The impedance of a half space of the material `material` at a face, in the coordinates of the
/// face (an anisotropic material has to be rotated into them already, cf.
/// getRotatedMaterialCoefficients). A viscoelastic or viscoacoustic material has the impedance of
/// the elastic or acoustic part that it poses its Riemann problem with.
template <typename MaterialT>
InterfaceImpedance faceImpedance(const MaterialT& material) {
  InterfaceImpedance impedance = InterfaceImpedance::Zero();
  if constexpr (std::is_same_v<RiemannMaterialOf<MaterialT>, AnisotropicMaterial>) {
    impedance.topLeftCorner<3, 3>() =
        ImpedanceCompute<AnisotropicMaterial>::admittance(material).inverse();
  } else if constexpr (std::is_same_v<RiemannQuantitiesOf<MaterialT>, ElasticMaterial>) {
    const auto& elastic = static_cast<const ElasticMaterial&>(material);
    impedance(0, 0) = std::sqrt(elastic.rho * (elastic.lambda + 2 * elastic.mu));
    impedance(1, 1) = std::sqrt(elastic.rho * elastic.mu);
    impedance(2, 2) = impedance(1, 1);
  } else if constexpr (std::is_same_v<RiemannQuantitiesOf<MaterialT>, AcousticMaterial>) {
    const auto& acoustic = static_cast<const AcousticMaterial&>(material);
    const Eigen::Vector4d direction(1, 0, 0, 1);
    impedance = std::sqrt(acoustic.rho * acoustic.lambda) * direction * direction.transpose();
  } else {
    static_assert(std::is_same_v<RiemannQuantitiesOf<MaterialT>, PoroElasticMaterial>,
                  "the material is neither a solid, nor a fluid, nor a poroelastic solid");
    impedance = ImpedanceCompute<PoroElasticMaterial>::impedance(material);
  }
  return impedance;
}

/// The hydraulic conductance of one side of a face of a cell of the material `material`: zero for
/// a solid without pores, unbounded for a fluid; for a poroelastic solid, its interface conductance
/// if the face bounds its region of the mesh (`regionBoundary`), and unbounded within its region.
template <typename MaterialT>
double sideConductance(const MaterialT& material, bool regionBoundary) {
  constexpr auto Unbounded = std::numeric_limits<double>::infinity();
  if constexpr (std::is_same_v<RiemannQuantitiesOf<MaterialT>, ElasticMaterial>) {
    return 0;
  } else if constexpr (std::is_same_v<RiemannQuantitiesOf<MaterialT>, AcousticMaterial>) {
    return Unbounded;
  } else {
    static_assert(std::is_same_v<RiemannQuantitiesOf<MaterialT>, PoroElasticMaterial>);
    return regionBoundary ? material.interfaceConductance : Unbounded;
  }
}

/// The conductance of the skins of both sides of a face, `first` and `second`, in series.
inline double seriesConductance(double first, double second) {
  if (first == 0 || second == 0) {
    return 0;
  }
  if (std::isinf(first)) {
    return second;
  }
  if (std::isinf(second)) {
    return first;
  }
  return first * second / (first + second);
}

/// The rows that take the interface variables out of the quantities of the Riemann problem of the
/// material `MaterialT`, in the coordinates of the face. A solid has the stress tensor in Voigt
/// order (xx, yy, zz, xy, yz, xz), whose traction on a face with the normal x is (xx, xy, xz); a
/// fluid the mean normal stress, its pressure (see the top of this file).
template <typename MaterialT>
Eigen::MatrixXd interfaceIn() {
  constexpr auto Groups = RiemannQuantitiesOf<MaterialT>::PrimaryGroups;
  constexpr auto Traction = roleOffset(Groups, FaceRole::Traction);
  constexpr auto Velocity = roleOffset(Groups, FaceRole::Velocity);
  Eigen::MatrixXd in = Eigen::MatrixXd::Zero(8, totalExtent(Groups));
  if constexpr (roleKind(Groups, FaceRole::Traction) == QuantityKind::SymTensor2) {
    in(0, Traction + 0) = 1;
    in(1, Traction + 3) = 1;
    in(2, Traction + 5) = 1;
  } else {
    static_assert(roleKind(Groups, FaceRole::Traction) == QuantityKind::Scalar);
    in(0, Traction) = 1;
    in(3, Traction) = 1;
  }
  for (std::size_t i = 0; i < 3; ++i) {
    in(4 + i, Velocity + i) = 1;
  }
  if constexpr (roleExtent(Groups, FaceRole::ExtraTraction) > 0) {
    // the pore pressure and the filtration velocity
    in(3, roleOffset(Groups, FaceRole::ExtraTraction)) = -1;
    in(7, roleOffset(Groups, FaceRole::ExtraVelocity)) = 1;
  }
  return in;
}

/// The quantities of the Riemann problem of the material `MaterialT` that the interface variables
/// at a face give, in the coordinates of the face: for a fluid, its normal velocity is the one of
/// the frame plus the filtration (see the top of this file).
template <typename MaterialT>
Eigen::MatrixXd interfaceOut() {
  constexpr auto Groups = RiemannQuantitiesOf<MaterialT>::PrimaryGroups;
  if constexpr (roleKind(Groups, FaceRole::Traction) == QuantityKind::Scalar) {
    constexpr auto Traction = roleOffset(Groups, FaceRole::Traction);
    constexpr auto Velocity = roleOffset(Groups, FaceRole::Velocity);
    Eigen::MatrixXd out = Eigen::MatrixXd::Zero(totalExtent(Groups), 8);
    out(Traction, 0) = 1;
    for (std::size_t i = 0; i < 3; ++i) {
      out(Velocity + i, 4 + i) = 1;
    }
    out(Velocity, 7) = 1;
    return out;
  } else {
    // the traction and the velocity as they are, the pore pressure with its sign
    return interfaceIn<MaterialT>().transpose();
  }
}

/// The state at a face from the one of the cell (`local`) and the one of the neighbor
/// (`neighbor`), both in the interface variables: the traction on the side of the cell and the
/// velocity.
struct InterfaceGodunov {
  InterfaceMap local;
  InterfaceMap neighbor;
};

/// The state at a face between half spaces of the impedances `zLocal` (the cell) and `zNeighbor`,
/// with a skin of the conductance `conductance` between them (see the top of this file).
inline InterfaceGodunov interfaceGodunov(const InterfaceImpedance& zLocal,
                                         const InterfaceImpedance& zNeighbor,
                                         double conductance) {
  const InterfaceImpedance sum = zLocal + zNeighbor;
  InterfaceImpedance w = InterfaceImpedance::Zero();
  if (conductance == 0) {
    // sealed: no filtration, and the pore pressures of both sides apart
    w.topLeftCorner<3, 3>() = sum.topLeftCorner<3, 3>().inverse();
  } else if (std::isinf(conductance)) {
    w = sum.inverse();
  } else {
    InterfaceImpedance resistance = sum;
    resistance(3, 3) += 1 / conductance;
    w = resistance.inverse();
  }
  const InterfaceImpedance identity = InterfaceImpedance::Identity();
  InterfaceGodunov godunov{.local = InterfaceMap::Zero(), .neighbor = InterfaceMap::Zero()};
  godunov.local.block<4, 4>(0, 0) = identity - zLocal * w;
  godunov.local.block<4, 4>(0, 4) = zLocal * w * zLocal - zLocal;
  godunov.local.block<4, 4>(4, 0) = -w;
  godunov.local.block<4, 4>(4, 4) = w * zLocal;
  godunov.neighbor.block<4, 4>(0, 0) = zLocal * w;
  godunov.neighbor.block<4, 4>(0, 4) = zLocal * w * zNeighbor;
  godunov.neighbor.block<4, 4>(4, 0) = w;
  godunov.neighbor.block<4, 4>(4, 4) = w * zNeighbor;
  return godunov;
}

/// The Godunov state at a face, in the quantities of the Riemann problems of the cell and its
/// neighbor and in the coordinates of the face: `local` (cell x cell) from the state of the cell,
/// `neighbor` (cell x neighbor) from the one of the neighbor.
struct CoupledGodunov {
  Eigen::MatrixXd local;
  Eigen::MatrixXd neighbor;
};

/// The Godunov state at a face between the cell of the material `local` and the neighbor
/// `neighbor`, both given in the coordinates of the face (cf. getRotatedMaterialCoefficients),
/// with a skin of the conductance `conductance` between them.
template <typename MaterialT, typename NeighborT>
CoupledGodunov
    coupledGodunov(const MaterialT& local, const NeighborT& neighbor, double conductance) {
  const auto godunov = interfaceGodunov(faceImpedance(local), faceImpedance(neighbor), conductance);
  const Eigen::MatrixXd in = interfaceIn<MaterialT>();
  const Eigen::MatrixXd out = interfaceOut<MaterialT>();
  // the quantities of the cell besides the interface variables stay as they are
  const Eigen::MatrixXd kept = Eigen::MatrixXd::Identity(in.cols(), in.cols()) - out * in;
  return {out * godunov.local * in + kept, out * godunov.neighbor * interfaceIn<NeighborT>()};
}

/// The rotation of the quantities of the groups `groups` from the coordinates of the face with the
/// normal `normal` and the tangents `tangent1` and `tangent2` into the global ones, or back if
/// `Inverse`: T and Tinv of getFaceRotationMatrix, as matrices that act on a column of quantities.
template <bool Inverse, std::size_t N>
Eigen::MatrixXd quantityRotation(const std::array<QuantityGroup, N>& groups,
                                 const CoordinateT& normal,
                                 const CoordinateT& tangent1,
                                 const CoordinateT& tangent2) {
  const auto count = static_cast<unsigned>(totalExtent(groups));
  std::vector<double> data(static_cast<std::size_t>(count) * count);
  yateto::DenseTensorView<2, double, unsigned> view(data.data(), {count, count});
  detail::writeRotationBlocks<Inverse>(groups, normal, tangent1, tangent2, view);
  return Eigen::Map<const Eigen::MatrixXd>(data.data(), count, count);
}

/// The contribution of the neighbor to the Godunov state of the cell at a face between different
/// families, in the global coordinates: the quantities of the Riemann problem of the neighbor
/// (rows) to those of the cell (columns), as the time integral of a neighbor, a (basis x quantity)
/// matrix, is multiplied with it from the right. The materials are given in the coordinates of the
/// face, which has the normal `normal` and the tangents `tangent1` and `tangent2`.
template <typename MaterialT, typename NeighborT>
Eigen::MatrixXd coupledNeighborMap(const CoupledGodunov& godunov,
                                   const CoordinateT& normal,
                                   const CoordinateT& tangent1,
                                   const CoordinateT& tangent2) {
  const auto toGlobal = quantityRotation<false>(
      RiemannQuantitiesOf<MaterialT>::PrimaryGroups, normal, tangent1, tangent2);
  const auto toFace = quantityRotation<true>(
      RiemannQuantitiesOf<NeighborT>::PrimaryGroups, normal, tangent1, tangent2);
  return (toGlobal * godunov.neighbor * toFace).transpose();
}

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_FACECOUPLING_H_
