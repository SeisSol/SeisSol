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
#include "Geometry/MeshDefinition.h"
#include "Model/Common.h"
#include "Model/Quantities.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>
#include <yateto.h>

/**
 * The Riemann problem at a face between a solid and a fluid, which pose it in different
 * quantities.
 *
 * Both sides meet in the traction on the face and in the velocity, in the coordinates of the face
 * (the normal first, pointing from the cell to its neighbor). A wave that leaves the face into a
 * half space carries the traction t and the velocity v of the half space to the state at the face
 * along its impedance Z, a symmetric 3x3 matrix: towards the cell, t* - t = Z (v* - v); towards the
 * neighbor, t* - t = -Z (v* - v). With W = (Z_cell + Z_neighbor)^-1, the state at the face is
 *
 *   v* = W (Z_cell v_cell - t_cell) + W (Z_neighbor v_neighbor + t_neighbor)
 *   t* = Z_neighbor W t_cell - Z_cell W Z_neighbor v_cell
 *        + Z_cell W t_neighbor + Z_cell W Z_neighbor v_neighbor.
 *
 * A fluid has no shear impedance, so that the shear traction at the face vanishes; the fluid only
 * shares the normal velocity with a solid. The flux of a cell reads the state at the face through
 * the traction and the velocity alone, so that the other stress components of a solid in the
 * Godunov state are those of the cell.
 */
namespace seissol::model {

/// The traction on a face and the velocity there, in the coordinates of the face.
using InterfaceMap = Eigen::Matrix<double, 6, 6>;

/// The impedance of a half space of the material `material` at a face, in the coordinates of the
/// face (for an anisotropic material, `material` has to be rotated into them already, cf.
/// getRotatedMaterialCoefficients). A viscoelastic or viscoacoustic material has the impedance of
/// the elastic or acoustic part that it poses its Riemann problem with.
template <typename MaterialT>
Eigen::Matrix3d faceImpedance(const MaterialT& material) {
  if constexpr (std::is_same_v<RiemannMaterialOf<MaterialT>, AnisotropicMaterial>) {
    return ImpedanceCompute<AnisotropicMaterial>::admittance(material).inverse();
  } else if constexpr (std::is_same_v<RiemannQuantitiesOf<MaterialT>, ElasticMaterial>) {
    const auto& elastic = static_cast<const ElasticMaterial&>(material);
    const double zp = std::sqrt(elastic.rho * (elastic.lambda + 2 * elastic.mu));
    const double zs = std::sqrt(elastic.rho * elastic.mu);
    return Eigen::Vector3d(zp, zs, zs).asDiagonal();
  } else {
    static_assert(std::is_same_v<RiemannQuantitiesOf<MaterialT>, AcousticMaterial>,
                  "the material is neither a solid nor a fluid");
    const auto& acoustic = static_cast<const AcousticMaterial&>(material);
    return Eigen::Vector3d(std::sqrt(acoustic.rho * acoustic.lambda), 0, 0).asDiagonal();
  }
}

/// The rows that take the traction on a face and the velocity out of the quantities of the Riemann
/// problem of the material `MaterialT`, in the coordinates of the face. A solid has the stress
/// tensor in Voigt order (xx, yy, zz, xy, yz, xz), whose traction on a face with the normal x is
/// (xx, xy, xz); a fluid the mean normal stress, whose traction is normal.
template <typename MaterialT>
Eigen::MatrixXd interfaceRows() {
  constexpr auto Groups = RiemannQuantitiesOf<MaterialT>::PrimaryGroups;
  constexpr auto Traction = roleOffset(Groups, FaceRole::Traction);
  constexpr auto Velocity = roleOffset(Groups, FaceRole::Velocity);
  Eigen::MatrixXd rows = Eigen::MatrixXd::Zero(6, totalExtent(Groups));
  if constexpr (roleKind(Groups, FaceRole::Traction) == QuantityKind::SymTensor2) {
    rows(0, Traction + 0) = 1;
    rows(1, Traction + 3) = 1;
    rows(2, Traction + 5) = 1;
  } else {
    static_assert(roleKind(Groups, FaceRole::Traction) == QuantityKind::Scalar);
    rows(0, Traction) = 1;
  }
  for (std::size_t i = 0; i < 3; ++i) {
    rows(3 + i, Velocity + i) = 1;
  }
  return rows;
}

/// The state at a face from the one of the cell (`local`) and the one of the neighbor
/// (`neighbor`), both as traction and velocity in the coordinates of the face.
struct InterfaceGodunov {
  InterfaceMap local;
  InterfaceMap neighbor;
};

/// The state at a face between half spaces of the impedances `zLocal` (the cell) and `zNeighbor`
/// (see the top of this file).
inline InterfaceGodunov interfaceGodunov(const Eigen::Matrix3d& zLocal,
                                         const Eigen::Matrix3d& zNeighbor) {
  const Eigen::Matrix3d w = (zLocal + zNeighbor).inverse();
  InterfaceGodunov godunov{InterfaceMap::Zero(), InterfaceMap::Zero()};
  godunov.local.block<3, 3>(0, 0) = zNeighbor * w;
  godunov.local.block<3, 3>(0, 3) = -zLocal * w * zNeighbor;
  godunov.local.block<3, 3>(3, 0) = -w;
  godunov.local.block<3, 3>(3, 3) = w * zLocal;
  godunov.neighbor.block<3, 3>(0, 0) = zLocal * w;
  godunov.neighbor.block<3, 3>(0, 3) = zLocal * w * zNeighbor;
  godunov.neighbor.block<3, 3>(3, 0) = w;
  godunov.neighbor.block<3, 3>(3, 3) = w * zNeighbor;
  return godunov;
}

/// The Godunov state at a face between a solid and a fluid, in the quantities of the Riemann
/// problems of both and in the coordinates of the face: `local` (cell x cell) from the state of the
/// cell, `neighbor` (cell x neighbor) from the one of the neighbor.
struct CoupledGodunov {
  Eigen::MatrixXd local;
  Eigen::MatrixXd neighbor;
};

/// The Godunov state at a face between the cell of the material `local` and the neighbor
/// `neighbor`, both given in the coordinates of the face (cf. getRotatedMaterialCoefficients).
template <typename MaterialT, typename NeighborT>
CoupledGodunov coupledGodunov(const MaterialT& local, const NeighborT& neighbor) {
  const auto godunov = interfaceGodunov(faceImpedance(local), faceImpedance(neighbor));
  const Eigen::MatrixXd localRows = interfaceRows<MaterialT>();
  const Eigen::MatrixXd neighborRows = interfaceRows<NeighborT>();
  const Eigen::MatrixXd embed = localRows.transpose();
  // the quantities of the cell besides the traction and the velocity stay as they are
  const Eigen::MatrixXd kept =
      Eigen::MatrixXd::Identity(localRows.cols(), localRows.cols()) - embed * localRows;
  return {embed * godunov.local * localRows + kept, embed * godunov.neighbor * neighborRows};
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

/// The contribution of the neighbor to the Godunov state of the cell at a face between a solid and
/// a fluid, in the global coordinates: the quantities of the Riemann problem of the neighbor (rows)
/// to those of the cell (columns), as the time integral of a neighbor, a (basis x quantity)
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
