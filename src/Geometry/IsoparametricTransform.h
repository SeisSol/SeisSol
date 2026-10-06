// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_ISOPARAMETRICTRANSFORM_H_
#define SEISSOL_SRC_GEOMETRY_ISOPARAMETRICTRANSFORM_H_

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <vector>

namespace seissol::geometry {

/**
 * A cell given by the Lagrange interpolant through its nodes: x(xi) = sum_i x_i L_i(xi).
 *
 * The nodes sit on the equispaced lattice of the reference cell, in the order latticeNodes()
 * returns. Order one is the affine map through the vertices; order two adds the six edge
 * midpoints, which is the isoparametric P2 cell a mesher hands out. The map is a polynomial of
 * the given order, so its Jacobian is one order lower and varies inside the cell from order two
 * on; the Jacobian determinant has the degree 3 * (order - 1), and the cofactor matrix
 * det(J) J^-1 the degree 2 * (order - 1).
 */
class IsoparametricTransform : public CellTransform {
  public:
  /// `nodes` in space coordinates, in the order of latticeNodes(order)
  IsoparametricTransform(std::size_t order, const std::vector<VectorEigenT>& nodes);

  /**
   * The lattice the nodes of a cell of `order` sit on, in reference coordinates: the four
   * vertices first, in the order of the reference cell, then everything else ordered by the
   * position on the lattice.
   */
  static auto latticeNodes(std::size_t order) -> std::vector<VectorEigenT>;

  /// a straight-sided cell, written with the nodes of `order`
  static auto fromVertices(const std::array<VectorEigenT, Cell::NumVertices>& vertices,
                           std::size_t order) -> IsoparametricTransform;

  /**
   * The P2 cell through the vertices and the midpoints of its edges, the edges ordered as
   * (0,1), (0,2), (0,3), (1,2), (1,3), (2,3) by their vertices.
   */
  static auto fromEdgeMidpoints(const std::array<VectorEigenT, Cell::NumVertices>& vertices,
                                const std::array<VectorEigenT, 6>& midpoints)
      -> IsoparametricTransform;

  [[nodiscard]] auto refToSpace(const VectorEigenT& input) const -> VectorEigenT override;
  [[nodiscard]] auto refToSpaceJacobian(const VectorEigenT& input) const -> MatrixEigenT override;
  using CellTransform::refToSpace;

  [[nodiscard]] auto order() const -> std::size_t { return order_; }

  /// the nodes the cell was given, in space coordinates
  [[nodiscard]] auto nodes() const -> const std::vector<VectorEigenT>& { return nodes_; }

  private:
  std::size_t order_;
  std::vector<VectorEigenT> nodes_;
  std::vector<std::array<int, Cell::Dim>> exponents_;
  // x(xi) = origin_ + coefficients_ * m(xi), m the monomials of exponents_; the first node is
  // taken out before the fit, so that a cell far away from the origin keeps its digits
  VectorEigenT origin_;
  Eigen::Matrix<double, Cell::Dim, Eigen::Dynamic> coefficients_;
};

/**
 * A face of an isoparametric cell. It is the cell's map restricted to the face, so it is curved
 * where the cell is, and two cells sharing a face give the same surface as long as they agree on
 * the nodes on it -- which a conforming mesh guarantees.
 */
class IsoparametricFaceTransform : public FaceTransform {
  public:
  IsoparametricFaceTransform(const IsoparametricTransform& cell, const ReferenceFaceMap& embedding);

  [[nodiscard]] auto refToSpace(const FaceVectorT& input) const -> VectorT override;
  [[nodiscard]] auto refToCell(const FaceVectorT& input) const -> VectorT override;
  [[nodiscard]] auto refToSpaceJacobian(const FaceVectorT& input) const -> JacobianT override;

  [[nodiscard]] auto map() const -> const ReferenceFaceMap& { return map_; }

  private:
  ReferenceFaceMap map_;
  IsoparametricTransform cell_;
};

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_ISOPARAMETRICTRANSFORM_H_
