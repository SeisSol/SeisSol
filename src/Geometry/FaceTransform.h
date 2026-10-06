// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_FACETRANSFORM_H_
#define SEISSOL_SRC_GEOMETRY_FACETRANSFORM_H_

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/MeshDefinition.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol::geometry {

class MeshReader;

/**
 * How a face is oriented with respect to the cell it is approached from.
 *
 * A face is shared by two cells which generally enumerate its vertices differently. Local denotes
 * the enumeration of the cell owning the face; the rotations denote the enumeration seen from the
 * neighbor, and correspond to the side orientation stored per element and face.
 */
enum class FaceOrientation : std::int8_t { Local = -1, Rotate0 = 0, Rotate1 = 1, Rotate2 = 2 };

/**
 * Embeds the reference face into the reference cell.
 *
 * The map is affine and depends on nothing but the side and the orientation, i.e. it carries no
 * mesh information whatsoever. There are Cell::NumFaces * (Face::NumOrientations + 1) instances in
 * total, all of which are known at compile time.
 */
class ReferenceFaceMap {
  public:
  using FaceVectorT = Eigen::Vector<double, Face::Dim>;
  using CellVectorT = Eigen::Vector<double, Cell::Dim>;
  using EmbeddingT = Eigen::Matrix<double, Cell::Dim, Face::Dim>;
  using ProjectionT = Eigen::Matrix<double, Face::Dim, Cell::Dim>;

  explicit ReferenceFaceMap(std::size_t side, FaceOrientation orientation = FaceOrientation::Local);

  /// maps a face coordinate (chi, tau) to a cell coordinate (xi, eta, zeta)
  [[nodiscard]] auto faceToCell(const FaceVectorT& input) const -> CellVectorT;

  /// maps a cell coordinate on this face back to a face coordinate
  [[nodiscard]] auto cellToFace(const CellVectorT& input) const -> FaceVectorT;

  /// the constant Jacobian d(cell)/d(face)
  [[nodiscard]] auto embedding() const -> const EmbeddingT& { return embedding_; }

  [[nodiscard]] auto side() const -> std::size_t { return side_; }

  [[nodiscard]] auto orientation() const -> FaceOrientation { return orientation_; }

  private:
  EmbeddingT embedding_;
  CellVectorT offset_;
  ProjectionT projection_;
  FaceVectorT projectionOffset_;
  std::size_t side_;
  FaceOrientation orientation_;
};

/**
 * Maps the reference face onto a face in physical space.
 *
 * A face transform is the composition of a ReferenceFaceMap with the transform of one of the two
 * adjacent cells. Deriving it that way keeps cell and face geometry consistent by construction; in
 * particular, a curved cell yields a curved face rather than a flat one inscribed into it.
 */
class FaceTransform {
  public:
  using FaceVectorT = Eigen::Vector<double, Face::Dim>;
  using VectorT = Eigen::Vector<double, Cell::Dim>;
  using JacobianT = Eigen::Matrix<double, Cell::Dim, Face::Dim>;

  FaceTransform() = default;
  virtual ~FaceTransform();

  protected:
  FaceTransform(const FaceTransform&) = default;
  FaceTransform(FaceTransform&&) = default;
  auto operator=(const FaceTransform&) -> FaceTransform& = default;
  auto operator=(FaceTransform&&) -> FaceTransform& = default;

  public:
  /// maps a face coordinate to a space coordinate
  [[nodiscard]] virtual auto refToSpace(const FaceVectorT& input) const -> VectorT = 0;

  /// maps a face coordinate to the reference coordinate of the adjacent cell
  [[nodiscard]] virtual auto refToCell(const FaceVectorT& input) const -> VectorT = 0;

  /// the Jacobian d(space)/d(face), evaluated at a **face** coordinate
  [[nodiscard]] virtual auto refToSpaceJacobian(const FaceVectorT& input) const -> JacobianT = 0;

  /**
   * The face normal, evaluated at a **face** coordinate. It points away from the cell the transform
   * was derived from, and is not normalized; its norm is the surface Jacobian.
   */
  [[nodiscard]] auto normal(const FaceVectorT& input) const -> VectorT;

  /// the norm of the (unnormalized) normal; for a straight-sided face, twice its area
  [[nodiscard]] auto surfaceJacobian(const FaceVectorT& input) const -> double;

  /// the area of the face; exact only for a straight-sided face
  [[nodiscard]] auto area() const -> double;

  /// the barycenter of the face, in space coordinates
  [[nodiscard]] auto center() const -> VectorT;

  /// the vertices of the face, in space coordinates and in the order of the underlying map
  [[nodiscard]] auto vertices() const -> std::array<VectorT, Face::NumVertices>;

  /**
   * The orthogonal basis (normal, tangent1, tangent2) used for face-aligned rotation matrices.
   *
   * This is deliberately not the parametric basis: tangent1 follows the first face edge and
   * tangent2 completes the right-handed system, which is what the rotation matrices of the
   * material models are defined against.
   */
  [[nodiscard]] auto faceAlignedBasis() const -> std::array<VectorT, Cell::Dim>;

  /**
   * The same basis at a point of the face, given as a **face** coordinate.
   *
   * On a curved face the normal turns from point to point, so the basis does too: the normal is
   * the one at the point, tangent1 is the first face edge projected onto the plane the face has
   * there, and tangent2 completes the right-handed system. Like the basis above it is not
   * normalized. Where the face is straight, it is that basis wherever it is evaluated.
   */
  [[nodiscard]] auto faceAlignedBasis(const FaceVectorT& input) const
      -> std::array<VectorT, Cell::Dim>;
};

/// A face transform derived from an AffineTransform; its Jacobian is constant.
class AffineFaceTransform : public FaceTransform {
  public:
  AffineFaceTransform(const AffineTransform& cell, const ReferenceFaceMap& embedding);

  /// the face of a mesh cell, as seen from that cell
  static auto fromMeshCell(std::size_t id,
                           std::size_t side,
                           const MeshReader& mesh,
                           FaceOrientation orientation = FaceOrientation::Local)
      -> AffineFaceTransform;

  /// the face carrying a fault, as seen from whichever adjacent cell is on this rank
  static auto fromMeshFault(std::size_t faultId, const MeshReader& mesh) -> AffineFaceTransform;

  [[nodiscard]] auto refToSpace(const FaceVectorT& input) const -> VectorT override;

  [[nodiscard]] auto refToCell(const FaceVectorT& input) const -> VectorT override;

  [[nodiscard]] auto refToSpaceJacobian(const FaceVectorT& /*input*/) const -> JacobianT override;

  [[nodiscard]] auto map() const -> const ReferenceFaceMap& { return map_; }

  private:
  ReferenceFaceMap map_;
  AffineTransform cell_;
  JacobianT transform_;
  VectorT offset_;
};

} // namespace seissol::geometry
#endif // SEISSOL_SRC_GEOMETRY_FACETRANSFORM_H_
