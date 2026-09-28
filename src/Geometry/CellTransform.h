// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_
#define SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_

#include "Common/Constants.h"
#include "Geometry/MeshDefinition.h"

#include <Eigen/Dense>
#include <Eigen/LU>
#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace seissol::geometry {

class MeshReader;

/**
 * Maps the reference cell onto a cell in physical space.
 *
 * Coordinates come in two flavours which must not be confused: reference coordinates live on the
 * reference cell, space coordinates live in the physical domain. Every member states which one it
 * expects; for a non-affine transform, passing the wrong one silently yields a wrong result.
 */
class CellTransform {
  public:
  using VectorEigenT = Eigen::Vector<double, Cell::Dim>;
  using MatrixEigenT = Eigen::Matrix<double, Cell::Dim, Cell::Dim>;

  using VectorT = std::array<double, Cell::Dim>;

  CellTransform() = default;
  virtual ~CellTransform();

  protected:
  // copying a transform through a base reference would slice it
  CellTransform(const CellTransform&) = default;
  CellTransform(CellTransform&&) = default;
  auto operator=(const CellTransform&) -> CellTransform& = default;
  auto operator=(CellTransform&&) -> CellTransform& = default;

  public:
  /// maps a reference coordinate to a space coordinate
  [[nodiscard]] virtual auto refToSpace(const VectorEigenT& input) const -> VectorEigenT = 0;

  /// the Jacobian d(space)/d(reference), evaluated at a **reference** coordinate
  [[nodiscard]] virtual auto refToSpaceJacobian(const VectorEigenT& input) const
      -> MatrixEigenT = 0;

  /// maps a space coordinate back to a reference coordinate
  [[nodiscard]] virtual auto spaceToRef(const VectorEigenT& input) const -> VectorEigenT;

  [[nodiscard]] auto refToSpace(const VectorT& input) const -> VectorT;

  [[nodiscard]] auto refToSpace(const std::vector<VectorT>& input) const -> std::vector<VectorT>;

  /**
   * Batched refToSpace writing into caller-provided storage; `output` must hold at least as many
   * entries as `input`. Preferable inside per-cell loops, where the allocating overload would
   * allocate once per cell.
   */
  void refToSpace(const VectorT* input, VectorT* output, std::size_t count) const;

  /**
   * The inverse Jacobian d(reference)/d(space), evaluated at a **reference** coordinate. Its rows
   * are grad(xi), grad(eta), grad(zeta).
   */
  [[nodiscard]] auto refToSpaceJacobianInverse(const VectorEigenT& input) const -> MatrixEigenT;

  /**
   * The inverse Jacobian d(reference)/d(space), evaluated at a **space** coordinate. Prefer
   * refToSpaceJacobianInverse where the reference coordinate is known: this overload has to invert
   * the transform first.
   */
  [[nodiscard]] auto spaceToRefJacobian(const VectorEigenT& input) const -> MatrixEigenT;
};

/**
 * A cell transform which is affine, i.e. its Jacobian does not depend on the evaluation point. This
 * is the transform of a straight-sided tetrahedron.
 *
 * The factorization needed by spaceToRef is computed on first use. Hence an instance must not be
 * shared between threads; construct one per thread instead (they are cheap to copy).
 */
class AffineTransform : public CellTransform {
  public:
  explicit AffineTransform(const std::array<CoordinateT, Cell::NumVertices>& vertices);

  explicit AffineTransform(const std::array<VectorEigenT, Cell::NumVertices>& vertices);

  static auto fromMeshCell(std::size_t id, const MeshReader& mesh) -> AffineTransform;

  [[nodiscard]] auto refToSpace(const VectorEigenT& input) const -> VectorEigenT override;

  [[nodiscard]] auto refToSpaceJacobian(const VectorEigenT& /*input*/) const
      -> MatrixEigenT override;

  [[nodiscard]] auto spaceToRef(const VectorEigenT& input) const -> VectorEigenT override;

  using CellTransform::refToSpace;

  /// the constant Jacobian d(space)/d(reference)
  [[nodiscard]] auto jacobian() const -> const MatrixEigenT& { return transform_; }

  /// the determinant of the Jacobian; six times the signed cell volume
  [[nodiscard]] auto determinant() const -> double { return determinant_; }

  private:
  [[nodiscard]] auto factorization() const -> const Eigen::PartialPivLU<MatrixEigenT>&;

  void setup(const std::array<VectorEigenT, Cell::NumVertices>& vertices);

  MatrixEigenT transform_;
  VectorEigenT offset_;
  double determinant_{};
  mutable std::optional<Eigen::PartialPivLU<MatrixEigenT>> factorization_;
};
} // namespace seissol::geometry
#endif // SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_
