// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_
#define SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_

#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
namespace seissol::geometry {
class CellTransform {
  public:
  using VectorEigenT = Eigen::Vector<double, Cell::Dim>;
  using VectorArrayEigenT = Eigen::Vector<double, Cell::Dim>;
  using MatrixEigenT = Eigen::Matrix<double, Cell::Dim, Cell::Dim>;

  using VectorT = std::array<double, Cell::Dim>;

  virtual ~CellTransform();

  [[nodiscard]] virtual auto refToSpaceImpl(const VectorEigenT& input) const -> VectorEigenT = 0;
  [[nodiscard]] virtual auto refToSpaceJacobianImpl(const VectorEigenT& input) const
      -> MatrixEigenT = 0;

  [[nodiscard]] auto refToSpace(const VectorEigenT& input) const -> VectorEigenT;
  [[nodiscard]] auto refToSpaceJacobian(const VectorEigenT& input) const -> MatrixEigenT;

  [[nodiscard]] auto refToSpace(const VectorT& input) const -> VectorT;

  [[nodiscard]] auto refToSpace(const std::vector<VectorT>& input) const -> std::vector<VectorT>;

  [[nodiscard]] auto spaceToRefJacobian(const VectorEigenT& input) const -> MatrixEigenT;

  [[nodiscard]] virtual auto spaceToRef(const VectorEigenT& input) const -> VectorEigenT;
};

class AffineTransform : public CellTransform {
  public:
  explicit AffineTransform(const std::array<CoordinateT, Cell::NumVertices>& vertices);

  explicit AffineTransform(const std::array<VectorEigenT, Cell::NumVertices>& vertices);

  [[nodiscard]] auto refToSpaceImpl(const VectorEigenT& input) const -> VectorEigenT override;

  [[nodiscard]] auto refToSpaceJacobianImpl(const VectorEigenT& /*input*/) const
      -> MatrixEigenT override;

  static auto fromMeshCell(std::size_t id, const MeshReader& mesh) -> AffineTransform;

  [[nodiscard]] auto spaceToRef(const VectorEigenT& input) const -> VectorEigenT override;

  private:
  MatrixEigenT transform_;
  Eigen::PartialPivLU<MatrixEigenT> itransform_;
  VectorEigenT offset_;
  double determinant_{};
};
} // namespace seissol::geometry
#endif // SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_
