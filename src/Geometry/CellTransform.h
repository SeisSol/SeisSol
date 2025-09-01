// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_
#define SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_

#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshReader.h>
namespace seissol::geometry {
class CellTransform {
  public:
  using VectorEigenT = Eigen::Vector<double, Cell::Dim>;
  using VectorArrayEigenT = Eigen::Vector<double, Cell::Dim>;
  using MatrixEigenT = Eigen::Matrix<double, Cell::Dim, Cell::Dim>;

  using VectorT = std::array<double, Cell::Dim>;

  virtual ~CellTransform() = default;

  [[nodiscard]] virtual VectorEigenT refToSpaceImpl(const VectorEigenT& input) const = 0;
  [[nodiscard]] virtual MatrixEigenT refToSpaceJacobianImpl(const VectorEigenT& input) const = 0;

  [[nodiscard]] VectorEigenT refToSpace(const VectorEigenT& input) const {
    return refToSpaceImpl(input);
  }
  [[nodiscard]] MatrixEigenT refToSpaceJacobian(const VectorEigenT& input) const {
    return refToSpaceJacobianImpl(input);
  }

  [[nodiscard]] VectorT refToSpace(const VectorT& input) const {
    const auto inputEigen = VectorEigenT(input.data());
    const auto outputEigen = refToSpaceImpl(inputEigen);
    VectorT output{};
    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      output[i] = outputEigen(i);
    }
    return output;
  }

  [[nodiscard]] std::vector<VectorT> refToSpace(const std::vector<VectorT>& input) const {
    std::vector<VectorT> output(input.size());
    for (std::size_t i = 0; i < output.size(); ++i) {
      output[i] = refToSpace(input[i]);
    }
    return output;
  }

  [[nodiscard]] MatrixEigenT spaceToRefJacobian(const VectorEigenT& input) const {
    const auto refToSpaceJ = refToSpaceJacobian(spaceToRef(input));
    return refToSpaceJ.inverse();
  }

  [[nodiscard]] virtual VectorEigenT spaceToRef(const VectorEigenT& input) const {
    // in the general case... we need to invert a function. So... Newton.
    // we want: f(y) = x; or: f(y) - x = 0

    auto iterate = VectorEigenT();
    constexpr double Eps = 1e-8;
    constexpr std::size_t Tries = 100000;
    for (std::size_t i = 0; i < Tries; ++i) {
      const auto inputProbe = refToSpace(iterate);
      if ((inputProbe - input).norm() < Eps) {
        return iterate;
      }
      const auto inputProbeDerivative = refToSpaceJacobian(iterate);
      iterate -= inputProbeDerivative.fullPivLu().solve(inputProbe);
    }

    logError() << "Root finding failed for" << input << "after" << Tries
               << "iterations. Last iterate:" << iterate;
    return iterate;
  }
};

class AffineTransform : public CellTransform {
  public:
  AffineTransform(const std::array<CoordinateT, Cell::NumVertices>& vertices) {
    offset = VectorEigenT(vertices[0].data());

    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      const auto v = VectorEigenT(vertices[i + 1].data()) - offset;
      for (std::size_t j = 0; j < Cell::Dim; ++j) {
        transform(j, i) = v(j);
      }
    }

    itransform = transform.inverse();
  }

  AffineTransform(const std::array<VectorEigenT, Cell::NumVertices>& vertices) {
    offset = VectorEigenT(vertices[0]);

    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      const auto v = VectorEigenT(vertices[i + 1]) - offset;
      for (std::size_t j = 0; j < Cell::Dim; ++j) {
        transform(j, i) = v(j);
      }
    }

    itransform = transform.inverse();
  }

  [[nodiscard]] VectorEigenT refToSpaceImpl(const VectorEigenT& input) const override {
    return transform * input + offset;
  }

  [[nodiscard]] MatrixEigenT refToSpaceJacobianImpl(const VectorEigenT& input) const override {
    // since we're linear—no dependency on the input vector here
    return transform;
  }

  static AffineTransform fromMeshCell(std::size_t id, const MeshReader& mesh) {
    const auto& vertexIndices = mesh.getElements()[id].vertices;
    std::array<CoordinateT, Cell::NumVertices> vertices;
    for (std::size_t i = 0; i < vertexIndices.size(); ++i) {
      vertices[i] = mesh.getVertices()[vertexIndices[i]].coords;
    }
    return AffineTransform(vertices);
  }

  [[nodiscard]] VectorEigenT spaceToRef(const VectorEigenT& input) const override {
    // we can invert pretty straight-forwardly
    return itransform * (input - offset);
  }

  private:
  MatrixEigenT transform;
  MatrixEigenT itransform;
  VectorEigenT offset;
};
} // namespace seissol::geometry
#endif // SEISSOL_SRC_GEOMETRY_CELLTRANSFORM_H_
