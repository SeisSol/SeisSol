// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CellTransform.h"

#include "Common/Constants.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"

#include <array>
#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::geometry {

CellTransform::~CellTransform() = default;

auto CellTransform::refToSpace(const VectorEigenT& input) const -> VectorEigenT {
  return refToSpaceImpl(input);
}
auto CellTransform::refToSpaceJacobian(const VectorEigenT& input) const -> MatrixEigenT {
  return refToSpaceJacobianImpl(input);
}

auto CellTransform::refToSpace(const VectorT& input) const -> VectorT {
  const auto inputEigen = VectorEigenT(input.data());
  const auto outputEigen = refToSpaceImpl(inputEigen);
  VectorT output{};
  for (std::size_t i = 0; i < Cell::Dim; ++i) {
    output[i] = outputEigen(i);
  }
  return output;
}

auto CellTransform::refToSpace(const std::vector<VectorT>& input) const -> std::vector<VectorT> {
  std::vector<VectorT> output(input.size());
  for (std::size_t i = 0; i < output.size(); ++i) {
    output[i] = refToSpace(input[i]);
  }
  return output;
}

auto CellTransform::spaceToRefJacobian(const CellTransform::VectorEigenT& input) const
    -> MatrixEigenT {
  const auto refToSpaceJ = refToSpaceJacobian(spaceToRef(input));
  return refToSpaceJ.inverse();
}

auto CellTransform::spaceToRef(const CellTransform::VectorEigenT& input) const -> VectorEigenT {
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

AffineTransform::AffineTransform(const std::array<CoordinateT, Cell::NumVertices>& vertices) {
  offset_ = VectorEigenT(vertices[0].data());

  for (std::size_t i = 0; i < Cell::Dim; ++i) {
    const auto v = VectorEigenT(vertices[i + 1].data()) - offset_;
    for (std::size_t j = 0; j < Cell::Dim; ++j) {
      transform_(j, i) = v(j);
    }
  }

  itransform_ = transform_.inverse();
}

AffineTransform::AffineTransform(const std::array<VectorEigenT, Cell::NumVertices>& vertices) {
  offset_ = VectorEigenT(vertices[0]);

  for (std::size_t i = 0; i < Cell::Dim; ++i) {
    const auto v = VectorEigenT(vertices[i + 1]) - offset_;
    for (std::size_t j = 0; j < Cell::Dim; ++j) {
      transform_(j, i) = v(j);
    }
  }

  itransform_ = transform_.inverse();
}

auto AffineTransform::refToSpaceImpl(const VectorEigenT& input) const -> VectorEigenT {
  return transform_ * input + offset_;
}

auto AffineTransform::refToSpaceJacobianImpl(const VectorEigenT& /*input*/) const -> MatrixEigenT {
  // since we're linear—no dependency on the input vector here
  return transform_;
}

auto AffineTransform::fromMeshCell(std::size_t id, const MeshReader& mesh) -> AffineTransform {
  const auto& vertexIndices = mesh.getElements()[id].vertices;
  std::array<CoordinateT, Cell::NumVertices> vertices{};
  for (std::size_t i = 0; i < vertexIndices.size(); ++i) {
    vertices[i] = mesh.getVertices()[vertexIndices[i]].coords;
  }
  return AffineTransform(vertices);
}

auto AffineTransform::spaceToRef(const VectorEigenT& input) const -> VectorEigenT {
  // we can invert pretty straight-forwardly
  return itransform_ * (input - offset_);
}

} // namespace seissol::geometry
